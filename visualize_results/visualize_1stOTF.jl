using DifferentialEquations
using FFTW
using LinearAlgebra
using CairoMakie
using Printf
using SpecialFunctions: erf  # Gaussian error function
# Krylov.jl is needed for the GMRES linear solver
using Krylov

# --- 1. Parameters and Grid Setup ---

# Physical Parameters
const L = 20.0       # Domain width in x [-L/2, L/2]
const H = 1.0        # Domain height in z [0, H]
const d_ekman = 0.1  # Ekman layer depth
Ro = 1.0           # Rossby number
ε = 0.05            # Perturbation parameter

# Numerical Parameters
Nx = 264       # Number of Fourier modes in x
Nz = 64        # Number of grid points in z
tmax = 20.0    # Maximum simulation time

# Parameters for de-aliasing (2/3 rule)
Nx_padded = Int(floor(3/2 * Nx))

# Grid Setup
x = L / Nx * (-Nx/2 : Nx/2 - 1); k = rfftfreq(Nx, Nx * 2π / L);
x_padded = L / Nx_padded * (-Nx_padded/2 : Nx_padded/2 - 1); k_padded = rfftfreq(Nx_padded, Nx_padded * 2π / L);
z = range(0, H, length=Nz); dz = z[2] - z[1];

# --- 2. Finite Difference Matrices (z-derivatives) ---

# First Derivative (D1) with one-sided boundary conditions
D1 = Tridiagonal(-ones(Nz-1), zeros(Nz), ones(Nz-1)) / (2*dz)
D1[1, 1:2] .= [-1, 1] / dz
D1[Nz, Nz-1:Nz] .= [-1, 1] / dz

# Second Derivative (D2) with Neumann BC at z=1, Dirichlet at z=0 for psi
# This matrix will be used to solve the Poisson equation: D2 * psi = omega
# BCs for psi: psi(z=0)=0, psi_z(z=1) approx 0 (from u1_z=0 -> psi_zz=0)
# We incorporate psi(0)=0 and psi(1)=0 directly.
D2 = Tridiagonal(ones(Nz-1), -2*ones(Nz), ones(Nz-1)) / dz^2
D2[1, :] .= 0; D2[1,1] = 1 # For psi(z=0)=0 -> omega is not defined here, but psi is.
D2[Nz, :] .= 0; D2[Nz,Nz] = 1 # For psi(z=1)=0
# Pre-calculate inverse for efficiency. Add a small regularization for stability.
Lz = inv(D2 + 1e-12 * I) 

# --- 3. Background State ---
println("Setting up exact background state...")

# Base functions
B0(x_val) = 0.5 * erf(x_val / sqrt(2))
B0_prime(x_val) = 1/sqrt(2*pi) * exp(-x_val^2 / 2)
B0_double_prime(x_val) = -x_val/sqrt(2*pi) * exp(-x_val^2 / 2)

# Create storage for all background fields
v0_exact = zeros(Nz, Nx)
b0_exact = zeros(Nz, Nx)
v0_x_exact = zeros(Nz, Nx)
v0_z_exact = zeros(Nz, Nx)
b0_x_exact = zeros(Nz, Nx)
# b0_z is also needed, though it will be small
b0_z_exact = zeros(Nz, Nx)
v0_zz_exact = zeros(Nz, Nx) # NOTE: d2v0/dz2 requires more complex implicit diff. We assume it's small.
b0_zz_exact = zeros(Nz, Nx) # NOTE: Assume small for simplicity of this script.

# Create storage for padded background fields
v0_padded = zeros(Nz, Nx_padded)
v0_x_padded = zeros(Nz, Nx_padded)
v0_z_padded = zeros(Nz, Nx_padded)

# Newton's method on the padded grid
for i in 1:Nx_padded, j in 1:Nz
    xi, zj = x_padded[i], z[j]
    f(v) = v - B0_prime(xi + v) * (zj - 0.5); f_prime(v) = 1 - B0_double_prime(xi + v) * (zj - 0.5);
    v_guess = 0.0
    for iter in 1:30; f_prime_val = f_prime(v_guess); if abs(f(v_guess) / f_prime(v_guess)) > 1e-10 v_guess -= f(v_guess) / f_prime_val end end
    v0_padded[j, i] = v_guess
    
    X = xi + v_guess
    J = 1 - B0_double_prime(X) * (zj - 0.5); if abs(J) < 1e-6 J = 1e-6 end
    dv0_dx = (B0_double_prime(X) * (zj - 0.5)) / J
    dv0_dz = B0_prime(X) / J
    v0_x_padded[j, i] = dv0_dx
    v0_z_padded[j, i] = dv0_dz
end

b0_x_padded = B0_prime.(x_padded' .+ v0_padded) .* (1 .+ v0_x_padded)
b0_z_padded = B0_prime.(x_padded' .+ v0_padded) .* v0_z_padded
# Assume v0_zz and b0_zz are negligible
v0_zz_padded = zeros(Nz, Nx_padded)
b0_zz_padded = zeros(Nz, Nx_padded)

# # First, v0_zz
# dfdz = B0_double_prime(X) * dv0_dz
# dgdz = - (B0_triple_prime(X) * dv0_dz * (zj - 0.5) + B0_double_prime(X))
# dv0_dzz = (dfdz * J - B0_prime(X) * dgdz) / J^2
# v0_zz_exact[j, i] = dv0_dzz

# # Second, b0_zz (using the result for v0_zz)
# db0_dzz = B0_double_prime(X) * (dv0_dz^2) + B0_prime(X) * dv0_dzz
# b0_zz_exact[j, i] = db0_dzz

println("Background state setup complete.")

# Ekman flow
UE = @. -0*exp((z-1)/d_ekman) * sin((z-1)/d_ekman - π/4)

# --- 4. ODE Functions for IMEX Scheme ---

# REVISION: Add a hyperviscosity parameter
const nu_hyper = 1e-5 # This is a tuning parameter. Start small.

# f1: Stiff, linear part to be handled implicitly
function imex_f_implicit!(du_flat, u_flat, p, t)
    # REVISION: Reshape the incoming flat vectors into our logical 3D array structure
    dims = p.implicit_params.dims
    du_hat = reshape(du_flat, dims..., 3)
    u_hat = reshape(u_flat, dims..., 3)

    k, D1, Lz, nu_hyper = p.implicit_params.operators

    omega_hat = @view u_hat[:, :, 1]; v1_hat = @view u_hat[:, :, 2]; b1_hat = @view u_hat[:, :, 3]
    psi_hat = Lz * omega_hat; u1_hat = D1 * psi_hat; 
    d_omega_hat_res = D1 * v1_hat; d_v1_hat_res = -u1_hat; d_b1_hat_res = zeros(ComplexF64, size(b1_hat))
    
    for ki in 1:dims[2]
        k_val = k[ki]
        # Linear dynamics
        @. d_omega_hat_res[:, ki] -= im * k_val * b1_hat[:, ki]
        
        # REVISION: Add hyperviscosity/hyperdiffusion to the stiff part
        # This term is -nu_hyper * k^4 * field
        dissipation_filter = -nu_hyper * k_val^4
        @. d_omega_hat_res[:, ki] += dissipation_filter * omega_hat[:, ki]
        @. d_v1_hat_res[:, ki]    += dissipation_filter * v1_hat[:, ki]
        @. d_b1_hat_res[:, ki]    += dissipation_filter * b1_hat[:, ki]
    end
    
    du_hat[:, :, 1] .= d_omega_hat_res
    du_hat[:, :, 2] .= d_v1_hat_res
    du_hat[:, :, 3] .= d_b1_hat_res
    du_hat[1, :, 1] .= 0; du_hat[Nz, :, 1] .= 0
    return nothing
end

# f2: Non-stiff, forcing part to be handled explicitly
function imex_f_explicit!(du_flat, u_flat, p, t)
    # REVISION: Reshape the incoming flat vectors
    dims = p.implicit_params.dims
    du_hat = reshape(du_flat, dims..., 3)
    u_hat = reshape(u_flat, dims..., 3)

    # Unpack full parameters for forcing calculation
    k, D1, Lz, x_padded, z, UE, v0_pad, v0_x_pad, v0_z_pad, v0_zz_pad, b0_x_pad, b0_z_pad, b0_zz_pad, Nx, Nx_padded, k_padded = p.explicit_params
    
    omega_hat = @view u_hat[:, :, 1]
    psi_hat   = Lz * omega_hat
    
    # --- DE-ALIASING PROCEDURE (still good practice) ---
    psi_hat_padded = zeros(ComplexF64, Nz, length(k_padded))
    psi_hat_padded[:, 1:size(psi_hat, 2)] .= psi_hat
    u1_phys_padded = D1 * irfft(psi_hat_padded, Nx_padded, 2)
    w1_phys_padded = -irfft(im * k_padded' .* psi_hat_padded, Nx_padded, 2)
    Fv_phys_padded = -v0_pad .+ v0_zz_pad .- (u1_phys_padded .- x_padded' .+ UE) .* v0_x_pad .+ w1_phys_padded .* v0_z_pad
    Fb_phys_padded = b0_zz_pad .- (u1_phys_padded .- x_padded' .+ UE) .* b0_x_pad .+ w1_phys_padded .* b0_z_pad
    Fv_hat_padded = rfft(Fv_phys_padded, 2)
    Fb_hat_padded = rfft(Fb_phys_padded, 2)
    Fv_hat = Fv_hat_padded[:, 1:size(u_hat, 2)]
    Fb_hat = Fb_hat_padded[:, 1:size(u_hat, 2)]
    # --- END DE-ALIASING ---

    d_omega_hat = @view du_hat[:, :, 1]
    d_v1_hat    = @view du_hat[:, :, 2]
    d_b1_hat    = @view du_hat[:, :, 3]
    
    # Non-stiff forcing
    d_omega_hat .= 0
    d_v1_hat    .= Fv_hat
    d_b1_hat    .= Fb_hat

    return nothing
end

# --- 5. Solve the ODE System using an IMEX Solver ---

# The initial condition must be a flat vector
u0_3D = zeros(ComplexF64, Nz, length(k), 3)
u0 = vec(u0_3D)

# Pass dimensions into the parameter tuple
params = (
    implicit_params = (operators = (k, D1, Lz, nu_hyper), dims=(Nz, length(k))),
    explicit_params = (k, D1, Lz, x_padded, z, UE, v0_padded, v0_x_padded, v0_z_padded, v0_zz_padded, b0_x_padded, b0_z_padded, b0_zz_padded, Nx, Nx_padded, k_padded)
);
tspan = (0.0, tmax);

function analytical_jvp!(Jv_flat, v_flat, u_flat, p, t)
    imex_f_implicit!(Jv_flat, v_flat, p, t)
end

f_implicit = ODEFunction(imex_f_implicit!; jvp=analytical_jvp!)
f_explicit = ODEFunction(imex_f_explicit!)
prob = SplitODEProblem(f_implicit, f_explicit, u0, tspan, params)

condition(u, t, integrator) = norm(u) > 1e8
action!(integrator) = (println("Solution norm > 1e8 at t=$(integrator.t). Terminating."); terminate!(integrator))
cb = DiscreteCallback(condition, action!)

println("Solving ODE system with IMEX + Hyperviscosity...")
sol = solve(prob, KenCarp4(linsolve=KrylovJL_GMRES()), reltol=1e-7, abstol=1e-7, saveat=0.2, callback=cb)
println("Solution complete. Final time: $(sol.t[end])")

# --- 6. Post-processing and Visualization ---
println("Reconstructing and animating solution...")

# Helper function to reconstruct all fields at a given time index
function reconstruct_fields(sol_u_flat, p)
    dims = p.implicit_params.dims
    operators = p.implicit_params.operators
    _, _, Lz = operators
    Nx = p.explicit_params[end-2]

    # REVISION: Reshape the solution vector before processing
    sol_u = reshape(sol_u_flat, dims..., 3)

    
    omega_hat = sol_u[:, :, 1]
    v1_hat = sol_u[:, :, 2]
    b1_hat = sol_u[:, :, 3]
    
    psi_hat = Lz * omega_hat
    
    psi1_phys = irfft(psi_hat, Nx, 2)
    v1_phys = irfft(v1_hat, Nx, 2)
    b1_phys = irfft(b1_hat, Nx, 2)
    
    return psi1_phys, v1_phys, b1_phys
end

b⁰ = Matrix{Float64}(undef, Nz, Nx)
for i in 1:Nx, j in 1:Nz
    xi, zj = x[i], z[j]
    f(v) = v - B0_prime(xi + v) * (zj - 0.5); f_prime(v) = 1 - B0_double_prime(xi + v) * (zj - 0.5);
    v_guess = 0.0
    for iter in 1:30; f_prime_val = f_prime(v_guess); if abs(f(v_guess) / f_prime(v_guess)) > 1e-10 v_guess -= f(v_guess) / f_prime_val end end 
    X = xi + v_guess
    b⁰[j,i] = B0(X) 
end

# Create figure and axis
fig = Figure(size = (640, 640))
ax_b = Axis(fig[1, 1]; ylabel = L"z", titlealign = :left, title = L"\text{(a) 1st-order buoyancy}~b^1")
ax_v = Axis(fig[2, 1]; ylabel = L"z", titlealign = :left, title = L"\text{(b) 1st-order velocity}~v^1")
ax_ψ = Axis(fig[3, 1]; xlabel = L"x", ylabel = L"z", titlealign = :left, title = L"\text{(c) 1st-order streamfunction}~\psi^1")
hidexdecorations!(ax_b, ticks = false)
hidexdecorations!(ax_v, ticks = false)

# Initialize the heatmap with the first frame
t = sol.t[1]
ψ1, v1, b1 = reconstruct_fields(sol.u[1], params)

# Create the heatmap
hm_b = heatmap!(ax_b, x, z, b1', colorrange = (-5,5), colormap = :diff)
hm_v = heatmap!(ax_v, x, z, v1', colorrange = (-2,2), colormap = :delta)
hm_ψ = heatmap!(ax_ψ, x, z, ψ1', colorrange = (-0.2, 0.2), colormap = :PuOr)
ct_b = contour!(ax_b, x, z, ε * b1' .+ b⁰'; levels = 10, linewidth = 1, color = :black)
ct_v = contour!(ax_v, x, z, ε * b1' .+ b⁰'; levels = 10, linewidth = 1, color = :black)
ct_ψ = contour!(ax_ψ, x, z, ε * b1' .+ b⁰'; levels = 10, linewidth = 1, color = :black)
# Add colorbar
Colorbar(fig[1, 2], hm_b)
Colorbar(fig[2, 2], hm_v)
Colorbar(fig[3, 2], hm_ψ)

# Add title
title_obs = Observable(@sprintf("t=%.2f", t))
fig[0, :] = Label(fig, title_obs, fontsize = 20)
colgap!(fig.layout, 1, Relative(0))
for i = 1:3
    rowgap!(fig.layout, i, Relative(0))
end
resize_to_layout!(fig)

# Create animation
record(fig, "b1v1psi1_evolution.gif", 1:length(sol.t); framerate = 10) do i
    t = sol.t[i]
    ψ1, v1, b1 = reconstruct_fields(sol.u[i], params)
    
    # Update the heatmap data
    hm_b[3] = b1'
    hm_v[3] = v1'
    hm_ψ[3] = ψ1'
    ct_b[3] = ε * b1' .+ b⁰'
    ct_v[3] = ε * b1' .+ b⁰'
    ct_ψ[3] = ε * b1' .+ b⁰'
    # Update the title
    title_obs[] = @sprintf("t=%.2f", t)
end
println("Animation saved to b1v1psi1_evolution.gif")