using DifferentialEquations
using FFTW
using LinearAlgebra
using CairoMakie
using Printf

# --- 1. Parameters and Grid Setup ---

# Physical Parameters
const L = 20.0
const H = 1.0
const d_ekman = 0.1

# Numerical Parameters
const Nx = 128
const Nz = 64
const tmax = 20.0 

# REVISION: Parameters for de-aliasing (2/3 rule)
const Nx_padded = Int(floor(3/2 * Nx))
const k_padded = rfftfreq(Nx_padded, Nx_padded * 2π / L)

# Grid Setup
x = L / Nx * (-Nx/2 : Nx/2 - 1)
k = rfftfreq(Nx, Nx * 2π / L)
z = range(0, H, length=Nz)
dz = z[2] - z[1]
# REVISION: Padded grid for de-aliasing
x_padded = L / Nx_padded * (-Nx_padded/2 : Nx_padded/2 - 1)

# --- 2. Finite Difference Matrices ---
# (No changes here)
D1 = Tridiagonal(-ones(Nz-1), zeros(Nz), ones(Nz-1)) / (2*dz); D1[1, 1:2] .= [-1, 1] / dz; D1[Nz, Nz-1:Nz] .= [-1, 1] / dz;
D2 = Tridiagonal(ones(Nz-1), -2*ones(Nz), ones(Nz-1)) / dz^2; D2[1, :] .= 0; D2[1,1] = 1.0; D2[Nz, :] .= 0; D2[Nz,Nz] = 1.0;
Lz = inv(D2 + 1e-12 * I);

# --- 3. Setup of Exact Background State ---
println("Setting up exact background state...")
# (No changes to the derivation logic)
B0(x_val) = 0.5 * erf(x_val / sqrt(2)); B0_prime(x_val) = 1/sqrt(2*pi) * exp(-x_val^2 / 2); B0_double_prime(x_val) = -x_val/sqrt(2*pi) * exp(-x_val^2 / 2); B0_triple_prime(x_val) = (x_val^2 - 1)/sqrt(2*pi) * exp(-x_val^2 / 2);

# REVISION: Create storage for padded background fields
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

println("Background state setup complete.")

# --- 4. ODE Function `rhs!` ---
function rhs!(du_hat, u_hat, p, t)
    k, D1, Lz, x_padded, _, v0_pad, v0_x_pad, v0_z_pad, b0_x_pad, b0_z_pad, _, Nx_padded = p
    
    omega_hat = @view u_hat[:, :, 1]
    v1_hat = @view u_hat[:, :, 2]
    b1_hat = @view u_hat[:, :, 3]

    psi_hat = Lz * omega_hat
    
    # --- DE-ALIASING PROCEDURE ---
    # 1. Pad spectral data with zeros
    psi_hat_padded = zeros(ComplexF64, Nz, length(k_padded))
    psi_hat_padded[:, 1:size(psi_hat, 2)] .= psi_hat

    # 2. Transform to padded physical grid
    u1_phys_padded = D1 * irfft(psi_hat_padded, Nx_padded, 2)
    w1_phys_padded = -irfft(im * k_padded' .* psi_hat_padded, Nx_padded, 2)
    
    # 3. Compute products on the padded grid
    Fv_phys_padded = -v0_pad .- (u1_phys_padded .- x_padded') .* v0_x_pad .+ w1_phys_padded .* v0_z_pad
    Fb_phys_padded = 0 .- (u1_phys_padded .- x_padded') .* b0_x_pad .+ w1_phys_padded .* b0_z_pad
    
    # 4. Transform back to padded spectral space
    Fv_hat_padded = rfft(Fv_phys_padded, 2)
    Fb_hat_padded = rfft(Fb_phys_padded, 2)
    
    # 5. Truncate to get de-aliased forcing
    Fv_hat = Fv_hat_padded[:, 1:size(u_hat, 2)]
    Fb_hat = Fb_hat_padded[:, 1:size(u_hat, 2)]
    # --- END DE-ALIASING ---

    d_omega_hat = @view du_hat[:, :, 1]
    d_v1_hat    = @view du_hat[:, :, 2]
    d_b1_hat    = @view du_hat[:, :, 3]

    d_omega_hat .= D1 * v1_hat - im * k' .* b1_hat
    d_v1_hat .= -(D1 * psi_hat) .+ Fv_hat
    d_b1_hat .= Fb_hat

    d_omega_hat[1, :] .= 0
    d_omega_hat[Nz, :] .= 0
end

# --- 5. Solve the ODE System ---
u0 = zeros(ComplexF64, Nz, length(k), 3)
params = (k, D1, Lz, x_padded, z, v0_padded, v0_x_padded, v0_z_padded, b0_x_padded, b0_z_padded, Nx, Nx_padded)
tspan = (0.0, tmax)

prob = ODEProblem(rhs!, u0, tspan, params)
println("Solving ODE system with de-aliasing...")
sol = solve(prob, Vern7(), reltol=1e-8, abstol=1e-8, saveat=0.2)
println("Solution complete.")

# --- 6. Post-processing and Visualization ---
println("Reconstructing and animating solution...")

# Helper function to reconstruct all fields at a given time index
function reconstruct_fields(sol_u, p)
    Lz, Nx = p[3], p[end-1]
    
    omega_hat = sol_u[:, :, 1]
    v1_hat = sol_u[:, :, 2]
    b1_hat = sol_u[:, :, 3]
    
    psi_hat = Lz * omega_hat
    
    psi1_phys = irfft(psi_hat, Nx, 2)
    v1_phys = irfft(v1_hat, Nx, 2)
    b1_phys = irfft(b1_hat, Nx, 2)
    
    return psi1_phys, v1_phys, b1_phys
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
record(fig, "b1v1psi1_evolutionv0.gif", 1:length(sol.t); framerate = 10) do i
    t = sol.t[i]
    ψ1, v1, b1 = reconstruct_fields(sol.u[i], params)
    
    # Update the heatmap data
    hm_b[3] = b1'
    hm_v[3] = v1'
    hm_ψ[3] = ψ1'
    
    # Update the title
    title_obs[] = @sprintf("t=%.2f", t)
end
println("Animation saved to b1v1psi1_evolutionv0.gif")