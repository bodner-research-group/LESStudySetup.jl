#=
File: compute_psi1.jl
Author: Gemini @ MIT
Date: 2025-07-07
Description: Computes and visualizes the 2D function ψ¹(x,z) based on the
             provided semi-analytic solution.
=#

using CairoMakie
using SixelTerm
using SpecialFunctions # For erf()
set_theme!(theme_latexfonts(), fontsize=12,figure_padding = 10) 
# ------------------------------------------------------------------
# 1.  Setup and Grid Definition
# ------------------------------------------------------------------
println("Setting up grid and parameters...")

# Physical domain
x_min, x_max = -4.0, 4.0
z_min, z_max = 0.0, 1.0
Ro = 1.0  # Rossby number

# Computational grids (moderate resolution; raise if desired)
nx, nz = 201, 101
x = collect(range(x_min, x_max; length = nx))   # <- make it a Vector!
z = collect(range(z_min, z_max; length = nz))   # (z isn’t mutated, but keep symmetrical)
const dx = x[2] - x[1]
const dz = z[2] - z[1]

# ------------------------------------------------------------------
# 2.  Background Buoyancy Profile and Derivatives
#
# B₀(X) is defined using the error function, representing a shear layer.
# We need up to the fourth derivative for the Aᵢ coefficients.
# ------------------------------------------------------------------
println("Defining buoyancy profiles...")
const invs2π = 1 / √(2π)

# --- Case 1: Error Function Profile (Original) ---
B0_erf(X)   = 0.5 * erf(X / √2)
dB0_erf(X)  = @. exp(-0.5 * X^2) * invs2π
d2B0_erf(X) = @. -X * exp(-0.5 * X^2) * invs2π
d3B0_erf(X) = @. (X^2 - 1) * exp(-0.5 * X^2) * invs2π
d4B0_erf(X) = @. (3X - X^3) * exp(-0.5 * X^2) * invs2π

const B_funcs_erf = (dB0_erf, d2B0_erf, d3B0_erf, d4B0_erf)

# --- Case 2: Gaussian Profile (New) ---
B0_gauss(X)   = @. -0.5 * exp(-0.5 * X^2)
dB0_gauss(X)  = @. 0.5 * X * exp(-0.5 * X^2)
d2B0_gauss(X) = @. -0.5 * (X^2 - 1) * exp(-0.5 * X^2)
d3B0_gauss(X) = @. 0.5 * X * (X^2 - 3) * exp(-0.5 * X^2)
d4B0_gauss(X) = @. -0.5 * (X^4 - 6*X^2 + 3) * exp(-0.5 * X^2) # Note: 2*1.5=3, 2*3=6

const B_funcs_gauss = (dB0_gauss, d2B0_gauss, d3B0_gauss, d4B0_gauss)

# ------------------------------------------------------------------
# 3.  Generic ψ¹ Computation
#
# These functions are now generic and take the required buoyancy
# derivatives as arguments, allowing us to reuse the same logic
# for different physical cases.
# ------------------------------------------------------------------

"""
    solve_X0_vec(x_vec, z_val, Ro, dB0, d2B0)

Generic Newton solver for X₀ that accepts buoyancy derivatives.
"""
function solve_X0_vec(x_vec, z_val::Float64, Ro::Float64, dB0::Function, d2B0::Function; maxiter=50, tol=1e-12)
    X = collect(x_vec) # Initial guess for X₀ is x
    Ro_sq_z = Ro^2 * (z_val - 0.5)

    for _ in 1:maxiter
        Bp  = dB0(X)
        Bpp = d2B0(X)
        f   = X .- (x_vec .+ Ro_sq_z .* Bp)
        df  = 1.0 .- Ro_sq_z .* Bpp
        dX  = f ./ df
        X .-= dX
        if maximum(abs, dX) < tol; break; end
    end
    return X
end

"""
    compute_coefficients(X0, Ro, B_funcs; use_new_A_coeffs)

Generic coefficient calculator. A boolean flag switches between Aᵢ sets.
"""
function compute_coefficients(X0::Float64, Ro::Float64, B_funcs::Tuple; use_new_A_coeffs::Bool = false)
    dB0, d2B0, d3B0, d4B0 = B_funcs
    
    B_p = dB0(X0); B_pp = d2B0(X0); B_ppp = d3B0(X0); B_pppp = d4B0(X0)

    Ro2 = Ro^2; Ro4 = Ro^4; Bp_sq = B_p^2; Bpp_sq = B_pp^2

    local A0, A1, A2
    if !use_new_A_coeffs
        # Original Aᵢ coefficient definitions
        A0 = 3 * B_p * (2 * Bpp_sq + B_ppp * B_p)
        A1 = -Ro2 * B_p * (12 * B_pp^3 - B_pppp * Bp_sq - 3 * B_ppp * B_p * B_pp)
        A2 = Ro4 * B_p * (6 * Bpp_sq * (Bpp_sq - B_ppp * B_p) - Bp_sq * (B_pppp * B_pp - 3 * B_ppp^2))
    else
        # New Aᵢ coefficient definitions
        A0 = B_ppp
        A1 = Ro2 * (B_ppp * B_pp + B_pppp * B_p)
        A2 = -Ro4 * (2 * B_ppp * Bpp_sq + B_p * (B_pppp * B_pp - 3 * B_ppp^2))
    end

    Ro2Bp = Ro2 * B_p
    Psi(Z) = 2 + Ro2Bp * (1 - 2Z)
    denom = (Ro2Bp)^4
    
    if abs(denom) < 1e-20; return ntuple(_ -> 0.0, 5); end

    C3 = A2 / (2 * denom)
    C4 = (A1 * Ro2Bp + 2 * A2) / denom
    C5 = (A0 * Ro2Bp^2 + A1 * Ro2Bp + A2) / denom

    function f(Psi_val)
        if Psi_val <= 0; return 0.0; end
        return C3 * Psi_val * (log(Psi_val / 2) - 1) + C4 * log(Psi_val) + C5 / Psi_val
    end

    C2 = -f(Psi(0.0))
    C1 = -(C2 + f(Psi(1.0)))
    
    return C1, C2, C3, C4, C5
end

"""
    compute_psi1_field(x, z, Ro, B_funcs)

Top-level function to compute the entire ψ¹ field for a given set of
buoyancy functions.
"""
function compute_psi1_field(x, z, Ro, B_funcs; use_new_A_coeffs::Bool = false)
    nx, nz = length(x), length(z)
    psi1_field = zeros(Float64, nx, nz)
    dB0_func, d2B0_func = B_funcs[1], B_funcs[2]

    for j in 1:nz
        z_val = z[j]
        X0_row = solve_X0_vec(x, z_val, Ro, dB0_func, d2B0_func)
        
        for i in 1:nx
            X0_val = X0_row[i]
            C1, C2, C3, C4, C5 = compute_coefficients(X0_val, Ro, B_funcs; use_new_A_coeffs=use_new_A_coeffs)
            
            Ro2Bp = Ro^2 * dB0_func(X0_val)
            Psi_val = 2 + Ro2Bp * (1 - 2 * z_val)
            
            if Psi_val <= 0; psi1_field[i, j] = NaN; continue; end

            term3 = C3 * Psi_val * (log(Psi_val / 2) - 1)
            term4 = C4 * log(Psi_val)
            term5 = C5 / Psi_val
            
            psi1_field[i, j] = C1 * z_val + C2 + term3 + term4 + term5
        end
    end
    return psi1_field
end

# ------------------------------------------------------------------
# 4.  Main Execution
# ------------------------------------------------------------------
println("Computing all four cases...")
# Case 1: Erf profile, Original A coeffs
psi1_erf_orig = compute_psi1_field(x, z, Ro, B_funcs_erf)
# Case 2: Gaussian profile, Original A coeffs
psi1_gauss_orig = compute_psi1_field(x, z, Ro, B_funcs_gauss)
# Case 3: Erf profile, New A coeffs
psi1_erf_new = compute_psi1_field(x, z, Ro, B_funcs_erf; use_new_A_coeffs=true)
# Case 4: Gaussian profile, New A coeffs
psi1_gauss_new = compute_psi1_field(x, z, Ro, B_funcs_gauss; use_new_A_coeffs=true)
println("All computations complete.")

# ------------------------------------------------------------------
# 5.  Visualization with CairoMakie
# ------------------------------------------------------------------
println("Generating comparison plot...")

# Create a figure and axis
fig1 = Figure(size = (640, 450));
gab = fig1[1, 1] = GridLayout()
ax_a  = Axis(gab[1,1], ylabel = L"z",
             titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
             title = L"\text{(a) Front}~\psi^1_{V,p}")
clims = (-maximum(abs, psi1_erf_orig), maximum(abs, psi1_erf_orig))
hm_a = heatmap!(ax_a, x, z, -psi1_erf_orig;
        colormap = :PuOr, rasterize = true, colorrange = clims)
Colorbar(gab[1,2], hm_a)
ax_b  = Axis(gab[1,3],  
             titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
             title = L"\text{(b) Filament}~\psi^1_{V,p}")
hm_b = heatmap!(ax_b, x, z, -psi1_gauss_orig;
        colormap = :PuOr, rasterize = true, 
        colorrange = (-maximum(abs, psi1_gauss_orig), maximum(abs, psi1_gauss_orig)))
Colorbar(gab[1,4], hm_b)
ax_c = Axis(gab[2,1], xlabel = L"x", ylabel = L"z", 
            titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
            title = L"\text{(c) Front}~\psi^1_{H,p}")
clims = (-maximum(abs, psi1_erf_new), maximum(abs, psi1_erf_new))
hm_c = heatmap!(ax_c, x, z, -psi1_erf_new;
colormap = :PuOr, rasterize = true, colorrange = clims)
Colorbar(gab[2,2], hm_c)
ax_d  = Axis(gab[2,3], xlabel = L"x", 
            titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
            title = L"\text{(d) Filament}~\psi^1_{H,p}")
hm_d = heatmap!(ax_d, x, z, -psi1_gauss_new;
colormap = :PuOr, rasterize = true, 
colorrange = (-maximum(abs, psi1_gauss_new), maximum(abs, psi1_gauss_new)))
Colorbar(gab[2,4], hm_d)
hideydecorations!(ax_b, ticks=false)
hideydecorations!(ax_d, ticks=false)
hidexdecorations!(ax_a, ticks=false)
hidexdecorations!(ax_b, ticks=false)
rowgap!(gab, 5)
colgap!(gab, 1, 3)
colgap!(gab, 2, 10)
colgap!(gab, 3, 3)
resize_to_layout!(fig1)

# ------------------------------------------------------------------
# 6.  Residual Analysis 
# ------------------------------------------------------------------
"""
    compute_residual_field(x, z, Ro, t, B_funcs, dx, dz)

Calculates the residual of the governing equation for the given approximate
solutions. The residual is defined as:
R = ∂/∂t(∂²ψ¹/∂z²) - ∂v¹/∂z + ∂b¹/∂x
Derivatives are computed using second-order finite differences.
"""
function compute_residual_field(x, z, Ro, t, B_funcs, dx, dz)
    nx, nz = length(x), length(z)
    dB0, d2B0 = B_funcs[1], B_funcs[2]

    # --- Step 1: Compute X₀ and buoyancy derivatives on the full grid ---
    X0_grid = zeros(nx, nz)
    for j in 1:nz
        X0_grid[:, j] = solve_X0_vec(x, z[j], Ro, dB0, d2B0)
    end
    B_p_grid = dB0(X0_grid)
    B_pp_grid = d2B0(X0_grid)

    # --- Step 2: Compute the physical fields v¹ and b¹ at time t ---
    v1 = zeros(nx, nz)
    b1 = zeros(nx, nz)
    
    # This is the spatial part of ∂ψ¹/∂t
    psi1_t_spatial = @. -B_p_grid * z' * (z' - 1)

    for j in 1:nz, i in 1:nx
        X0 = X0_grid[i,j]; Bp = B_p_grid[i,j]; Bpp = B_pp_grid[i,j]
        
        N = t * (Bpp * X0 + Bp) - 2 * Bp * sin(t)
        D = 1 - Ro^2 * Bpp * (z[j] - 0.5)
        
        v1[i,j] = (z[j] - 0.5) * N / D
        b1[i,j] = Bp * (t * X0 + v1[i,j])
    end

    # --- Step 3: Compute residual using finite differences ---
    residual = fill(NaN, nx, nz) # Use NaN for boundaries
    
    # We compute the time derivative of ∂²ψ¹/∂z² as ∂²/∂z²(∂ψ¹/∂t)
    # The ∂ψ¹/∂t term is psi1_t_spatial * sin(t)
    term1_spatial_part = zeros(nx, nz)
    for j in 2:nz-1, i in 1:nx # No x-derivatives needed here
        term1_spatial_part[i,j] = (psi1_t_spatial[i,j+1] - 2*psi1_t_spatial[i,j] + psi1_t_spatial[i,j-1]) / dz^2
    end
    
    term1 = term1_spatial_part .* sin(t)

    for j in 2:nz-1, i in 2:nx-1
        # ∂v¹/∂z
        dv1_dz = (v1[i, j+1] - v1[i, j-1]) / (2 * dz)
        
        # ∂b¹/∂x
        db1_dx = (b1[i+1, j] - b1[i-1, j]) / (2 * dx)
        
        residual[i,j] = (term1[i,j] - dv1_dz + db1_dx)
    end
    
    return residual/maximum(abs, term1_spatial_part)
end

println("\nStarting residual analysis for approximate solutions...")
t_eval = π / 2
residual_erf = compute_residual_field(x, z, Ro, t_eval, B_funcs_erf, dx, dz)
residual_gauss = compute_residual_field(x, z, Ro, t_eval, B_funcs_gauss, dx, dz)
println("Residual analysis complete.")

# ------------------------------------------------------------------
# 7.  Identity Verification for Stationary Solutions
# ------------------------------------------------------------------

"""
    compute_identity_error(psi1_field, x, z, Ro, B_funcs, dx, dz)

Computes the error in the identity J₀(∂/∂z - Ro²B'₀∂/∂x)²ψ¹ - ∂²ψ¹/∂z² = 0.
Derivatives are computed with second-order finite differences.
"""
function compute_identity_error(psi1_field, x, z, Ro, B_funcs, dx, dz)
    nx, nz = length(x), length(z)
    dB0, d2B0 = B_funcs[1], B_funcs[2]

    # Step 1: Compute X₀ and required buoyancy derivatives on the full grid
    X0_grid = zeros(nx, nz)
    for j in 1:nz
        X0_grid[:, j] = solve_X0_vec(x, z[j], Ro, dB0, d2B0)
    end
    Bp_grid = dB0(X0_grid); Bpp_grid = d2B0(X0_grid)

    # Step 2: Compute intermediate field ϕ = (∂/∂z - Ro²B'₀∂/∂x)ψ¹
    phi = fill(NaN, nx, nz)
    dpsi1_dx = fill(NaN, nx, nz); dpsi1_dz = fill(NaN, nx, nz)
    for j in 2:nz-1, i in 2:nx-1
        dpsi1_dx[i,j] = (psi1_field[i+1, j] - psi1_field[i-1, j]) / (2 * dx)
        dpsi1_dz[i,j] = (psi1_field[i, j+1] - psi1_field[i, j-1]) / (2 * dz)
        phi[i,j] = dpsi1_dz[i,j] - Ro^2 * Bp_grid[i,j] * dpsi1_dx[i,j]
    end

    # Step 3: Compute L²ψ¹ = (∂/∂z - Ro²B'₀∂/∂x)ϕ
    L2_psi1 = fill(NaN, nx, nz)
    for j in 3:nz-2, i in 3:nx-2
        dphi_dx = (phi[i+1, j] - phi[i-1, j]) / (2 * dx)
        dphi_dz = (phi[i, j+1] - phi[i, j-1]) / (2 * dz)
        L2_psi1[i,j] = dphi_dz - Ro^2 * Bp_grid[i,j] * dphi_dx
    end

    # Step 4: Compute ∂²ψ¹/∂z²
    d2psi1_dz2 = fill(NaN, nx, nz)
    for j in 2:nz-1, i in 1:nx
        d2psi1_dz2[i,j] = (psi1_field[i, j+1] - 2*psi1_field[i,j] + psi1_field[i, j-1]) / dz^2
    end

    # Step 5: Compute J₀ and final error
    z_grid = z'
    J0_grid = @. 1 / (1 - Ro^2 * Bpp_grid * (z_grid - 0.5))
    error_field = @. J0_grid * L2_psi1 - d2psi1_dz2
    
    return error_field./maximum(abs, filter(!isnan, d2psi1_dz2))
end

println("\nStarting identity verification analysis...")
identity_error_erf_orig = compute_identity_error(psi1_erf_orig, x, z, Ro, B_funcs_erf, dx, dz)
identity_error_gauss_orig = compute_identity_error(psi1_gauss_orig, x, z, Ro, B_funcs_gauss, dx, dz)
identity_error_erf_new = compute_identity_error(psi1_erf_new, x, z, Ro, B_funcs_erf, dx, dz)
identity_error_gauss_new = compute_identity_error(psi1_gauss_new, x, z, Ro, B_funcs_gauss, dx, dz)
println("Identity verification complete.")

# --- Figure 2: Appendix plot for the residual analysis ---
fig2 = Figure(size = (640, 250))
ax_r1 = Axis(fig2[1, 1], titlealign = :left, 
    title = L"\text{(a) Front}~\psi^1_{\gamma}~/text{error}",
    xlabel = L"x", ylabel = L"z", 
)
max_abs_resid = maximum(abs, filter(!isnan, residual_erf))
hm_r1 = heatmap!(ax_r1, x, z, residual_erf; 
                 rasterize = true,colormap=:balance,colorrange=(-max_abs_resid, max_abs_resid))

ax_r2 = Axis(fig2[1, 3], titlealign = :left, 
    title = L"\text{(b) Filament}~\psi^1_{\gamma}~/text{error}",
    xlabel = L"x", 
)
max_abs_resid = maximum(abs, filter(!isnan, residual_gauss))
hm_r2 = heatmap!(ax_r2, x, z, residual_gauss; 
                 rasterize = true,colormap=:balance,colorrange=(-max_abs_resid, max_abs_resid))

linkyaxes!(ax_r1, ax_r2)
hideydecorations!(ax_r2, grid=false)
Colorbar(fig2[1, 2], hm_r1)
Colorbar(fig2[1, 4], hm_r1)
colgap!(fig2.layout, 5)
resize_to_layout!(fig2)

# --- Figure 3: Appendix plot for the stationary identity verification ---


fig3 = Figure(size = (640, 450));
gab = fig3[1, 1] = GridLayout()
ax_a  = Axis(gab[1,1], ylabel = L"z",
             titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
             title = L"\text{(a) Front}~\psi^1_{V}~\text{error}")
max_abs_id_error = maximum(abs, filter(!isnan,identity_error_erf_orig))
id_range = (-max_abs_id_error, max_abs_id_error)
hm_a = heatmap!(ax_a, x, z, -identity_error_erf_orig;
                colormap = :balance, rasterize = true, colorrange = id_range)
Colorbar(gab[1,2], hm_a)
ax_b  = Axis(gab[1,3],  
             titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
             title = L"\text{(b) Filament}~\psi^1_{V}~\text{error}")
max_abs_id_error = maximum(abs, filter(!isnan,identity_error_gauss_orig))
id_range = (-max_abs_id_error, max_abs_id_error)
hm_b = heatmap!(ax_b, x, z, -identity_error_gauss_orig;
                colormap = :balance, rasterize = true, colorrange = id_range)
Colorbar(gab[1,4], hm_b)
ax_c = Axis(gab[2,1], xlabel = L"x", ylabel = L"z", 
            titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
            title = L"\text{(c) Front}~\psi^1_{H}~\text{error}")
max_abs_id_error = maximum(abs, filter(!isnan,identity_error_erf_new))
id_range = (-max_abs_id_error, max_abs_id_error)
hm_c = heatmap!(ax_c, x, z, -identity_error_erf_new;
                colormap = :balance, rasterize = true, colorrange = id_range)
Colorbar(gab[2,2], hm_c)
ax_d  = Axis(gab[2,3], xlabel = L"x", 
            titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
            title = L"\text{(d) Filament}~\psi^1_{H}~\text{error}")
max_abs_id_error = maximum(abs, filter(!isnan,identity_error_gauss_new))
id_range = (-max_abs_id_error, max_abs_id_error)
hm_d = heatmap!(ax_d, x, z, -identity_error_gauss_new;
                colormap = :balance, rasterize = true, colorrange = id_range)
Colorbar(gab[2,4], hm_d)
hideydecorations!(ax_b, ticks=false)
hideydecorations!(ax_d, ticks=false)
hidexdecorations!(ax_a, ticks=false)
hidexdecorations!(ax_b, ticks=false)
rowgap!(gab, 5)
colgap!(gab, 1, 3)
colgap!(gab, 2, 10)
colgap!(gab, 3, 3)
resize_to_layout!(fig3)

# --- Display Both Figures ---
display(fig1)
display(fig2)
display(fig3)

println("Script finished.")