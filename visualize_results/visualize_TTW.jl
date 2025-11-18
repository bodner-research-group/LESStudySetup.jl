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
γ = 0.03 # Strain parameter

# Computational grids (moderate resolution; raise if desired)
nx, nz = 401, 101
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
function compute_coefficients(X0, Ro, B_funcs; use_new_A_coeffs = false, return_A = false)
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

    if return_A
        return A0, A1, A2
    end

    # Compute the Cᵢ coefficients

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
function compute_psi1_field(x, z, Ro, B_funcs; use_new_A_coeffs::Bool = false, strain=false)
    nx, nz = length(x), length(z)
    psi1_field = zeros(Float64, nx, nz)
    dB0_func, d2B0_func = B_funcs[1], B_funcs[2]

    for j in 1:nz
        z_val = z[j]
        X0_row = solve_X0_vec(x, z_val, Ro, dB0_func, d2B0_func)
        
        for i in 1:nx
            X0_val = X0_row[i]
            if strain 
                psi1_field[i, j] = -Ro * dB0_func(X0_val) * z_val * (z_val-1)
            else
                C1, C2, C3, C4, C5 = compute_coefficients(X0_val, Ro, B_funcs; use_new_A_coeffs)
                
                Ro2Bp = Ro^2 * dB0_func(X0_val)
                Psi_val = 2 + Ro2Bp * (1 - 2 * z_val)
                
                if Psi_val <= 0; psi1_field[i, j] = NaN; continue; end

                term3 = C3 * Psi_val * (log(Psi_val / 2) - 1)
                term4 = C4 * log(Psi_val)
                term5 = C5 / Psi_val
                
                psi1_field[i, j] = C1 * z_val + C2 + term3 + term4 + term5
            end
        end
    end
    return psi1_field
end

function compute_blt_fields(x, z, Ro, B_funcs)
    nx, nz = length(x), length(z)
    ϕzzz_field = zeros(Float64, nx, nz)
    ϕxzz_field = zeros(Float64, nx, nz)
    ϕxxz_field = zeros(Float64, nx, nz)
    ϕxxx_field = zeros(Float64, nx, nz)
    ϕγv_field = zeros(Float64, nx, nz)
    ϕγb_field = zeros(Float64, nx, nz)
    dB0_func, d2B0_func, d3B0_func = B_funcs[1], B_funcs[2], B_funcs[3]

    for j in 1:nz
        z_val = z[j]
        X0_row = solve_X0_vec(x, z_val, Ro, dB0_func, d2B0_func)
        
        for i in 1:nx
            X0_val = X0_row[i]
            dB0_val, d2B0_val, d3B0_val = dB0_func(X0_val), d2B0_func(X0_val), d3B0_func(X0_val)
            J0_val = 1 / (1 - Ro^2 * d2B0_val * (z_val - 0.5))
            ϕzzz_field[i, j] = Ro^4 * dB0_val^2 * J0_val^2 * (3d2B0_val+Ro^2*dB0_val*d3B0_val*(z_val-0.5)*J0_val)
            ϕxzz_field[i, j] = Ro^2 * dB0_val * J0_val^2 * (2d2B0_val+Ro^2*dB0_val*d3B0_val*(z_val-0.5)*J0_val)
            ϕxxz_field[i, j] = J0_val^2 * (d2B0_val + Ro^2*dB0_val*d3B0_val*(z_val - 0.5)*J0_val)
            ϕxxx_field[i, j] = d3B0_val * J0_val^3 * (z_val - 0.5)
            ϕγv_field[i, j] = - dB0_val*(z_val - 0.5) + x[i]*J0_val*d2B0_val*(z_val - 0.5)
            ϕγb_field[i, j] = x[i]*J0_val*dB0_val
        end
    end

    return ϕzzz_field, ϕxzz_field, ϕxxz_field, ϕxxx_field, ϕγv_field, ϕγb_field
end

function compute_b0v0_field(x, z, Ro, B_funcs; gauss = false)
    nx, nz = length(x), length(z)
    b_field = zeros(Float64, nx, nz)
    v_field = zeros(Float64, nx, nz)
    dB0_func, d2B0_func = B_funcs[1], B_funcs[2]

    for j in 1:nz
        z_val = z[j]
        X0_row = solve_X0_vec(x, z_val, Ro, dB0_func, d2B0_func)
        
        for i in 1:nx
            X0_val = X0_row[i]
            if gauss
                b_field[i,j] = B0_gauss(X0_val)
            else
                b_field[i,j] = B0_erf(X0_val)
            end
            v_field[i, j] = dB0_func(X0_val) * (z_val - 0.5)
        end
    end

    return b_field, v_field
end

function compute_ϕ12_field(x, z, Ro, B_funcs)
    nx, nz = length(x), length(z)
    ϕxx_field = zeros(Float64, nx, nz)
    ϕxz_field = zeros(Float64, nx, nz)
    ϕzz_field = zeros(Float64, nx, nz)
    dB0_func, d2B0_func = B_funcs[1], B_funcs[2]

    for j in 1:nz
        z_val = z[j]
        X0_row = solve_X0_vec(x, z_val, Ro, dB0_func, d2B0_func)
        
        for i in 1:nx
            X0_val = X0_row[i]
            dB0_val, d2B0_val = dB0_func(X0_val), d2B0_func(X0_val)
            J0_val = 1 / (1 - Ro^2 * d2B0_val * (z_val - 0.5))
            ϕxx_field[i, j] = J0_val*d2B0_val*(z_val - 0.5)
            ϕxz_field[i, j] = J0_val*dB0_val
            ϕzz_field[i, j] = J0_val*dB0_val^2
        end
    end

    return ϕxx_field, ϕxz_field, ϕzz_field
end

function central_diff(A, ds, dim)
    out = fill(NaN, size(A))
    if dim == 1 # d/dx for (nx, nz) matrix
        C = 2:size(A,1)-1; K = 1:size(A,2)
        @views @. out[C,K] = (A[C.+1,K] - A[C.-1,K]) / (2*ds)
        out[1, :] = (A[2, :] - A[1, :]) / ds
        out[end, :] = (A[end, :] - A[end-1, :]) / ds
    elseif dim == 2 # d/dz for (nx, nz) matrix
        C = 1:size(A,1); K = 2:size(A,2)-1
        @views @. out[C,K] = (A[C,K.+1] - A[C,K.-1]) / (2*ds)
        out[:, 1] = (A[:, 2] - A[:, 1]) / ds
        out[:, end] = (A[:, end] - A[:, end-1]) / ds
    end
    return out
end

function compute_b1v1_field(x, z, Ro, Fb, Fv, ψ¹ₚ, ϕ12s, t)
    dx, dz = x[2] - x[1], z[2] - z[1]
    ϕxx, ϕxz, ϕzz = ϕ12s[1], ϕ12s[2], ϕ12s[3]
    ψ¹ₚz = central_diff(ψ¹ₚ, dz, 2)
    ψ¹ₚx = central_diff(ψ¹ₚ, dx, 1)

    b_field = Ro * Fb * t
    v_field = Ro * Fv * t
    b_field .-= Ro * (ϕxz.* ψ¹ₚz .- ϕzz.* ψ¹ₚx) * (t - sin(t))
    v_field .-= Ro * ((1/Ro^2 .+ ϕxx).* ψ¹ₚz .- ϕxz.* ψ¹ₚx) * (t - sin(t))
    return b_field, v_field
end

function compute_T1b_field(x, z, Ro, Fb, Fv, ψ¹ₚ, ϕ12s, t; γ̃ = 1, dt = 0.001)
    dx, dz = x[2] - x[1], z[2] - z[1]
    _, ϕxz, _ = ϕ12s[1], ϕ12s[2], ϕ12s[3]
    ψ¹ₚz = central_diff(ψ¹ₚ, dz, 2)
    ψ¹ₚx = central_diff(ψ¹ₚ, dx, 1)
    b1⁺, _ = compute_b1v1_field(x, z, Ro, Fb, Fv, ψ¹ₚ, ϕ12s, t+dt)
    b1⁻, _ = compute_b1v1_field(x, z, Ro, Fb, Fv, ψ¹ₚ, ϕ12s, t-dt)
    b1⁺dx = central_diff(b1⁺, dx, 1)
    b1⁻dx = central_diff(b1⁻, dx, 1)

    T¹b = ϕxz .* (b1⁺dx .- b1⁻dx) / (2*dt) / Ro
    T¹b .+= 0.5*((ψ¹ₚz .- γ̃ / Ro * reshape(x,:,1)) .* central_diff(ϕxz.^2, dx, 1) .- ψ¹ₚx .* central_diff(ϕxz.^2, dz, 2))

    return T¹b
end

function compute_bVbH_fields(x, z, Ro, b0, ε, Fbs, Fvs, ψ¹ₚs, ϕ12s, t; ψ¹ₚγ = nothing, γ = 0.03)
    if isnothing(ψ¹ₚγ)
        b1, _ = compute_b1v1_field(x, z, Ro, Fbs[1], 0.0*Fvs[1], ψ¹ₚs[1], ϕ12s, t)
        bV = b0 .+ ε * b1
        b1, _ = compute_b1v1_field(x, z, Ro, Fbs[2], 0.0*Fvs[2], ψ¹ₚs[2], ϕ12s, t)
        bH = b0 .+ ε * b1

        b1, _ = compute_b1v1_field(x, z, Ro, 0.0*Fbs[1], Fvs[1], -ψ¹ₚs[1], ϕ12s, t)
        bVv = b0 .+ ε * b1
        b1, _ = compute_b1v1_field(x, z, Ro, 0.0*Fbs[2], Fvs[2], -ψ¹ₚs[2], ϕ12s, t)
        bHv = b0 .+ ε * b1
    else
        b1, _ = compute_b1v1_field(x, z, Ro, Fbs[1].+ε/γ*Fbs[2], Fvs[1], ψ¹ₚγ .+ ε/γ * ψ¹ₚs[1], ϕ12s, t)
        bV = b0 .+ γ * b1
        b1, _ = compute_b1v1_field(x, z, Ro, Fbs[1].+ε/γ*Fbs[3], Fvs[1], ψ¹ₚγ .+ ε/γ * ψ¹ₚs[2], ϕ12s, t)
        bH = b0 .+ γ * b1

        b1, _ = compute_b1v1_field(x, z, Ro, Fbs[1], Fvs[1].+ε/γ*Fvs[2], ψ¹ₚγ .- ε/γ * ψ¹ₚs[1], ϕ12s, t)
        bVv = b0 .+ γ * b1
        b1, _ = compute_b1v1_field(x, z, Ro, Fbs[1], Fvs[2].+ε/γ*Fvs[3], ψ¹ₚγ .- ε/γ * ψ¹ₚs[2], ϕ12s, t)
        bHv = b0 .+ γ * b1
    end
    return bV, bH, bVv, bHv
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
# Case 5: Erf profile, strain
psi1_erf_strain = compute_psi1_field(x, z, Ro, B_funcs_erf; strain=true)
# Case 6: Gaussian profile, strain
psi1_gauss_strain = compute_psi1_field(x, z, Ro, B_funcs_gauss; strain=true)
println("All computations complete.")

b_erf, v_erf = compute_b0v0_field(x, z, Ro, B_funcs_erf)
b_gauss, v_gauss = compute_b0v0_field(x, z, Ro, B_funcs_gauss; gauss=true)
# ------------------------------------------------------------------
# 5.  Visualization with CairoMakie
# ------------------------------------------------------------------
println("Generating comparison plot...")

# Create a figure and axis
fig0 = Figure(size = (640, 450));
gab = fig0[1, 1] = GridLayout()
ax_a  = Axis(gab[1,1], ylabel = L"z",
             titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
             title = L"\text{(a) Front}~b^0")
clims = (-maximum(abs, b_erf), maximum(abs, b_erf))
hm_a = heatmap!(ax_a, x, z, b_erf;
        colormap = :diff, rasterize = true, colorrange = clims)
Colorbar(gab[1,2], hm_a)
ax_b  = Axis(gab[1,3],  
             titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
             title = L"\text{(b) Filament}~b^0")
hm_b = heatmap!(ax_b, x, z, b_gauss;
        colormap = :diff, rasterize = true, 
        colorrange = (-maximum(abs, b_gauss), maximum(abs, b_gauss)))
Colorbar(gab[1,4], hm_b)
ax_c = Axis(gab[2,1], xlabel = L"x", ylabel = L"z", 
            titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
            title = L"\text{(c) Front}~v^0")
clims = (-maximum(abs, v_erf), maximum(abs, v_erf))
hm_c = heatmap!(ax_c, x, z, v_erf;
colormap = :delta, rasterize = true, colorrange = clims)
Colorbar(gab[2,2], hm_c)
contour!(ax_c, x, z, b_erf; levels = 10, linewidth = 1, color = :black)
ax_d  = Axis(gab[2,3], xlabel = L"x", 
            titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
            title = L"\text{(d) Filament}~v^0")
hm_d = heatmap!(ax_d, x, z, v_gauss;
colormap = :delta, rasterize = true, 
colorrange = (-maximum(abs, v_gauss), maximum(abs, v_gauss)))
Colorbar(gab[2,4], hm_d)
contour!(ax_d, x, z, b_gauss; levels = 10, linewidth = 1, color = :black)
hideydecorations!(ax_b, ticks=false)
hideydecorations!(ax_d, ticks=false)
hidexdecorations!(ax_a, ticks=false)
hidexdecorations!(ax_b, ticks=false)
rowgap!(gab, 5)
colgap!(gab, 1, 3)
colgap!(gab, 2, 10)
colgap!(gab, 3, 3)
resize_to_layout!(fig0)

# Create a figure and axis
fig1 = Figure(size = (640, 450));
gab = fig1[1, 1] = GridLayout()
ax_a  = Axis(gab[1,1], ylabel = L"z",
             titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
             title = L"\text{(a) Front}~\psi^1_{p,V}")
clims = (-maximum(abs, psi1_erf_orig), maximum(abs, psi1_erf_orig))
hm_a = heatmap!(ax_a, x, z, -psi1_erf_orig;
        colormap = :PuOr, rasterize = true, colorrange = clims)
Colorbar(gab[1,2], hm_a)
ax_b  = Axis(gab[1,3],  
             titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
             title = L"\text{(b) Filament}~\psi^1_{p,V}")
hm_b = heatmap!(ax_b, x, z, -psi1_gauss_orig;
        colormap = :PuOr, rasterize = true, 
        colorrange = (-maximum(abs, psi1_gauss_orig), maximum(abs, psi1_gauss_orig)))
Colorbar(gab[1,4], hm_b)
ax_c = Axis(gab[2,1], xlabel = L"x", ylabel = L"z", 
            titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
            title = L"\text{(c) Front}~\psi^1_{p,H}")
clims = (-maximum(abs, psi1_erf_new), maximum(abs, psi1_erf_new))
hm_c = heatmap!(ax_c, x, z, -psi1_erf_new;
colormap = :PuOr, rasterize = true, colorrange = clims)
Colorbar(gab[2,2], hm_c)
ax_d  = Axis(gab[2,3], xlabel = L"x", 
            titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
            title = L"\text{(d) Filament}~\psi^1_{p,H}")
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

ϕxx_erf, ϕxz_erf, ϕzz_erf = compute_ϕ12_field(x, z, Ro, B_funcs_erf)
ϕxx_gauss, ϕxz_gauss, ϕzz_gauss = compute_ϕ12_field(x, z, Ro, B_funcs_gauss)
ϕzzz_erf, ϕxzz_erf, ϕzxx_erf, ϕxxx_erf, Fγv_erf, Fγb_erf = compute_blt_fields(x, z, Ro, B_funcs_erf)
ϕzzz_gauss, ϕxzz_gauss, ϕzxx_gauss, ϕxxx_gauss, Fγv_gauss, Fγb_gauss = compute_blt_fields(x, z, Ro, B_funcs_gauss)

t = 10
b1_erf, v1_erf = compute_b1v1_field(x, z, Ro, Fγb_erf, Fγv_erf, psi1_erf_strain, (ϕxx_erf, ϕxz_erf, ϕzz_erf), t)
b1_gauss, v1_gauss = compute_b1v1_field(x, z, Ro, Fγb_gauss, Fγv_gauss, psi1_gauss_strain, (ϕxx_gauss, ϕxz_gauss, ϕzz_gauss), t)
bγ_erf, vγ_erf = b_erf .+ γ * b1_erf, v_erf .+ γ * v1_erf
bγ_gauss, vγ_gauss = b_gauss .+ γ * b1_gauss, v_gauss .+ γ * v1_gauss
T¹b_erf = compute_T1b_field(x, z, Ro, Fγb_erf, Fγv_erf, psi1_erf_strain, (ϕxx_erf, ϕxz_erf, ϕzz_erf), t)
T¹b_gauss = compute_T1b_field(x, z, Ro, Fγb_gauss, Fγv_gauss, psi1_gauss_strain, (ϕxx_gauss, ϕxz_gauss, ϕzz_gauss), t)
figγ = Figure(size = (640, 650));
gab = figγ[1, 1] = GridLayout()
ax_a  = Axis(gab[1,1], ylabel = L"z",
             titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
             title = L"\text{(a) Front}~\psi^1_{p,\gamma}")
clims = (-maximum(abs, psi1_erf_strain), maximum(abs, psi1_erf_strain))
hm_a = heatmap!(ax_a, x, z, psi1_erf_strain;
        colormap = :PuOr, rasterize = true, colorrange = clims)
Colorbar(gab[1,2], hm_a)
ax_b  = Axis(gab[1,3],  
             titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
             title = L"\text{(b) Filament}~\psi^1_{p,\gamma}")
hm_b = heatmap!(ax_b, x, z, psi1_gauss_strain;
        colormap = :PuOr, rasterize = true, 
        colorrange = (-maximum(abs, psi1_gauss_strain), maximum(abs, psi1_gauss_strain)))
Colorbar(gab[1,4], hm_b)
ax_c  = Axis(gab[2,1], ylabel = L"z",
             titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
             title = L"\text{(c) Front} v^0+\gamma v^1_\gamma(t=10)")
clims = (-maximum(abs, vγ_erf), maximum(abs, vγ_erf))
hm_c = heatmap!(ax_c, x, z, vγ_erf;
        colormap = :delta, rasterize = true, colorrange = clims)
contour!(ax_c, x, z, bγ_erf; levels = 10, linewidth = 1, color = :black)
Colorbar(gab[2,2], hm_c)
ax_d  = Axis(gab[2,3], ylabel = L"z",
             titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
             title = L"\text{(d) Filament} v^0+\gamma v^1_\gamma(t=10)")
hm_d = heatmap!(ax_d, x, z, vγ_gauss;
        colormap = :delta, rasterize = true, 
        colorrange = (-maximum(abs, vγ_gauss), maximum(abs, vγ_gauss)))
contour!(ax_d, x, z, bγ_gauss; levels = 10, linewidth = 1, color = :black)
Colorbar(gab[2,4], hm_d)
ax_e = Axis(gab[3,1], xlabel = L"x", ylabel = L"z", 
            titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
            title = L"\text{(e) Front} T^1_{b,\gamma}(t=10)")
clims = (-maximum(abs, T¹b_erf), maximum(abs, T¹b_erf))
hm_e = heatmap!(ax_e, x, z, T¹b_erf;
colormap = :balance, rasterize = true, colorrange = clims)
contour!(ax_e, x, z, bγ_erf; levels = 10, linewidth = 1, color = :black)
Colorbar(gab[3,2], hm_e)
ax_f  = Axis(gab[3,3], xlabel = L"x", 
            titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
            title = L"\text{(f) Filament} T^1_{b,\gamma}(t=10)")
hm_f = heatmap!(ax_f, x, z, T¹b_gauss;
colormap = :balance, rasterize = true, 
colorrange = (-maximum(abs, T¹b_gauss), maximum(abs, T¹b_gauss)))
contour!(ax_f, x, z, bγ_gauss; levels = 10, linewidth = 1, color = :black)
Colorbar(gab[3,4], hm_f)
hideydecorations!(ax_b, ticks=false)
hideydecorations!(ax_d, ticks=false)
hideydecorations!(ax_f, ticks=false)
hidexdecorations!(ax_a, ticks=false)
hidexdecorations!(ax_b, ticks=false)
hidexdecorations!(ax_c, ticks=false)
hidexdecorations!(ax_d, ticks=false)
rowgap!(gab, 5)
colgap!(gab, 1, 3)
colgap!(gab, 2, 10)
colgap!(gab, 3, 3)
resize_to_layout!(figγ)

ε = 0.03
b1_erf, v1_erf = compute_b1v1_field(x, z, Ro, ϕzzz_erf, 0.0*ϕxzz_erf, psi1_erf_orig, (ϕxx_erf, ϕxz_erf, ϕzz_erf), t)
b1_gauss, v1_gauss = compute_b1v1_field(x, z, Ro, ϕzzz_gauss, 0.0*ϕxzz_gauss, psi1_gauss_orig, (ϕxx_gauss, ϕxz_gauss, ϕzz_gauss), t)
bV_erf, vV_erf = b_erf .+ ε * b1_erf, v_erf .+ ε * v1_erf
bV_gauss, vV_gauss = b_gauss .+ ε * b1_gauss, v_gauss .+ ε * v1_gauss
b1_erf, v1_erf = compute_b1v1_field(x, z, Ro, ϕzxx_erf, 0.0*ϕxxx_erf, psi1_erf_new, (ϕxx_erf, ϕxz_erf, ϕzz_erf), t)
b1_gauss, v1_gauss = compute_b1v1_field(x, z, Ro, ϕzxx_gauss, 0.0*ϕxxx_gauss, psi1_gauss_new, (ϕxx_gauss, ϕxz_gauss, ϕzz_gauss), t)
bH_erf, vH_erf = b_erf .+ ε * b1_erf, v_erf .+ ε * v1_erf
bH_gauss, vH_gauss = b_gauss .+ ε * b1_gauss, v_gauss .+ ε * v1_gauss

b1_erf, v1_erf = compute_b1v1_field(x, z, Ro, 0.0*ϕzzz_erf, ϕxzz_erf, -psi1_erf_orig, (ϕxx_erf, ϕxz_erf, ϕzz_erf), t)
b1_gauss, v1_gauss = compute_b1v1_field(x, z, Ro, 0.0*ϕzzz_gauss, ϕxzz_gauss, -psi1_gauss_orig, (ϕxx_gauss, ϕxz_gauss, ϕzz_gauss), t)
bVv_erf, vVv_erf = b_erf .+ ε * b1_erf, v_erf .+ ε * v1_erf
bVv_gauss, vVv_gauss = b_gauss .+ ε * b1_gauss, v_gauss .+ ε * v1_gauss
b1_erf, v1_erf = compute_b1v1_field(x, z, Ro, 0.0*ϕzxx_erf, ϕxxx_erf, -psi1_erf_new, (ϕxx_erf, ϕxz_erf, ϕzz_erf), t)
b1_gauss, v1_gauss = compute_b1v1_field(x, z, Ro, 0.0*ϕzxx_gauss, ϕxxx_gauss, -psi1_gauss_new, (ϕxx_gauss, ϕxz_gauss, ϕzz_gauss), t)
bHv_erf, vHv_erf = b_erf .+ ε * b1_erf, v_erf .+ ε * v1_erf
bHv_gauss, vHv_gauss = b_gauss .+ ε * b1_gauss, v_gauss .+ ε * v1_gauss

figv = Figure(size = (640, 850));
gab = figv[1, 1] = GridLayout()
ax_a  = Axis(gab[1,1], ylabel = L"z",
             titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
             title = L"\text{(a) Front vertical diffusivity}")
clims = (-maximum(abs, vV_erf), maximum(abs, vV_erf))
hm_a = heatmap!(ax_a, x, z, vV_erf;
        colormap = :delta, rasterize = true, colorrange = clims)
Colorbar(gab[1,2], hm_a)
contour!(ax_a, x, z, bV_erf; levels = 10, linewidth = 1, color = :black)
ax_b  = Axis(gab[1,3],  
             titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
             title = L"\text{(b) Filament vertical diffusivity}")
hm_b = heatmap!(ax_b, x, z, vV_gauss;
        colormap = :delta, rasterize = true, 
        colorrange = (-maximum(abs, vV_gauss), maximum(abs, vV_gauss)))
Colorbar(gab[1,4], hm_b)
contour!(ax_b, x, z, bV_gauss; levels = 10, linewidth = 1, color = :black)
ax_c  = Axis(gab[2,1], ylabel = L"z",
             titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
             title = L"\text{(c) Front vertical viscosity}")
clims = (-maximum(abs, vVv_erf), maximum(abs, vVv_erf))
hm_c = heatmap!(ax_c, x, z, vVv_erf;
        colormap = :delta, rasterize = true, colorrange = clims)
Colorbar(gab[2,2], hm_c)
contour!(ax_c, x, z, bVv_erf; levels = 10, linewidth = 1, color = :black)
ax_d  = Axis(gab[2,3], 
             titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
             title = L"\text{(d) Filament vertical viscosity}")
hm_d = heatmap!(ax_d, x, z, vVv_gauss;
        colormap = :delta, rasterize = true, 
        colorrange = (-maximum(abs, vVv_gauss), maximum(abs, vVv_gauss)))
Colorbar(gab[2,4], hm_d)
contour!(ax_d, x, z, bVv_gauss; levels = 10, linewidth = 1, color = :black)
ax_e  = Axis(gab[3,1], ylabel = L"z",
            titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
            title = L"\text{(e) Front horizontal diffusivity}")
clims = (-maximum(abs, vH_erf), maximum(abs, vH_erf))
hm_e = heatmap!(ax_e, x, z, vH_erf;
        colormap = :delta, rasterize = true, colorrange = clims)
Colorbar(gab[3,2], hm_e)
contour!(ax_e, x, z, bH_erf; levels = 10, linewidth = 1, color = :black)
ax_f  = Axis(gab[3,3], 
             titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
             title = L"\text{(f) Filament horizontal diffusivity}")
hm_f = heatmap!(ax_f, x, z, vH_gauss;
        colormap = :delta, rasterize = true, 
        colorrange = (-maximum(abs, vH_gauss), maximum(abs, vH_gauss)))
Colorbar(gab[3,4], hm_f)
contour!(ax_f, x, z, bH_gauss; levels = 10, linewidth = 1, color = :black)
ax_g = Axis(gab[4,1], xlabel = L"x", ylabel = L"z", 
            titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
            title = L"\text{(g) Front horizontal viscosity}")
clims = (-maximum(abs, vHv_erf), maximum(abs, vHv_erf))
hm_g = heatmap!(ax_g, x, z, vHv_erf;
colormap = :delta, rasterize = true, colorrange = clims)
Colorbar(gab[4,2], hm_g)
contour!(ax_g, x, z, bHv_erf; levels = 10, linewidth = 1, color = :black)
ax_h  = Axis(gab[4,3], xlabel = L"x", 
            titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
            title = L"\text{(h) Filament horizontal viscosity}")
hm_h = heatmap!(ax_h, x, z, vHv_gauss;
colormap = :delta, rasterize = true, 
colorrange = (-maximum(abs, vHv_gauss), maximum(abs, vHv_gauss)))
Colorbar(gab[4,4], hm_h)
contour!(ax_h, x, z, bHv_gauss; levels = 10, linewidth = 1, color = :black)
hideydecorations!(ax_b, ticks=false)
hideydecorations!(ax_d, ticks=false)
hideydecorations!(ax_f, ticks=false)
hideydecorations!(ax_h, ticks=false)
hidexdecorations!(ax_a, ticks=false)
hidexdecorations!(ax_b, ticks=false)
hidexdecorations!(ax_c, ticks=false)
hidexdecorations!(ax_d, ticks=false)
hidexdecorations!(ax_e, ticks=false)
hidexdecorations!(ax_f, ticks=false)
rowgap!(gab, 5)
colgap!(gab, 1, 3)
colgap!(gab, 2, 10)
colgap!(gab, 3, 3)
resize_to_layout!(figv)

T¹bV_erf = compute_T1b_field(x, z, Ro, ϕzzz_erf, 0.0*ϕxzz_erf, psi1_erf_orig, (ϕxx_erf, ϕxz_erf, ϕzz_erf), t; γ̃ = 0.0)
T¹bV_gauss = compute_T1b_field(x, z, Ro, ϕzzz_gauss, 0.0*ϕxzz_gauss, psi1_gauss_orig, (ϕxx_gauss, ϕxz_gauss, ϕzz_gauss), t; γ̃ = 0.0)
T¹bH_erf = compute_T1b_field(x, z, Ro, ϕzxx_erf, 0.0*ϕxxx_erf, psi1_erf_new, (ϕxx_erf, ϕxz_erf, ϕzz_erf), t; γ̃ = 0.0)
T¹bH_gauss = compute_T1b_field(x, z, Ro, ϕzxx_gauss, 0.0*ϕxxx_gauss, psi1_gauss_new, (ϕxx_gauss, ϕxz_gauss, ϕzz_gauss), t; γ̃ = 0.0)
T¹bVv_erf = compute_T1b_field(x, z, Ro, 0.0*ϕzzz_erf, ϕxzz_erf, -psi1_erf_orig, (ϕxx_erf, ϕxz_erf, ϕzz_erf), t; γ̃ = 0.0)
T¹bVv_gauss = compute_T1b_field(x, z, Ro, 0.0*ϕzzz_gauss, ϕxzz_gauss, -psi1_gauss_orig, (ϕxx_gauss, ϕxz_gauss, ϕzz_gauss), t; γ̃ = 0.0)
T¹bHv_erf = compute_T1b_field(x, z, Ro, 0.0*ϕzxx_erf, ϕxxx_erf, -psi1_erf_new, (ϕxx_erf, ϕxz_erf, ϕzz_erf), t; γ̃ = 0.0)
T¹bHv_gauss = compute_T1b_field(x, z, Ro, 0.0*ϕzxx_gauss, ϕxxx_gauss, -psi1_gauss_new, (ϕxx_gauss, ϕxz_gauss, ϕzz_gauss), t; γ̃ = 0.0)
figTb = Figure(size = (640, 850));
gab = figTb[1, 1] = GridLayout()
ax_a  = Axis(gab[1,1], ylabel = L"z",
             titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
             title = L"\text{(a) Front vertical diffusivity}")
clims = (-maximum(abs, T¹bV_erf), maximum(abs, T¹bV_erf))
hm_a = heatmap!(ax_a, x, z, T¹bV_erf;
        colormap = :balance, rasterize = true, colorrange = clims)
Colorbar(gab[1,2], hm_a)
contour!(ax_a, x, z, bV_erf; levels = 10, linewidth = 1, color = :black)
ax_b  = Axis(gab[1,3],  
             titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
             title = L"\text{(b) Filament vertical diffusivity}")
hm_b = heatmap!(ax_b, x, z, T¹bV_gauss;
        colormap = :balance, rasterize = true, 
        colorrange = (-maximum(abs, T¹bV_gauss), maximum(abs, T¹bV_gauss)))
Colorbar(gab[1,4], hm_b)
contour!(ax_b, x, z, bV_gauss; levels = 10, linewidth = 1, color = :black)
ax_c = Axis(gab[2,1], xlabel = L"x", ylabel = L"z", 
            titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
            title = L"\text{(c) Front vertical viscosity}")
clims = (-maximum(abs, T¹bVv_erf), maximum(abs, T¹bVv_erf))
hm_c = heatmap!(ax_c, x, z, T¹bVv_erf;
colormap = :balance, rasterize = true, colorrange = clims)
Colorbar(gab[2,2], hm_c)
contour!(ax_c, x, z, bVv_erf; levels = 10, linewidth = 1, color = :black)
ax_d  = Axis(gab[2,3], xlabel = L"x", 
            titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
            title = L"\text{(d) Filament vertical viscosity}")
hm_d = heatmap!(ax_d, x, z, T¹bVv_gauss;
colormap = :balance, rasterize = true, 
colorrange = (-maximum(abs, T¹bVv_gauss), maximum(abs, T¹bVv_gauss)))
Colorbar(gab[2,4], hm_d)
contour!(ax_d, x, z, bVv_gauss; levels = 10, linewidth = 1, color = :black)
ax_e = Axis(gab[3,1], xlabel = L"x", ylabel = L"z", 
            titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
            title = L"\text{(e) Front horizontal diffusivity}")
clims = (-maximum(abs, T¹bH_erf), maximum(abs, T¹bH_erf))
hm_e = heatmap!(ax_e, x, z, T¹bH_erf; colormap = :balance, rasterize = true, colorrange = clims)
Colorbar(gab[3,2], hm_e)
contour!(ax_e, x, z, bH_erf; levels = 10, linewidth = 1, color = :black)
ax_f  = Axis(gab[3,3], xlabel = L"x", 
            titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
            title = L"\text{(f) Filament horizontal diffusivity}")
hm_f = heatmap!(ax_f, x, z, T¹bH_gauss; colormap = :balance, rasterize = true, 
colorrange = (-maximum(abs, T¹bH_gauss), maximum(abs, T¹bH_gauss)))
Colorbar(gab[3,4], hm_f)
contour!(ax_f, x, z, bH_gauss; levels = 10, linewidth = 1, color = :black)
ax_g = Axis(gab[4,1], xlabel = L"x", ylabel = L"z", 
            titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
            title = L"\text{(g) Front horizontal viscosity}")
clims = (-maximum(abs, T¹bHv_erf), maximum(abs, T¹bHv_erf))
hm_g = heatmap!(ax_g, x, z, T¹bHv_erf;
colormap = :balance, rasterize = true, colorrange = clims)
Colorbar(gab[4,2], hm_g)
contour!(ax_g, x, z, bHv_erf; levels = 10, linewidth = 1, color = :black)
ax_h  = Axis(gab[4,3], xlabel = L"x", 
            titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
            title = L"\text{(h) Filament horizontal viscosity}")
hm_h = heatmap!(ax_h, x, z, T¹bHv_gauss;
colormap = :balance, rasterize = true, 
colorrange = (-maximum(abs, T¹bHv_gauss), maximum(abs, T¹bHv_gauss)))
Colorbar(gab[4,4], hm_h)
contour!(ax_h, x, z, bHv_gauss; levels = 10, linewidth = 1, color = :black)
hideydecorations!(ax_b, ticks=false)
hideydecorations!(ax_d, ticks=false)
hideydecorations!(ax_f, ticks=false)
hideydecorations!(ax_h, ticks=false)
hidexdecorations!(ax_a, ticks=false)
hidexdecorations!(ax_b, ticks=false)
hidexdecorations!(ax_c, ticks=false)
hidexdecorations!(ax_d, ticks=false)
hidexdecorations!(ax_e, ticks=false)
hidexdecorations!(ax_f, ticks=false)
rowgap!(gab, 5)
colgap!(gab, 1, 3)
colgap!(gab, 2, 10)
colgap!(gab, 3, 3)
resize_to_layout!(figTb)

dt, t_max = 0.1, 15.0
t_all = 0:dt:t_max
nt = length(t_all)
dvd, dvv, dhd, dhv = zeros(nt,2), zeros(nt,2), zeros(nt,2), zeros(nt,2)
dvdg, dvvg, dhdg, dhvg = zeros(nt,2), zeros(nt,2), zeros(nt,2), zeros(nt,2)

ε_all = [0.0, 0.01, 0.03, 0.05, 0.07, 0.09]
nε = length(ε_all)
dγvd, dγvv, dγhd, dγhv = zeros(nt,nε), zeros(nt,nε), zeros(nt,nε), zeros(nt,nε)
dγvdg, dγvvg, dγhdg, dγhvg = zeros(nt,nε), zeros(nt,nε), zeros(nt,nε), zeros(nt,nε)
psi1s_erf = (psi1_erf_orig, psi1_erf_new)
psi1s_gauss = (psi1_gauss_orig, psi1_gauss_new)
Fbs_erf = (Fγb_erf, ϕzzz_erf, ϕzxx_erf)
Fvs_erf = (Fγv_erf, ϕxzz_erf, ϕxxx_erf)
Fbs_gauss = (Fγb_gauss, ϕzzz_gauss, ϕzxx_gauss)
Fvs_gauss = (Fγv_gauss, ϕxzz_gauss, ϕxxx_gauss)
for (i, t) in enumerate(t_all)
    bVi_erf, bHi_erf, bVvi_erf, bHvi_erf = compute_bVbH_fields(x, z, Ro, b_erf, ε, (ϕzzz_erf, ϕzxx_erf), (ϕxzz_erf, ϕxxx_erf), psi1s_erf, (ϕxx_erf, ϕxz_erf, ϕzz_erf), t)
    bVi_gauss, bHi_gauss, bVvi_gauss, bHvi_gauss = compute_bVbH_fields(x, z, Ro, b_gauss, ε, (ϕzzz_gauss, ϕzxx_gauss), (ϕxzz_gauss, ϕxxx_gauss), psi1s_gauss, (ϕxx_gauss, ϕxz_gauss, ϕzz_gauss), t)

    dbvddx = central_diff(bVi_erf, dx, 1)
    dbvvdx = central_diff(bVvi_erf, dx, 1)
    dbhddx = central_diff(bHi_erf, dx, 1)
    dbhvdx = central_diff(bHvi_erf, dx, 1)
    dbvdgdx = central_diff(bVi_gauss, dx, 1)
    dbvvgdx = central_diff(bVvi_gauss, dx, 1)
    dbhdgdx = central_diff(bHi_gauss, dx, 1)
    dbhvgdx = central_diff(bHvi_gauss, dx, 1)
    if t == 0
        global dbdx0 = maximum(abs, dbvddx)
        global dbdx0_gauss = maximum(abs, dbvdgdx)
    end
    for (j, n) in enumerate([1, nz])
        dvd[i,j] = dbdx0/maximum(abs, dbvddx[:,n])
        dvv[i,j] = dbdx0/maximum(abs, dbvvdx[:,n])
        dhd[i,j] = dbdx0/maximum(abs, dbhddx[:,n])
        dhv[i,j] = dbdx0/maximum(abs, dbhvdx[:,n])

        dvdg[i,j] = dbdx0_gauss/maximum(abs, dbvdgdx[:,n])
        dvvg[i,j] = dbdx0_gauss/maximum(abs, dbvvgdx[:,n])
        dhdg[i,j] = dbdx0_gauss/maximum(abs, dbhdgdx[:,n])
        dhvg[i,j] = dbdx0_gauss/maximum(abs, dbhvgdx[:,n])
    end

    for (j, εj) in enumerate(ε_all)
        bVij_erf, bHij_erf, bVvij_erf, bHvij_erf = compute_bVbH_fields(x, z, Ro, b_erf, εj, Fbs_erf, Fvs_erf, psi1s_erf, (ϕxx_erf, ϕxz_erf, ϕzz_erf), t; ψ¹ₚγ=psi1_erf_strain)
        bVij_gauss, bHij_gauss, bVvij_gauss, bHvij_gauss = compute_bVbH_fields(x, z, Ro, b_gauss, εj, Fbs_gauss, Fvs_gauss, psi1s_gauss, (ϕxx_gauss, ϕxz_gauss, ϕzz_gauss), t; ψ¹ₚγ=psi1_gauss_strain)

        dbvddx = central_diff(bVij_erf, dx, 1)
        dbvvdx = central_diff(bVvij_erf, dx, 1)
        dbhddx = central_diff(bHij_erf, dx, 1)
        dbhvdx = central_diff(bHvij_erf, dx, 1)
        dbvdgdx = central_diff(bVij_gauss, dx, 1)
        dbvvgdx = central_diff(bVvij_gauss, dx, 1)
        dbhdgdx = central_diff(bHij_gauss, dx, 1)
        dbhvgdx = central_diff(bHvij_gauss, dx, 1)

        dγvd[i,j] = dbdx0/maximum(abs, dbvddx[:,nz])
        dγvv[i,j] = dbdx0/maximum(abs, dbvvdx[:,nz])
        dγhd[i,j] = dbdx0/maximum(abs, dbhddx[:,nz])
        dγhv[i,j] = dbdx0/maximum(abs, dbhvdx[:,nz])
        dγvdg[i,j] = dbdx0_gauss/maximum(abs, dbvdgdx[:,nz])
        dγvvg[i,j] = dbdx0_gauss/maximum(abs, dbvvgdx[:,nz])
        dγhdg[i,j] = dbdx0_gauss/maximum(abs, dbhdgdx[:,nz])
        dγhvg[i,j] = dbdx0_gauss/maximum(abs, dbhvgdx[:,nz])
    end
        
end

using Makie
wcolors = Makie.wong_colors()
bcolors = Makie.to_colormap(:Blues_9)
# Create figure and axis
figts = Figure(size = (640, 320))
ax_b = Axis(figts[1, 1]; ylabel = L"d", xlabel = L"t", limits = ((0,t_max),nothing), titlealign = :left, title = L"\text{(a) Front}")
ax_b2 = Axis(figts[1, 2]; xlabel = L"t", limits = ((0,t_max),nothing), titlealign = :left, title = L"\text{(b) Filament}")
lines!(ax_b, t_all, dvd[:,2]; label = "vertical diffusivity")
lines!(ax_b, t_all, dvv[:,2]; label = "vertical viscosity")
lines!(ax_b, t_all, dhd[:,2]; label = "horizontal diffusivity")
lines!(ax_b, t_all, dhv[:,2]; label = "horizontal viscosity")
lines!(ax_b2, t_all, dvdg[:,2])
lines!(ax_b2, t_all, dvvg[:,2])
lines!(ax_b2, t_all, dhdg[:,2])
lines!(ax_b2, t_all, dhvg[:,2])
axislegend(ax_b, framevisible = false, position = :lt, 
           padding = (0f0, 0f0, 0f0, 0f0), patchlabelgap = 3)
resize_to_layout!(figts)

labels = [L"\epsilon/\lambda=0", L"\epsilon/\lambda=1/3", L"\epsilon/\lambda=3/3", L"\epsilon/\lambda=5/3", L"\epsilon/\lambda=7/3", L"\epsilon/\lambda=9/3"]
figts2 = Figure(size = (640, 900))
ax_a = Axis(figts2[1, 1]; ylabel = L"d", limits = ((0,t_max),nothing), titlealign = :left, title = L"\text{(a) Front vertical diffusivity}")
ax_b = Axis(figts2[1, 2]; limits = ((0,t_max),nothing), titlealign = :left, title = L"\text{(b) Filament vertical diffusivity}")
lines!(ax_a, t_all, dγvd[:,1], color=:black, label=L"\epsilon/\lambda=0")
lines!(ax_b, t_all, dγvdg[:,1], color=:black)
for i in 2:nε
    lines!(ax_a, t_all, dγvd[:,i], color=bcolors[end+2-i], label=labels[i])
    lines!(ax_b, t_all, dγvdg[:,i], color=bcolors[end+2-i])
end
axislegend(ax_a, framevisible = false, position = :rt, nbanks=2,
           padding = (0f0, 0f0, 0f0, 0f0), patchlabelgap = 3)
ax_c = Axis(figts2[2, 1]; ylabel = L"d", limits = ((0,t_max),nothing), titlealign = :left, title = L"\text{(c) Front vertical viscosity}")
ax_d = Axis(figts2[2, 2]; limits = ((0,t_max),nothing), titlealign = :left, title = L"\text{(d) Filament vertical viscosity}")
lines!(ax_c, t_all, dγvv[:,1], color=:black)
lines!(ax_d, t_all, dγvvg[:,1], color=:black)
for i in 2:nε
    lines!(ax_c, t_all, dγvv[:,i], color=bcolors[end+2-i])
    lines!(ax_d, t_all, dγvvg[:,i], color=bcolors[end+2-i])
end
ax_e = Axis(figts2[3, 1]; ylabel = L"d", limits = ((0,t_max),nothing), titlealign = :left, title = L"\text{(e) Front horizontal diffusivity}")
ax_f = Axis(figts2[3, 2]; limits = ((0,t_max),nothing), titlealign = :left, title = L"\text{(f) Filament horizontal diffusivity}")
lines!(ax_e, t_all, dγhd[:,1], color=:black)
lines!(ax_f, t_all, dγhdg[:,1], color=:black)
for i in 2:nε
    lines!(ax_e, t_all, dγhd[:,i], color=bcolors[end+2-i])
    lines!(ax_f, t_all, dγhdg[:,i], color=bcolors[end+2-i])
end
ax_g = Axis(figts2[4, 1]; xlabel = L"t", ylabel = L"d", limits = ((0,t_max),nothing), titlealign = :left, title = L"\text{(g) Front horizontal viscosity}")
ax_h = Axis(figts2[4, 2]; xlabel = L"t", limits = ((0,t_max),nothing), titlealign = :left, title = L"\text{(h) Filament horizontal viscosity}")
lines!(ax_g, t_all, dγhv[:,1], color=:black)
lines!(ax_h, t_all, dγhvg[:,1], color=:black)
for i in 2:nε
    lines!(ax_g, t_all, dγhv[:,i], color=bcolors[end+2-i])
    lines!(ax_h, t_all, dγhvg[:,i], color=bcolors[end+2-i])
end
hidexdecorations!(ax_a, grid = false)
hidexdecorations!(ax_b, grid = false)
hidexdecorations!(ax_c, grid = false)
hidexdecorations!(ax_d, grid = false)
hidexdecorations!(ax_e, grid = false)
hidexdecorations!(ax_f, grid = false)
resize_to_layout!(figts2)

Fblt = false
if Fblt
    figblt = Figure(size = (640, 800));
    gab = figblt[1, 1] = GridLayout()
    ax_a  = Axis(gab[1,1], ylabel = L"z",
                titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
                title = L"\text{(a) Front horizontal viscosity}")
    clims = (-maximum(abs, ϕxxx_erf), maximum(abs, ϕxxx_erf))
    hm_a = heatmap!(ax_a, x, z, ϕxxx_erf;
            colormap = :delta, rasterize = true, colorrange = clims)
    Colorbar(gab[1,2], hm_a)
    ax_b  = Axis(gab[1,3],  
                titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
                title = L"\text{(b) Filament horizontal viscosity}")
    hm_b = heatmap!(ax_b, x, z, ϕxxx_gauss;
            colormap = :delta, rasterize = true, 
            colorrange = (-maximum(abs, ϕxxx_gauss), maximum(abs, ϕxxx_gauss)))
    Colorbar(gab[1,4], hm_b)
    ax_c = Axis(gab[2,1], xlabel = L"x", ylabel = L"z", 
                titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
                title = L"\text{(c) Front vertical viscosity}")
    clims = (-maximum(abs, ϕxzz_erf), maximum(abs, ϕxzz_erf))
    hm_c = heatmap!(ax_c, x, z, ϕxzz_erf; colormap = :delta, rasterize = true, colorrange = clims)
    Colorbar(gab[2,2], hm_c)
    ax_d  = Axis(gab[2,3], xlabel = L"x", 
                titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
                title = L"\text{(d) Filament vertical viscosity}")
    hm_d = heatmap!(ax_d, x, z, ϕxzz_gauss; colormap = :delta, rasterize = true, 
    colorrange = (-maximum(abs, ϕxzz_gauss), maximum(abs, ϕxzz_gauss)))
    Colorbar(gab[2,4], hm_d)
    ax_e = Axis(gab[3,1], xlabel = L"x", ylabel = L"z", 
                titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
                title = L"\text{(e) Front horizontal diffusivity}")
    clims = (-maximum(abs, ϕzxx_erf), maximum(abs, ϕzxx_erf))
    hm_e = heatmap!(ax_e, x, z, ϕzxx_erf; colormap = :diff, rasterize = true, colorrange = clims)
    Colorbar(gab[3,2], hm_e)
    ax_f = Axis(gab[3,3], xlabel = L"x", 
                titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
                title = L"\text{(f) Filament horizontal diffusivity}")
    hm_f = heatmap!(ax_f, x, z, ϕzxx_gauss; colormap = :diff, rasterize = true, 
    colorrange = (-maximum(abs, ϕzxx_gauss), maximum(abs, ϕzxx_gauss)))
    Colorbar(gab[3,4], hm_f)
    ax_g = Axis(gab[4,1], xlabel = L"x", ylabel = L"z", 
                titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
                title = L"\text{(g) Front vertical diffusivity}")
    clims = (-maximum(abs, ϕzzz_erf), maximum(abs, ϕzzz_erf))
    hm_g = heatmap!(ax_g, x, z, ϕzzz_erf; colormap = :diff, rasterize = true, colorrange = clims)
    Colorbar(gab[4,2], hm_g)
    ax_h  = Axis(gab[4,3], xlabel = L"x", 
                titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
                title = L"\text{(h) Filament vertical diffusivity}")
    hm_h = heatmap!(ax_h, x, z, ϕzzz_gauss; colormap = :diff, rasterize = true, 
    colorrange = (-maximum(abs, ϕzzz_gauss), maximum(abs, ϕzzz_gauss)))
    Colorbar(gab[4,4], hm_h)
    hideydecorations!(ax_b, ticks=false)
    hideydecorations!(ax_d, ticks=false)
    hideydecorations!(ax_f, ticks=false)
    hideydecorations!(ax_h, ticks=false)
    hidexdecorations!(ax_a, ticks=false)
    hidexdecorations!(ax_b, ticks=false)
    hidexdecorations!(ax_c, ticks=false)
    hidexdecorations!(ax_d, ticks=false)
    hidexdecorations!(ax_e, ticks=false)
    hidexdecorations!(ax_f, ticks=false)
    rowgap!(gab, 5)
    colgap!(gab, 1, 3)
    colgap!(gab, 2, 10)
    colgap!(gab, 3, 3)
    resize_to_layout!(figblt)
end

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
    title = L"\text{(a) Front}~\psi^1_{\gamma}~\text{error}",
    xlabel = L"x", ylabel = L"z", 
)
max_abs_resid = maximum(abs, filter(!isnan, residual_erf))
hm_r1 = heatmap!(ax_r1, x, z, residual_erf; 
                 rasterize = true,colormap=:balance,colorrange=(-max_abs_resid, max_abs_resid))

ax_r2 = Axis(fig2[1, 3], titlealign = :left, 
    title = L"\text{(b) Filament}~\psi^1_{\gamma}~\text{error}",
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
display(fig0)
display(fig1)
display(fig2)
display(fig3)
display(figγ)
display(figv)
display(figTb)
display(figts)
display(figts2)
# display(figblt)

# save("b0v0_ff.pdf", fig0; pt_per_unit = 1)
# save("pVpH_ff.pdf", fig1; pt_per_unit = 1)
# save("psigamma_error.pdf", fig2; pt_per_unit = 1)
# save("psiHV_error.pdf", fig3; pt_per_unit = 1)
save("strain_ff.pdf", figγ; pt_per_unit = 1)
# save("blt_ff.pdf", figblt; pt_per_unit = 1)
save("vVvH_ff.pdf", figv; pt_per_unit = 1)
save("Tb_ff.pdf", figTb; pt_per_unit = 1)
save("ts_ff.pdf", figts; pt_per_unit = 1)
save("ts2_ff.pdf", figts2; pt_per_unit = 1)

checkF = false
if checkF
    function compute_F1_field(x, z, Ro, B_funcs; use_new_A_coeffs::Bool = false)
        nx, nz = length(x), length(z)
        F1_field = zeros(Float64, nx, nz)
        dB0_func, d2B0_func = B_funcs[1], B_funcs[2]
    
        for j in 1:nz
            z_val = z[j]
            X0_row = solve_X0_vec(x, z_val, Ro, dB0_func, d2B0_func)
            
            for i in 1:nx
                X0_val = X0_row[i]
                A0, A1, A2 = compute_coefficients(X0_val, Ro, B_funcs; use_new_A_coeffs, return_A=true)
                J0_val = 1 / (1 - Ro^2 * d2B0_func(X0_val) * (z_val - 0.5))
                F1_field[i, j] = J0_val^4 * (A0 + A1 * (z_val - 1/2) + A2 * (z_val - 1/2)^2)
            end
        end
    
        return F1_field
    end

    function central_diff(A, ds, dim)
        out = fill(NaN, size(A))
        if dim == 1 # d/dx for (nx, nz) matrix
            R = 2:size(A,1)-1; C = 1:size(A,2)
            @views @. out[R,C] = (A[R.+1,C] - A[R.-1,C]) / (2*ds)
            out[1,:] = (A[2,1:nz] - A[1,1:nz]) / ds
            out[end,:] = (A[end,1:nz] - A[end-1,1:nz]) / ds
        elseif dim == 2 # d/dz for (nx, nz) matrix
            R = 1:size(A,1); C = 2:size(A,2)-1
            @views @. out[R,C] = (A[R,C.+1] - A[R,C.-1]) / (2*ds)
            out[:,1] = (A[1:nx,2] - A[1:nx,1]) / ds
            out[:,end] = (A[1:nx,end] - A[1:nx,end-1]) / ds
        end
        return out
    end

    function compute_F1_cdiff_method(x, z, Ro, B0_func, B_funcs, dx, dz; F_type::Symbol)
        nx, nz = length(x), length(z)
        dB0, d2B0 = B_funcs[1], B_funcs[2]
        X0_grid = zeros(nx, nz)
        for j in 1:nz; X0_grid[:, j] = solve_X0_vec(x, z[j], Ro, dB0, d2B0); end
        B0_grid = B0_func.(X0_grid)

        F = if F_type == :orig # F = ∂³B₀/∂x∂z²
            central_diff(central_diff(central_diff(B0_grid, dz, 2), dz, 2), dx, 1)
        else # F_type == :new, F = ∂³B₀/∂x³
            central_diff(central_diff(central_diff(B0_grid, dx, 1), dx, 1), dx, 1)
        end
        F[isnan.(F)] .= 0.0 # Set boundary NaNs from differentiation to zero

        return F
    end

    F1_erf_orig = compute_F1_field(x, z, Ro, B_funcs_erf)
    F1_gauss_orig = compute_F1_field(x, z, Ro, B_funcs_gauss)
    F1_erf_new = compute_F1_field(x, z, Ro, B_funcs_erf; use_new_A_coeffs=true)
    F1_gauss_new = compute_F1_field(x, z, Ro, B_funcs_gauss; use_new_A_coeffs=true)
    figF1 = Figure(size = (640, 450));
    gab = figF1[1, 1] = GridLayout()
    ax_a  = Axis(gab[1,1], ylabel = L"z",
                titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
                title = L"\text{(a) Front}~F_{V}")
    clims = (-maximum(abs, F1_erf_orig), maximum(abs, F1_erf_orig))
    hm_a = heatmap!(ax_a, x, z, F1_erf_orig;
            colormap = :PuOr, rasterize = true, colorrange = clims)
    Colorbar(gab[1,2], hm_a)
    ax_b  = Axis(gab[1,3],  
                titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
                title = L"\text{(b) Filament}~F_{V}")
    hm_b = heatmap!(ax_b, x, z, F1_gauss_orig;
            colormap = :PuOr, rasterize = true, 
            colorrange = (-maximum(abs, F1_gauss_orig), maximum(abs, F1_gauss_orig)))
    Colorbar(gab[1,4], hm_b)
    ax_c = Axis(gab[2,1], xlabel = L"x", ylabel = L"z", 
                titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
                title = L"\text{(c) Front}~F_{H}")
    clims = (-maximum(abs, F1_erf_new), maximum(abs, F1_erf_new))
    hm_c = heatmap!(ax_c, x, z, F1_erf_new; colormap = :PuOr, rasterize = true, colorrange = clims)
    Colorbar(gab[2,2], hm_c)
    ax_d  = Axis(gab[2,3], xlabel = L"x", 
                titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
                title = L"\text{(d) Filament}~F_{H}")
    hm_d = heatmap!(ax_d, x, z, F1_gauss_new; colormap = :PuOr, rasterize = true, 
    colorrange = (-maximum(abs, F1_gauss_new), maximum(abs, F1_gauss_new)))
    Colorbar(gab[2,4], hm_d)
    hideydecorations!(ax_b, ticks=false)
    hideydecorations!(ax_d, ticks=false)
    hidexdecorations!(ax_a, ticks=false)
    hidexdecorations!(ax_b, ticks=false)
    rowgap!(gab, 5)
    colgap!(gab, 1, 3)
    colgap!(gab, 2, 10)
    colgap!(gab, 3, 3)
    resize_to_layout!(figF1)

    F1_erf_orig_cdiff = compute_F1_cdiff_method(x, z, Ro, B0_erf, B_funcs_erf, dx, dz, F_type=:orig)
    F1_gauss_orig_cdiff = compute_F1_cdiff_method(x, z, Ro, B0_gauss, B_funcs_gauss, dx, dz, F_type=:orig)
    F1_erf_new_cdiff = compute_F1_cdiff_method(x, z, Ro, B0_erf, B_funcs_erf, dx, dz, F_type=:new)
    F1_gauss_new_cdiff = compute_F1_cdiff_method(x, z, Ro, B0_gauss, B_funcs_gauss, dx, dz, F_type=:new)

    figF2 = Figure(size = (640, 450));
    gab = figF2[1, 1] = GridLayout()
    ax_a  = Axis(gab[1,1], ylabel = L"z",
                titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
                title = L"\text{(a) Front}~F^1_{V,p}")
    clims = (-maximum(abs, F1_erf_orig_cdiff), maximum(abs, F1_erf_orig_cdiff))
    hm_a = heatmap!(ax_a, x, z, F1_erf_orig_cdiff;
            colormap = :PuOr, rasterize = true, colorrange = clims)
    Colorbar(gab[1,2], hm_a)
    ax_b  = Axis(gab[1,3],  
                titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
                title = L"\text{(b) Filament}~F^1_{V,p}")
    hm_b = heatmap!(ax_b, x, z, F1_gauss_orig_cdiff;
            colormap = :PuOr, rasterize = true, 
            colorrange = (-maximum(abs, F1_gauss_orig_cdiff), maximum(abs, F1_gauss_orig_cdiff)))
    Colorbar(gab[1,4], hm_b)
    ax_c = Axis(gab[2,1], xlabel = L"x", ylabel = L"z", 
                titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
                title = L"\text{(c) Front}~F^1_{H,p}")
    clims = (-maximum(abs, F1_erf_new_cdiff), maximum(abs, F1_erf_new_cdiff))
    hm_c = heatmap!(ax_c, x, z, F1_erf_new_cdiff; colormap = :PuOr, rasterize = true, colorrange = clims)
    Colorbar(gab[2,2], hm_c)
    ax_d  = Axis(gab[2,3], xlabel = L"x", 
                titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
                title = L"\text{(d) Filament}~F^1_{H,p}")
    hm_d = heatmap!(ax_d, x, z, F1_gauss_new_cdiff; colormap = :PuOr, rasterize = true, 
    colorrange = (-maximum(abs, F1_gauss_new_cdiff), maximum(abs, F1_gauss_new_cdiff)))
    Colorbar(gab[2,4], hm_d)
    hideydecorations!(ax_b, ticks=false)
    hideydecorations!(ax_d, ticks=false)
    hidexdecorations!(ax_a, ticks=false)
    hidexdecorations!(ax_b, ticks=false)
    rowgap!(gab, 5)
    colgap!(gab, 1, 3)
    colgap!(gab, 2, 10)
    colgap!(gab, 3, 3)
    resize_to_layout!(figF2)
    display(figF1)
    display(figF2)
end

println("Script finished.")