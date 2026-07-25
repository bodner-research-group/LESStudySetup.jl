################################################################################
#  FIGURE 1  –  Similarity solution of surface-gradient frontogenesis
#  Julia 1.10+, CairoMakie ≥ 0.11, SpecialFunctions ≥ 2.3
################################################################################
using CairoMakie, SixelTerm  # high-quality vector / raster plotting
using SpecialFunctions: erf  # Gaussian error function
using Printf
using MathTeXEngine
set_theme!(theme_latexfonts(), fontsize=12,figure_padding = 10)

# ------------------------------------------------------------------
# 1.  Non-dimensional parameters   (Appendix A.1 and γ=0.1)
# ------------------------------------------------------------------
const Ro   = 1.0
const Bu   = 0.1
const Fr   = Ro / Bu
γ          = 0.03           # growth-rate parameter requested
ϵ          = 0.00           # small parameter for geostrophic departure
t_max      = 10.0           # blow-up time
t_early    = 0.4            # diagnostic snapshot
x_min,x_max = -4, 4
z_min,z_max =  0.0, 1.0

# computational grids (moderate resolution; raise if desired)
nx, nz = 301, 101
x = collect(range(x_min, x_max; length = nx))   # <- make it a Vector!
z = collect(range(z_min, z_max; length = nz))   # (z isn’t mutated, but keep symmetrical)
dx = x[2] - x[1]
dz = z[2] - z[1]

# ------------------------------------------------------------------
# 2.  Background buoyancy and its derivatives
# ------------------------------------------------------------------
invs2π = 1 / √(2π)

B0(X)   = 0.5 .* erf.(X ./ √2)
dB0(X)  = @. exp(-0.5 * X^2) * invs2π
d2B0(X) = @. -X * exp(-0.5 * X^2) * invs2π

# B0(X)   = -0.5 .* exp.(-(X ./ √2).^2)
# dB0(X)  = @. -X .* B0(X)
# d2B0(X) = @. (X^2 .- 1) .* B0(X)

# ------------------------------------------------------------------
# 3.  Vectorised Newton solve for the mapping X(x,z,T)
# ------------------------------------------------------------------
function solve_X_vec(x_vec, z_val, T; γ = 0, ϵ = 0, maxiter = 50, tol = 1e-12)
    eT = exp(γ*T)
    X  = eT .* x_vec                     # this is a Vector (mutable)
    A = (eT-ϵ*cos(sqrt(1-γ^2)*T)+γ*(ϵ-2)*sin(sqrt(1-γ^2)*T)/sqrt(1-γ^2))
    for _ in 1:maxiter
        Bp  = dB0(X)
        Bpp = d2B0(X)
        f   = X .- eT .* (x_vec .+ Ro^2 * (z_val - 0.5) * A * Bp)
        df  = 1 .- eT * Ro^2 * (z_val - 0.5) * A * Bpp 
        dX  = f ./ df
        X  .-= dX
        maximum(abs, dX) < tol && break
    end
    return X
end

# ------------------------------------------------------------------
# 4.  Inviscid fields at any given time T
# ------------------------------------------------------------------
function fields_at_time(T; γ = 0, ϵ = 0, sg = false)
    eT= exp(γ*T)

    b = Matrix{Float64}(undef, nz, nx)
    v = similar(b)
    w = similar(b)
    u = similar(b)
    ψ = similar(b)

    if sg
        A = eT
        B = γ * eT
        C = 2γ * eT
    else
        A = (eT-ϵ*cos(sqrt(1-γ^2)*T)+γ*(ϵ-2)*sin(sqrt(1-γ^2)*T)/sqrt(1-γ^2))
        B = (γ * (eT-cos(sqrt(1-γ^2)*T))+(ϵ-2γ^2)*sin(sqrt(1-γ^2)*T)/2/sqrt(1-γ^2))
        C = (2γ * (eT-cos(sqrt(1-γ^2)*T))+(ϵ-2γ^2)*sin(sqrt(1-γ^2)*T)/sqrt(1-γ^2))
    end
        
    for (k, zz) in enumerate(z)
        Xvals = solve_X_vec(x, zz, T; γ, ϵ)

        Bp  = dB0(Xvals)
        Bpp = d2B0(Xvals)

        dvdX = Bpp .* (zz - 0.5) * A
        dvdZ = Bp * A
        v[k,:] =  Bp  .* (zz - 0.5) * A                            # (b)
        w[k,:] =  Ro * B .* Bpp .* eT .* zz*(zz-1) ./ (1 .- Ro^2 * eT .* dvdX) # (c)
        u[k,:] = -Ro * (C * Bp .* (zz - 0.5) .+ Ro * w[k,:] .* dvdZ)   # (d)
        ψ[k,:] = -Ro * B .* Bp .* eT .* zz*(zz-1)

        #intZ   = 0.5*zz^2 - 0.5*zz
        #Δb     = -(Ro^2)/(Fr^2) * A .* Bpp .* intZ
        b[k,:] =  B0(Xvals) #.+ 1/Fr^2*zz .+ Δb                               # (a)
    end

    #ψ = -cumsum(w0; dims = 2) * dx      # integrate w = −∂ψ/∂x, ψ=0 at x_min
    return b', v', ψ', w', u'
end


function firstfs_at_time(T)

    b = Matrix{Float64}(undef, nz, nx)
    v = similar(b)
    w = similar(b)
    u = similar(b)
    ψ = similar(b)

    for (k, zz) in enumerate(z)
        Xvals = solve_X_vec(x, zz, T)

        Bp  = dB0(Xvals)
        Bpp = d2B0(Xvals)

        dvdX = Bpp .* (zz - 0.5)
        dvdZ = Bp 
        v[k,:] =  (zz - 0.5) * (T*(Bpp .* Xvals .+ Bp) .- 2*sin(T)*Bp) ./ (1 .- Ro^2 * dvdX)
        w[k,:] =  Ro * Bpp .* zz*(zz-1)*(1-cos(T)) ./ (1 .- Ro^2 * dvdX) 
        ψ[k,:] = -Ro * Bp .* zz*(zz-1)*(1-cos(T))
        u[k,:] = -Ro * Bp .* (2*zz-1)*(1-cos(T)) .- Ro^2 * w[k,:] .* dvdZ   

        b[k,:] =  Bp .* (T*Xvals .+ Ro^2 * v[k,:])
    end
    return b', v', ψ', w', u'
end

function central_diff(A, ds, dim)
    out = fill(NaN, size(A))
    if dim == 1 # d/dt for (nt, nx, nz) matrix
        R = 2:size(A,1)-1; C = 1:size(A,2); K = 1:size(A,3)
        @views @. out[R,C,K] = (A[R.+1,C,K] - A[R.-1,C,K]) / (2*ds)
        out[1, :, :] = (A[2, :, :] - A[1, :, :]) / ds
        out[end, :, :] = (A[end, :, :] - A[end-1, :, :]) / ds
    elseif dim == 2 # d/dx for (nt, nx, nz) matrix
        R = 1:size(A,1); C = 2:size(A,2)-1; K = 1:size(A,3)
        @views @. out[R,C,K] = (A[R,C.+1,K] - A[R,C.-1,K]) / (2*ds)
        out[:, 1, :] = (A[:, 2, :] - A[:, 1, :]) / ds
        out[:, end, :] = (A[:, end, :] - A[:, end-1, :]) / ds
    elseif dim == 3 # d/dz for (nt, nx, nz) matrix
        R = 1:size(A,1); C = 1:size(A,2); K = 2:size(A,3)-1
        @views @. out[R,C,K] = (A[R,C,K.+1] - A[R,C,K.-1]) / (2*ds)
        out[:, :, 1] = (A[:, :, 2] - A[:, :, 1]) / ds
        out[:, :, end] = (A[:, :, end] - A[:, :, end-1]) / ds
    end
    return out
end

# -------------------- early-time fields for panels (b) & (c) -------------
b0, _, _, _, _ = fields_at_time(t_early)
b0 = reshape(b0, (1, nx, nz))

# ------------------------------------------------------------------
# 5.  Plotting with CairoMakie
# ------------------------------------------------------------------
set_theme!(theme_latexfonts(), fontsize=12,figure_padding = 10)            # optional stylistic baseline

dt = 0.1
t_all = Array(0:dt:t_max)
nt = length(t_all)
ψf, vf, bf, uf, wf = zeros(nt,nx,nz),zeros(nt,nx,nz),zeros(nt,nx,nz),zeros(nt,nx,nz),zeros(nt,nx,nz)
ψ1, v1, b1, u1, w1 = zeros(nt,nx,nz),zeros(nt,nx,nz),zeros(nt,nx,nz),zeros(nt,nx,nz),zeros(nt,nx,nz)
for i = 1:nt
    bf[i,:,:], vf[i,:,:], ψf[i,:,:], wf[i,:,:], uf[i,:,:] = fields_at_time(t_all[i]; γ=γ)
    b1[i,:,:], v1[i,:,:], ψ1[i,:,:], w1[i,:,:], u1[i,:,:] = firstfs_at_time(t_all[i])
end

dbfdx = central_diff(bf, dx, 2)
db0dx = central_diff(b0, dx, 2)
db1dx = central_diff(b1, dx, 2)
Tbf = 0.5 * central_diff(dbfdx.^2, dt, 1) 
Tbfmax = vec(maximum(abs, Tbf[:, :, end], dims = 2))
T¹b = db0dx .* central_diff(db1dx, dt, 1) 
T¹b .+= 0.5*((u1 .- reshape(x, (1,nx,1))) .* central_diff(db0dx.^2, dx, 2) .+ w1 .* central_diff(db0dx.^2, dz, 3))
T¹bmax = vec(maximum(abs, T¹b[:, :, end], dims = 2))

fig = Figure(size = (640, 450));
g4 = fig[1, 1] = GridLayout()
axa = Axis(g4[1, 1], ylabel = L"\text{Frontal tendency}", titlealign = :left, title = L"(a)")
axb = Axis(g4[1, 2], titlealign = :left, title = L"(b)")
axc = Axis(g4[2, 1], xlabel = L"\text{Time}", ylabel = L"\text{Frontal tendency}", titlealign = :left, title = L"(c)")
axd = Axis(g4[2, 2], xlabel = L"\text{Time}", titlealign = :left, title = L"(d)")

lines!(axa, t_all, γ*T¹bmax)
lines!(axb, t_all, γ*T¹bmax)
lines!(axc, t_all, γ*T¹bmax)
lines!(axd, t_all, γ*T¹bmax)
lines!(axa, t_all, Tbfmax)
lines!(axb, t_all, Tbfmax)
lines!(axc, t_all, Tbfmax)
lines!(axd, t_all, Tbfmax)

fig

d1 = invs2π ./ maximum(abs, γ*db1dx .+ db0dx, dims=(2,3))