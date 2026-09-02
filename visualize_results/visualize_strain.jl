################################################################################
#  FIGURE 1  –  Similarity solution of surface-gradient frontogenesis
#  Julia 1.10+, CairoMakie ≥ 0.11, SpecialFunctions ≥ 2.3
################################################################################
using CairoMakie             # high-quality vector / raster plotting
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
γ₁         = 0.03           # strain rate in X
γ₂         = 0.00           # strain rate in solution
ϵ          = 0.00           # small parameter for geostrophic departure
t_max      = 20.0           # blow-up time
t_early    = 0.4            # diagnostic snapshot
x_min,x_max = -4, 4
z_min,z_max =  0.0, 1.0

# computational grids (moderate resolution; raise if desired)
nx, nz = 201, 101
x = collect(range(x_min, x_max; length = nx))   # <- make it a Vector!
z = collect(range(z_min, z_max; length = nz))   # (z isn’t mutated, but keep symmetrical)
dx = x[2] - x[1]

# ------------------------------------------------------------------
# 2.  Background buoyancy and its derivatives
# ------------------------------------------------------------------
invs2π = 1 / √(2π)

B0(X)   = 0.5 .* erf.(X ./ √2)
dB0(X)  = @. exp(-0.5 * X^2) * invs2π
d2B0(X) = @. -X * exp(-0.5 * X^2) * invs2π

# ------------------------------------------------------------------
# 3.  Vectorised Newton solve for the mapping X(x,z,T)
# ------------------------------------------------------------------
function solve_X_vec(x_vec, z_val, T; γ₁ = 0, γ₂ = 0, ϵ = 0, maxiter = 50, tol = 1e-12)
    eT = exp(γ₂*T)
    X  = exp(γ₁*T) .* x_vec                     # this is a Vector (mutable)
    A = (eT-ϵ*cos(sqrt(1-γ₂^2)*T)+γ₂*(ϵ-2)*sin(sqrt(1-γ₂^2)*T)/sqrt(1-γ₂^2))
    for _ in 1:maxiter
        Bp  = dB0(X)
        Bpp = d2B0(X)
        f   = X .- exp(γ₁*T) .* (x_vec .+ Ro^2 * (z_val - 0.5) * A * Bp)
        df  = 1 .- exp(γ₁*T) * Ro^2 * (z_val - 0.5) * A * Bpp 
        dX  = f ./ df
        X  .-= dX
        maximum(abs, dX) < tol && break
    end
    return X
end

# ------------------------------------------------------------------
# 4.  Inviscid fields at any given time T
# ------------------------------------------------------------------
function fields_at_time(T; γ₁ = 0, γ₂ = 0, ϵ = 0, sg = false)
    eT= exp(γ₂*T)

    b = Matrix{Float64}(undef, nz, nx)
    v = similar(b)
    w = similar(b)
    u = similar(b)
    ψ = similar(b)

    if sg
        A = eT
        B = γ₂ * eT
        C = 2γ₂ * eT
    else
        A = (eT-ϵ*cos(sqrt(1-γ₂^2)*T)+γ₂*(ϵ-2)*sin(sqrt(1-γ₂^2)*T)/sqrt(1-γ₂^2))
        B = (γ₂ * (eT-cos(sqrt(1-γ₂^2)*T))+(ϵ-2γ₂^2)*sin(sqrt(1-γ₂^2)*T)/2/sqrt(1-γ₂^2))
        C = (2γ₂ * (eT-cos(sqrt(1-γ₂^2)*T))+(ϵ-2γ₂^2)*sin(sqrt(1-γ₂^2)*T)/sqrt(1-γ₂^2))
    end
        
    for (k, zz) in enumerate(z)
        Xvals = solve_X_vec(x, zz, T; γ₁, γ₂, ϵ)

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

function sgfields_at_time(T; γ₁ = 0, γ₂ = 0, ϵ = 0)
    eT= exp(γ₂*T)

    b = Matrix{Float64}(undef, nz, nx)
    v = similar(b)
    w = similar(b)
    u = similar(b)
    ψ = similar(b)

    A = (eT)
    B = (γ₂ * (eT-cos(sqrt(1-γ₂^2)*T))+(ϵ-2γ₂^2)*sin(sqrt(1-γ₂^2)*T)/2/sqrt(1-γ₂^2))
    C = (2γ₂ * (eT-cos(sqrt(1-γ₂^2)*T))+(ϵ-2γ₂^2)*sin(sqrt(1-γ₂^2)*T)/sqrt(1-γ₂^2))
    for (k, zz) in enumerate(z)
        Xvals = solve_X_vec(x, zz, T; γ₁, γ₂, ϵ)

        Bp  = dB0(Xvals)
        Bpp = d2B0(Xvals)

        dvdX = Bpp .* (zz - 0.5) * A
        dvdZ = Bp * A
        v[k,:] =  Bp  .* (zz - 0.5) * A                            # (b)
        w[k,:] =  Ro * B .* Bpp .* eT .* zz*(zz-1) ./ (1 .- Ro^2 * eT .* dvdX) # (c)
        u[k,:] = -Ro * (C * Bp .* (zz - 0.5) .+ Ro * w[k,:] .* dvdZ)   # (d)
        ψ[k,:] = -Ro * B .* Bp .* eT .* zz*(zz-1)

        intZ   = 0.5*zz^2 - 0.5*zz
        Δb     = -(Ro^2)/(Fr^2) * A .* Bpp .* intZ
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

# -------------------- early-time fields for panels (b) & (c) -------------
b_e, v_e, ψ_e, _, _ = fields_at_time(t_early)

# ------------------------------------------------------------------
# 5.  Plotting with CairoMakie
# ------------------------------------------------------------------
set_theme!(theme_latexfonts(), fontsize=12,figure_padding = 10)            # optional stylistic baseline

dt = 0.2
t_all = 0:dt:t_max
nt = length(t_all)
ψf, vf, bf, uf, wf = zeros(nt,nx,nz),zeros(nt,nx,nz),zeros(nt,nx,nz),zeros(nt,nx,nz),zeros(nt,nx,nz)
ψ1, v1, b1, u1, w1 = zeros(nt,nx,nz),zeros(nt,nx,nz),zeros(nt,nx,nz),zeros(nt,nx,nz),zeros(nt,nx,nz)
df, d1, uff, u1f, vff, v1f, wff, w1f = zeros(nt), zeros(nt), zeros(nt), zeros(nt), zeros(nt), zeros(nt), zeros(nt), zeros(nt)
dsg, usg, vsg, wsg = zeros(nt), zeros(nt), zeros(nt), zeros(nt)
for i = 1:nt
    bf[i,:,:], vf[i,:,:], ψf[i,:,:], wf[i,:,:], uf[i,:,:] = fields_at_time(t_all[i]; γ₁, γ₂=γ₁)
    bsi, vsi, ψsi, wsi, usi = fields_at_time(t_all[i]; γ₁, γ₂=γ₁, sg = true)
    b1[i,:,:], v1[i,:,:], ψ1[i,:,:], w1[i,:,:], u1[i,:,:] = fields_at_time(t_all[i]; γ₁)
    dbfdx = (bf[i,3:end,:].-bf[i,1:end-2,:])./(2*dx)
    df[i] = invs2π/maximum(abs, dbfdx)
    argmax_dbfdx = argmax(abs.(dbfdx[:,1]))
    uff[i] = uf[i,argmax_dbfdx,1]
    vff[i] = vf[i,argmax_dbfdx,1]
    argmax_dbfdx = argmax(abs.(dbfdx[:,51]))
    wff[i] = wf[i,argmax_dbfdx,51]
    dbsdx = (bsi[3:end,:].-bsi[1:end-2,:])./(2*dx)
    dsg[i] = invs2π/maximum(abs, dbsdx)
    argmax_dbsdx = argmax(abs.(dbsdx[:,1]))
    usg[i] = usi[argmax_dbsdx,1]
    vsg[i] = vsi[argmax_dbsdx,1]
    argmax_dbsdx = argmax(abs.(dbsdx[:,51]))
    wsg[i] = wsi[argmax_dbsdx,51]
    db1dx = (b1[i,3:end,:].-b1[i,1:end-2,:])./(2*dx)
    d1[i] = invs2π/maximum(abs, db1dx)
    argmax_db1dx = argmax(abs.(db1dx[:,1]))
    u1f[i] = u1[i,argmax_db1dx,1]
    v1f[i] = v1[i,argmax_db1dx,1]
    argmax_db1dx = argmax(abs.(db1dx[:,51]))
    w1f[i] = w1[i,argmax_db1dx,51]
end

# Create figure and axis
fig = Figure(size = (640, 640))
ax_b = Axis(fig[1, 1]; ylabel = L"d", limits = ((0,20),nothing), titlealign = :left, title = L"\text{(a)}")
ax_u = Axis(fig[2, 1]; ylabel = L"u", limits = ((0,20),nothing), titlealign = :left, title = L"\text{(b)}")
ax_v = Axis(fig[3, 1]; ylabel = L"v", limits = ((0,20),nothing), titlealign = :left, title = L"\text{(c)}")
ax_w = Axis(fig[4, 1]; xlabel = L"t", ylabel = L"w", limits = ((0,20),nothing), titlealign = :left, title = L"\text{(d)}")
lines!(ax_b, t_all, d1; label = "advection")
lines!(ax_b, t_all, df; label = "full")
lines!(ax_b, t_all, dsg; label = "SG")
lines!(ax_u, t_all, u1f)
lines!(ax_u, t_all, uff)
lines!(ax_u, t_all, usg)
lines!(ax_v, t_all, v1f)
lines!(ax_v, t_all, vff)
lines!(ax_v, t_all, vsg)
lines!(ax_w, t_all, w1f)
lines!(ax_w, t_all, wff)
lines!(ax_w, t_all, wsg)
axislegend(ax_b, framevisible = false, position = :ct, orientation = :horizontal,
            padding = (0f0, 0f0, 0f0, 0f0), patchlabelgap = 3)
hidexdecorations!(ax_b, ticks=false)
hidexdecorations!(ax_u, ticks=false)
hidexdecorations!(ax_v, ticks=false)
for i = 1:3
    rowgap!(fig.layout, i, Relative(0))
end
resize_to_layout!(fig)
save("st_duvw.pdf", fig; pt_per_unit = 1)