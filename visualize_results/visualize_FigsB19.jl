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
γ          = 0.03           # growth-rate parameter requested
ϵ          = 0.00           # small parameter for geostrophic departure
t_max      = 15.0           # blow-up time
t_early    = 0.4            # diagnostic snapshot
x_min,x_max = -4, 4
z_min,z_max =  0.0, 1.0

# computational grids (moderate resolution; raise if desired)
nx, nz = 301, 101
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

# B0(X)   = -0.5 .* exp.(-(X ./ √2).^2)
# dB0(X)  = @. -X .* B0(X)
# d2B0(X) = @. (X^2 .- 1) .* B0(X)

# ------------------------------------------------------------------
# 3.  Vectorised Newton solve for the mapping X(x,z,T)
# ------------------------------------------------------------------
function solve_X_vec(x_vec, z_val, T; γ = 0, ϵ = 0, maxiter = 500, tol = 1e-12)
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

function sgfields_at_time(T; γ = 0, ϵ = 0)
    eT= exp(γ*T)

    b = Matrix{Float64}(undef, nz, nx)
    v = similar(b)
    w = similar(b)
    u = similar(b)
    ψ = similar(b)

    A = (eT)
    B = (γ * (eT-cos(sqrt(1-γ^2)*T))+(ϵ-2γ^2)*sin(sqrt(1-γ^2)*T)/2/sqrt(1-γ^2))
    C = (2γ * (eT-cos(sqrt(1-γ^2)*T))+(ϵ-2γ^2)*sin(sqrt(1-γ^2)*T)/sqrt(1-γ^2))
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

# Figure 2 – b-v geostrophic balance
if ϵ<1e-6
    if γ < 1e-6
        x_min,x_max = -4, 4
        fig1 = Figure(size = (640, 240));
        gab = fig1[1, 1] = GridLayout()
        ax_a  = Axis(gab[1,1], xlabel = L"x", ylabel = L"z",
                     titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
                     title = L"\text{(a) Buoyancy}~b^0")
        clims = (-maximum(abs, b_e), maximum(abs, b_e))#1.1 .* extrema(b_e)
        hm_a = heatmap!(ax_a, x, z, b_e;
                colormap = :diff, rasterize = true, colorrange = clims)
        Colorbar(gab[1,2], hm_a)
        ax_b  = Axis(gab[1,3], xlabel = L"x", 
                     titlealign = :left, limits = ((x_min, x_max), (z_min, z_max)),
                     title = L"\text{(b) Along-front velocity}~v^0=v_g")
        hm_b = heatmap!(ax_b, x, z, v_e;
                colormap = :delta, rasterize = true, 
                colorrange = (-maximum(abs, v_e), maximum(abs, v_e)))
        Colorbar(gab[1,4], hm_b)
        contour!(ax_b, x, z, b_e; levels = 10, linewidth = 1, color = :black)
        hideydecorations!(ax_b, ticks=false)
        colgap!(gab, 1, 3)
        colgap!(gab, 2, 10)
        colgap!(gab, 3, 3)
        resize_to_layout!(fig1)
        save("figure3_bgvg.pdf", fig1; pt_per_unit = 1)
        #save("figure2_bgvg.png", fig1)
    else
        dt = 0.1
        t_all = 0:dt:t_max
        nt = length(t_all)
        ψf, vf, bf, uf, wf = zeros(nt,nx,nz),zeros(nt,nx,nz),zeros(nt,nx,nz),zeros(nt,nx,nz),zeros(nt,nx,nz)
        ψ1, v1, b1, u1, w1 = zeros(nt,nx,nz),zeros(nt,nx,nz),zeros(nt,nx,nz),zeros(nt,nx,nz),zeros(nt,nx,nz)
        df, d1, uff, u1f, vff, v1f, wff, w1f = zeros(nt), zeros(nt), zeros(nt), zeros(nt), zeros(nt), zeros(nt), zeros(nt), zeros(nt)
        dsg, usg, vsg, wsg = zeros(nt), zeros(nt), zeros(nt), zeros(nt)
        for i = 1:nt
            bf[i,:,:], vf[i,:,:], ψf[i,:,:], wf[i,:,:], uf[i,:,:] = fields_at_time(t_all[i]; γ)
            bsi, vsi, ψsi, wsi, usi = fields_at_time(t_all[i]; γ, sg = true)
            b1[i,:,:], v1[i,:,:], ψ1[i,:,:], w1[i,:,:], u1[i,:,:] = firstfs_at_time(t_all[i])
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
            db1dx = (γ*(b1[i,3:end,:].-b1[i,1:end-2,:]).+(b_e[3:end,:].-b_e[1:end-2,:]))./(2*dx)
            d1[i] = invs2π/maximum(abs, db1dx)
            argmax_db1dx = argmax(abs.(db1dx[:,1]))
            u1f[i] = γ * u1[i,argmax_db1dx,1]
            v1f[i] = γ * v1[i,argmax_db1dx,1] + v_e[argmax_db1dx,1]
            v1[i,:,:] = γ * v1[i,:,:] .+ v_e
            argmax_db1dx = argmax(abs.(db1dx[:,51]))
            w1f[i] = γ * w1[i,argmax_db1dx,51]
        end

        # Create figure and axis
        fig = Figure(size = (640, 640))
        ax_b = Axis(fig[1, 1]; ylabel = L"d", limits = ((0,t_max),nothing), titlealign = :left, title = L"\text{(a)}")
        ax_u = Axis(fig[2, 1]; ylabel = L"u", limits = ((0,t_max),nothing), titlealign = :left, title = L"\text{(b)}")
        ax_v = Axis(fig[3, 1]; ylabel = L"v", limits = ((0,t_max),nothing), titlealign = :left, title = L"\text{(c)}")
        ax_w = Axis(fig[4, 1]; xlabel = L"t", ylabel = L"w", limits = ((0,t_max),nothing), titlealign = :left, title = L"\text{(d)}")
        lines!(ax_b, t_all, d1; label = "linear")
        lines!(ax_b, t_all, df; label = "ST13")
        lines!(ax_b, t_all, dsg; label = "HB72")
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
        rowgap!(fig.layout, 0.2)
        resize_to_layout!(fig)
        save("duvw_front.pdf", fig; pt_per_unit = 1)
        println("Time series saved to front_duvw.pdf")

        # alphabet = [letter for letter in 'a':'z'];
        # fig = Figure(size = (640, 750))
        # axis_kwargs = (titlealign = :left, titlefont=texfont(), limits = ((-4,4),(0,1)))
        # for i = 1:4
        #     ax_1 = Axis(fig[i, 1]; ylabel = L"z", title = "("*alphabet[2*(i-1)+1]*") ", axis_kwargs...)
        #     ax_2 = Axis(fig[i, 2]; title = "("*alphabet[2*(i-1)+2]*") ", axis_kwargs...)
        #     if i < 4
        #         hidexdecorations!(ax_1, ticks=false)
        #         hidexdecorations!(ax_2, ticks=false)
        #     else
        #         ax_1.xlabel = L"x"
        #         ax_2.xlabel = L"x"
        #     end
        #     hideydecorations!(ax_2, ticks=false) 
        #     idx = 1+Int((i-1)*5÷dt)
        #     cphi = (-0.02, 0.02)
        #     hm_1 = heatmap!(ax_1, x, z, γ * ψ1[idx,:,:], colorrange = cphi, colormap = :PuOr)
        #     hm_2 = heatmap!(ax_2, x, z, ψf[idx,:,:], colorrange = cphi, colormap = :PuOr)
        #     Colorbar(fig[i, 3], hm_1)
        #     clevels=-0.4:0.1:0.4
        #     ct_1 = contour!(ax_1, x, z, b_e .+ γ * b1[idx,:,:]; levels = clevels, linewidth = 1, color = :black)
        #     ct_2 = contour!(ax_2, x, z, bf[idx,:,:];levels = clevels, linewidth = 1, color = :black)
        # end
        # for i = 1:3
        #     rowgap!(fig.layout, i, Relative(0))
        # end
        # colgap!(fig.layout, 2, Relative(0))
        # resize_to_layout!(fig)
        # save("front_phib.pdf", fig; pt_per_unit = 1)

        # # Create figure and axis
        # fig = Figure(size = (640, 640))
        # ax_b = Axis(fig[1, 1]; ylabel = L"z", titlealign = :left, title = L"\text{(a) 1st-order buoyancy}~b^1")
        # ax_v = Axis(fig[2, 1]; ylabel = L"z", titlealign = :left, title = L"\text{(b) 1st-order velocity}~v^1")
        # ax_ψ = Axis(fig[3, 1]; xlabel = L"x", ylabel = L"z", titlealign = :left, title = L"\text{(c) 1st-order streamfunction}~\psi^1")
        # hidexdecorations!(ax_b, ticks = false)
        # hidexdecorations!(ax_v, ticks = false)

        # # Create the heatmap
        # hm_b = heatmap!(ax_b, x, z, b1[1,:,:], colorrange = (-5,5), colormap = :diff)
        # hm_v = heatmap!(ax_v, x, z, v1[1,:,:], colorrange = (-2,2), colormap = :delta)
        # hm_ψ = heatmap!(ax_ψ, x, z, ψ1[1,:,:], colorrange = (-0.2, 0.2), colormap = :PuOr)
        # levels=-0.4:0.1:0.4
        # ct_b = contour!(ax_b, x, z, b_e .+ γ * b1[1,:,:]; levels = levels, linewidth = 1, color = :black)
        # ct_v = contour!(ax_v, x, z, b_e .+ γ * b1[1,:,:]; levels = levels, linewidth = 1, color = :black)
        # ct_ψ = contour!(ax_ψ, x, z, b_e .+ γ * b1[1,:,:]; levels = levels, linewidth = 1, color = :black)

        # # Add colorbar
        # Colorbar(fig[1, 2], hm_b)
        # Colorbar(fig[2, 2], hm_v)
        # Colorbar(fig[3, 2], hm_ψ)

        # # Add title
        # title_obs = Observable(@sprintf("t=%.2f", t_all[1]))
        # fig[0, :] = Label(fig, title_obs, fontsize = 20)
        # colgap!(fig.layout, 1, Relative(0))
        # for i = 1:3
        #     rowgap!(fig.layout, i, Relative(0))
        # end
        # resize_to_layout!(fig)

        # # Create animation
        # record(fig, "b1v1psi1_evolution_linear.gif", 1:nt; framerate = 10) do i
        #     t = t_all[i]
            
        #     # Update the heatmap data
        #     hm_b[3] = b1[i,:,:]
        #     hm_v[3] = v1[i,:,:]
        #     hm_ψ[3] = ψ1[i,:,:]
        #     ct_b[3] = b_e .+ γ * b1[i,:,:]
        #     ct_v[3] = b_e .+ γ * b1[i,:,:]
        #     ct_ψ[3] = b_e .+ γ * b1[i,:,:]
            
        #     # Update the title
        #     title_obs[] = @sprintf("t=%.2f", t)
        # end
        # println("Animation saved to b1v1psi1_evolution_linear.gif")

    end
else

    # Figure 3,4,5 – geostrophic adjustment

    # ----------------------------- helper: tiny meshgrid clone -------------------
    """
        meshgrid(x, y) → X, Y   (MATLAB-style)

    Give 2-D matrices whose shape is (length(x), length(y)):
        X = [x  …  x]   (columns   = copies of x)
        Y = [yᵀ …  yᵀ]  (rows      = copies of yᵀ)
    """
    meshgrid(x::AbstractVector, y::AbstractVector) =
        repeat(x, 1, length(y)), repeat(y', length(x), 1)

    # ------------------------------ problem geometry -----------------------------
    xgrid  = x                  # cross-stream coordinate       (Nx =  100)
    zgrid  = z                    # vertical coordinate           (Nz =  101)

    T      = 0:0.01:t_max                 # (MATLAB x)  → physical time
    Nx, Nz = length(xgrid), length(zgrid)
    i1, i2 = 30, 130                          # “front” and “aft” time indices
    Trng   = T[i1:i2]                         # 30 ≤ t ≤ 13 (length 111)

    # reverse x once – MATLAB used  y_dir = fliplr(y)
    revx!(A) = (reverse!(A, dims = 1); A)     # in-place to avoid temp allocations

    # ------------------------- collect the fields we need ------------------------
    bottom_b = Matrix{Float32}(undef, Nx, length(Trng))
    top_b    = similar(bottom_b)

    for (k, t) in enumerate(Trng)
        b, _, _ = fields_at_time(t)
        bottom_b[:, k] .= revx!(b[:,   1])    # z = 0 slice, flipped in x
        top_b[:,    k] .= revx!(b[:, end])    # z = 1 slice
    end

    b_i1, v0_i1, ψ_i1 = fields_at_time(T[i1])
    b_i2, v0_i2, ψ_i2 = fields_at_time(T[i2])
    revx!.((b_i1, v0_i1, ψ_i1, b_i2, v0_i2, ψ_i2))   # flip all once, MATLAB style

    # --------------------------- build the slice grids ---------------------------
    # 1) constant-z slices (bottom & top)
    t_bt, _        = meshgrid(Trng, xgrid)            # Ntʹ × Nx  (→ transpose later)
    x_bt, _        = meshgrid(xgrid, Trng)            # same size
    z0  = zeros(size(x_bt));     z1 = ones(size(x_bt))

    # 2) constant-t slices (front & aft)
    x_xt, z_xt     = meshgrid(xgrid, zgrid)           # Nx × Nz
    t_front        = fill(T[i1], size(x_xt))
    t_aft          = fill(T[i2], size(x_xt))

    # ----------------------------- figure & axis ---------------------------------
    fig = Figure(size = (640, 320))
    ax  = Axis3(fig[2, 1]; 
            xlabel = L"\text{Time}~t ", ylabel = L"\text{Cross-front}~x", zlabel = L"\text{Vertical}~z",
            aspect = (1.1, 1, 0.6),
            limits = ((T[i1], T[i2]),
                        (minimum(xgrid), maximum(xgrid)),
                        (minimum(zgrid), maximum(zgrid))),
            elevation = 0.3, azimuth = 0.24π,
            xspinesvisible = false, yspinesvisible = false, zspinesvisible = false,
            xgridvisible = false, ygridvisible = false, zgridvisible = false,
            perspectiveness = 0.7,protrusions = (40,10,0,0))
    ax.xreversed = true

    # common keyword bundle -------------------------------------------------------
    cm_buoy  = :diff
    cm_vel   = :delta
    cm_psi   = :PuOr

    clims = 1.1 .* extrema(top_b)
    kw_buoy  = (colorrange = clims, colormap = cm_buoy,
                rasterize = true, shading=NoShading)
    kw_vel   = (colormap = cm_vel, rasterize = true, shading=NoShading)
    kw_psi   = (colormap   = cm_psi, rasterize = true, shading=NoShading)

    # ----------------------------- draw the surfaces -----------------------------                      # velocity
    sf = surface!(ax, t_aft,    x_xt,  z_xt;     color = b_i2,      kw_buoy...)   # t = t₂ 
    surface!(ax, t_bt, x_bt', z0';   color = bottom_b', kw_buoy...)   # bottom 
    surface!(ax, t_bt, x_bt', z1';         color = top_b',    kw_buoy...)   # top
    surface!(ax, t_front,  x_xt,  z_xt;     color = b_i1,      kw_buoy...)   # t = t₁ 
    Colorbar(fig[2, 2], sf, height = Relative(0.3), tellheight=false)
    title = L"\text{Buoyancy evolution}~b^0(x,z,t)"
    fig[1, 1:2] = Label(fig, title; tellwidth = false, padding = (0, 0, -120, 0))
    rowgap!(fig.layout, 1, Relative(-0.3))
    colgap!(fig.layout, 1, Relative(0))    
    resize_to_layout!(fig)
    save("figure3_B00.pdf", fig; pt_per_unit = 1)
    println("✓  saved → figure3_B00.pdf")

    fig4 = Figure(size = (640, 300));
    gab = fig4[1, 1] = GridLayout()
    ax_a  = Axis(gab[1,1], xlabel = L"x", ylabel = L"z", titlealign = :left, 
                title = L"\text{(a) Early frontogenesis}~v^0")
    hm_a = heatmap!(ax_a, x, z, v0_i1;
            colormap = cm_vel, rasterize = true,
            colorrange = (-maximum(abs, v0_i1), maximum(abs, v0_i1)))
    Colorbar(gab[1,2], hm_a)
    contour!(ax_a, x, z, b_i1; levels = 10, linewidth = 1, color = :black)
    ax_b  = Axis(gab[1,3], xlabel = L"x", titlealign = :left, 
                title = L"\text{(b) Late frontogenesis}~v^0")
    hm_b = heatmap!(ax_b, x, z, v0_i2;
            colormap = cm_vel, rasterize = true, 
            colorrange = (-maximum(abs, v0_i2), maximum(abs, v0_i2)))
    Colorbar(gab[1,4], hm_b)
    contour!(ax_b, x, z, b_i2; levels = 10, linewidth = 1, color = :black)
    hideydecorations!(ax_b, ticks=false)
    colgap!(gab, 1, 3)
    colgap!(gab, 2, 10)
    colgap!(gab, 3, 3)
    resize_to_layout!(fig4)
    save("figure4_v1v2.pdf", fig4; pt_per_unit = 1)

    fig5 = Figure(size = (640, 300));
    gab = fig5[1, 1] = GridLayout()
    ax_a  = Axis(gab[1,1], xlabel = L"x", ylabel = L"z", titlealign = :left, 
                title = L"\text{(a) Early frontogenesis}~\psi^0")
    hm_a = heatmap!(ax_a, x, z, ψ_i1;
            colormap = cm_psi, rasterize = true,
            colorrange = (-maximum(abs, ψ_i1), maximum(abs, ψ_i1)))
    Colorbar(gab[1,2], hm_a)
    contour!(ax_a, x, z, b_i1; levels = 10, linewidth = 1, color = :black)
    ax_b  = Axis(gab[1,3], xlabel = L"x", titlealign = :left, 
                title = L"\text{(b) Late frontogenesis}~\psi^0")
    hm_b = heatmap!(ax_b, x, z, ψ_i2;
            colormap = cm_psi, rasterize = true, 
            colorrange = (-maximum(abs, ψ_i2), maximum(abs, ψ_i2)))
    Colorbar(gab[1,4], hm_b)
    contour!(ax_b, x, z, b_i2; levels = 10, linewidth = 1, color = :black)
    hideydecorations!(ax_b, ticks=false)
    colgap!(gab, 1, 3)
    colgap!(gab, 2, 10)
    colgap!(gab, 3, 3)
    resize_to_layout!(fig5)
    save("figure5_psi1psi2.pdf", fig5; pt_per_unit = 1)
end