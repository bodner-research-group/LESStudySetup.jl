################################################################################
#  FIGURE 1  –  Similarity solution of surface-gradient frontogenesis
#  Julia 1.10+, CairoMakie ≥ 0.11, SpecialFunctions ≥ 2.3
################################################################################
using CairoMakie             # high-quality vector / raster plotting
using SpecialFunctions: erf  # Gaussian error function
using Printf

# ------------------------------------------------------------------
# 1.  Non-dimensional parameters   (Appendix A.1 and γ=0.1)
# ------------------------------------------------------------------
const Ro   = 1.0
const Bu   = 0.1
const Fr   = Ro / Bu
γ          = 0.05           # growth-rate parameter requested
ϵ          = 0.00           # small parameter for geostrophic departure
t_max      = 20.0           # blow-up time
t_early    = 0.4            # diagnostic snapshot
x_min,x_max = -10, 10
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
# 4.  Zeroth-order fields at any given time T
# ------------------------------------------------------------------
function fields_at_time(T; γ = 0, ϵ = 0)
    eT= exp(γ*T)

    b   = Matrix{Float64}(undef, nz, nx)
    v0  = similar(b)
    w0  = similar(b)

    A = (eT-ϵ*cos(sqrt(1-γ^2)*T)+γ*(ϵ-2)*sin(sqrt(1-γ^2)*T)/sqrt(1-γ^2))
    B = (γ * (eT-cos(sqrt(1-γ^2)*T))+(ϵ-2γ^2)*sin(sqrt(1-γ^2)*T)/2/sqrt(1-γ^2))
    for (k, zz) in enumerate(z)
        Xvals = solve_X_vec(x, zz, T; γ, ϵ)

        Bp  = dB0(Xvals)
        Bpp = d2B0(Xvals)

        dv0dX = Bpp .* (zz - 0.5) * A
        v0[k,:] =  Bp  .* (zz - 0.5) * A                            # (b)
        w0[k,:] =  Ro * B .* Bpp .* eT .* zz*(zz-1) ./ (1 .- Ro^2 * eT .* dv0dX) # (c)

        intZ   = 0.5*zz^2 - 0.5*zz
        Δb     = -(Ro^2)/(Fr^2) * A .* Bpp .* intZ
        b[k,:] =  B0(Xvals) #.+ 1/Fr^2*zz .+ Δb                               # (a)
    end

    ψ = -cumsum(w0; dims = 2) * dx      # integrate w = −∂ψ/∂x, ψ=0 at x_min
    return b', v0', ψ'
end

# -------------------- early-time fields for panels (b) & (c) -------------
b_e, v_e, ψ_e = fields_at_time(t_early)

# ------------------------------------------------------------------
# 5.  Plotting with CairoMakie
# ------------------------------------------------------------------
set_theme!(theme_latexfonts(), fontsize=12,figure_padding = 10)            # optional stylistic baseline

# Figure 2 – b-v geostrophic balance
if ϵ<1e-6
    if γ < 1e-6
        fig1 = Figure(size = (640, 300));
        gab = fig1[1, 1] = GridLayout()
        ax_a  = Axis(gab[1,1], xlabel = L"x", ylabel = L"z", titlealign = :left, 
                    title = L"\text{(a) Buoyancy}~b^0")
        clims = 1.1 .* extrema(b_e)
        hm_a = heatmap!(ax_a, x, z, b_e;
                colormap = :diff, rasterize = true, colorrange = clims)
        Colorbar(gab[1,2], hm_a)
        ax_b  = Axis(gab[1,3], xlabel = L"x", titlealign = :left, 
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
        save("figure2_bgvg.pdf", fig1; pt_per_unit = 1)
    else
        dt = 0.2
        t_all = 0:dt:t_max
        nt = length(t_all)
        ψ1, v1, b1 = zeros(nt,nx,nz),zeros(nt,nx,nz),zeros(nt,nx,nz)
        bfull = zeros(nt,nx,nz)
        for i = 1:nt
            bfulli, vfulli, ψfulli = fields_at_time(t_all[i]; γ, ϵ)
            bfull[i,:,:] = bfulli
            ψ1[i,:,:], v1[i,:,:], b1[i,:,:] = ψfulli .- ψ_e, vfulli .- v_e, bfulli .- b_e
        end
        ψ1, v1, b1 = ψ1./γ, v1./γ, b1./γ

        # Create figure and axis
        fig = Figure(size = (640, 640))
        ax_b = Axis(fig[1, 1]; ylabel = L"z", titlealign = :left, title = L"\text{(a) 1st-order buoyancy}~b^1")
        ax_v = Axis(fig[2, 1]; ylabel = L"z", titlealign = :left, title = L"\text{(b) 1st-order velocity}~v^1")
        ax_ψ = Axis(fig[3, 1]; xlabel = L"x", ylabel = L"z", titlealign = :left, title = L"\text{(c) 1st-order streamfunction}~\psi^1")
        hidexdecorations!(ax_b, ticks = false)
        hidexdecorations!(ax_v, ticks = false)

        # Create the heatmap
        hm_b = heatmap!(ax_b, x, z, b1[1,:,:], colorrange = (-5,5), colormap = :diff)
        hm_v = heatmap!(ax_v, x, z, v1[1,:,:], colorrange = (-2,2), colormap = :delta)
        hm_ψ = heatmap!(ax_ψ, x, z, ψ1[1,:,:], colorrange = (-0.2, 0.2), colormap = :PuOr)
        ct_b = contour!(ax_b, x, z, bfull[1,:,:]; levels = 10, linewidth = 1, color = :black)
        ct_v = contour!(ax_v, x, z, bfull[1,:,:]; levels = 10, linewidth = 1, color = :black)
        ct_ψ = contour!(ax_ψ, x, z, bfull[1,:,:]; levels = 10, linewidth = 1, color = :black)

        # Add colorbar
        Colorbar(fig[1, 2], hm_b)
        Colorbar(fig[2, 2], hm_v)
        Colorbar(fig[3, 2], hm_ψ)

        # Add title
        title_obs = Observable(@sprintf("t=%.2f", t_all[1]))
        fig[0, :] = Label(fig, title_obs, fontsize = 20)
        colgap!(fig.layout, 1, Relative(0))
        for i = 1:3
            rowgap!(fig.layout, i, Relative(0))
        end
        resize_to_layout!(fig)

        # Create animation
        record(fig, "b1v1psi1_evolutionv2.gif", 1:nt; framerate = 10) do i
            t = t_all[i]
            
            # Update the heatmap data
            hm_b[3] = b1[i,:,:]
            hm_v[3] = v1[i,:,:]
            hm_ψ[3] = ψ1[i,:,:]
            ct_b[3] = bfull[i,:,:]
            ct_v[3] = bfull[i,:,:]
            ct_ψ[3] = bfull[i,:,:]
            
            # Update the title
            title_obs[] = @sprintf("t=%.2f", t)
        end
        println("Animation saved to b1v1psi1_evolutionv2.gif")

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