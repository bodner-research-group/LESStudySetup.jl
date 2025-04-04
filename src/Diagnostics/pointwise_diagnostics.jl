using Oceananigans.Operators: div_xyᶜᶜᶜ
using Oceananigans.Operators
using Oceananigans.Utils: launch!
using KernelAbstractions: @kernel, @index
using Statistics: mean, var

""" propagate a diagnostic over the timeseries found in snapshots 
    and save the results in a FieldTimeSeries object that is backed up in 
    `filename`.
"""
function propagate_function(func, snapshots; filename = "temp.jld2")
    first_operation = func(snapshots, 1)
    grid = first_operation.grid
    Nx, Ny, Nz = size(grid)
    field = Field(first_operation; indices = (2:Nx-1, 2:Ny-1, 2:Nz-1))
    compute!(field)

    loc  = location(field)
    grid = field.grid

    saved_times = times(snapshots)
    func_name   = String(Symbol(func))
    Nt          = length(saved_times)

    field_time_series = FieldTimeSeries{loc...}(grid, saved_times; 
                                                backend = OnDisk(),
                                                path = filename,
                                                name = func_name)

    @info "calculating $func_name over the timeseries..."

    for i in 1:Nt
        set!(field, func(snapshots, i))
        set!(field_time_series, field, i)
        @info "calculated $func_name at iteration $i of $Nt"
    end

    return field_time_series
end

""" x-z momentum flux """
function uw(snapshots, i)
    u = snapshots[:u][i]
    w = snapshots[:w][i]

    return u * w
end

""" y-z momentum flux """
function vw(snapshots, i)
    v = snapshots[:u][i]
    w = snapshots[:w][i]

    return v * w
end

""" filtered x-z momentum flux """
function uw_filtered(snapshots, i; cutoff = 20kilometer)
    u = snapshots[:u][i]
    w = snapshots[:w][i]
    ul, uh = symmetric_filtering(u; cutoff)
    wl, wh = symmetric_filtering(w; cutoff)

    uw = (u - mean(u,dims=(1,2))) * (w - mean(w,dims=(1,2)))
    ulwl = (ul - mean(ul,dims=(1,2))) * (wl - mean(wl,dims=(1,2)))
    uhwh = (uh - mean(uh,dims=(1,2))) * (wh - mean(wh,dims=(1,2)))

    return uw, ulwl, uhwh
end

""" filtered y-z momentum flux """
function vw_filtered(snapshots, i; cutoff = 20kilometer)
    v = snapshots[:v][i]
    w = snapshots[:w][i]
    vl, vh = symmetric_filtering(v; cutoff)
    wl, wh = symmetric_filtering(w; cutoff)

    vw = (v - mean(v,dims=(1,2))) * (w - mean(w,dims=(1,2)))
    vlwl = (vl - mean(vl,dims=(1,2))) * (wl - mean(wl,dims=(1,2)))
    vhwh = (vh - mean(vh,dims=(1,2))) * (wh - mean(wh,dims=(1,2)))

    return vw, vlwl, vhwh
end

""" filtered x-y momentum flux """
function uv_filtered(snapshots, i; cutoff = 20kilometer)
    u = snapshots[:u][i]
    v = snapshots[:v][i]
    ul, uh = symmetric_filtering(u; cutoff)
    vl, vh = symmetric_filtering(v; cutoff)

    uv = (u - mean(u,dims=(1,2))) * (v - mean(v,dims=(1,2)))
    ulvl = (ul - mean(ul,dims=(1,2))) * (vl - mean(vl,dims=(1,2)))
    uhvh = (uh - mean(uh,dims=(1,2))) * (vh - mean(vh,dims=(1,2)))

    return uv, ulvl, uhvh
end

""" filtered x-x momentum flux """
function u²_filtered(snapshots, i; cutoff = 20kilometer)
    u = snapshots[:u][i]
    ul, uh = symmetric_filtering(u; cutoff)

    u² = (u - mean(u,dims=(1,2)))^2
    ul² = (ul - mean(ul,dims=(1,2)))^2
    uh² = (uh - mean(uh,dims=(1,2)))^2

    return u², ul², uh²
end

""" filtered y-y momentum flux """
function v²_filtered(snapshots, i; cutoff = 20kilometer)
    v = snapshots[:v][i]
    vl, vh = symmetric_filtering(v; cutoff)

    v² = (v - mean(v,dims=(1,2)))^2
    vl² = (vl - mean(vl,dims=(1,2)))^2
    vh² = (vh - mean(vh,dims=(1,2)))^2

    return v², vl², vh²
end

""" filtered z-z momentum flux """
function w²_filtered(snapshots, i; cutoff = 20kilometer)
    w = snapshots[:w][i]
    wl, wh = symmetric_filtering(w; cutoff)

    w²  = (w - mean(w,dims=(1,2)))^2
    wl² = (wl - mean(wl,dims=(1,2)))^2
    wh² = (wh - mean(wh,dims=(1,2)))^2

    return w², wl², wh²
end

""" zonal buoyancy flux """
function ub(snapshots, i)
    α = parameters.α
    g = parameters.g

    u = snapshots[:u][i]
    T = snapshots[:T][i]
    
    return α * g * T * u
end

""" meridional buoyancy flux """
function vb(snapshots, i)
    α = parameters.α
    g = parameters.g

    v = snapshots[:v][i]
    T = snapshots[:T][i]
    
    return α * g * T * v
end

""" vertical buoyancy flux """
function wb(snapshots, i=nothing)
    α = parameters.α
    g = parameters.g

    w = snapshots[:w]
    T = snapshots[:T]

    if !isnothing(i)
      w = w[i]
      T = T[i]
    end

    return α * g * T * w
end

""" filtered zonal buoyancy flux """
function ub_filtered(snapshots, i=nothing; cutoff = 20kilometer)
    α = parameters.α
    g = parameters.g

    u = snapshots[:u]
    b = α * g * snapshots[:T]

    if !isnothing(i)
      u = u[i]
      b = b[i]
    end

    ul, uh = symmetric_filtering(u; cutoff)
    bl, bh = symmetric_filtering(b; cutoff)

    ub = (u - mean(u,dims=(1,2))) * (b - mean(b,dims=(1,2)))
    ulbl = (ul - mean(ul,dims=(1,2))) * (bl - mean(bl,dims=(1,2)))
    uhwh = (uh - mean(uh,dims=(1,2))) * (bh - mean(bh,dims=(1,2)))

    return ub, ulbl, uhwh
end

""" filtered meridional buoyancy flux """
function vb_filtered(snapshots, i=nothing; cutoff = 20kilometer)
    α = parameters.α
    g = parameters.g

    v = snapshots[:v]
    b = α * g * snapshots[:T]

    if !isnothing(i)
      v = v[i]
      b = b[i]
    end

    vl, vh = symmetric_filtering(v; cutoff)
    bl, bh = symmetric_filtering(b; cutoff)

    vb = (v - mean(v,dims=(1,2))) * (b - mean(b,dims=(1,2)))
    vlbl = (vl - mean(vl,dims=(1,2))) * (bl - mean(bl,dims=(1,2)))
    vhwh = (vh - mean(vh,dims=(1,2))) * (bh - mean(bh,dims=(1,2)))

    return vb, vlbl, vhwh
end

""" filtered vertical buoyancy flux """
function wb_filtered(snapshots, i=nothing; cutoff = 20kilometer)
    α = parameters.α
    g = parameters.g

    w = snapshots[:w]
    b = α * g * snapshots[:T]
    
    if !isnothing(i)
      w = w[i]
      b = b[i]
    end

    wl, wh = symmetric_filtering(w; cutoff)
    bl, bh = symmetric_filtering(b; cutoff)

    wb = (w - mean(w,dims=(1,2))) * (b - mean(b,dims=(1,2)))
    wlbl = (wl - mean(wl,dims=(1,2))) * (bl - mean(bl,dims=(1,2)))
    whwh = (wh - mean(wh,dims=(1,2))) * (bh - mean(bh,dims=(1,2)))

    return wb, wlbl, whwh
end

""" horizontal kinetic energy """
function KE(snapshots, i=nothing)
    u = snapshots[:u]
    v = snapshots[:v]
    
    if !isnothing(i)
      u = u[i]
      v = v[i]
    end

    return 0.5 * (u^2 + v^2)
end

""" mixed layer depth """
function MLD(snapshots, i=nothing; threshold = 0.03, surface = false)
    α  = parameters.α
    ρ₀ = parameters.ρ₀
    T  = snapshots[:T][i]
    T  = isnothing(i) ? T : T[i]
    grid = T.grid
    h    = MixedLayerDepth(grid, (; T); ΔT = abs(threshold / ρ₀ / α), surface)
    return h
end

""" boundary layer depth """
function BLD1D(snapshots, i; threshold = 0.3)
    α  = parameters.α
    g  = parameters.g
    T  = snapshots[:T][i]
    B  = mean(compute!(Field(α * g * T)), dims = (1, 2))
    uc = snapshots[:u][i]
    vc = snapshots[:v][i]
    wc = compute!(Field(@at (Center, Center, Center) snapshots[:w][i]))
    Uₜ² = (var(uc,dims=(1, 2))+var(vc,dims=(1, 2))+var(wc,dims=(1, 2)))/3
    uc = mean(uc, dims = (1, 2))
    vc = mean(vc, dims = (1, 2))
    wc = mean(wc, dims = (1, 2))
    
    grid = B.grid
    h    = BoundaryLayerDepth(grid, (; B, uc, vc, wc, Uₜ²); Ric = threshold)
    return h
end

""" stratification """
function N²(snapshots, i)
    α = parameters.α
    g = parameters.g
    
    T = snapshots[:T][i]
    B = α * g * T
    
    return ∂z(B)
end

""" horizontalstratification """
function M²(snapshots, i)
    α = parameters.α
    g = parameters.g
    
    T = snapshots[:T][i]
    B = α * g * T
    
    return ∂x(B), ∂y(B)
end

""" vertical vorticity """
function ζ(snapshots, i=nothing; u0=0, v0=0)
    ui = isnothing(i) ? snapshots[:u] : snapshots[:u][i]
    vi = isnothing(i) ? snapshots[:v] : snapshots[:v][i]

    u = ui + u0
    v = vi + v0

    return ∂x(v) - ∂y(u)
end

""" horizontal divergence """
function δ(snapshots, i; U=0, V=0)
    u0 = compute!(Field(snapshots[:u][i] + U))
    v0 = compute!(Field(snapshots[:v][i] + V))
    u = XFaceField(u0.grid,Float32)
    v = YFaceField(v0.grid,Float32)
    set!(u, interior(u0))
    set!(v, interior(v0))
    fill_halo_regions!(u)
    fill_halo_regions!(v)
    
    grid = u.grid
    return KernelFunctionOperation{Center, Center, Center}(div_xyᶜᶜᶜ, grid, u, v)
end

""" potential vorticity """
function PV(snapshots, i)
    u = snapshots[:u][i]
    v = snapshots[:v][i]
    w = snapshots[:w][i]
    T = snapshots[:T][i]

    f = parameters.f
    α = parameters.α
    g = parameters.g

    ωz = ζ(snapshots, i) + f
    ωx = ∂z(v) - ∂y(w)
    ωy = ∂x(w) - ∂z(u)

    bx = α * g * ∂x(T)
    by = α * g * ∂y(T)
    bz = α * g * ∂z(T)

    return ωx * bx + ωy * by + ωz * bz
end

""" horizontal frontogenesis """
function Bₕ(snapshots, i)
    u = snapshots[:u][i]
    v = snapshots[:v][i]
    T = snapshots[:T][i]

    α = parameters.α
    g = parameters.g

    b = compute!(Field(α * g * T))

    return -(∂x(b)^2 * ∂x(u) + ∂y(b)^2 * ∂y(v))-∂x(b)*∂y(b)*(∂x(v) + ∂y(u))
end

""" vertical frontogenesis """
function Bᵥ(snapshots, i)
    w = snapshots[:w][i]
    T = snapshots[:T][i]

    α = parameters.α
    g = parameters.g

    b = compute!(Field(α * g * T))

    return -(∂z(b) * (∂x(w) * ∂x(b) + ∂y(w) * ∂y(b)))
end

""" spectral horizontal advection """
function Ah(snapshots, i; U=0, V=0)
    u = snapshots[:u][i] 
    v = snapshots[:v][i]
    #KE = 0.5 * (u^2 + v^2)
    #ζ = ∂x(v) - ∂y(u)
    xu, yu, zu = nodes(u)
    xv, yv, _ = nodes(v)
    Nz = length(zu)

    #Au = compute!(Field(- ∂x(KE) + v * ζ))
    #Av = compute!(Field(- ∂y(KE) - u * ζ))
    Au = compute!(Field(- (u + U) * ∂x(u) - (v + V) * ∂y(u)))
    Av = compute!(Field(- (u + U) * ∂x(v) - (v + V) * ∂y(v)))

    # Fourier transform
    Au1 = isotropic_powerspectrum(interior(u, :, :, 1), interior(Au, :, :, 1), xu, yu)
    Av1 = isotropic_powerspectrum(interior(v, :, :, 1), interior(Av, :, :, 1), xv, yv)
    A = zeros(Nz,length(Au1.spec))
    A[1,:] = real.(Au1.spec + Av1.spec)
    for k = 2:Nz
        Auk = isotropic_powerspectrum(interior(Au, :, :, k), interior(u, :, :, k), xu, yu)
        Avk = isotropic_powerspectrum(interior(Av, :, :, k), interior(v, :, :, k), xv, yv)
        A[k,:] = real.(Auk.spec + Avk.spec)
    end

    return A, zu, Au1.freq
end

""" spectral vertical advection """
function Av(snapshots, i)
    u = snapshots[:u][i]
    v = snapshots[:v][i]
    w = compute!(Field(@at (Center, Center, Center) snapshots[:w][i]))

    xu, yu, zu = nodes(u)
    xv, yv, _ = nodes(v)
    Nz = length(zu)

    Au = compute!(Field(- w * ∂z(u)))
    Av = compute!(Field(- w * ∂z(v)))

    Au1 = isotropic_powerspectrum(interior(Au, :, :, 1), interior(u, :, :, 1), xu, yu)
    A = zeros(Nz,length(Au1.spec))
    for k = 2:Nz
        Auk = isotropic_powerspectrum(interior(Au, :, :, k), interior(u, :, :, k), xu, yu)
        Avk = isotropic_powerspectrum(interior(Av, :, :, k), interior(v, :, :, k), xv, yv)
        A[k,:] = real.(Auk.spec .+ Avk.spec)
    end

    return A, zu, Au1.freq

end

""" spectral convertion of potential to KE """
function Ac(snapshots, i)
    w = compute!(Field(@at (Center, Center, Center) snapshots[:w][i]))
    T = snapshots[:T][i]

    xT, yT, zT = nodes(T)
    Nz = length(zT)

    α = parameters.α
    g = parameters.g
    B = compute!(Field(α * g * T))
    C1 = isotropic_powerspectrum(interior(B, :, :, 1), interior(w, :, :, 1), xT, yT)
    C = zeros(Nz,length(C1.spec))
    C[1,:] = real.(C1.spec)
    for k = 2:Nz
        Ck = isotropic_powerspectrum(interior(B, :, :, k), interior(w, :, :, k), xT, yT)
        C[k,:] = real.(Ck.spec)
    end

    return C, zT, C1.freq
end

""" spectral 3D pressure work """
function Ap(snapshots, i)
    u = snapshots[:u][i]
    v = snapshots[:v][i]
    w = compute!(Field(@at (Center, Center, Center) snapshots[:w][i]))
    p = snapshots[:p][i]
    T = snapshots[:T][i]

    α = parameters.α
    g = parameters.g
    ρ₀ = parameters.ρ₀
    px = compute!(Field(∂x(p)))
    py = compute!(Field(∂y(p)))
    pz = compute!(Field(α * g * T))

    xu, yu, zu = nodes(u)
    xv, yv, _ = nodes(v)
    xw, yw, _ = nodes(w)

    Au1 = isotropic_powerspectrum(interior(px, :, :, 1), interior(u, :, :, 1), xu, yu)
    Av1 = isotropic_powerspectrum(interior(py, :, :, 1), interior(v, :, :, 1), xv, yv)
    Aw1 = isotropic_powerspectrum(interior(pz, :, :, 1), interior(w, :, :, 1), xw, yw)
    A = zeros(length(zu),length(Au1.spec))
    A[1,:] = - real.(Au1.spec + Av1.spec + Aw1.spec)
    for k = 2:length(zu)
        Auk = isotropic_powerspectrum(interior(px, :, :, k), interior(u, :, :, k), xu, yu)
        Avk = isotropic_powerspectrum(interior(py, :, :, k), interior(v, :, :, k), xv, yv)
        Awk = isotropic_powerspectrum(interior(pz, :, :, k), interior(w, :, :, k), xw, yw)
        A[k,:] = - real.(Auk.spec + Avk.spec + Awk.spec)
    end

    return A, zu, Au1.freq
end

"""
    subfilter_stress!(τ, u, v, u̅, v̅; kernel=:tophat, cutoff=4kilometer)

Compute the subfilter stress (residual stress) field τ from fields `u` and `v` 
given their coarse-grained versions `u̅` and `v̅`. The subfilter stress is defined as 
τ = coarse_graining(u*v) - u̅ * v̅.
"""
function subfilter_stress!(τ, u, v, u̅, v̅; kernel=:tophat, cutoff=4kilometer)
    # Compute the product field (u*v) on the fly.
    product_field = compute!(Field(u * v))
    # Coarse-grain the product field and store the result in τ.
    coarse_graining!(product_field, τ; kernel=kernel, cutoff=cutoff)

    # Compute the product of the filtered fields and subtract.
    # We assume that the multiplication is elementwise.
    stress = compute!(Field(τ - (u̅ * v̅)))
    # Replace the interior of τ with the computed stress.
    set!(τ, interior(stress))
    fill_halo_regions!(τ)
    return nothing
end

"""
    coarse_grained_fluxes(snapshots, i; kernel=:tophat, cutoff=4kilometer)

Compute the coarse-grained velocities and cross-scale fluxes (i.e. residual stresses and
transfer terms) from a snapshot indexed by `i` in `snapshots`. The function uses a
coarse-graining (filtering) procedure based on a specified kernel and cutoff. It assumes
periodic boundary conditions (so that FFTs can be used) and the use of Oceananigans Field
types.
"""
function coarse_grained_fluxes(snapshots, i; kernel=:tophat, cutoff=4kilometer)
    t0 = time()

    # Extract snapshot fields.
    u0 = snapshots[:u][i]
    v0 = snapshots[:v][i]
    w0 = snapshots[:w][i]
    T0 = snapshots[:T][i]
    u = XFaceField(u0.grid,Float32)
    v = YFaceField(v0.grid,Float32)
    w = ZFaceField(w0.grid,Float32)
    T = CenterField(T0.grid,Float32)
    set!(u, interior(u0))
    set!(v, interior(v0))
    set!(w, interior(w0))
    set!(T, interior(T0))
    fill_halo_regions!(u)
    fill_halo_regions!(v)
    fill_halo_regions!(w)
    fill_halo_regions!(T)
    _, _, zT = nodes(T)
    kT = length(zT)
    k0 = findlast(zT .< -60)

    # Retrieve physical parameters.
    α = parameters.α
    g = parameters.g
    f = parameters.f

    # Compute the buoyancy field B = α*g*T.
    B = compute!(Field(α * g * T))

    # Allocate filtered velocity and buoyancy fields.
    u̅ = XFaceField(u.grid,Float32)
    v̅ = YFaceField(v.grid,Float32)
    w̅ = ZFaceField(w.grid,Float32)
    B̅ = CenterField(B.grid,Float32)

    # --- Coarse-grain (filter) the primary fields. ---
    for (orig, filt) in ((u, u̅), (v, v̅), (w, w̅), (B, B̅))
        coarse_graining!(orig, filt; kernel=kernel, cutoff=cutoff)
    end
    println("Computed filtered fields at $(time() - t0)s...")

    # --- Compute subfilter stresses ---
    # First set: fluxes associated with the u-momentum.
    τuu = XFaceField(u.grid,Float32)
    τuv = YFaceField(v.grid,Float32)
    τuw = ZFaceField(w.grid,Float32)
    subfilter_stress!(τuu, u, u, u̅, u̅; kernel=kernel, cutoff=cutoff)
    subfilter_stress!(τuv, v, u, v̅, u̅; kernel=kernel, cutoff=cutoff)
    subfilter_stress!(τuw, w, u, w̅, u̅; kernel=kernel, cutoff=cutoff)

    # Second set: fluxes associated with the v-momentum.
    τvv = YFaceField(v.grid,Float32)
    τvu = XFaceField(u.grid,Float32)
    τvw = ZFaceField(w.grid,Float32)
    subfilter_stress!(τvv, v, v, v̅, v̅; kernel=kernel, cutoff=cutoff)
    subfilter_stress!(τvu, u, v, u̅, v̅; kernel=kernel, cutoff=cutoff)
    subfilter_stress!(τvw, w, v, w̅, v̅; kernel=kernel, cutoff=cutoff)
    println("Mean surface τvv is $(mean(interior(τvv, :, :, kT)))")

    # Third set: fluxes associated with the w-momentum.
    τww = ZFaceField(w.grid,Float32)
    subfilter_stress!(τww, w, w, w̅, w̅; kernel=kernel, cutoff=cutoff)
    println("Computed residual stress terms at $(time() - t0)s...")

    # --- Compute transfer (flux) terms using derivatives ---
    # Note: The derivative operators (∂x, ∂y, ∂z) are assumed to be available.
    # Πₕ: Horizontal transfer term.
    Πₕ = compute!(Field(-(τuu * ∂x(u̅) +
                           τvv * ∂y(v̅) +
                           τuv * ∂y(u̅) +
                           τvu * ∂x(v̅))))
    println("Transfer term Πₕ done at $(time() - t0)s, ML-average $(mean(interior(Πₕ, :, :,k0:kT)))")

    # Πδ: Diagonal (divergence-related) transfer term.
    Πδ = compute!(Field(-(τuu + τvv) * (∂x(u̅) + ∂y(v̅)) / 2))
    println("Transfer term Πδ done at $(time() - t0)s")

    # Πᵥ: Vertical transfer term.
    Πᵥ = compute!(Field(@at (Center, Center, Center) -(τuw * ∂z(u̅) +
                                                         τvw * ∂z(v̅))))
    println("Transfer term Πᵥ done at $(time() - t0)s")

    # Πvgl: Flux transfer term associated with buoyancy gradients, scaled by 1/f.
    Πvgl = compute!(Field(@at (Center, Center, Center) -(τuw * ∂z(-∂y(B̅)) +
                                                          τvw * ∂z(∂x(B̅))) / f))
    println("Transfer term Πvgl done at $(time() - t0)s")

    return u̅, v̅, w̅, τuu, τvv, τww, Πₕ, Πδ, Πᵥ, Πvgl
end


""" filtered energy budgets """
function filtered_budgets(snapshots, i; budget=true, U=0, V=0, kernel = :tophat, cutoff = 25kilometer)
    t0 = time()
    u = snapshots[:u][i]
    v = snapshots[:v][i]
    w = snapshots[:w][i]
    T = snapshots[:T][i]
    p = snapshots[:pHY′][i]
    κu = snapshots[:κu][i]
    κc = snapshots[:κc][i]
    grid = u.grid
    _, _, zT = nodes(T)
    α = parameters.α
    g = parameters.g
    Au = compute!(Field(U * ∂x(u) + V * ∂y(u)))
    Av = compute!(Field(V * ∂y(v) + U * ∂x(v)))
    b = compute!(Field(α * g * T))
    Ab = compute!(Field(U * ∂x(b) + V * ∂y(b)))
    Fx, Fy = compute!(Field(∂z(κu * ∂z(u)))), compute!(Field(∂z(κu * ∂z(v))))
    Q = compute!(Field(∂z(κc * ∂z(b))))

    u̅ = XFaceField(u.grid)
    v̅ = YFaceField(v.grid)
    w̅ = ZFaceField(w.grid)
    b̅ = CenterField(b.grid)
    E̅ₚ = CenterField(b.grid)
    # U̅ = XFaceField(u.grid)
    # V̅ = YFaceField(v.grid)
    p̅ = CenterField(p.grid)
    F̅x, F̅y = XFaceField(Fx.grid), YFaceField(Fy.grid)
    A̅u, A̅v = XFaceField(Au.grid), YFaceField(Av.grid)
    A̅b = XFaceField(Ab.grid)
    Q̅ = CenterField(Q.grid)

    println("computing filtered fields at $(time()-t0)s...")
    coarse_graining!(u, u̅; kernel, cutoff)
    coarse_graining!(v, v̅; kernel, cutoff)
    coarse_graining!(w, w̅; kernel, cutoff)
    coarse_graining!(b, b̅; kernel, cutoff)
    # coarse_graining!(U, U̅; kernel, cutoff)
    # coarse_graining!(V, V̅; kernel, cutoff)
    coarse_graining!(p, p̅; kernel, cutoff)
    coarse_graining!(Fx, F̅x; kernel, cutoff)
    coarse_graining!(Fy, F̅y; kernel, cutoff)
    coarse_graining!(Au, A̅u; kernel, cutoff)
    coarse_graining!(Av, A̅v; kernel, cutoff)
    coarse_graining!(Ab, A̅b; kernel, cutoff)
    coarse_graining!(Q, Q̅; kernel, cutoff)
    set!(E̅ₚ, interior(b̅,:,:,:).*reshape(-zT,(1,1,:)))

    println("computing finescale fields at $(time()-t0)s...")
    uᵖ = compute!(Field(u - u̅))
    vᵖ = compute!(Field(v - v̅))
    # Uᵖ = compute!(Field(U - U̅))
    # Vᵖ = compute!(Field(V - V̅))
    wᵖ = compute!(Field(w - w̅))
    bᵖ = compute!(Field(b - b̅))
    pᵖ = compute!(Field(p - p̅))
    Fxᵖ = compute!(Field(Fx - F̅x))
    Fyᵖ = compute!(Field(Fy - F̅y))
    Auᵖ = compute!(Field(Au - A̅u))
    Avᵖ = compute!(Field(Av - A̅v))
    Abᵖ = compute!(Field(Ab - A̅b))
    Qᵖ = compute!(Field(Q - Q̅))
    Eₚᵖ = CenterField(b.grid)
    set!(Eₚᵖ, interior(bᵖ,:,:,:).*reshape(-zT,(1,1,:)))

    E̅ₖ = compute!(Field((u̅^2 + v̅^2)/2))
    Eₖᵖ = compute!(Field((uᵖ^2 + vᵖ^2)/2))
    Σ̅b = compute!(Field(b̅^2 / 2))
    Σbᵖ = compute!(Field(bᵖ^2 / 2))

    havg(x) = vec(mean(x, dims = (1,2)))
    if !budget
        return havg(E̅ₚ), havg(E̅ₖ), havg(Eₚᵖ), havg(Eₖᵖ), havg(Σ̅b), havg(Σbᵖ)
    else

        println("computing residual stress terms at $(time()-t0)s...")
        τuu = XFaceField(u.grid)
        τuv = YFaceField(v.grid)
        τuw = ZFaceField(w.grid)
        subfilter_stress!(τuu,u,u,u̅,u̅; kernel, cutoff)
        subfilter_stress!(τuv,v,u,v̅,u̅; kernel, cutoff)
        subfilter_stress!(τuw,w,u,w̅,u̅; kernel, cutoff)
        
        τvv = YFaceField(v.grid)
        τvu = XFaceField(u.grid)
        τvw = ZFaceField(w.grid)
        subfilter_stress!(τvv,v,v,v̅,v̅; kernel, cutoff)
        subfilter_stress!(τvu,u,v,u̅,v̅; kernel, cutoff)
        subfilter_stress!(τvw,w,v,w̅,v̅; kernel, cutoff)

        println("computing KE terms at $(time()-t0)s...")
        ue, ve = compute!(Field(u̅ * E̅ₖ)), compute!(Field(v̅ * E̅ₖ))
        D = KernelFunctionOperation{Center, Center, Center}(Oceananigans.Operators.div_xyᶜᶜᶜ, grid, ue,ve)
        T1 = compute!(Field(D + ∂z(w̅ * E̅ₖ)))
        println("mean KE term 1 done at $(time()-t0)s")
        up, vp = compute!(Field(u̅ * p̅)), compute!(Field(v̅ * p̅))
        D = KernelFunctionOperation{Center, Center, Center}(Oceananigans.Operators.div_xyᶜᶜᶜ, grid, up,vp)
        T2 = compute!(Field(D + ∂z(w̅ * p̅)))
        println("mean KE term 2 done at $(time()-t0)s")
        uτx = compute!(Field(τuu * u̅ + τvu * v̅))
        uτy = compute!(Field(τuv * u̅ + τvv * v̅))
        uτz = compute!(Field(τuw * u̅ + τvw * v̅))
        D = KernelFunctionOperation{Center, Center, Center}(Oceananigans.Operators.div_xyᶜᶜᶜ, grid, uτx,uτy)
        T3 = compute!(Field(D + ∂z(uτz)))
        println("mean KE term 3 done at $(time()-t0)s")

        B̅ = compute!(Field(b̅ * w̅)) 
        println("mean KE term B̅ done at $(time()-t0)s")
        uF̅ = compute!(Field(u̅ * F̅x + v̅ * F̅y))
        println("mean KE term uF̅ done at $(time()-t0)s")
        uA̅ = compute!(Field(u̅ * A̅u + v̅ * A̅v))
        println("mean KE term uA̅ done at $(time()-t0)s")

        Πₕ = compute!(Field(-(τuu * ∂x(u̅) + τvv * ∂y(v̅) + τuv * ∂y(u̅) + τvu * ∂x(v̅))))
        println("mean KE term Πₕ done at $(time()-t0)s")
        #Πδ = compute!(Field(-(τuU + τvV) * (∂x(u̅) + ∂y(v̅))/2))
        #println("mean KE term Πδ done at $(time()-t0)s")
        Πᵥ = compute!(Field(@at (Center, Center, Center) -(τuw * ∂z(u̅) + τvw * ∂z(v̅))))
        println("mean KE term Πᵥ done at $(time()-t0)s")
        mke = [havg(T1), havg(T2), havg(T3), havg(uA̅), havg(B̅), havg(uF̅), havg(Πₕ), havg(Πᵥ)]

        ue, ve = compute!(Field(u * Eₖᵖ)), compute!(Field(v * Eₖᵖ))
        D = KernelFunctionOperation{Center, Center, Center}(Oceananigans.Operators.div_xyᶜᶜᶜ, grid, ue,ve)
        T1 = compute!(Field(D + ∂z(w * Eₖᵖ)))
        println("finescale KE term 1 done at $(time()-t0)s")
        up, vp = compute!(Field(uᵖ * pᵖ)), compute!(Field(vᵖ * pᵖ))
        D = KernelFunctionOperation{Center, Center, Center}(Oceananigans.Operators.div_xyᶜᶜᶜ, grid, up,vp)
        T2 = compute!(Field(D + ∂z(wᵖ * pᵖ)))
        println("finescale KE term 2 done at $(time()-t0)s")
        uτx = compute!(Field(τuu * uᵖ + τvu * vᵖ))
        uτy = compute!(Field(τuv * uᵖ + τvv * vᵖ))
        uτz = compute!(Field(τuw * uᵖ + τvw * vᵖ))
        D = KernelFunctionOperation{Center, Center, Center}(Oceananigans.Operators.div_xyᶜᶜᶜ, grid, uτx,uτy)
        T3 = compute!(Field(D + ∂z(uτz)))
        println("finescale KE term 3 done at $(time()-t0)s")

        Bᵖ  = compute!(Field(bᵖ * wᵖ))
        println("finescale KE term Bᵖ done at $(time()-t0)s")
        uFᵖ = compute!(Field(uᵖ * Fxᵖ + vᵖ * Fyᵖ))
        println("finescale KE term uFᵖ done at $(time()-t0)s")
        uAᵖ = compute!(Field(uᵖ * Auᵖ + vᵖ * Avᵖ))
        println("finescale KE term uAᵖ done at $(time()-t0)s")

        Trₕ = compute!(Field(-(uᵖ * (uᵖ * ∂x(u̅) + vᵖ * ∂y(u̅)) + vᵖ * (uᵖ * ∂x(v̅) + vᵖ * ∂y(v̅)))))
        println("finescale KE term Trₕ done at $(time()-t0)s")
        Trᵥ = compute!(Field(-(uᵖ * ∂z(u̅) + vᵖ * ∂z(v̅)) * wᵖ))
        println("finescale KE term Trᵥ done at $(time()-t0)s")

        Pₕ = compute!(Field(-(τuu * ∂x(uᵖ) + τvv * ∂y(vᵖ) + τuv * ∂y(uᵖ) + τvu * ∂x(vᵖ))))
        println("finescale KE term Pₕ done at $(time()-t0)s")
        Pᵥ = compute!(Field(@at (Center, Center, Center) -(τuw * ∂z(uᵖ) + τvw * ∂z(vᵖ))))
        println("finescale KE term Pᵥ done at $(time()-t0)s")
        fke = [havg(T1), havg(T2), havg(T3), havg(uAᵖ), havg(uFᵖ), havg(Bᵖ), havg(Trₕ), havg(Trᵥ), havg(Pₕ), havg(Pᵥ)]

        println("computing residual heat fluxes at $(time()-t0)s...")
        qub = XFaceField(u.grid)
        qvb = YFaceField(v.grid)
        qwb = ZFaceField(w.grid)
        subfilter_stress!(qub,u,b,u̅,b̅; kernel, cutoff)
        subfilter_stress!(qvb,v,b,v̅,b̅; kernel, cutoff)
        subfilter_stress!(qwb,w,b,w̅,b̅; kernel, cutoff)
    
        println("computing PE terms at $(time()-t0)s...")
        ue, ve = compute!(Field(u̅ * E̅ₚ)), compute!(Field(v̅ * E̅ₚ))
        D = KernelFunctionOperation{Center, Center, Center}(Oceananigans.Operators.div_xyᶜᶜᶜ, grid, ue,ve)
        T1 = compute!(Field(D + ∂z(w̅ * E̅ₚ)))
        println("mean PE term 1 done at $(time()-t0)s")
        D = KernelFunctionOperation{Center, Center, Center}(Oceananigans.Operators.div_xyᶜᶜᶜ, grid, qub,qvb)
        Dq = compute!(Field(D + ∂z(qwb)))
        T2 = CenterField(b.grid)
        set!(T2, interior(Dq,:,:,:).*reshape(zT,(1,1,:))) 
        println("mean PE term 2 done at $(time()-t0)s")
        T3 = CenterField(b.grid)
        set!(T3, interior(Q̅, :,:,:).*reshape(-zT,(1,1,:)))
        println("mean PE term 3 done at $(time()-t0)s")
        bA̅ = XFaceField(u.grid)
        set!(bA̅, interior(A̅b, :,:,:).*reshape(-zT,(1,1,:)))
        println("mean PE term bA̅ done at $(time()-t0)s")
        mpe = [havg(T1), havg(bA̅), havg(T2), havg(T3)]

        ue, ve = compute!(Field(u * Eₚᵖ)), compute!(Field(v * Eₚᵖ))
        D = KernelFunctionOperation{Center, Center, Center}(Oceananigans.Operators.div_xyᶜᶜᶜ, grid, ue,ve)
        T1 = compute!(Field(D + ∂z(w * Eₚᵖ)))
        println("finescale PE term 1 done at $(time()-t0)s")
        ub, vb = compute!(Field(uᵖ * b̅)), compute!(Field(vᵖ * b̅))
        D = KernelFunctionOperation{Center, Center, Center}(Oceananigans.Operators.div_xyᶜᶜᶜ, grid, ub,vb)
        Dub = compute!(Field(D + ∂z(wᵖ * b̅)))
        T2 = CenterField(b.grid)
        set!(T2, interior(Dub,:,:,:).*reshape(zT,(1,1,:)))
        println("finescale PE term 2 done at $(time()-t0)s")
        T3 = CenterField(b.grid)
        set!(T3, interior(Qᵖ, :,:,:).*reshape(-zT,(1,1,:)))
        println("finescale PE term 3 done at $(time()-t0)s")
        bAᵖ = XFaceField(u.grid)
        set!(bAᵖ, interior(Abᵖ, :,:,:).*reshape(-zT,(1,1,:)))
        println("finescale PE term bAᵖ done at $(time()-t0)s")
        fpe = [havg(T1), havg(bAᵖ), havg(T2), havg(T3)]

        ue, ve = compute!(Field(u̅ * Σ̅b)), compute!(Field(v̅ * Σ̅b))
        D = KernelFunctionOperation{Center, Center, Center}(Oceananigans.Operators.div_xyᶜᶜᶜ, grid, ue,ve)
        T1 = compute!(Field(D + ∂z(w̅ * Σ̅b)))
        println("mean BV term 1 done at $(time()-t0)s")
        D = KernelFunctionOperation{Center, Center, Center}(Oceananigans.Operators.div_xyᶜᶜᶜ, grid, qub,qvb)
        T2 = compute!(Field(-Dq * b̅))
        println("mean BV term 2 done at $(time()-t0)s")
        T3 = compute!(Field(Q̅ * b̅))
        println("mean BV term 3 done at $(time()-t0)s")
        bA̅ = compute!(Field(A̅b * b̅))
        println("mean BV term bA̅ done at $(time()-t0)s")
        mbv = [havg(T1), havg(bA̅), havg(T2), havg(T3)]

        ue, ve = compute!(Field(u * Σbᵖ)), compute!(Field(v * Σbᵖ))
        D = KernelFunctionOperation{Center, Center, Center}(Oceananigans.Operators.div_xyᶜᶜᶜ, grid, ue,ve)
        T1 = compute!(Field(D + ∂z(w * Σbᵖ)))
        println("finescale BV term 1 done at $(time()-t0)s")
        T2 = compute!(Field(-Dub * bᵖ))
        println("finescale BV term 2 done at $(time()-t0)s")
        T3 = compute!(Field(-Dq * bᵖ))
        println("finescale BV term 3 done at $(time()-t0)s")
        T4 = compute!(Field(Qᵖ * bᵖ))
        println("finescale BV term 4 done at $(time()-t0)s")
        bAᵖ = compute!(Field(Abᵖ * bᵖ))
        println("finescale BV term bAᵖ done at $(time()-t0)s")
        fbv = [havg(T1), havg(bAᵖ), havg(T2), havg(T3), havg(T4)]
             
        return E̅ₖ, Eₖᵖ, mke, fke, E̅ₚ, Eₚᵖ, mpe, fpe, Σ̅b, Σbᵖ, mbv, fbv
    end
end

""" streamfunction computation kernel """
@kernel function _streamfunction!(ψ, u, grid)
    i, k = @index(Global, NTuple)

    @inbounds ψ[i, 1, k] = 0
    for j in 2:grid.Ny+1
        @inbounds ψ[i, j, k] = ψ[i, j-1, k] - u[i, j-1, k] * Δyᶠᶜᶜ(i, j, k, grid)
    end
end

""" streamfunction computation kernel """
@kernel function _zstreamfunction!(ψ, v, grid)
    i, j = @index(Global, NTuple)

    @inbounds ψ[i, j, 1] = 0
    for k in 2:grid.Nz+1
        @inbounds ψ[i, j, k] = ψ[i, j-1, k] - v[i, j-1, k] * Δzᶜᶠᶜ(i, j, k, grid)
    end
end

""" barotropic streamfunction """
function Ψ(snapshots, i)
    u = snapshots[:u][i]
    grid = u.grid
    ψ = Field{Face, Face, Center}(grid)
    arch = architecture(grid)
    launch!(arch, grid, :xz, _streamfunction!, ψ, u, grid)
    return ψ
end

""" vertical streamfunction """
function Ψz(snapshots, i)
    v = snapshots[:v][i]
    grid = v.grid
    ψ = Field{Center, Face, Face}(grid)
    arch = architecture(grid)
    launch!(arch, grid, :xy, _zstreamfunction!, ψ, v, grid)
    return ψ
end

""" mixed layer average kernel """
@kernel function _zMLaverage!(ψ, v, h, z, grid)
    i, j = @index(Global, NTuple)   

    @inbounds ψ[i, j, 1] = 0
    k0 = findfirst(z .> -h[i,j,1] + Δzᶜᶠᶜ(i, j, 1, grid)/2)
    for k in k0-1:grid.Nz
        if k == k0-1
            @inbounds ψ[i, j, 1] = v[i, j, k] * (Δzᶜᶠᶜ(i, j, k, grid)/2 + z[k] + h[i,j,1])
        else
            @inbounds ψ[i, j, 1] = ψ[i, j, 1] + v[i, j, k] * Δzᶜᶠᶜ(i, j, k, grid)
        end
    end
    @inbounds ψ[i, j, 1] = h[i,j,1] > 0 ? ψ[i, j, 1]/h[i,j,1] : 0
end

""" mixed layer average function """
function MLaverage(snapshots, i, v)
    _, _, z = nodes(v)
    h = compute!(MLD(snapshots,i; threshold = 0.03))
    grid = v.grid
    ψ = Field{Center, Center, Nothing}(grid)
    arch = architecture(grid)
    launch!(arch, grid, :xy, _zMLaverage!, ψ, v, h, z, grid)
    return ψ
end
