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

    w² = (w - mean(w,dims=(1,2)))^2
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
function wb(snapshots, i)
    α = parameters.α
    g = parameters.g

    w = snapshots[:w][i]
    T = snapshots[:T][i]

    return α * g * T * w
end

""" filtered zonal buoyancy flux """
function ub_filtered(snapshots, i; cutoff = 20kilometer)
    α = parameters.α
    g = parameters.g

    u = snapshots[:u][i]
    b = α * g * snapshots[:T][i]
    ul, uh = symmetric_filtering(u; cutoff)
    bl, bh = symmetric_filtering(b; cutoff)

    ub = (u - mean(u,dims=(1,2))) * (b - mean(b,dims=(1,2)))
    ulbl = (ul - mean(ul,dims=(1,2))) * (bl - mean(bl,dims=(1,2)))
    uhwh = (uh - mean(uh,dims=(1,2))) * (bh - mean(bh,dims=(1,2)))

    return ub, ulbl, uhwh
end

""" filtered meridional buoyancy flux """
function vb_filtered(snapshots, i; cutoff = 20kilometer)
    α = parameters.α
    g = parameters.g

    v = snapshots[:v][i]
    b = α * g * snapshots[:T][i]
    vl, vh = symmetric_filtering(v; cutoff)
    bl, bh = symmetric_filtering(b; cutoff)

    vb = (v - mean(v,dims=(1,2))) * (b - mean(b,dims=(1,2)))
    vlbl = (vl - mean(vl,dims=(1,2))) * (bl - mean(bl,dims=(1,2)))
    vhwh = (vh - mean(vh,dims=(1,2))) * (bh - mean(bh,dims=(1,2)))

    return vb, vlbl, vhwh
end

""" filtered vertical buoyancy flux """
function wb_filtered(snapshots, i; cutoff = 20kilometer)
    α = parameters.α
    g = parameters.g

    w = snapshots[:w][i]
    b = α * g * snapshots[:T][i]
    wl, wh = symmetric_filtering(w; cutoff)
    bl, bh = symmetric_filtering(b; cutoff)

    wb = (w - mean(w,dims=(1,2))) * (b - mean(b,dims=(1,2)))
    wlbl = (wl - mean(wl,dims=(1,2))) * (bl - mean(bl,dims=(1,2)))
    whwh = (wh - mean(wh,dims=(1,2))) * (bh - mean(bh,dims=(1,2)))

    return wb, wlbl, whwh
end

""" horizontal kinetic energy """
function KE(snapshots, i)
    u = snapshots[:u][i]
    v = snapshots[:v][i]

    return 0.5 * (u^2 + v^2)
end

""" mixed layer depth """
function MLD(snapshots, i; threshold = 0.03)
    α  = parameters.α
    ρ₀ = parameters.ρ₀
    T    = snapshots[:T][i]
    grid = T.grid
    h    = MixedLayerDepth(grid, (; T); ΔT = abs(threshold / ρ₀ / α))
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
function ζ(snapshots, i)
    u = snapshots[:u][i]
    v = snapshots[:v][i]

    return ∂x(v) - ∂y(u)
end

""" horizontal divergence """
function δ(snapshots, i)
    u = snapshots[:u][i]
    v = snapshots[:v][i]

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

""" spectral horizontal advection """
function Ah(snapshots, i; nfactor=100)
    u = snapshots[:u][i]
    v = snapshots[:v][i]
    KE = 0.5 * (u^2 + v^2)
    ζ = ∂x(v) - ∂y(u)
    xu, yu, zu = nodes(u)
    xv, yv, _ = nodes(v)
    Nz = length(zu)

    Au = compute!(Field(- ∂x(KE) + v * ζ))
    Av = compute!(Field(- ∂y(KE) - u * ζ))

    # Fourier transform
    Au1 = isotropic_powerspectrum(interior(u, :, :, 1), interior(Au, :, :, 1), xu, yu; nfactor)
    Av1 = isotropic_powerspectrum(interior(v, :, :, 1), interior(Av, :, :, 1), xv, yv; nfactor)
    A = zeros(Nz,length(Au1.spec))
    A[1,:] = real.(Au1.spec + Av1.spec)
    for k = 2:Nz
        Auk = isotropic_powerspectrum(interior(Au, :, :, k), interior(u, :, :, k), xu, yu; nfactor)
        Avk = isotropic_powerspectrum(interior(Av, :, :, k), interior(v, :, :, k), xv, yv; nfactor)
        A[k,:] = real.(Auk.spec + Avk.spec)
    end

    return A, zu, Au1.freq
end

""" spectral vertical advection """
function Av(snapshots, i; nfactor=100)
    u = compute!(Field(@at (Center, Center, Face) snapshots[:u][i]))
    v = compute!(Field(@at (Center, Center, Face) snapshots[:v][i]))
    w = snapshots[:w][i]

    xu, yu, zu = nodes(u)
    xv, yv, _ = nodes(v)
    Nz = length(zu)

    Au = compute!(Field(- w * ∂z(u)))
    Av = compute!(Field(- w * ∂z(v)))

    Au1 = isotropic_powerspectrum(interior(Au, :, :, 1), interior(u, :, :, 1), xu, yu; nfactor)
    A = zeros(Nz,length(Au1.spec))
    for k = 2:Nz
        Auk = isotropic_powerspectrum(interior(Au, :, :, k), interior(u, :, :, k), xu, yu; nfactor)
        Avk = isotropic_powerspectrum(interior(Av, :, :, k), interior(v, :, :, k), xv, yv; nfactor)
        A[k,:] = real.(Auk.spec .+ Avk.spec)
    end

    return A, zu, Au1.freq

end

""" spectral convertion of potential to kinetic energy """
function Ac(snapshots, i; nfactor=100)
    w = compute!(Field(@at (Center, Center, Center) snapshots[:w][i]))
    T = snapshots[:T][i]

    xT, yT, zT = nodes(T)
    Nz = length(zT)

    α = parameters.α
    g = parameters.g
    B = compute!(Field(α * g * T))
    C1 = isotropic_powerspectrum(interior(B, :, :, 1), interior(w, :, :, 1), xT, yT; nfactor)
    C = zeros(Nz,length(C1.spec))
    C[1,:] = real.(C1.spec)
    for k = 2:Nz
        Ck = isotropic_powerspectrum(interior(B, :, :, k), interior(w, :, :, k), xT, yT; nfactor)
        C[k,:] = real.(Ck.spec)
    end

    return C, zT, C1.freq
end

""" spectral 3D pressure work """
function Ap(snapshots, i; nfactor=100)
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

    Au1 = isotropic_powerspectrum(interior(px, :, :, 1), interior(u, :, :, 1), xu, yu; nfactor)
    Av1 = isotropic_powerspectrum(interior(py, :, :, 1), interior(v, :, :, 1), xv, yv; nfactor)
    Aw1 = isotropic_powerspectrum(interior(pz, :, :, 1), interior(w, :, :, 1), xw, yw; nfactor)
    A = zeros(length(zu),length(Au1.spec))
    A[1,:] = - real.(Au1.spec + Av1.spec + Aw1.spec)
    for k = 2:length(zu)
        Auk = isotropic_powerspectrum(interior(px, :, :, k), interior(u, :, :, k), xu, yu; nfactor)
        Avk = isotropic_powerspectrum(interior(py, :, :, k), interior(v, :, :, k), xv, yv; nfactor)
        Awk = isotropic_powerspectrum(interior(pz, :, :, k), interior(w, :, :, k), xw, yw; nfactor)
        A[k,:] = - real.(Auk.spec + Avk.spec + Awk.spec)
    end

    return A, zu, Au1.freq
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