using Oceananigans.Operators: div_xyᶜᶜᶜ
using Oceananigans.Operators
using Oceananigans.Utils: launch!
using Oceananigans.Architectures: on_architecture
using Oceananigans.Fields: interpolate!, interpolate
using KernelAbstractions: @kernel, @index
using Statistics: mean, var
using Interpolations
using Oceananigans.Operators                 

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
function BLD(snapshots, i; threshold = 1e-5)
    κ  = snapshots[:κu][i]
    
    grid = κ.grid
    h    = BoundaryLayerDepth(grid, (; κ); κc = threshold)
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
function subfilter_stress!(τ, u, v, u̅, v̅; kernel=:tophat, cutoff=4kilometer, border=:circular,method=:physical,use_gpu=false,plans=nothing)
    # Compute the product field (u*v) on the fly.
    product_field = compute!(Field(u * v))
    # Coarse-grain the product field and store the result in τ.
    coarse_graining!(product_field, τ; kernel, cutoff, border, method, use_gpu, plans)

    # Compute the product of the filtered fields and subtract.
    # We assume that the multiplication is elementwise.
    stress = compute!(Field(τ - (u̅ * v̅)))
    # Replace the interior of τ with the computed stress.
    set!(τ, interior(stress))
    fill_halo_regions!(τ)
    return nothing
end

"""
    TKE(snapshot; kernel=:lanczos, cutoff=300, border=:circular)

Compute the turbulent kinetic energy (TKE) from a snapshot. The function uses a
coarse-graining (filtering) procedure based on a specified kernel, cutoff and border 
for (non)periodic boundary conditions. It assumes the use of Oceananigans Field types.
"""
function TKE(snapshot; kernel=:gaussian, cutoff=100, Δh = 4.8828125, Lx = 1e5, Ly = 1e5, 
                       border=:circular, method=:physical,use_gpu=false,plans=nothing)
    t0 = time()

    set_value!(; Δh, Lx, Ly)

    # Extract snapshot fields.
    u0 = snapshot[:u]
    v0 = snapshot[:v]
    w0 = snapshot[:w]

    u = XFaceField(u0.grid,Float32)
    v = YFaceField(v0.grid,Float32)
    w = ZFaceField(w0.grid,Float32)
    set!(u, interior(u0))
    set!(v, interior(v0))
    set!(w, interior(w0))
    fill_halo_regions!(u)
    fill_halo_regions!(v)
    fill_halo_regions!(w)

    # Allocate filtered velocity and buoyancy fields.
    u̅ = XFaceField(u.grid,Float32)
    v̅ = YFaceField(v.grid,Float32)
    w̅ = ZFaceField(w.grid,Float32)

    # --- Coarse-grain (filter) the primary fields. ---
    for (orig, filt) in ((u, u̅), (v, v̅), (w, w̅))
        coarse_graining!(orig, filt; kernel, cutoff, border, method,use_gpu,plans)
    end
    println("Computed filtered velocity fields at $(time() - t0)s...")

    # --- Compute subfilter TKE ---
    τuu = XFaceField(u.grid,Float32)
    τvv = YFaceField(v.grid,Float32)
    τww = ZFaceField(w.grid,Float32)
    subfilter_stress!(τww, w, w, w̅, w̅; kernel, cutoff, border, method,use_gpu,plans)
    subfilter_stress!(τuu, u, u, u̅, u̅; kernel, cutoff, border, method,use_gpu,plans)
    subfilter_stress!(τvv, v, v, v̅, v̅; kernel, cutoff, border, method,use_gpu,plans)

    can_use_gpu = use_gpu && CUDA.functional()
    if can_use_gpu
        # Move the filtered fields to the GPU.
        gpu_field(cpu_field) = on_architecture(GPU(), cpu_field)
        return gpu_field(u̅), gpu_field(v̅), gpu_field(w̅), gpu_field(τuu), gpu_field(τvv), gpu_field(τww)
    else
        return u̅, v̅, w̅, τuu, τvv, τww
    end
end

function along_front_averages(snapshots; i=0, kernel=:gaussian, cutoff=100, Δh = 4.8828125, Lx = 1e5, Ly = 1e5, 
    border=:circular, method=:physical, use_gpu=false,plans=nothing, m₀ = 60,to_grid=nothing)
    t0 = time()

    set_value!(; Δh, Lx, Ly, m₀)
    @info "to grid $(!isnothing(to_grid)): $(to_grid)"

    # Extract snapshot fields.
    u0 = i >= 1 ? snapshots[:u][i] : snapshots[:u]
    v0 = i >= 1 ? snapshots[:v][i] : snapshots[:v]
    w0 = i >= 1 ? snapshots[:w][i] : snapshots[:w]
    T0 = i >= 1 ? snapshots[:T][i] : snapshots[:T]
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

    # Retrieve physical parameters.
    α = parameters.α
    g = parameters.g

    # Compute the buoyancy field B = α*g*T.
    B = compute!(Field(α * g * T))

    # Allocate filtered velocity and buoyancy fields.
    u̅ = XFaceField(u.grid,Float32)
    v̅ = YFaceField(v.grid,Float32)
    w̅ = ZFaceField(w.grid,Float32)
    B̅ = CenterField(B.grid,Float32)

    # --- Coarse-grain (filter) the primary fields. ---
    for (orig, filt) in ((u, u̅), (v, v̅), (w, w̅), (B, B̅))
    coarse_graining!(orig, filt; kernel, cutoff, border, method,use_gpu,plans)
    end
    @info "Computed filtered fields at $(time() - t0)s..."

    if !isnothing(to_grid)
        itp_u̅ = XFaceField(to_grid,Float32)
        itp_v̅ = YFaceField(to_grid,Float32)
        itp_w̅ = ZFaceField(to_grid,Float32)
        itp_B̅ = CenterField(to_grid,Float32)
        interpolate!(itp_u̅, u̅)
        interpolate!(itp_v̅, v̅)
        interpolate!(itp_w̅, w̅)
        interpolate!(itp_B̅, B̅)
        uᵃ = mean(itp_u̅, dims=2)
        vᵃ = mean(itp_v̅, dims=2)
        wᵃ = mean(itp_w̅, dims=2)
        Bᵃ = mean(itp_B̅, dims=2)
    else
        uᵃ = mean(u̅, dims=2)
        vᵃ = mean(v̅, dims=2)
        wᵃ = mean(w̅, dims=2)
        Bᵃ = mean(B̅, dims=2)
    end
    return uᵃ, vᵃ, wᵃ, Bᵃ

end


"""
    coarse_grained_fluxes(snapshots, i; kernel=:tophat, cutoff=4kilometer)

Compute the coarse-grained velocities and cross-scale fluxes (i.e. residual stresses and
transfer terms) from a snapshot indexed by `i` in `snapshots`. The function uses a
coarse-graining (filtering) procedure based on a specified kernel and cutoff. It assumes
periodic boundary conditions (so that FFTs can be used) and the use of Oceananigans Field
types.
"""
function coarse_grained_fluxes(snapshots, iU, iV, varᵃ; i=0, hy=false, kernel=:gaussian, cutoff=100, Δh = 4.8828125, Lx = 1e5, Ly = 1e5, 
                                          border=:circular, method=:physical,use_gpu=false,plans=nothing, m₀ = 60,to_grid=nothing)
    t0 = time()

    set_value!(; Δh, Lx, Ly, m₀)
    @info "to grid $(!isnothing(to_grid)): $(to_grid)"
    can_use_gpu = CUDA.functional()
    @info "Using GPU: $(can_use_gpu)"

    # Extract snapshot fields.
    u0 = i >= 1 ? snapshots[:u][i] : snapshots[:u]
    v0 = i >= 1 ? snapshots[:v][i] : snapshots[:v]
    w0 = i >= 1 ? snapshots[:w][i] : snapshots[:w]
    T0 = i >= 1 ? snapshots[:T][i] : snapshots[:T]
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
    xT, yT, zT = nodes(T)
    kT = length(zT)
    k0 = findlast(zT .< -parameters.m₀)

    # Construct background (U,V) fields.
    U = XFaceField(u.grid,Float32)
    V = YFaceField(v.grid,Float32)
    interpolate!(U, iU)
    interpolate!(V, iV)
    @info "Interpolated U size $(size(interior(U)))"

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
    U̅ = XFaceField(u.grid,Float32)
    V̅ = YFaceField(v.grid,Float32)

    # --- Coarse-grain (filter) the primary fields. ---
    for (orig, filt) in ((u, u̅), (v, v̅), (w, w̅), (B, B̅), (U, U̅), (V, V̅))
        coarse_graining!(orig, filt; kernel, cutoff, border, method,use_gpu,plans)
    end
    @info "Computed filtered fields at $(time() - t0)s..."

    # Move the filtered fields to the GPU.
    # Pre-transfer all needed fields once
    if can_use_gpu
        if !isnothing(to_grid)
            itp_u̅ = XFaceField(to_grid,Float32)
            itp_v̅ = YFaceField(to_grid,Float32)
            itp_w̅ = ZFaceField(to_grid,Float32)
            itp_B̅ = CenterField(to_grid,Float32)
            itp_U̅ = XFaceField(to_grid,Float32)
            itp_V̅ = YFaceField(to_grid,Float32)
            itp_τuuₜ = CenterField(to_grid,Float32)
            itp_τvvₜ = CenterField(to_grid,Float32)
            itp_τuvₜ = CenterField(to_grid,Float32)
            itp_τvuₜ = CenterField(to_grid,Float32)
            itp_τuw = CenterField(to_grid,Float32)
            itp_τvw = CenterField(to_grid,Float32)
            interpolate!(itp_u̅, u̅)
            interpolate!(itp_v̅, v̅)
            interpolate!(itp_w̅, w̅)
            interpolate!(itp_B̅, B̅)
            interpolate!(itp_U̅, U̅)
            interpolate!(itp_V̅, V̅)
            xi, _, _ = nodes(itp_B̅)
        end
        u̅_gpu = on_architecture(GPU(), !isnothing(to_grid) ? itp_u̅ : u̅)
        v̅_gpu = on_architecture(GPU(), !isnothing(to_grid) ? itp_v̅ : v̅)
        w̅_gpu = on_architecture(GPU(), !isnothing(to_grid) ? itp_w̅ : w̅)
        B̅_gpu = on_architecture(GPU(), !isnothing(to_grid) ? itp_B̅ : B̅)
        U̅_gpu = on_architecture(GPU(), !isnothing(to_grid) ? itp_U̅ : U̅)
        V̅_gpu = on_architecture(GPU(), !isnothing(to_grid) ? itp_V̅ : V̅)
    else
        u̅_gpu = u̅
        v̅_gpu = v̅
        w̅_gpu = w̅
        U̅_gpu = U̅
        V̅_gpu = V̅
    end
    @info "Interpolated filtered fields at $(time() - t0)s..."
    @info "Interpolated u̅_gpu size $(size(interior(u̅_gpu)))" 
    
    du_dy = compute!(Field(∂y(u̅_gpu))) 
    dv_dy = compute!(Field(∂y(v̅_gpu)))
    dw_dy = compute!(Field(∂y(w̅_gpu))) 
    CUDA.pool_status()

    if !hy
        # Compute area averages and broadcast them to full fields
        # uᵃ = mean(u̅_gpu, dims=2)
        # vᵃ = mean(v̅_gpu, dims=2)
        # wᵃ = mean(w̅_gpu, dims=2)
        # Bᵃ = mean(B̅_gpu, dims=2)

        # Create full fields from the averaged values for GPU kernel compatibility
        grid = u̅_gpu.grid
        uᵃ = Field{Face, Nothing, Center}(grid); set!(uᵃ, varᵃ[1])
        vᵃ = Field{Center, Nothing, Center}(grid); set!(vᵃ, varᵃ[2])
        wᵃ = Field{Center, Nothing, Face}(grid);   set!(wᵃ, varᵃ[3])
        Bᵃ = Field{Center, Nothing, Center}(grid); set!(Bᵃ, varᵃ[4])
        @info "average velocities and buoyancy done at $(time() - t0)s"

        # Compute deviations in-place
        uˢ = compute!(Field(u̅_gpu - uᵃ))
        vˢ = compute!(Field(v̅_gpu - vᵃ))
        wˢ = compute!(Field(w̅_gpu - wᵃ))
        Bˢ = compute!(Field(B̅_gpu - Bᵃ))
        @info "submeso velocities and buoyancy done at $(time() - t0)s"
        CUDA.pool_status()

        uᵃrep = Field{Face, Center, Center}(grid); set!(uᵃrep, repeat(varᵃ[1],1,512*2,1))
        vᵃrep = Field{Center, Center, Center}(grid); set!(vᵃrep, repeat(varᵃ[2],1,512*2,1))
        wᵃrep = Field{Center, Center, Face}(grid);   set!(wᵃrep, repeat(varᵃ[3],1,512*2,1))
        Bᵃrep = Field{Center, Center, Center}(grid); set!(Bᵃrep, repeat(varᵃ[4],1,512*2,1))
        CUDA.pool_status()

        # Compute Pᵃ
        Pᵃₕ = compute!(Field(- uˢ * (uˢ+U̅_gpu) * ∂x(uᵃrep) - vˢ * (uˢ+U̅_gpu) * ∂x(vᵃrep)))
        @info "Transfer term Pᵃₕ done at $(time() - t0)s"
        Pᵃᵥg = compute!(Field(- vˢ * wˢ * ∂x(Bᵃrep) / f))
        Pᵃᵥ = compute!(Field(- (uˢ * wˢ * ∂z(uᵃrep) + vˢ * wˢ * ∂z(vᵃrep) +
                                wˢ * (uˢ+U̅_gpu) * ∂x(wᵃrep) + wˢ^2 * ∂z(wᵃrep))))
        @info "Transfer term Pᵃᵥ done at $(time() - t0)s"
        Pᵃ = compute!(Field(Pᵃₕ + Pᵃᵥ))
        @info "Transfer term Pᵃ done at $(time() - t0)s"
        CUDA.pool_status()

        # Compute averages of submeso velocities
        du_dx = compute!(Field(∂x(uˢ)))  
        dv_dx = compute!(Field(∂x(vˢ))) 
        dw_dx = compute!(Field(∂x(wˢ))) 
        du_dz = compute!(Field(∂z(uˢ)))
        dv_dz = compute!(Field(∂z(vˢ))) 
        dw_dz = compute!(Field(∂z(wˢ))) 
        uˢuˢ_avg = mean(uˢ * (uˢ+U̅_gpu), dims=2) #Field{Center, Center, Center}(grid)
        uˢvˢ_avg = mean(uˢ * (vˢ+V̅_gpu), dims=2) #Field{Center, Center, Center}(grid)
        uˢwˢ_avg = mean(uˢ * wˢ, dims=2)#Field{Center, Center, Center}(grid)
        vˢuˢ_avg = mean(vˢ * (uˢ+U̅_gpu), dims=2) #Field{Center, Center, Center}(grid)
        vˢvˢ_avg = mean(vˢ * (vˢ+V̅_gpu), dims=2) #Field{Center, Center, Center}(grid)
        vˢwˢ_avg = mean(vˢ * wˢ, dims=2)#Field{Center, Center, Center}(grid)
        wˢuˢ_avg = mean(wˢ * (uˢ+U̅_gpu), dims=2) #Field{Center, Center, Center}(grid)
        wˢvˢ_avg = mean(wˢ * (vˢ+V̅_gpu), dims=2) #Field{Center, Center, Center}(grid)
        wˢwˢ_avg = mean(wˢ^2, dims=2) #Field{Center, Center, Face}(grid)
        # set!(uˢuˢ_avg, mean(uˢ^2, dims=2))
        # set!(uˢvˢ_avg, mean(uˢ*vˢ, dims=2))
        # set!(uˢwˢ_avg, mean(uˢ*wˢ, dims=2))
        # set!(vˢvˢ_avg, mean(vˢ^2, dims=2))
        # set!(vˢwˢ_avg, mean(vˢ*wˢ, dims=2))
        # set!(wˢwˢ_avg, mean(wˢ^2, dims=2))
        @info "average submeso velocities done at $(time() - t0)s"
        CUDA.pool_status()

        # Compute Pˢ
        Pˢₕ = compute!(Field(-(uˢuˢ_avg * du_dx + uˢvˢ_avg * du_dy  +
                                vˢuˢ_avg * dv_dx + vˢvˢ_avg * dv_dy)))
        Pˢᵥg = compute!(Field(-(-uˢwˢ_avg * ∂y(Bˢ) + vˢwˢ_avg * ∂x(Bˢ))/f))
        Pˢᵥ = compute!(Field(-(uˢwˢ_avg * du_dz + vˢwˢ_avg * dv_dz + 
                                wˢuˢ_avg * dw_dx + wˢvˢ_avg * dw_dy + wˢwˢ_avg * dw_dz)))
        Pˢ = compute!(Field(Pˢₕ + Pˢᵥ))
        @info "Transfer term Pˢ done at $(time() - t0)s"
        CUDA.pool_status()
    end

    # --- Compute subfilter stresses ---
    # First set: fluxes associated with the u-momentum.
    τuuₜ = XFaceField(u.grid,Float32)
    τuvₜ = YFaceField(v.grid,Float32)
    τuw = ZFaceField(w.grid,Float32)
    subfilter_stress!(τuuₜ, u, u+U, u̅, u̅+U̅; kernel, cutoff, border, method,use_gpu,plans)
    subfilter_stress!(τuvₜ, v+V, u, v̅+V̅, u̅; kernel, cutoff, border, method,use_gpu,plans)
    subfilter_stress!(τuw, w, u, w̅, u̅; kernel, cutoff, border, method,use_gpu,plans)

    # Second set: fluxes associated with the v-momentum.
    τvvₜ = YFaceField(v.grid,Float32)
    τvuₜ = XFaceField(u.grid,Float32)
    τvw = ZFaceField(w.grid,Float32)
    subfilter_stress!(τvvₜ, v, v+V, v̅, v̅+V̅; kernel, cutoff, border, method,use_gpu,plans)
    subfilter_stress!(τvuₜ, u+U, v, u̅+U̅, v̅; kernel, cutoff, border, method,use_gpu,plans)
    subfilter_stress!(τvw, w, v, w̅, v̅; kernel, cutoff, border, method,use_gpu,plans)

    # Third set: fluxes associated with the w-momentum.
    τww = ZFaceField(w.grid,Float32)
    subfilter_stress!(τww, w, w, w̅, w̅; kernel, cutoff, border, method,use_gpu,plans)
    if !hy
        τwuₜ = XFaceField(w.grid,Float32)
        τwvₜ = YFaceField(w.grid,Float32)
        subfilter_stress!(τwuₜ, u+U, w, u̅+U̅, w̅; kernel, cutoff, border, method,use_gpu,plans)
        subfilter_stress!(τwvₜ, v+V, w, v̅+V̅, w̅; kernel, cutoff, border, method,use_gpu,plans)
    end

    τuu = XFaceField(u.grid,Float32)
    τvv = YFaceField(v.grid,Float32)
    subfilter_stress!(τuu, u, u, u̅, u̅; kernel, cutoff, border, method,use_gpu,plans)
    subfilter_stress!(τvv, v, v, v̅, v̅; kernel, cutoff, border, method,use_gpu,plans)
    @info "Computed residual stress terms at $(time() - t0)s..."

    @info "τuuₜ extrema: $(extrema(interior(τuuₜ)))"
    @info "τvvₜ extrema: $(extrema(interior(τvvₜ)))"
    @info "τuu extrema: $(extrema(interior(τuu)))"
    @info "τvv extrema: $(extrema(interior(τvv)))"
    @info "τww extrema: $(extrema(interior(τww)))"

    if can_use_gpu
        if !isnothing(to_grid)
            interpolate!(itp_τuuₜ, τuuₜ)
            interpolate!(itp_τvvₜ, τvvₜ)
            interpolate!(itp_τuvₜ, τuvₜ)
            interpolate!(itp_τvuₜ, τvuₜ)
            interpolate!(itp_τuw, τuw)
            interpolate!(itp_τvw, τvw)
            if !hy  
                itp_τwuₜ = CenterField(to_grid,Float32)
                itp_τwvₜ = CenterField(to_grid,Float32)
                itp_τww = CenterField(to_grid,Float32)
                interpolate!(itp_τwuₜ, τwuₜ)
                interpolate!(itp_τwvₜ, τwvₜ)
                interpolate!(itp_τww, τww)
            end
        end
        τuuₜ_gpu = on_architecture(GPU(), !isnothing(to_grid) ? itp_τuuₜ : τuuₜ)
        τvvₜ_gpu = on_architecture(GPU(), !isnothing(to_grid) ? itp_τvvₜ : τvvₜ)
        τuvₜ_gpu = on_architecture(GPU(), !isnothing(to_grid) ? itp_τuvₜ : τuvₜ)
        τvuₜ_gpu = on_architecture(GPU(), !isnothing(to_grid) ? itp_τvuₜ : τvuₜ)
        τuw_gpu = on_architecture(GPU(), !isnothing(to_grid) ? itp_τuw : τuw)
        τvw_gpu = on_architecture(GPU(), !isnothing(to_grid) ? itp_τvw : τvw)
        if !hy  
            τwuₜ_gpu = on_architecture(GPU(), !isnothing(to_grid) ? itp_τwuₜ : τwuₜ)
            τwvₜ_gpu = on_architecture(GPU(), !isnothing(to_grid) ? itp_τwvₜ : τwvₜ)
            τww_gpu = on_architecture(GPU(), !isnothing(to_grid) ? itp_τww : τww)
        end
    else
        τuuₜ_gpu = τuuₜ
        τvvₜ_gpu = τvvₜ
        τuvₜ_gpu = τuvₜ
        τvuₜ_gpu = τvuₜ
        τuw_gpu = τuw
        τvw_gpu = τvw
        if !hy  
            τwuₜ_gpu = τwuₜ
            τwvₜ_gpu = τwvₜ
            τww_gpu = τww
        end
    end

    # --- Compute transfer (flux) terms using derivatives ---
    # Note: The derivative operators (∂x, ∂y, ∂z) are assumed to be available.
    # Πₕ: Horizontal transfer term.
    Πₕ = compute!(Field(-(τuuₜ_gpu * ∂x(u̅_gpu) + τvvₜ_gpu * dv_dy +
                          τuvₜ_gpu * du_dy + τvuₜ_gpu * ∂x(v̅_gpu))))

    @info "Transfer term Πₕ done at $(time() - t0)s"
    CUDA.pool_status()

    # Πᵥ: Vertical transfer term.
    if hy
        Πᵥ = compute!(Field(-(τuw_gpu * ∂z(u̅_gpu) + τvw_gpu * ∂z(v̅_gpu))))
    else
        Πᵥ = compute!(Field(-(τuw_gpu * ∂z(u̅_gpu) + τvw_gpu * ∂z(v̅_gpu) +
                              τww_gpu * ∂z(w̅_gpu) + τwuₜ_gpu * ∂x(w̅_gpu) + τwvₜ_gpu * dw_dy)))
        Πᵥg = compute!(Field( -(-τuw_gpu * ∂y(B̅_gpu) + τvw_gpu * ∂x(B̅_gpu)) / f))
    end
    @info "Transfer term Πᵥ done at $(time() - t0)s"
    CUDA.pool_status()

    if hy
        # Πδ: Diagonal (divergence-related) transfer term.
        Πδ = compute!(Field(-(τuu + τvv) * (∂x(u̅) + ∂y(v̅)) / 2))
        println("Transfer term Πδ done at $(time() - t0)s")

        # Πvgl: Flux transfer term associated with buoyancy gradients, scaled by 1/f.
        Πvgl = compute!(Field( -(-τuw * ∂y(B̅) + τvw * ∂x(B̅)) / f))
        println("Transfer term Πvgl done at $(time() - t0)s")

        return Πₕ, Πδ, Πᵥ, Πvgl
    else
        τwb = CenterField(B.grid,Float32)
        subfilter_stress!(τwb, B, w, B̅, w̅; kernel, cutoff, border, method,use_gpu,plans)
        @info "Transfer term τwb done at $(time() - t0)s"


        # Compute Pᵀ
        τuuₜˢ = compute!(Field(τuuₜ_gpu - mean(τuuₜ_gpu, dims=2)))
        τvvₜˢ = compute!(Field(τvvₜ_gpu - mean(τvvₜ_gpu, dims=2)))
        τuvₜˢ = compute!(Field(τuvₜ_gpu - mean(τuvₜ_gpu, dims=2)))
        τvuₜˢ = compute!(Field(τvuₜ_gpu - mean(τvuₜ_gpu, dims=2)))
        τuwˢ = compute!(Field(τuw_gpu - mean(τuw_gpu, dims=2)))
        τvwˢ = compute!(Field(τvw_gpu - mean(τvw_gpu, dims=2)))
        τwwˢ = compute!(Field(τww_gpu - mean(τww_gpu, dims=2)))
        τwuₜˢ = compute!(Field(τwuₜ_gpu - mean(τwuₜ_gpu, dims=2)))
        τwvₜˢ = compute!(Field(τwvₜ_gpu - mean(τwvₜ_gpu, dims=2)))
        Pᵀₕ = compute!(Field(τuuₜˢ * du_dx + τvvₜˢ * dv_dy +
                             τuvₜˢ * du_dy + τvuₜˢ * dv_dx))
        Pᵀᵥ = compute!(Field(τuwˢ * du_dz + τvwˢ * dv_dz +
                             τwwˢ * dw_dz + τwuₜˢ * dw_dx + τwvₜˢ * dw_dy))
        Pᵀ = compute!(Field(Pᵀₕ + Pᵀᵥ))
        @info "Transfer term Pᵀ done at $(time() - t0)s"
        CUDA.pool_status()

        wˢbˢ = compute!(Field(Bˢ*wˢ))
        @info "Transfer term wˢbˢ done at $(time() - t0)s"
        CUDA.pool_status()
        
        g2c(gpu_field) = on_architecture(CPU(), gpu_field)
        yc = 0.5 * (yT[640] + yT[641])
        itp2nodes(field) = Array([interpolate((x, yc, z), field) for x in xi, z in zT'])
        τ̅uu = mean(τuu, dims=2)
        τ̅ww = mean(τww, dims=2)
        τ̅vv = mean(τvv, dims=2)
        τ̅wb = mean(τwb, dims=2)
        τuu = itp2nodes(τuu)
        τww = itp2nodes(τww)
        τvv = itp2nodes(τvv)
        τwb = itp2nodes(τwb)
        u̅ʸ = mean(itp_u̅, dims=2)
        v̅ʸ = mean(itp_v̅, dims=2)
        w̅ʸ = mean(itp_w̅, dims=2)
        B̅ʸ = mean(itp_B̅, dims=2)
        u̅ = itp2nodes(u̅)
        v̅ = itp2nodes(v̅)
        w̅ = itp2nodes(w̅)
        B̅ = itp2nodes(B̅)
        println(τ̅uu)
      #  τ̅uu = itp2nodes(τ̅uu)
      #  τ̅ww = itp2nodes(τ̅ww)
      #  τ̅vv = itp2nodes(τ̅vv)
      #  τ̅wb = itp2nodes(τ̅wb)
      #  u̅ʸ = itp2nodes(u̅ʸ)
      #  v̅ʸ = itp2nodes(v̅ʸ)
      #  w̅ʸ = itp2nodes(w̅ʸ)
      #  B̅ʸ = itp2nodes(B̅ʸ)

        return (u̅,u̅ʸ), (v̅,v̅ʸ), (w̅,w̅ʸ), (B̅,B̅ʸ), g2c(uˢ), g2c(vˢ), g2c(wˢ), g2c(Bˢ), (τuu,τ̅uu), (τvv,τ̅vv), (τww,τ̅ww), (τwb,τ̅wb), g2c(Πₕ), g2c(Πᵥ), g2c(Πᵥg), g2c(Pᵃ), g2c(Pᵃᵥg), g2c(Pˢ), g2c(Pˢᵥg), g2c(Pᵀ), g2c(wˢbˢ)
    end
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
function MLaverage(snapshots, i, v; kernel=:tophat, scale=20kilometer)
    _, _, z = nodes(v)
    h = MLD(snapshots,i; threshold = 0.03)
    grid = v.grid
    H = Field{Center, Center, Nothing}(grid)
    coarse_graining!(h, H; kernel, cutoff = scale)
    ψ = Field{Center, Center, Nothing}(grid)
    arch = architecture(grid)
    launch!(arch, grid, :xy, _zMLaverage!, ψ, v, H, z, grid)
    return ψ
end

""" mixed layer instability """
function MLI(snapshots, i; kernel=:tophat, scale=20kilometer)
    α = parameters.α
    g = parameters.g
    f = parameters.f
    Ti = snapshots[:T][i]
    h = MLD(snapshots,i; threshold = 0.03)

    H = Field{Center, Center, Nothing}(Ti.grid)
    T = CenterField(Ti.grid)
    coarse_graining!(h, H; kernel, cutoff = scale)
    coarse_graining!(Ti, T; kernel, cutoff = scale)

    ∇b = compute!(Field(α * g * (∂x(T)^2 + ∂y(T)^2)^0.5))

    μ = CenterField(Ti.grid)
    _,_,z = nodes(Ti)
    data = 2 * reshape(z, (1, 1, :)) ./ interior(H, :,:,1)
    set!(μ, (1 .- (data .+ 1).^2) .* (1 .+ 5/21 * (data .+ 1).^2))
    fill_halo_regions!(μ)

    return MLaverage(snapshots,i,compute!(Field(μ * ∇b^2)); kernel, scale) * H^2 / f
end

# Define instability category constants as per your request
const STABLE = 0
const I_SI   = 1  # Inertial/symmetric instability
const SI     = 2  # Symmetric instability
const SI_G   = 3  # Symmetric/gravitational instability
const G      = 4  # Gravitational Instability

"""
    _classify_instability_kernel!(categories, f, q, ωz, Ri)

An Oceananigans.jl kernel function to classify instability at each grid point (i, j, k).

Arguments:
- `categories`: The output 3D field where integer category codes will be stored.
- `f`: Coriolis frequency (scalar).
- `q`: Potential vorticity (3D field).
- `ωz`: Vertical component of vorticity (3D field).
- `Ri`: Richardson number (3D field).
"""
@kernel function _classify_instability_kernel!(categories, f, q, ωz, Ri)
    i, j, k = @index(Global, NTuple)

    # Default to Stable
    category_value = STABLE

    ϕRi = atan(-1/Ri[i, j, k])
    ϕRo = atan(-ωz[i, j, k]/f)
    if q[i, j, k] < 0  # Check for instability prerequisite
        # Anticyclonic Vorticity Condition: ωz < f
        if ωz[i, j, k] < f
            if -π/4 < ϕRi && ϕRi <= ϕRo
                category_value = I_SI
            elseif -π/2 < ϕRi && ϕRi <= -π/4
                category_value = SI
            elseif -3π/4 < ϕRi && ϕRi <= -π/2
                category_value = SI_G
            elseif -π <= ϕRi && ϕRi <= -3π/4
                category_value = G
            elseif ϕRo < ϕRi && ϕRi <= 0
                category_value = STABLE # Stable case within q < 0
            end
        # Cyclonic Vorticity Condition: ωz > f
        elseif ωz[i, j, k] > f
            # Note: No I/SI for cyclonic case in the table
            if -π/2 < ϕRi && ϕRi <= ϕRo # SI condition uses ϕRo here
                category_value = SI
            elseif -3π/4 < ϕRi && ϕRi <= -π/2
                category_value = SI_G
            elseif -π <= ϕRi && ϕRi <= -3π/4
                category_value = G
            elseif ϕRo < ϕRi && ϕRi <= 0
                category_value = STABLE # Stable case within q < 0
            end
        # If ωz[i, j, k] == f, it remains STABLE (0) if q < 0 but no other conditions met,
        # as category_value was initialized to STABLE.
        end
    else # q[i, j, k] >= 0
        category_value = STABLE # Stable if potential vorticity is not negative
    end
    categories[i, j, k] = category_value
end

"""
    classify_instability(f, q, ωz, Ri)

Returns a field of instability categories based on input oceanographic fields.

Arguments:
- `f`: Coriolis frequency (scalar, s^{-1}).
- `q`: Potential vorticity (Oceananigans.jl field, s^{-3}).
- `ωz`: Vertical component of relative vorticity (Oceananigans.jl field, s^{-1}).
         The table specifies `ωz > f` for cyclonic and `ωz < f` for anticyclonic.
- `Ri`: the Richardson number (Oceananigans.jl field).

Returns:
- `categories`: An Oceananigans.jl `Field` containing integer codes:
    - 0: Stable (S)
    - 1: Inertial/symmetric instability (I/SI)
    - 2: Symmetric instability (SI)
    - 3: Symmetric/gravitational instability (SI/G)
    - 4: Gravitational Instability (G)

Assumes `q`, `ωz`, `Ri` are all defined on the same `grid` and at the same cell locations.
"""
function classify_instability(f, q, ωz, Ri)
    
    grid = q.grid
    arch = architecture(grid)
    # Assumes q, ωz, ϕRi, ϕRo are all at the same location (e.g., Center, Center, Center)
    # Create an output field for the categories. It must store integers.
    # Initialize with STABLE (0). The kernel also sets a default for each point.
    categories = Field{Center, Center, Center}(grid)
    # fill!(categories, STABLE) # Optional: kernel initializes each point anyway

    # Launch the kernel function to populate the categories field
    launch!(arch, grid, :xyz, 
            _classify_instability_kernel!, 
            categories, f, q, ωz, Ri)

    return categories
end

# --- Helper function to calculate categories for 2D averaged data ---
# This function replicates the logic from `_classify_instability_kernel!`
# but operates on 2D arrays directly.
function calculate_instability_categories_2d(f, q_2d, ωz_2d, Ri_2d)
    Nx, Nz = size(q_2d)
    categories_2d = zeros(Int, Nx, Nz) # Initialize with STABLE

    for k_idx in 1:Nz # z-dimension
        for i_idx in 1:Nx # x-dimension
            q_val = q_2d[i_idx, k_idx]
            ωz_val = ωz_2d[i_idx, k_idx]
            ϕRi_val = atan(-1/Ri_2d[i_idx, k_idx])
            ϕRo_val = atan(-ωz_2d[i_idx, k_idx]/f)
            
            category_value = STABLE # Default for current point

            if q_val < 0
                if ωz_val < f # Anticyclonic
                    if -π/4 < ϕRi_val && ϕRi_val <= ϕRo_val
                        category_value = I_SI
                    elseif -π/2 < ϕRi_val && ϕRi_val <= -π/4
                        category_value = SI
                    elseif -3π/4 < ϕRi_val && ϕRi_val <= -π/2
                        category_value = SI_G
                    elseif -π <= ϕRi_val && ϕRi_val <= -3π/4
                        category_value = G
                    elseif ϕRo_val < ϕRi_val && ϕRi_val <= 0 
                        category_value = STABLE
                    # If none of the above, it remains STABLE due to initialization if q_val < 0 and no specific instability or the explicit stable condition is met.
                    # However, for clarity and to ensure it's STABLE if no other instability type within q<0, ωz<f is found:
                    # else category_value = STABLE; (already set, but good to keep in mind)
                    end
                elseif ωz_val > f # Cyclonic
                    if -π/2 < ϕRi_val && ϕRi_val <= ϕRo_val
                        category_value = SI
                    elseif -3π/4 < ϕRi_val && ϕRi_val <= -π/2
                        category_value = SI_G
                    elseif -π <= ϕRi_val && ϕRi_val <= -3π/4
                        category_value = G
                    elseif ϕRo_val < ϕRi_val && ϕRi_val <= 0
                        category_value = STABLE
                    # else category_value = STABLE; (as above)
                    end
                # If ωz_val == f, it remains STABLE (as initialized)
                end
            else # q_val >= 0
                category_value = STABLE
            end
            categories_2d[i_idx, k_idx] = category_value
        end
    end
    return categories_2d
end
