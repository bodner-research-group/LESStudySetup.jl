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