using Oceananigans.Utils
using Oceananigans.Grids: node, architecture
using Oceananigans.BoundaryConditions
using Oceananigans.Advection: TracerAdvection
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEMixingLength, CATKEVerticalDiffusivity
using Statistics: mean
using KernelAbstractions: @index, @kernel

model_type(::Val{true})  = HydrostaticFreeSurfaceModel
model_type(::Val{false}) = NonhydrostaticModel

isforced(model::HydrostaticFreeSurfaceModel) = model.advection.momentum isa ForcedAdvection
isforced(model::NonhydrostaticModel) = model.advection isa ForcedAdvection

@kernel function _compute_v_from_continuity!(v, grid, u)
    i, k = @index(Global, NTuple)

    @inbounds v[i, 1, k] = 0
    for j in 2:size(grid, 2) + 1
        @inbounds v[i, j, k] = v[i, j-1, k] - Δyᶜᶜᶜ(i, j-1, k, grid) * ∂xᶜᶜᶜ(i, j-1, k, grid, u)
    end
end

function model_settings(model_type, grid; background_forcing = false)
    
    advection = WENO(; order = 9)

    if background_forcing
        u_background = XFaceField(grid)
        v_background = YFaceField(grid)

        set!(u_background, uᵢ)
        fill_halo_regions!(u_background)

        launch!(architecture(grid), grid, :xz, _compute_v_from_continuity!, v_background, grid, u_background)
        fill_halo_regions!(v_background)

        fill_halo_regions!(u_background)
        fill_halo_regions!(v_background)

        advection = ForcedAdvection(; scheme = advection,
                                      u_background,
                                      v_background)
    end

    if model_type == HydrostaticFreeSurfaceModel # Additional stuff to add if 
        mixing_length = CATKEMixingLength(Cᵇ = 0.01)
        closure = CATKEVerticalDiffusivity(; mixing_length)
        tracers = (:T, :e)

        free_surface = SplitExplicitFreeSurface(grid; substeps = 75, gravitational_acceleration = parameters.g)
        @info "running with $(length(free_surface.settings.substepping.averaging_weights)) substeps"
       
        return (; closure, 
                  tracers, 
                  free_surface, 
                  momentum_advection = advection, 
                  tracer_advection = advection)
    else
        return (; tracers = :T, 
                  timestepper = :RungeKutta3,
                  hydrostatic_pressure_anomaly = CenterField(grid),
                  advection)
    end
end

function progress(sim) 
    u, v, w = sim.model.velocities
    T = sim.model.tracers.T

    ui = interior(u)
    vi = interior(v)
    wi = interior(w)

    Ti = interior(T)

    msg0 = @sprintf("Time: %s, iteration: %d, Δt: %s ", prettytime(sim.model.clock.time), 
                                                        sim.model.clock.iteration,
                                                        prettytime(sim.Δt))
    msg1 = @sprintf("(u, v, w): %.2e %.2e %.2e ", maximum(ui), maximum(vi), maximum(wi))
    msg2 = @sprintf("T: %.2e %.2e ", minimum(Ti), maximum(Ti))

    @info msg0 * msg1 * msg2 
end