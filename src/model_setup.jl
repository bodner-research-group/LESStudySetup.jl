using Oceananigans.Utils
using Oceananigans.Grids: node, architecture
using Oceananigans.DistributedComputations
using Oceananigans.BoundaryConditions
using Oceananigans.Advection: FluxFormAdvection
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEMixingLength, CATKEVerticalDiffusivity
using Statistics: mean
using KernelAbstractions: @index, @kernel
using MPI

model_type(::Val{true})  = HydrostaticFreeSurfaceModel
model_type(::Val{false}) = NonhydrostaticModel

isforced(model::HydrostaticFreeSurfaceModel) = model.advection.momentum isa ForcedAdvection
isforced(model::NonhydrostaticModel) = model.advection isa ForcedAdvection

function compute_v_from_continuity!(v_background, arch, grid, u_background)
    launch!(architecture(grid), grid, :xz, _compute_v_from_continuity!, v_background, grid, u_background)
    fill_halo_regions!(v_background)
    return nothing
end

function compute_v_from_continuity!(v_background, arch::Distributed, grid, u_background)
    ranks = DistributedComputations.ranks(arch)[2]

    for jᵣ in 1:ranks
        fill_halo_regions!(v_background)

        if jᵣ == arch.local_index[2]
            launch!(arch, grid, :xz, _compute_v_from_continuity!, v_background, grid, u_background)
        end

        MPI.Barrier(MPI.COMM_WORLD)
    end

    fill_halo_regions!(v_background)
        
    return nothing
end

@kernel function _compute_v_from_continuity!(v, grid, u)
    i, k = @index(Global, NTuple)

    @inbounds v[i, 1, k] = v[i, 0, k] - Δyᶜᶜᶜ(i, 0, k, grid) * ∂xᶜᶜᶜ(i, 0, k, grid, u)
    for j in 2:size(grid, 2)
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

        compute_v_from_continuity!(v_background, architecture(grid), grid, u_background)

        advection = ForcedAdvection(; scheme = advection,
                                      u_background,
                                      v_background)
    end

    if model_type == HydrostaticFreeSurfaceModel # Additional stuff to add if 
        mixing_length = CATKEMixingLength(Cᵇ = 0.01)
        closure = CATKEVerticalDiffusivity(; mixing_length)
        tracers = (:T, :e)

        free_surface = SplitExplicitFreeSurface(grid; substeps = 75, gravitational_acceleration = parameters.g)
        #@info "running with $(length(free_surface.settings.substepping.averaging_weights)) substeps"
       
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