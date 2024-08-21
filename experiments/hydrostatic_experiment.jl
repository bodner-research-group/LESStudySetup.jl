using LESStudySetup
using LESStudySetup.Oceananigans.Units
using LESStudySetup.Oceananigans.Utils: ConsecutiveIterations
using JLD2 

# Architecture (CPU, GPU, or Distributed)
architecture = GPU()

# Setting some initial values (Q = heat flux in W/m², Δz = vertical spacing)
set_value!(Δh = 250, Δz = 2)

function run_experiment(experiment; 
                        Q   = 0.0,  # Cooling heat flux in W/m²
                        τw  = 0.0,  # Wind stress in N/m²
                        θ   = 30.0, # Wind stress angle in degrees (0 correspond to zonal wind stress)
                        Δh  = 250,  # Horizontal resolution [m]
                        ΔTᵉ = 0.5,  # Eddy temperature difference
                        ΔTᶠ = 0.5,  # Meridional temperature difference
                        Lf  = 0.9,  # Size of temperature front (large numbers correspond to steeper fronts)
                        σ²  = 0.15, # Initial spread of the barotropic eddy
                        N²  = 2e-6, # Initial stratification below the thermocline
                        output_frequency = 3hours,
                        checkpoint_frequency = 3hours,
                        stop_time = 20days,
                        restoring = false,
                        restart_file = nothing)
    
    set_value!(; Q, τw, θ, ΔTᵉ, ΔTᶠ, Lf, N², σ², Δh)

    @info "Simulation parameters: " parameters

    # Let's start with an hydrostatic setup running for 20 days
    simulation = idealized_setup(architecture; stop_time, hydrostatic_approximation = true, restoring)

    jldsave("experiment_$(experiment)_metadata.jld2", parameters = parameters)

    # Show the configuration of the simulation
    @info simulation

    # Let's attach some outputs
    model         = simulation.model
    output_fields = merge(model.velocities, model.tracers)

    simulation.output_writers[:checkpoint] = Chekpointer(model;
                                                         schedule = TimeInterval(checkpoint_frequency),
                                                         prefix = "hydrostatic_",
                                                         overwrite_existing = true)

    simulation.output_writers[:snapshots] = JLD2OutputWriter(model, output_fields;
                                                             schedule = ConsecutiveIterations(TimeInterval(output_frequency)),
                                                             overwrite_existing = true,
                                                             array_type = Array{Float32},
                                                             with_halos = true,
                                                             filename = "hydrostatic_snapshots_$(experiment)")


    simulation.output_writers[:free_surface] = JLD2OutputWriter(model, (; η = model.free_surface.η);
                                                                schedule = ConsecutiveIterations(TimeInterval(output_frequency)),
                                                                overwrite_existing = true,
                                                                array_type = Array{Float32},
                                                                with_halos = true,
                                                                filename = "hydrostatic_free_surface_$(experiment)")

    #####
    ##### Let's run!!!!
    #####

    run!(simulation; pickup = restart_file)

    return simulation
end
