using LESStudySetup
using LESStudySetup.Oceananigans.Units
using LESStudySetup.Oceananigans.Utils: ConsecutiveIterations
using JLD2 

# Architecture (CPU, GPU, or Distributed)
architecture = GPU()

function generate_initial_turbulence(τw  = 0.0,  # Wind stress in N/m²
                                     θ   = 30.0, # Wind stress angle in degrees (0 correspond to zonal wind stress)
                                     Δh  = 2,    # Horizontal resolution [m]
                                     a   = 1,    # Eddy temperature amplitude
                                     Δz  = 1,    # Vertical resolution [m]
                                     ΔTᵉ = 1,    # Eddy temperature difference
                                     output_frequency = 1hours,
                                     checkpoint_frequency = 10hours,
                                     stop_time = 10hours)
    
    set_value!(; τw, θ, ΔTᵉ, a, Δz, Δh)

    @info "Simulation parameters: " parameters

    # Let's start with an hydrostatic setup running for 20 days
    simulation = turbulence_generator_setup(architecture; stop_time)

    jldsave("turbulence_generator_metadata.jld2", parameters = parameters)

    # Show the configuration of the simulation
    @info simulation

    # Let's attach some outputs
    model         = simulation.model
    output_fields = merge(model.velocities, model.tracers)

    simulation.output_writers[:checkpoint] = Chekpointer(model;
                                                         schedule = TimeInterval(checkpoint_frequency),
                                                         prefix = "turbulence_generator_checkpoint",
                                                         overwrite_existing = true)

    simulation.output_writers[:snapshots] = JLD2OutputWriter(model, output_fields;
                                                             schedule = TimeInterval(output_frequency),
                                                             overwrite_existing = true,
                                                             array_type = Array{Float32},
                                                             with_halos = true,
                                                             filename = "turbulence_generator_output")

    #####
    ##### Let's run!!!!
    #####

    run!(simulation)

    return simulation
end
