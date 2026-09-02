using LESStudySetup
using LESStudySetup.Oceananigans.Units
using LESStudySetup.Oceananigans.Utils: ConsecutiveIterations
using LESStudySetup.Oceananigans.OutputWriters: Checkpointer
using JLD2 

# Architecture (CPU, GPU, or Distributed)
architecture = GPU()

function generate_initial_turbulence(;τw  = 0.1,  # Wind stress in N/m²
                                     θ   = 30.0, # Wind stress angle in degrees (0 correspond to zonal wind stress)
                                     Δh  = 4.8828125,    # Horizontal resolution [m]
                                     a   = 0,    # Eddy temperature amplitude
                                     Δz  = 1.125,    # Vertical resolution [m]
                                     ΔTᵉ = 0,    # Eddy temperature difference
                                     Lz  = 252,
                                     Q   = 0,
                                     m₀ = 60,
                                     T₀ = 20,
                                     output_frequency = 1hours,
                                     checkpoint_frequency = 10hours,
                                     stop_time = 10hours)
    
    set_value!(; τw, θ, ΔTᵉ, a, Δz, Δh, Lz, Q, m₀, T₀)

    # Let's start with an hydrostatic setup running for 20 days
    simulation = turbulence_generator_setup(architecture; stop_time)
    
    @info "Simulation parameters: " parameters

    filehead = "/orcd/data/abodner/002/shirui/LESStudySetup.jl/"

    jldsave(filehead * "turbulence_generator_metadata_10h.jld2", parameters = parameters)

    # Show the configuration of the simulation
    @info simulation

    # Let's attach some outputs
    model         = simulation.model
    output_fields = merge(model.velocities, model.tracers)

    simulation.output_writers[:checkpoint] = Checkpointer(model;
                                                         schedule = TimeInterval(checkpoint_frequency),
                                                         prefix = filehead * "turbulence_generator_checkpoint_10h",
                                                         overwrite_existing = true)

    simulation.output_writers[:snapshots] = JLD2OutputWriter(model, output_fields;
                                                             schedule = TimeInterval(output_frequency),
                                                             overwrite_existing = true,
                                                             array_type = Array{Float32},
                                                             with_halos = true,
                                                             filename = filehead * "turbulence_generator_output_10h")

    #####
    ##### Let's run!!!!
    #####

    run!(simulation)

    return simulation
end
