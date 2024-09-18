using Preferences
@debug "Preloading GTL library" iscray
import Libdl
Libdl.dlopen_e("libmpi_gtl_cuda", Libdl.RTLD_LAZY | Libdl.RTLD_GLOBAL)

using MPI
MPI.Init()
using LESStudySetup
using LESStudySetup.Oceananigans.Units
using JLD2

# Architecture (CPU, GPU, or Distributed)
arch = Distributed(GPU(), partition = Partition(32, 32))

# Domain size is 100km x 100km x 250m, the mesh size 100000 / Δh, 100000 / Δh and 250 / Δz
# Setting some initial values (Q = heat flux in W/m², Δz = vertical spacing)
set_value!(; # Forcing
             Q = 50.0, 
            τw = 0.1, 
             # Spacing -> BEWARE: this leads to a grid which is 50000 × 50000 × 250 in size!!!!
            Δh = 4.01875, 
            Δz = 1.28, 
             # Initial condition
           ΔTᶠ = 1.0, 
           ΔTᵉ = 0.5, 
             Φ = 0.025, 
             a = 1, 
            σ² = 0.15, 
           Lf  = 0.9)

# Show all the parameters we are using
@info "Simulation parameters: " parameters

# Output writer details
output_frequency = 10minutes
checkpoint_frequency = 10minutes
stop_time = 1hour

background_forcing = true
restart_file = false

# Let's start with an nonhydrostatic setup running for 30 days
simulation = idealized_setup(arch; stop_time, background_forcing)

if arch.local_rank == 0
    jldsave("nonhydrostatic_experiment_metadata.jld2", parameters = parameters)
end

# Show the configuration of the simulation
@info simulation

# Let's attach some outputs
model         = simulation.model
output_fields = merge(model.velocities, model.tracers)

simulation.output_writers[:snapshots] = JLD2OutputWriter(model, output_fields;
                                                         schedule = ConsecutiveIterations(TimeInterval(output_frequency)),
                                                         overwrite_existing = true,
                                                         filename = "nonhydrostatic_snapshots_$(arch.local_rank)")

simulation.output_writers[:checkpoint] = Checkpointer(model;
                                                        schedule = TimeInterval(checkpoint_frequency),
                                                        prefix = "nonhydrostatic_checkpoint_$(arch.local_rank)",
                                                        overwrite_existing = true,
                                                        cleanup = true)

#####
##### Let's run!!!!
#####

run!(simulation; pickup = restart_file)
