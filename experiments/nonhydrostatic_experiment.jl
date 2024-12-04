using Preferences
@debug "Preloading GTL library" iscray
import Libdl
Libdl.dlopen_e("libmpi_gtl_cuda", Libdl.RTLD_LAZY | Libdl.RTLD_GLOBAL)

using MPI
MPI.Init()
using LESStudySetup
using LESStudySetup.Oceananigans.Units
using LESStudySetup.Oceananigans.Utils: ConsecutiveIterations
using LESStudySetup.Oceananigans.OutputWriters: Checkpointer
using JLD2

proc_max = 32 # Maximum number of GPUs (32 * 32 GPUs == 1024 GPUs) for the large experiment
proc = 32 # Actually used number of GPUs in each direction

scaling = Int(proc_max / proc)

# Architecture (CPU, GPU, or Distributed)
arch = Distributed(GPU(), partition = Partition(proc, proc))

# Domain size is 100km x 100km x 250m, 
# the mesh size 100000 / Δh, 100000 / Δh and 252 / Δz

# Grid size
Δh = 4.8828125 * scaling
Δz = 1.125 

LESStudySetup.default_experimental_setup!(; Δh, Δz)
                  
# Show all the parameters we are using
if arch.local_rank == 0
    @info "Simulation parameters: " parameters
end

# Output writer details
output_frequency = 10minutes
checkpoint_frequency = 2hours
stop_time = 10days

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
							 array_type = Array{Float32},
                                                         filename = "nonhydrostatic_snapshots_$(arch.local_rank)")

simulation.output_writers[:checkpoint] = Checkpointer(model;
                                                      schedule = TimeInterval(checkpoint_frequency),
                                                      prefix = "nonhydrostatic_checkpoint_$(arch.local_rank)",
                                                      overwrite_existing = true)

#####
##### Let's run!!!!
#####

run!(simulation; pickup = restart_file)
