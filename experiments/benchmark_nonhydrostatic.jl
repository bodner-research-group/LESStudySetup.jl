using Preferences
using MPI
MPI.Init()

using LESStudySetup
using LESStudySetup.Oceananigans.Units
using JLD2

# The architecture specifies the partitioning of the domain (i.e., how many GPUs we use)
# We partition only in x and y, there are limitations on the possible configuration:
# -> Nx should be a multiple integer of the partitions in x
# -> Ny should be a multiple integer of the partitions in y and x
# -> Nz should be a multiple integer of the partitions in y
arch = Distributed(GPU(), partition = Partition(x = 2, y = 2))

# The connectivity with ranks for the west, east, north, south, 
# northwest, northeast, southwest and southeast ranks is in 
@show arch.connectivity

# Domain size is 100km x 100km x 250m, the mesh size is
# Nx = Lx / Δh
# Ny = Ly / Δh
# Nz = Lz / Δz
# Set the grid size parameters
set_value!(; Δh = 2, Δz = 1, Lx = 10000, Ly = 10000, Lz = 250)

# Show all the parameters we are using
@info "Simulation parameters: " parameters

# Output writer details
stop_iteration = 1000

background_forcing = true
restart_file = false

# Let's start with an nonhydrostatic setup running for 30 days
simulation = idealized_setup(arch; stop_time, background_forcing, stop_iteration)

# Show the configuration of the simulation
@info simulation

#####
##### Let's run!!!!
#####

run!(simulation; pickup = restart_file)
