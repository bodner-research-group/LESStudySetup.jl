# include("hydrostatic_experiment.jl")

# using Statistics: mean

# run_experiment("free_surface_short_test_050_wind_01_dTf_2"; Q = 50.0, τw = 0.1, Δh = 100)
# run_experiment("free_surface_short_test_050_wind_02_dTf_2"; Q = 50.0, τw = 0.2, Δh = 100)
# run_experiment("free_surface_short_test_000_wind_00_dTf_2"; Q = 00.0, τw = 0.0, Δh = 100)
# run_experiment("free_surface_short_test_000_wind_00_dTf_5"; Q = 00.0, τw = 0.0, Δh = 100, ΔTᶠ = 5.0)

# Final configuration to use!!! (twin experiment for the nonhydrostatic one)
#run_experiment("hydrostatic_twin_simulation"; Q = 50.0, τw = 0.1, Δh = 200, ΔTᶠ = 1.0, a = 1, stop_time = 10days)
proc_max = 32 # Maximum number of GPUs (32 * 32 GPUs == 1024 GPUs) for the large experiment
proc = 2 # Actually used number of GPUs in each direction

scaling = Int(proc_max / proc)

# Grid size
Δh = 4.8828125 * scaling
Δz = 1.125 

run_experiment("hydrostatic_twin_simulation"; Δh, Δz)

#include("initial_turbulence_generator.jl")
#generate_initial_turbulence(;τw  = 0.1, ΔTᵉ = 0.5, output_frequency = 10minutes, stop_time = 10hours)