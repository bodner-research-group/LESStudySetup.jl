# include("hydrostatic_experiment.jl")

# using Statistics: mean

# run_experiment("free_surface_short_test_050_wind_01_dTf_2"; Q = 50.0, τw = 0.1, Δh = 100)
# run_experiment("free_surface_short_test_050_wind_02_dTf_2"; Q = 50.0, τw = 0.2, Δh = 100)
# run_experiment("free_surface_short_test_000_wind_00_dTf_2"; Q = 00.0, τw = 0.0, Δh = 100)
# run_experiment("free_surface_short_test_000_wind_00_dTf_5"; Q = 00.0, τw = 0.0, Δh = 100, ΔTᶠ = 5.0)

# Final configuration to use!!! (twin experiment for the nonhydrostatic one)
#run_experiment("hydrostatic_twin_simulation"; Q = 50.0, τw = 0.1, Δh = 200, ΔTᶠ = 1.0, a = 1, stop_time = 10days)

include("initial_turbulence_generator.jl")
generate_initial_turbulence(;τw  = 0.1, ΔTᵉ = 0.5, output_frequency = 10minutes, stop_time = 10hours)
