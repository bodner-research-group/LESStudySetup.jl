include("hydrostatic_experiment.jl")

using Statistics: mean

#run_experiment!("free_surface_short_test_050_wind_01_dTf_1_a_08"; Q = 50.0, τw = 0.1, Δh = 200, ΔTᶠ = 1.0, a = 0.8, stop_time = 10days)
#run_experiment!("free_surface_short_test_050_wind_01_dTf_1_a_09"; Q = 50.0, τw = 0.1, Δh = 200, ΔTᶠ = 1.0, a = 0.9, stop_time = 10days)
run_experiment!("free_surface_short_test_020_wind_01_dTf_1_a_10"; Q = 20.0, τw = 0.1, Δh = 200, ΔTᶠ = 1.0, a = 1, stop_time = 10days)
