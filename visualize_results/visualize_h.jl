using LESStudySetup
using CairoMakie
using SixelTerm
using Printf, Dates
using Oceananigans: compute!
using Oceananigans.Grids: xnodes, ynodes, znodes
using Statistics: mean
using LESStudySetup.Diagnostics
using LESStudySetup.Diagnostics: load_snapshots, MLD 
using JLD2
set_theme!(Theme(fontsize = 12))

# Read in the data
wind,dTf,a = 0.1,1,1.0
wind = replace("$(wind)","." => "" )
a = replace("$(a)","." => "" )

# Plot the results
fig = Figure(size = (900, 300))
ax1 = Axis(fig[1,1]; ylabel = "mean mixed layer depth (m)",xlabel="time (days)",limits = ((0, 20), (50, 150)))
for Q in [50,75,100]
    cooling = @sprintf("%03d", Q)
    filename = "results/h_$(cooling)_wind_$(wind)_dTf_$(dTf)_a_$(a).jld2"
    t = load(filename, "t")/60^2/24
    h = replace(vec(mean(load(filename, "h"), dims=(2,3))), -Inf => NaN)
    lines!(ax1, t, h; label = "cooling Q = $(Q) W/m²")
end
axislegend(ax1; position = :ct, orientation = :horizontal)
save("results/h_wind_$(wind)_dTf_$(dTf)_a_$(a).pdf", fig)