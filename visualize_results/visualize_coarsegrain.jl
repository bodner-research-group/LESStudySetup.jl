using LESStudySetup
using CairoMakie
using SixelTerm
using Printf, Dates
using Oceananigans: compute!
using Oceananigans.Grids: xnodes, ynodes, znodes
using Oceananigans.Units: kilometer
using Statistics: mean
using LESStudySetup.Diagnostics
using LESStudySetup.Diagnostics: load_snapshots, symmetric_filtering 
using LESStudySetup.Diagnostics: coarse_grained_fluxes, δ
set_theme!(Theme(fontsize = 12))

cooling, wind, dTf,a = 50, 0.1, -1,1.0

# Examples! (fill in the correct filename and metadata filename)
# cooling, wind, dTf = 25, 0.02, -1
cooling = @sprintf("%03d", cooling)
wind = replace("$(wind)","." => "" )
a = replace("$(a)","." => "" )
if dTf < 0
    fileparams = "hydrostatic_twin_simulation"
else
    if length(wind) < 2
        wind = "0" * wind
    end
    dTf = @sprintf("%1d", dTf)
    fileparams = "free_surface_short_test_$(cooling)_wind_$(wind)_dTf_$(dTf)_a_$(a)"
end
filehead = "./"
filename = filehead * "hydrostatic_snapshots_" * fileparams * ".jld2"
metadata = filehead * "experiment_" * fileparams * "_metadata.jld2"
filesave = "./results/"
initfile = filehead * "hydrostatic_snapshots_hydrostatic_background.jld2"

# load all the data!!
println("Loading data from $filename...")
snapshots = load_snapshots(filename; metadata)
initsnaps = load_snapshots(initfile)
u0 = initsnaps[:u][1]
v0 = initsnaps[:v][1]

# Let's pick the last snapshot!
times = snapshots[:T].times
snapshot_number = length(times)÷2 + 48
nday = @sprintf("%2.0f", (times[snapshot_number])/60^2/24)
println("Plotting snapshot $snapshot_number on day $(nday)...")

### Plot the coarse-grained cross-scale fluxes 
u̅l, v̅l, w̅l, Πhl, Πvl = coarse_grained_fluxes(snapshots, snapshot_number; u0, v0, cutoff = 4kilometer)
ωz = ζ(snapshots, snapshot_number; u0, v0)
δz = δ(snapshots, snapshot_number; u0, v0)
T = snapshots[:T][snapshot_number]
x, y, _ = nodes(T)
α = parameters.α
g = parameters.g
f = parameters.f
∇b = α * g * (∂x(T)^2 + ∂y(T)^2 + ∂z(T)^2)^0.5

Πl = compute!(Field(Πhl))
ωz = compute!(Field(ωz))
δz = compute!(Field(δz))
∇b = compute!(Field(∇b))

# Plotting!
fig = Figure(size = (900, 750))
axis_kwargs1 = (xlabel = "x (km)", ylabel = "y (km)",
            limits = ((0, 100), (0, 100)))
axis_kwargs2 = NamedTuple{(:xlabel,:limits)}(axis_kwargs1)
axis_kwargs3 = NamedTuple{(:ylabel,:limits)}(axis_kwargs1)
ax_Π = Axis(fig[1, 1][1,1]; title=L"\Pi_h^4", axis_kwargs3...)
ax_Ro = Axis(fig[1, 2][1,1]; title=L"\zeta/f", limits = ((0, 100), (0, 100)))
ax_δ = Axis(fig[2, 1][1,1]; title=L"\delta/f", axis_kwargs1...)
ax_b = Axis(fig[2, 2][1,1]; title=L"|\nabla b|", axis_kwargs2...)

k = 123
hm_Π = heatmap!(ax_Π, 1e-3x, 1e-3y, interior(Πl,:,:,k); rasterize = true, colormap = :balance, colorrange = (-1e-6, 1e-6))
hm_Ro = heatmap!(ax_Ro, 1e-3x, 1e-3y, interior(ωz,:,:,k)/f; rasterize = true, colormap = :balance, colorrange = (-2, 2))
hm_δ = heatmap!(ax_δ, 1e-3x, 1e-3y, interior(δz,:,:,k)/f; rasterize = true, colormap = :balance, colorrange = (-1, 1))
hm_b = heatmap!(ax_b, 1e-3x, 1e-3y, interior(∇b,:,:,k); rasterize = true, colormap = :binary, colorrange = (0, 1e-5))

Colorbar(fig[1, 1][1, 2], hm_Π)
Colorbar(fig[1, 2][1, 2], hm_Ro)
Colorbar(fig[2, 1][1, 2], hm_δ)
Colorbar(fig[2, 2][1, 2], hm_b)

hidexdecorations!(ax_Π, ticks = false)
hidexdecorations!(ax_Ro, ticks = false)
hideydecorations!(ax_Ro, ticks = false)
hideydecorations!(ax_b, ticks = false)

save(filesave * "coarse_grained_" * fileparams * "_d$(nday).pdf", fig)