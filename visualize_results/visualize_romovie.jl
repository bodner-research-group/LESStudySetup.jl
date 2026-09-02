using LESStudySetup
using CairoMakie
using SixelTerm
using Printf, Dates
using Oceananigans: compute!
using Oceananigans.Grids: xnodes, ynodes, znodes
using Statistics: mean
using LESStudySetup.Diagnostics
using LESStudySetup.Diagnostics: load_snapshots,ζ
set_theme!(theme_latexfonts(), fontsize=12,figure_padding = 10)
shift(x) = [x[3size(x,1)÷4+1:end, :]; x[1:3size(x,1)÷4, :]]

# Examples! (fill in the correct filename and metadata filename)
fileparams = "hydrostatic_twin_simulation"
filehead = "/orcd/data/abodner/002/shirui/LESStudySetup.jl/"
filename = filehead * "hydrostatic_snapshots_" * fileparams * ".jld2"
metadata = filehead * "experiment_" * fileparams * "_metadata.jld2"
filesave = filehead * "results/"

# load all the data!!
println("Loading data from $filename...")
snapshots = load_snapshots(filename; metadata)

# Let's pick the last snapshot!
times = snapshots[:T].times

t0 = now()

# Coordinate arrays
xT, yT, _ = nodes(snapshots[:T][end])

# Plot the fields
k = 202 # index of the vertical level to plot
vbnd = 3
f = parameters.f

fig = Figure(size = (500, 450))
gl = fig[1, 1] = GridLayout()
axis_kwargs = (xlabel = "x (km)", ylabel = "y (km)",
            limits = ((0, 100), (0, 100)), aspect = 1)  

n = Observable(1)
ax_v = Axis(gl[1,1]; 
            title = @lift("t = " * string(round(times[$n]/3600/24, digits=3)) * " days"),
            subtitle="ζ/f, z=-25.3m", axis_kwargs...)  
v = @lift interior(compute!(Field(ζ(snapshots, $n)/f)),:,:,k)

hm_v = heatmap!(ax_v, 1e-3xT, 1e-3yT, v; rasterize = true, colormap = :balance, colorrange = (-vbnd, vbnd))
Colorbar(gl[1, 2], hm_v)
colgap!(gl, 1, 1)
frames = 1:2:length(times)
@info "Making a neat animation of Ro..."
record(fig, filesave * "ro_" * fileparams * ".mp4", frames, framerate=4) do i
    println("Loading fields $i wall time: $((now() - t0).value/1e3) seconds.")
    n[] = i
end
