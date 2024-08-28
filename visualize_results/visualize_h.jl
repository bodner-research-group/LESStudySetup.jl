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

function visualize(cooling, wind, dTf, a)
    # Examples! (fill in the correct filename and metadata filename)
    # cooling, wind, dTf = 25, 0.02, -1
    cooling = @sprintf("%03d", cooling)
    wind = replace("$(wind)","." => "" )
    a = replace("$(a)","." => "" )
    if dTf < 0
        fileparams = "four_vortices_cooling_$(cooling)_wind_$(wind)"
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

    # load all the data!!
    println("Loading data from $filename...")
    snapshots = load_snapshots(filename; metadata)

    # Let's pick the last snapshot!
    times = snapshots[:T].times
    x,y,z = nodes(snapshots[:T][1])
    snapshot_number = length(times)
    nday = @sprintf("%2.0f", (times[snapshot_number])/60^2/24)

    t0 = now()
    
    idx = 1:2:snapshot_number
    hs = zeros(length(idx),length(x),length(y))
    for (i,n) = enumerate(idx)
        hs[i,:,:] = interior(compute!(MLD(snapshots,n; threshold = 0.03)))
        println("Loading fields $n day $(times[n]/60^2/24) wall time: $((now() - t0).value/1e3) seconds.")
        
        # save to file
        jldsave("results/h_$(cooling)_wind_$(wind)_dTf_$(dTf)_a_$(a).jld2"; t=times[idx], x=x, y=y, h=hs)
    end

    return 
end

coolings = [100]
winds = 0.1*ones(1)
dTfs = 1*ones(1)
as = 1.2*ones(1)

for i in 1:length(coolings)
    visualize(coolings[i], winds[i], dTfs[i], as[i])
end

wind,dTf,a = 0.1,1,1.2
wind = replace("$(wind)","." => "" )
a = replace("$(a)","." => "" )
fig = Figure(size = (900, 300))
ax1 = Axis(fig[1,1]; ylabel = "mean mixed layer depth (m)",xlabel="time (days)",limits = ((0, 20), (50, 100)))
for Q in [50,75,100]
    cooling = @sprintf("%03d", Q)
    filename = "results/h_$(cooling)_wind_$(wind)_dTf_$(dTf)_a_$(a).jld2"
    t = load(filename, "t")/60^2/24
    h = replace(vec(mean(load(filename, "h"), dims=(2,3))), -Inf => NaN)
    lines!(ax1, t, h; label = "cooling = $(Q) W/m²")
end
axislegend(ax1; position = :ct, orientation = :horizontal)
save("results/h_wind_$(wind)_dTf_$(dTf)_a_$(a).pdf", fig)

# Read in the data
cooling,wind,dTf,a = 50,0.1,1,1.1
cooling = @sprintf("%03d", cooling)
wind = replace("$(wind)","." => "" )
a = replace("$(a)","." => "" )
filename = "results/h_$(cooling)_wind_$(wind)_dTf_$(dTf)_a_$(a).jld2"
t = load(filename, "t")/60^2/24
x = load(filename, "x")
y = load(filename, "y")
h = load(filename, "h")

# Plot the results
fig = Figure(size = (1200, 900))
ax1 = Axis(fig[1,1]; ylabel = "mixed layer depth (m)",title="Four eddy corners")
ax2 = Axis(fig[2,1]; ylabel = "mixed layer depth (m)",title="Two eddy corners no fronts")
ax3 = Axis(fig[3,1]; ylabel = "mixed layer depth (m)",xlabel = "time (days)",title="Two eddy corners at fronts")
hidexdecorations!(ax1, ticks = false)
hidexdecorations!(ax2, ticks = false)
bounds = 1e3*[12, 38, 62, 88]

xidx1 = (x .<= bounds[1]) .| (x .>= bounds[end])
yidx1 = (y .<= bounds[1]) .| (y .>= bounds[end])
xidx2 = bounds[2] .<= x .<= bounds[3]
yidx2 = bounds[2] .<= y .<= bounds[3]
lines!(ax1, t, vec(mean(h[:,xidx1,yidx1], dims = (2,3))), label = "edge")
lines!(ax1, t, vec(mean(h[:,xidx1,yidx2], dims = (2,3))), label = "mid y")
lines!(ax1, t, vec(mean(h[:,xidx2,yidx1], dims = (2,3))), label = "mid x")
lines!(ax1, t, vec(mean(h[:,xidx2,yidx2], dims = (2,3))), label = "center")

xidx1 = bounds[1] .<= x .<= bounds[2]
xidx2 = bounds[3] .<= x .<= bounds[end]
lines!(ax2, t, vec(mean(h[:,xidx1,yidx1], dims = (2,3))), label = "left x, edge y")
lines!(ax2, t, vec(mean(h[:,xidx1,yidx2], dims = (2,3))), label = "left x, mid y")
lines!(ax2, t, vec(mean(h[:,xidx2,yidx1], dims = (2,3))), label = "right x, edge y")
lines!(ax2, t, vec(mean(h[:,xidx2,yidx2], dims = (2,3))), label = "right x, mid y")

xidx1 = (x .<= bounds[1]) .| (x .>= bounds[end])
xidx2 = bounds[2] .<= x .<= bounds[3]
yidx1 = bounds[1] .<= y .<= bounds[2]
yidx2 = bounds[3] .<= y .<= bounds[end]
lines!(ax3, t, vec(mean(h[:,xidx1,yidx1], dims = (2,3))), label = "edge x, down y")
lines!(ax3, t, vec(mean(h[:,xidx1,yidx2], dims = (2,3))), label = "edge x, up y")
lines!(ax3, t, vec(mean(h[:,xidx2,yidx1], dims = (2,3))), label = "mid x, down y")
lines!(ax3, t, vec(mean(h[:,xidx2,yidx2], dims = (2,3))), label = "mid x, up y")

axislegend(ax1, position = :ct, orientation = :horizontal)
axislegend(ax2, position = :ct, orientation = :horizontal)
axislegend(ax3, position = :ct, orientation = :horizontal)

save("results/h_$(cooling)_wind_$(wind)_dTf_$(dTf)_a_$(a).pdf", fig)