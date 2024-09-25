using LESStudySetup
using CairoMakie
using SixelTerm
using Printf, Dates
using Oceananigans: compute!
using Oceananigans.Grids: xnodes, ynodes, znodes
using Statistics: mean
using LESStudySetup.Diagnostics
using LESStudySetup.Diagnostics: load_snapshots,ζ
set_theme!(Theme(fontsize = 12))

function visualize(cooling, wind, dTf, a)
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

    # load all the data!!
    println("Loading data from $filename...")
    snapshots = load_snapshots(filename; metadata)

    # Let's pick the last snapshot!
    times = snapshots[:T].times
    
    t0 = now()

    # Coordinate arrays
    xT, yT, _ = nodes(snapshots[:T][end])

    # Plot the fields
    k = 113 # index of the vertical level to plot
    vbnd = 7
    f = parameters.f

    fig = Figure(size = (500, 400))
    axis_kwargs = (xlabel = "x (km)", ylabel = "y (km)",
                limits = ((0, 100), (0, 100)), aspect = 1)  

    n = Observable(1)
    ax_v = Axis(fig[1, 1][1,1]; 
                title = @lift("t = " * string(round(times[$n]/3600/24, digits=3)) * " days"),
                subtitle="ζ/f, z=-25m", axis_kwargs...)  
    v = @lift interior(compute!(Field(ζ(snapshots, $n)/f)),:,:,k)

    hm_v = heatmap!(ax_v, 1e-3xT, 1e-3yT, v; rasterize = true, colormap = :balance, colorrange = (-vbnd, vbnd))
    Colorbar(fig[1, 1][1, 2], hm_v)
    #Label(fig[1, 1][1,1:2], title, fontsize=15, tellwidth=false)
    frames = 1:2:length(times)
    @info "Making a neat animation of Ro..."
    record(fig, filesave * "ro_" * fileparams * ".mp4", frames, framerate=4) do i
        println("Loading fields $i wall time: $((now() - t0).value/1e3) seconds.")
        n[] = i
    end

    return
end

coolings = [50]
winds = [0.1]
dTfs = [-1]
as = [1.0]
i = 1
visualize(coolings[i], winds[i], dTfs[i], as[i])