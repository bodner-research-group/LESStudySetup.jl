#visualize_wmovie.jl
using CairoMakie
using SixelTerm
using Printf, Dates
using Oceananigans.Fields: interior
using Oceananigans: compute!, Field
using Oceananigans.Grids: nodes
using Statistics: mean
# LESStudySetup.Diagnostics
using LESStudySetup.Diagnostics: load_snapshots
set_theme!(Theme(fontsize = 12))

function visualize(cooling, wind, dTf)
    # Examples! (fill in the correct filename and metadata filename)
    # cooling, wind, dTf = 25, 0.02, -1
    cooling = @sprintf("%03d", cooling)
    wind = replace("$(wind)","." => "" )
    if dTf < 0
        fileparams = "four_vortices_cooling_$(cooling)_wind_$(wind)"
    else
        if length(wind) < 2
            wind = "0" * wind
        end
        fileparams = "free_surface_short_test_$(cooling)_wind_$(wind)_dTf_$(dTf)"
    end
    filehead = "/orcd/data/abodner/001/shirui/LESStudySetup.jl/experiments/"
    filename = filehead * "hydrostatic_snapshots_" * fileparams * ".jld2"
    metadata = filehead * "experiment_" * fileparams * "_metadata.jld2"
    filesave = "/orcd/data/abodner/001/shirui/LESStudySetup.jl/results/"

    # load all the data!!
    println("Loading data from $filename...")
    snapshots = load_snapshots(filename; metadata)

    # Let's pick the last snapshot!
    times = snapshots[:T].times
    
    t0 = now()

    # Coordinate arrays
    xw, yw, zw = nodes(snapshots[:w][end])

    # Plot the fields
    kw = 113 # index of the vertical level to plot
    wbnd = 0.005

    fig = Figure(size = (500, 400))
    axis_kwargs = (xlabel = "x (km)", ylabel = "y (km)",
                limits = ((0, 100), (0, 100)), aspect = 1)  

    n = Observable(1)
    ax_w = Axis(fig[1, 1][1,1]; 
                title = @lift("t = " * string(round(times[$n]/3600/24, digits=2)) * " days"),
                subtitle="w (m s⁻¹), z=-26m", axis_kwargs...)  
    w = @lift interior(snapshots[:w][$n],:,:,kw)

    hm_w = heatmap!(ax_w, 1e-3xw, 1e-3yw, w; rasterize = true, colormap = :balance, colorrange = (-wbnd, wbnd))
    Colorbar(fig[1, 1][1, 2], hm_w)
    #Label(fig[1, 1][1,1:2], title, fontsize=15, tellwidth=false)
    frames = 1:2:length(times)
    @info "Making a neat animation of w..."
    record(fig, filesave * "w_" * fileparams * ".mp4", frames, framerate=4) do i
        println("Loading fields $i wall time: $((now() - t0).value/1e3) seconds.")
        n[] = i
    end

    return
end

coolings = [50]
winds = [0.1]
dTfs = [2]
for i in 1:length(coolings)
    visualize(coolings[i], winds[i], dTfs[i])
end