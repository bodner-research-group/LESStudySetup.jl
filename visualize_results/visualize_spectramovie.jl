using LESStudySetup
using CairoMakie
using SixelTerm
using Printf, Dates
using Oceananigans: compute!
using Oceananigans.Grids: xnodes, ynodes, znodes
using LESStudySetup.Diagnostics
using LESStudySetup.Diagnostics: load_snapshots, isotropic_powerspectrum
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
    snapshot_number = length(times)
    nday = @sprintf("%02.0f", (times[snapshot_number])/60^2/24)
    println("Plotting snapshot $snapshot_number on day $(nday)...")
    t0 = now()

    T = snapshots[:T][snapshot_number]
    u = snapshots[:u][snapshot_number]
    v = snapshots[:v][snapshot_number]
    w = snapshots[:w][snapshot_number]
    println("Loading fields wall time: $((now() - t0).value/1e3) seconds.")

    # Coordinate arrays
    xu, yu, zu = nodes(u)
    xv, yv, zv = nodes(v)
    xw, yw, zw = nodes(w)
    xT, yT, zT = nodes(T)

    # Compute the horizontal spectrum of T, u, v, w 
    axis_kwargs = (xlabel = "Wavenumber (rad⋅m⁻¹)", 
                   ylabel = L"E_i(k)/E_{T,v}(k_{min},z=-3~m)",
                   xscale = log10, yscale = log10,
                   limits = ((8e-5, 4e-2), (1e-11,7)))
    fig = Figure(size = (400, 300))
    klev = 113
    println("Plotting spectra at z = $(zT[klev])m...")
    ax = Axis(fig[1, 1][1,1]; 
              title = @lift("t = " * string(round(times[$n]/3600/24, digits=2)) * " days"),
              subtitle="z = $(zT[klev])m", axis_kwargs...) 

    Su = @lift isotropic_powerspectrum(interior(snapshots[:u][$n], :, :, klev), interior(snapshots[:u][$n], :, :, klev), xu, yu)
    Sv = @lift isotropic_powerspectrum(interior(snapshots[:v][$n], :, :, klev), interior(snapshots[:v][$n], :, :, klev), xv, yv)
    wk = @lift (interior(snapshots[:w][$n], :, :, klev)+interior(snapshots[:w][$n], :, :, klev+1))/2
    Sw = @lift isotropic_powerspectrum($wk, $wk, xw, yw)
    St = @lift isotropic_powerspectrum(interior(snapshots[:T][$n], :, :, klev), interior(snapshots[:T][$n], :, :, klev), xT, yT)

    lines!(ax, @lift($Su.freq), 1e-8@lift($Su.freq).^-2, linestyle = :dash, color = :black)
    text!(ax, 1e-3, 1e-2; text = L"k^{-2}")
    lines!(ax, @lift($Su.freq), 1e-15@lift($Su.freq).^-3, linestyle = :dash, color = :gray)
    text!(ax, 10^-3.5, 1e-6; text = L"k^{-3}")
    lines!(ax, @lift($St.freq), @lift(Real.($St.spec./$St.spec[1])), color = :red, label = L"E_T")
    lines!(ax, @lift($Su.freq), @lift(Real.($Su.spec./$Sv.spec[1])), color = :blue, label = L"E_u")
    lines!(ax, @lift($Sv.freq), @lift(Real.($Sv.spec./$Sv.spec[1])), color = :green, label = L"E_v")
    lines!(ax, @lift($Sw.freq), @lift(Real.($Sw.spec./$Sv.spec[1])), color = :black, label = L"E_w")
    axislegend(ax, labelsize=10, patchsize = (20, 5))

    frames = 1:2:length(times)
    @info "Making a neat animation of spectra..."
    record(fig, filesave * "spectra_" * fileparams * ".mp4", frames, framerate=4) do i
        println("Loading fields $i day $(round(times[i]/3600/24, digits=2)) wall time: $((now() - t0).value/1e3) seconds.")
        n[] = i
    end

    return
end

coolings = [0]
winds = [0]
dTfs = [4]
as = [1.0]
for i in 1:length(coolings)
    visualize(coolings[i], winds[i], dTfs[i], as[i])
end