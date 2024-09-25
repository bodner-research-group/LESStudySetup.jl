using LESStudySetup
using CairoMakie
using SixelTerm
using Printf, Dates
using Oceananigans: compute!
using Oceananigans.Grids: xnodes, ynodes, znodes
using LESStudySetup.Diagnostics
using LESStudySetup.Diagnostics: load_snapshots, isotropic_powerspectrum, δ
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

# load all the data!!
println("Loading data from $filename...")
snapshots = load_snapshots(filename; metadata)

# Let's pick the last snapshot!
times = snapshots[:T].times
snapshot_number = length(times)÷2 + 48
nday = @sprintf("%02.0f", (times[snapshot_number])/60^2/24)
println("Plotting snapshot $snapshot_number on day $(nday)...")
t0 = now()

T = snapshots[:T][snapshot_number]
u = snapshots[:u][snapshot_number]
v = snapshots[:v][snapshot_number]
w = snapshots[:w][snapshot_number]
f = parameters.f 
ro = compute!(Field(ζ(snapshots, snapshot_number)/f))
rd = compute!(Field(δ(snapshots, snapshot_number)/f))

println("Loading fields wall time: $((now() - t0).value/1e3) seconds.")

# Coordinate arrays
xu, yu, zu = nodes(u)
xv, yv, zv = nodes(v)
xw, yw, zw = nodes(w)
xT, yT, zT = nodes(T)

# Compute the horizontal spectrum of T, u, v, w 
axis_kwargs1 = (xlabel = "Wavenumber (rad⋅m⁻¹)", 
                ylabel = L"E_i(k)/E_{T,v}(k_{min},z=-3~m)",
                xscale = log10, yscale = log10,
                limits = ((8e-5, 4e-2), (1e-9,1e1)))
axis_kwargs2 = NamedTuple{(:xlabel,:xscale,:yscale,:limits)}(axis_kwargs1)
fig = Figure(size = (800, 300))
for (i,klev) in enumerate([127, 116, 98])
    println("Plotting spectra at z = $(zT[klev])m...")
    if i == 1
        ax = Axis(fig[1, i]; title="z=$(zT[klev]) m", axis_kwargs1...)
    else
        ax = Axis(fig[1, i]; title="z=$(zT[klev]) m", axis_kwargs2...)
        hideydecorations!(ax, ticks = false)
    end

    Su = isotropic_powerspectrum(interior(u, :, :, klev), interior(u, :, :, klev), xu, yu)
    Sv = isotropic_powerspectrum(interior(v, :, :, klev), interior(v, :, :, klev), xv, yv)
    wk = (interior(w, :, :, klev)+interior(w, :, :, klev+1))/2
    Sw = isotropic_powerspectrum(wk, wk, xw, yw)
    St = isotropic_powerspectrum(interior(T, :, :, klev), interior(T, :, :, klev), xT, yT)

    if i == 1
        global Sv0,St0 = Sv,St
    end

    lines!(ax, Su.freq, 1e-8Su.freq.^-2, linestyle = :dash, color = :black)
    text!(ax, 1e-3, 1e-2; text = L"k^{-2}")
    lines!(ax, Su.freq, 1e-15Su.freq.^-3, linestyle = :dash, color = :gray)
    text!(ax, 10^-3.5, 1e-6; text = L"k^{-3}")
    lines!(ax, St.freq, Real.(St.spec./St0.spec[1]), color = :red, label = L"E_T")
    lines!(ax, Su.freq, Real.(Su.spec./Sv0.spec[1]), color = :blue, label = L"E_u")
    lines!(ax, Sv.freq, Real.(Sv.spec./Sv0.spec[1]), color = :green, label = L"E_v")
    lines!(ax, Sw.freq, Real.(Sw.spec./Sv0.spec[1]), color = :black, label = L"E_w")
    xlims!(ax, (8e-5, 4e-2))
    vlines!(ax, [2π/10^4]; color = :black, linewidth = 0.5)
    axislegend(ax, labelsize=10, patchsize = (20, 5))
end
save(filesave * "spectra_" * fileparams * "_d$(nday).pdf", fig)
println("Finished plotting spectra, wall time: $((now() - t0).value/1e3) seconds.")

axis_kwargs1 = (xlabel = "Wavenumber (rad⋅m⁻¹)", 
                ylabel = L"E_i(k)/E_{\zeta/f, \delta/f}(k_{min},z=-3~m)",
                xscale = log10, yscale = log10,
                limits = ((8e-5, 4e-2), (1e-1,1e4)))
axis_kwargs2 = NamedTuple{(:xlabel,:xscale,:yscale,:limits)}(axis_kwargs1)
fig = Figure(size = (800, 300))
for (i,klev) in enumerate([127, 116, 98])
    println("Plotting spectra at z = $(zT[klev])m...")
    if i == 1
        ax = Axis(fig[1, i]; title="z=$(zT[klev]) m", axis_kwargs1...)
    else
        ax = Axis(fig[1, i]; title="z=$(zT[klev]) m", axis_kwargs2...)
        hideydecorations!(ax, ticks = false)
    end

    So = isotropic_powerspectrum(interior(ro, :, :, klev), interior(ro, :, :, klev), xT, yT)
    Sd = isotropic_powerspectrum(interior(rd, :, :, klev), interior(rd, :, :, klev), xT, yT)
    if i == 1
        global So0,Sd0 = So,Sd
    end

    lines!(ax, So.freq, 10So.freq.^(1/3), linestyle = :dash, color = :black)
    text!(ax, 10^-3.5, 0.2; text = L"k^{1/3}")
    lines!(ax, So.freq, 1e3So.freq.^(2/3), linestyle = :dash, color = :gray)
    text!(ax, 10^-2, 20; text = L"k^{2/3}")
    lines!(ax, So.freq, Real.(So.spec./So0.spec[1]), color = :blue, label = L"E_{\zeta/f}")
    lines!(ax, Sd.freq, Real.(Sd.spec./Sd0.spec[1]), color = :green, label = L"E_{\delta/f}")
    xlims!(ax, (8e-5, 4e-2))
    vlines!(ax, [2π/10^4]; color = :black, linewidth = 0.5)
    axislegend(ax, labelsize=10, patchsize = (20, 5))
end
save(filesave * "spectra2_" * fileparams * "_d$(nday).pdf", fig)
println("Finished plotting spectra, wall time: $((now() - t0).value/1e3) seconds.")
