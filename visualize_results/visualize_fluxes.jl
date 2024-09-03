using LESStudySetup
using CairoMakie
using SixelTerm
using Printf, Dates
using Oceananigans: compute!
using Oceananigans.Grids: xnodes, ynodes, znodes
using Oceananigans.Units: kilometer
using Statistics: mean
using LESStudySetup.Diagnostics
using LESStudySetup.Diagnostics: load_snapshots, symmetric_filtering, u²_filtered, v²_filtered, w²_filtered
using LESStudySetup.Diagnostics: uv_filtered, uw_filtered, vw_filtered
using LESStudySetup.Diagnostics: ub_filtered, vb_filtered, wb_filtered
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
    snapshot_number = length(times)÷2 + 48
    nday = @sprintf("%2.0f", (times[snapshot_number])/60^2/24)
    println("Plotting snapshot $snapshot_number on day $(nday)...")

    axis_kwargs0 = (xlabel = L"\langle \phi^\prime \varphi^\prime \rangle_{xy}",
                    ylabel = "z (m)",
                    limits = (nothing, (-250, 0)))
    axis_kwargs1 = NamedTuple{(:xlabel,:limits)}(axis_kwargs0)
    axis_kwargs2 = NamedTuple{(:ylabel,:limits)}(axis_kwargs0)
    fig = Figure(size = (900, 900))

    cutoff = 10kilometer
    u², ul², uh² = u²_filtered(snapshots, snapshot_number; cutoff)
    v², vl², vh² = v²_filtered(snapshots, snapshot_number; cutoff)
    w², wl², wh² = w²_filtered(snapshots, snapshot_number; cutoff)
    uv, uvl, uvh = uv_filtered(snapshots, snapshot_number; cutoff)
    uw, uwl, uwh = uw_filtered(snapshots, snapshot_number; cutoff)
    vw, vwl, vwh = vw_filtered(snapshots, snapshot_number; cutoff)
    ub, ubl, ubh = ub_filtered(snapshots, snapshot_number; cutoff)
    vb, vbl, vbh = vb_filtered(snapshots, snapshot_number; cutoff)
    wb, wbl, wbh = wb_filtered(snapshots, snapshot_number; cutoff)
    
    _, _, zu = nodes(u²)
    _, _, zw = nodes(w²)
    zmin = minimum(zw)

    # Plot the fluxes
    ax_uu = Axis(fig[1, 1]; title=L"\langle u^\prime u^\prime \rangle_{xy}", axis_kwargs2...)
    ax_uv = Axis(fig[2, 1]; title=L"\langle u^\prime v^\prime \rangle_{xy}", axis_kwargs2...)
    ax_ub = Axis(fig[3, 1]; title=L"\langle u^\prime b^\prime \rangle_{xy}", axis_kwargs0...)
    ax_vv = Axis(fig[1, 2]; title=L"\langle v^\prime v^\prime \rangle_{xy}", limits = (nothing, (zmin, 0)))
    ax_uw = Axis(fig[2, 2]; title=L"\langle u^\prime w^\prime \rangle_{xy}", limits = (nothing, (zmin, 0)))
    ax_vb = Axis(fig[3, 2]; title=L"\langle v^\prime b^\prime \rangle_{xy}", axis_kwargs1...)
    ax_ww = Axis(fig[1, 3]; title=L"\langle w^\prime w^\prime \rangle_{xy}", limits = (nothing, (zmin, 0)))
    ax_vw = Axis(fig[2, 3]; title=L"\langle v^\prime w^\prime \rangle_{xy}", limits = (nothing, (zmin, 0)))
    ax_wb = Axis(fig[3, 3]; title=L"\langle w^\prime b^\prime \rangle_{xy}", axis_kwargs1...)

    lines!(ax_uu, vec(mean(compute!(Field(u²)),dims=(1,2))),zu; color = :black, label = "Total")
    lines!(ax_uu, vec(mean(compute!(Field(ul²)),dims=(1,2))),zu; color = :red, label = "Low-pass")
    lines!(ax_uu, vec(mean(compute!(Field(uh²)),dims=(1,2))),zu; color = :blue, label = "High-pass")
    axislegend(ax_uu, position = :rb)
    lines!(ax_vv, vec(mean(compute!(Field(v²)),dims=(1,2))),zu; color = :black)
    lines!(ax_vv, vec(mean(compute!(Field(vl²)),dims=(1,2))),zu; color = :red)
    lines!(ax_vv, vec(mean(compute!(Field(vh²)),dims=(1,2))),zu; color = :blue)
    lines!(ax_ww, vec(mean(compute!(Field(w²)),dims=(1,2))),zw; color = :black)
    lines!(ax_ww, vec(mean(compute!(Field(wl²)),dims=(1,2))),zw; color = :red)
    lines!(ax_ww, vec(mean(compute!(Field(wh²)),dims=(1,2))),zw; color = :blue)
    println("first row done...")
    lines!(ax_uv, vec(mean(compute!(Field(uv)),dims=(1,2))),zu; color = :black)
    lines!(ax_uv, vec(mean(compute!(Field(uvl)),dims=(1,2))),zu; color = :red)
    lines!(ax_uv, vec(mean(compute!(Field(uvh)),dims=(1,2))),zu; color = :blue)
    lines!(ax_uw, vec(mean(compute!(Field(uw)),dims=(1,2))),zu; color = :black)
    lines!(ax_uw, vec(mean(compute!(Field(uwl)),dims=(1,2))),zu; color = :red)
    lines!(ax_uw, vec(mean(compute!(Field(uwh)),dims=(1,2))),zu; color = :blue)
    lines!(ax_vw, vec(mean(compute!(Field(vw)),dims=(1,2))),zu; color = :black)
    lines!(ax_vw, vec(mean(compute!(Field(vwl)),dims=(1,2))),zu; color = :red)
    lines!(ax_vw, vec(mean(compute!(Field(vwh)),dims=(1,2))),zu; color = :blue)
    println("second row done...")
    lines!(ax_ub, vec(mean(compute!(Field(ub)),dims=(1,2))),zu; color = :black)
    lines!(ax_ub, vec(mean(compute!(Field(ubl)),dims=(1,2))),zu; color = :red)
    lines!(ax_ub, vec(mean(compute!(Field(ubh)),dims=(1,2))),zu; color = :blue)
    lines!(ax_vb, vec(mean(compute!(Field(vb)),dims=(1,2))),zu; color = :black)
    lines!(ax_vb, vec(mean(compute!(Field(vbl)),dims=(1,2))),zu; color = :red)
    lines!(ax_vb, vec(mean(compute!(Field(vbh)),dims=(1,2))),zu; color = :blue)
    lines!(ax_wb, vec(mean(compute!(Field(wb)),dims=(1,2))),zw; color = :black)
    lines!(ax_wb, vec(mean(compute!(Field(wbl)),dims=(1,2))),zw; color = :red)
    lines!(ax_wb, vec(mean(compute!(Field(wbh)),dims=(1,2))),zw; color = :blue)

    save(filesave * "10kmfilteredfluxes_" * fileparams * "_d$(nday).pdf", fig)

    return

end

cooling, wind, dTf,a = 50, 0.1, -1,1.0
visualize(cooling, wind, dTf, a)


cooling = @sprintf("%03d", cooling)
wind = replace("$(wind)","." => "" )
a = replace("$(a)","." => "" )
dTf = @sprintf("%1d", dTf)
fileparams = "free_surface_short_test_$(cooling)_wind_$(wind)_dTf_$(dTf)_a_$(a)"
fileparams = "hydrostatic_twin_simulation"
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
nday = @sprintf("%2.0f", (times[snapshot_number])/60^2/24)

ke = compute!(Field(KE(snapshots, snapshot_number)))
cutoff = 10kilometer
kel, keh = symmetric_filtering(ke;cutoff)
bnd, lbnd, hbnd = maximum(abs, interior(ke)), maximum(abs, interior(kel)), maximum(abs, interior(keh))

# Coordinate arrays
xu, yu, zu = nodes(ke)
ku = 113

fig = Figure(size = (800, 250))
axis_kwargs = (xlabel = "x (km)", ylabel = "y (km)",
            limits = ((0, 100), (0, 100)))
axis_kwargs2 = NamedTuple{(:xlabel,:limits)}(axis_kwargs)
ax_u = Axis(fig[1, 1][1,1]; title="Kinetic energy (m²/s²), z=$(zu[ku])m", axis_kwargs...)
ax_ul = Axis(fig[1, 2][1,1]; title="Low-passed KE (m²/s²), z=$(zu[ku])m", axis_kwargs2...)
ax_uh = Axis(fig[1, 3][1,1]; title="High-passed KE (m²/s²), z=$(zu[ku])m", axis_kwargs2...)    
hm_u =  heatmap!(ax_u,  1e-3xu, 1e-3yu, interior(ke,:,:,ku); rasterize = true, colormap = :amp, colorrange = (0,bnd))
hm_ul = heatmap!(ax_ul, 1e-3xu, 1e-3yu, interior(kel,:,:,ku); rasterize = true, colormap = :amp, colorrange = (0,lbnd))
hm_uh = heatmap!(ax_uh, 1e-3xu, 1e-3yu, interior(keh,:,:,ku); rasterize = true, colormap = :amp, colorrange = (0,hbnd))
hideydecorations!(ax_ul, ticks = false)
hideydecorations!(ax_uh, ticks = false)
Colorbar(fig[1, 1][1, 2], hm_u)
Colorbar(fig[1, 2][1, 2], hm_ul)
Colorbar(fig[1, 3][1, 2], hm_uh)
save(filesave * "10kmfilteredKE_" * fileparams * "_d$(nday).pdf", fig)