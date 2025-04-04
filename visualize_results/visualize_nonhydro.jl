using LESStudySetup
using CairoMakie
using Printf, Dates
using Statistics: mean, std
using LESStudySetup.Diagnostics
using LESStudySetup.Diagnostics: load_distributed_checkpoint
set_theme!(theme_latexfonts(), fontsize=12,figure_padding = 5)

filename = "../nonhydro1024gpus/nonhydrostatic_checkpoint_"
iteration = 520
level = 224

snapshot = load_distributed_checkpoint(filename, iteration; level)

rmin, rmax = -0.4,0.4
x, y, z = nodes(snapshot[:T])
fig = Figure(size = (750, 660))
gabc = fig[1, 1] = GridLayout()
axis_kwargs = (ylabel = "y (km)", aspect=1)
ax_a = Axis(gabc[1,1]; titlealign = :left, title=L"\text{(a)}~T~\text{({^\circ}C)}", axis_kwargs...)
ax_b = Axis(gabc[1,3]; titlealign = :left, title=L"\text{(b)}~u~\text{(m/s)}", aspect=1)
ax_c = Axis(gabc[2,1]; titlealign = :left, title=L"\text{(c)}~v~\text{(m/s)}", xlabel = "x (km)", ylabel = "y (km)", aspect=1) 
ax_d = Axis(gabc[2,3]; titlealign = :left, title=L"\text{(d)}~w~\text{(m/s)}", xlabel = "x (km)", aspect=1) 
hm_a = heatmap!(ax_a, 1e-3x, 1e-3y, interior(snapshot[:T],:,:,1); rasterize = true, colormap = :thermal, colorrange = (19.6, 20.2))
hm_b = heatmap!(ax_b, 1e-3x, 1e-3y, interior(snapshot[:u],:,:,1); rasterize = true, colormap = :balance, colorrange = (2rmin, 2rmax))
hm_c = heatmap!(ax_c, 1e-3x, 1e-3y, interior(snapshot[:v],:,:,1); rasterize = true, colormap = :balance, colorrange = (rmin, rmax))
hm_d = heatmap!(ax_d, 1e-3x, 1e-3y, interior(snapshot[:w],:,:,1); rasterize = true, colormap = :balance, colorrange = (rmin/40, rmax/40))
Colorbar(gabc[1,2], hm_a)
Colorbar(gabc[1,4], hm_b)
Colorbar(gabc[2,2], hm_c)
Colorbar(gabc[2,4], hm_d)
hidexdecorations!(ax_a, ticks = false)
hidexdecorations!(ax_b, ticks = false)
hideydecorations!(ax_b, ticks = false)
hideydecorations!(ax_d, ticks = false)
rowgap!(gabc, 3)
colgap!(gabc, 1, 3)
colgap!(gabc, 2, 10)
colgap!(gabc, 3, 3)
resize_to_layout!(fig)
save("Tuvw1h_nonhydro_1024gpus.pdf", fig; pt_per_unit = 1)

#############################
for rank in 272:299
    @info "loading rank $rank of $(Px * Py - 1)"

    file = jldopen(filename * "$(rank).jld2")

    Rx = file["grid/architecture/local_index/1"]
    Ry = file["grid/architecture/local_index/2"]

    udata = file["timeseries/u/" * "$(Ni)"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    vdata = file["timeseries/v/" * "$(Ni)"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    wdata = file["timeseries/w/" * "$(Ni)"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    Tdata = file["timeseries/T/" * "$(Ni)"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    
    irange = 1 + (Rx - 1) * nx : Rx * nx
    jrange = 1 + (Ry - 1) * ny : Ry * ny

    interior(u, irange, jrange, :) .= udata[:, :, indices[3]] 
    interior(v, irange, jrange, :) .= vdata[:, :, indices[3]]
    interior(w, irange, jrange, :) .= wdata[:, :, indices[3]]
    interior(T, irange, jrange, :) .= Tdata[:, :, indices[3]]
end