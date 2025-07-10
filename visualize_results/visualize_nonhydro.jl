using LESStudySetup
using CairoMakie
using Printf, Dates
using Statistics: mean, std
using LESStudySetup.Diagnostics
using LESStudySetup.Diagnostics: load_distributed_checkpoint,load_subdomain_snapshot
set_theme!(theme_latexfonts(), fontsize=12, figure_padding = 5)

# --- Simulation and Subdomain Parameters ---
filehead = "/orcd/data/abodner/002/nhyles_output/" 
iteration = 32207
fileparam = "subdomain1"
# 1. Define the filename of the saved snapshot
output_filename = filehead * "subdomains/" * fileparam * "_snapshot_iter$(iteration).jld2"

# 2. Load the snapshot using the new function
snapshot = load_subdomain_snapshot(output_filename)

x, y, z = nodes(snapshot[:T]);
k = 1
Tmap, wmap, vmap = :thermal,:delta,:balance
wmax,umax,vmax=1e-3,1e-2,1e-2#0.02,0.15,0.2
fig = Figure(size = (640, 480))
gabc = fig[1, 1] = GridLayout()
aspect = 2
axis_kwargs = (ylabel = L"y~\text{(km)}", aspect=aspect)
ax_a = Axis(gabc[1,1]; titlealign = :left, title=L"\text{(a)}~T~\text{({^\circ}C)}", axis_kwargs...)
ax_b = Axis(gabc[1,3]; titlealign = :left, title=L"\text{(b)}~w~\text{(m s^{-1})}", aspect=aspect)
ax_c = Axis(gabc[3,1]; titlealign = :left, title=L"\text{(c)}~u~\text{(m s^{-1})}", axis_kwargs...) 
ax_d = Axis(gabc[3,3]; titlealign = :left, title=L"\text{(d)}~v~\text{(m s^{-1})}", aspect=aspect) 
hm_a = heatmap!(ax_a, 1e-3x, 1e-3y, (interior(snapshot[:T],:,:,k)); rasterize = true, colormap = Tmap)
hm_b = heatmap!(ax_b, 1e-3x, 1e-3y, (interior(snapshot[:w],:,:,k)); rasterize = true, colormap = wmap, colorrange = (-wmax, wmax))
hm_c = heatmap!(ax_c, 1e-3x, 1e-3y, (interior(snapshot[:u],:,:,k)); rasterize = true, colormap = vmap, colorrange = (-umax, umax))
hm_d = heatmap!(ax_d, 1e-3x, 1e-3y, (interior(snapshot[:v],:,:,k)); rasterize = true, colormap = vmap, colorrange = (-vmax, vmax))
Colorbar(gabc[1,2], hm_a)
Colorbar(gabc[1,4], hm_b)
Colorbar(gabc[3,2], hm_c)
Colorbar(gabc[3,4], hm_d)
hidexdecorations!(ax_a, ticks = false)
hidexdecorations!(ax_b, ticks = false)
hidexdecorations!(ax_c, ticks = false)
hidexdecorations!(ax_d, ticks = false)
hideydecorations!(ax_b, ticks = false)
hideydecorations!(ax_d, ticks = false)

# T̄ = (mean(snapshot[:T], dims = 2));
# w̄ = (mean(snapshot[:w], dims = 2));
# ū = (mean(snapshot[:u], dims = 2));
# v̄ = (mean(snapshot[:v], dims = 2));
zmin = -252
kz = findfirst(z .≥ zmin)
Nz = length(z)
axis_kwargs0 = (xlabel = L"x~\text{(km)}", ylabel = L"z~\text{(m)}", limits = (nothing, (zmin, 0)))
axis_kwargs1 = NamedTuple{(:xlabel,:ylabel)}(axis_kwargs0)
ax_a = Axis(gabc[2,1]; titlealign = :left, axis_kwargs0...)
ax_b = Axis(gabc[2,3]; titlealign = :left, xlabel = L"x~\text{(km)}", limits = (nothing, (zmin, 0)))
ax_c = Axis(gabc[4,1]; titlealign = :left,  limits = (nothing, (zmin, 0)), axis_kwargs1...)
ax_d = Axis(gabc[4,3]; titlealign = :left, xlabel = L"x~\text{(km)}", limits = (nothing, (zmin, 0)))
wmax,umax,vmax=0.005,0.05,0.2
hm_a = heatmap!(ax_a, 1e-3x, z[kz:Nz], (interior(snapshot[:T],:,1,kz:Nz)); rasterize = true, colormap = Tmap)#, colorrange = (vmin, vmax))
hm_b = heatmap!(ax_b, 1e-3x, z[kz:Nz], (interior(snapshot[:w],:,1,kz:Nz)); rasterize = true, colormap = wmap, colorrange = (-wmax, wmax))
hm_c = heatmap!(ax_c, 1e-3x, z[kz:Nz], (interior(snapshot[:u],:,1,kz:Nz)); rasterize = true, colormap = vmap, colorrange = (-umax, umax))
hm_d = heatmap!(ax_d, 1e-3x, z[kz:Nz], (interior(snapshot[:v],:,1,kz:Nz)); rasterize = true, colormap = vmap, colorrange = (-vmax, vmax))
hideydecorations!(ax_b, ticks = false)
hideydecorations!(ax_d, ticks = false)
Colorbar(gabc[2,2], hm_a)
Colorbar(gabc[2,4], hm_b)
Colorbar(gabc[4,2], hm_c)
Colorbar(gabc[4,4], hm_d)
rowgap!(gabc, 4)
colgap!(gabc, 1, 1)
colgap!(gabc, 3, 1)
colgap!(gabc, 2, 5)
for row = [2,4]
    rowsize!(gabc, row, Relative(0.14))
end
resize_to_layout!(fig)
filesave = "results/"
save(filesave * "Twuv_" * fileparam * "_36h_iter$(iteration).pdf", fig; pt_per_unit = 1)
println("Finished plotting Twuv fields")

# filename = "../../nhyles_output/nonhydrostatic_checkpoint_"
# iteration = 41360
# level = 216
# snapshot = load_distributed_checkpoint(filename, iteration; level);

# rmin, rmax = -0.4,0.4
# x, y, z = nodes(snapshot[:T]);
# fig = Figure(size = (750, 660))
# gabc = fig[1, 1] = GridLayout()
# axis_kwargs = (ylabel = "y (km)", aspect=1)
# ax_a = Axis(gabc[1,1]; titlealign = :left, title=L"\text{(a)}~T~\text{({^\circ}C)}", axis_kwargs...)
# ax_b = Axis(gabc[1,3]; titlealign = :left, title=L"\text{(b)}~u~\text{(m/s)}", aspect=1)
# ax_c = Axis(gabc[2,1]; titlealign = :left, title=L"\text{(c)}~v~\text{(m/s)}", xlabel = "x (km)", ylabel = "y (km)", aspect=1) 
# ax_d = Axis(gabc[2,3]; titlealign = :left, title=L"\text{(d)}~w~\text{(m/s)}", xlabel = "x (km)", aspect=1) 
# hm_a = heatmap!(ax_a, 1e-3x, 1e-3y, interior(snapshot[:T],:,:,1); rasterize = true, colormap = :thermal, colorrange = (19.6, 20.2))
# hm_b = heatmap!(ax_b, 1e-3x, 1e-3y, interior(snapshot[:u],:,:,1); rasterize = true, colormap = :balance, colorrange = (2rmin, 2rmax))
# hm_c = heatmap!(ax_c, 1e-3x, 1e-3y, interior(snapshot[:v],:,:,1); rasterize = true, colormap = :balance, colorrange = (rmin, rmax))
# hm_d = heatmap!(ax_d, 1e-3x, 1e-3y, interior(snapshot[:w],:,:,1); rasterize = true, colormap = :balance, colorrange = (rmin/40, rmax/40))
# Colorbar(gabc[1,2], hm_a)
# Colorbar(gabc[1,4], hm_b)
# Colorbar(gabc[2,2], hm_c)
# Colorbar(gabc[2,4], hm_d)
# hidexdecorations!(ax_a, ticks = false)
# hidexdecorations!(ax_b, ticks = false)
# hideydecorations!(ax_b, ticks = false)
# hideydecorations!(ax_d, ticks = false)
# rowgap!(gabc, 3)
# colgap!(gabc, 1, 3)
# colgap!(gabc, 2, 10)
# colgap!(gabc, 3, 3)
# resize_to_layout!(fig)
# save("Tuvw35h_nonhydro_1024gpus_k216.pdf", fig; pt_per_unit = 1)