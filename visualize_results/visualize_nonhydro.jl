using LESStudySetup
using CairoMakie
using SixelTerm
using Printf, Dates
using Statistics: mean, std
set_theme!(theme_latexfonts(), fontsize=12,figure_padding = 5)

filename = "./1hnonhydrostatic_snapshots_0.jld2"
T1 = FieldTimeSeries(filename, string(:T); architecture=CPU(), backend = OnDisk())[5]
x1, y1, z1 = nodes(T1)

fig = Figure()
gl = fig[1, 1] = GridLayout()
ax = Axis(gl[1,1]; xlabel="x (m)", ylabel="y (m)", title = "1024 GPUs")
xmin, xmax = 50, 50
ymin, ymax = 50, 50
Nxmin, Nxmax = 320, 361
Tmin, Tmax = minimum(interior(T1, Nxmin:Nxmax, :, length(z1))), maximum(interior(T1, Nxmin:Nxmax, :, length(z1)))
for i = 512:575#; 576:591; 608:623; 640:655; 672:687; 704:719; 736:751; 768:783; 800:815; 832:847; 864:879; 896:911; 928:943; 960:975; 992:999]
    println(i)
    filename = "../nonhydro1024gpus/nonhydrostatic_snapshots_$(i).jld2"
    T = FieldTimeSeries(filename, string(:T); architecture=CPU(), backend = OnDisk())[1]
    x, y, z = nodes(T)
    xmin, xmax = min(xmin, minimum(x/1e3)), max(xmax, maximum(x/1e3))
    ymin, ymax = min(ymin, minimum(y/1e3)), max(ymax, maximum(y/1e3))
    hm = heatmap!(ax, x/1e3, y/1e3, interior(T, :, :, length(z)); rasterize = true, colormap = :thermal, colorrange = (Tmin, Tmax))
end
xlims!(ax, (xmin, xmax))
ylims!(ax, (ymin, ymax))

ax = Axis(gl[1,2]; xlabel="x (m)", ylabel="y (m)", limits = ((xmin, xmax), (ymin, ymax)), title = "16 GPUs")
for i = 8:11
    filename1 = "../nonhydro16gpus/nonhydrostatic_snapshots_$(i).jld2"
    T4 = FieldTimeSeries(filename1, string(:T); architecture=CPU(), backend = OnDisk())[3]
    x4, y4, z4 = nodes(T4)
    hm = heatmap!(ax, x4/1e3, y4/1e3, interior(T4, :, :, length(z4)); rasterize = true, colormap = :thermal, colorrange = (Tmin, Tmax))
end
hideydecorations!(ax, ticks = false)
ax = Axis(gl[1,3]; xlabel="x (m)", ylabel="y (m)", limits = ((xmin, xmax), (ymin, ymax)), title = "1 GPU")
hm = heatmap!(ax, x1/1e3, y1/1e3, interior(T1, :, :, length(z1)); rasterize = true, colormap = :thermal, colorrange = (Tmin, Tmax))
hideydecorations!(ax, ticks = false)
Colorbar(gl[1,4], hm, label = "Temperature (°C)")
colgap!(gl, 1, 10)
colgap!(gl, 2, 10)
colgap!(gl, 3, 3)
resize_to_layout!(fig)
save("T1h_nonhydro_1024vs16vs1gpus.pdf", fig; pt_per_unit = 1)