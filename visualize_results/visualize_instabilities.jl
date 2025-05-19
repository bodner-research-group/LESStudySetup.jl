using LESStudySetup
using CairoMakie
using Printf, Dates
using Oceananigans: compute!
using Oceananigans.BoundaryConditions
using Oceananigans.Grids: xnodes, ynodes, znodes
using Oceananigans.Operators: div_xyᶜᶜᶜ
using Statistics: mean, std, median
#using StatsBase: mode, Histogram, fit
using LESStudySetup.Diagnostics
using LESStudySetup.Oceananigans.Units
using LESStudySetup.Diagnostics: N², M², Bₕ, wb
using LESStudySetup.Diagnostics: load_snapshots, MLD, MLaverage, BLD
using LESStudySetup.Diagnostics: isotropic_powerspectrum, coarse_grained_fluxes, δ
using LESStudySetup.Diagnostics: MixedLayerN², MixedLayerDepth, MLI, spatial_filtering
using LESStudySetup.Diagnostics: subfilter_stress!, coarse_graining!, _horizontal_gauss_filter!
using LESStudySetup.Diagnostics: build_lanczos_kernel, build_tophat_kernel
using MathTeXEngine,GibbsSeaWater
set_theme!(theme_latexfonts(), fontsize=12,figure_padding = 10)

function cumtrapz(X, Y) 
    # Check matching vector length
    @assert length(X) == length(Y)
    # Initialize Output
    out = similar(X)
    out[1] = 0
    # Iterate over arrays
    for i in 2:length(X)
      out[i] = out[i-1] + 0.5*(X[i] - X[i-1])*(Y[i] + Y[i-1])
    end
    # Return output
    out
end

# Reduced Major Axis Regression (through the origin)
function rma_slope(X::AbstractVector, Y::AbstractVector)
    # Check that X and Y have the same length
    if length(X) != length(Y)
        error("Vectors X and Y must have the same length.")
    end

    # Compute sum of squares of X and Y
    sumX2 = sum(x -> x^2, X)
    sumY2 = sum(y -> y^2, Y)
    
    # Compute the dot product
    dotXY = sum(x * y for (x, y) in zip(X, Y))
    
    # Calculate the slope using the RMA formula for a line through the origin
    k = sign(dotXY) * sqrt(sumY2 / sumX2)
    
    return k
end
filehead = "/orcd/data/abodner/002/shirui/LESStudySetup.jl/"

fileparams = "hydrostatic_twin_simulation"
filename = filehead * "hydrostatic_snapshots_" * fileparams * ".jld2"
metadata = filehead * "experiment_" * fileparams * "_metadata.jld2"
freename = filehead * "hydrostatic_free_surface_" * fileparams * ".jld2"
#fileparams = "nonhydrostatic"
#filename = filehead * fileparams * "_snapshots.jld2"
#metadata = filehead * "nonhydrostatic_experiment_metadata.jld2"
filesave = filehead * "results/"

# load all the data!!
println("Loading data from $filename...")
snapshots = load_snapshots(filename; metadata,variables = (:u, :v, :w, :T, :κu))
times = snapshots[:T].times;
g = parameters.g;
α = parameters.α;
M²₀ = parameters.M²₀;
f = parameters.f;
N²s = parameters.N²s;

ls,lc = (1/4e-4), (1/4e-4);
kernel = :tophat;
cutoff = (1/4e-3)/2.4*2*π #785.40
kernel! = _horizontal_gauss_filter!;

i = 97
nday = @sprintf("%2.1f", (times[i])/60^2/24)
println("Plotting snapshot $i on day $(nday)...")
t0 = now();
wi,Ti = snapshots[:w][i],snapshots[:T][i];
ui,vi,κi = snapshots[:u][i],snapshots[:v][i],snapshots[:κu][i];
hi = compute!(MLD(snapshots,i; threshold = 0.09));
∇bi = compute!(Field(α * g * (∂x(Ti)^2 + ∂y(Ti)^2)^0.5));
h0 = compute!(MLD(snapshots,i; threshold = 0.09));
h1 = compute!(BLD(snapshots,i));

# Plot the fields
pB, pC = (-20, 0), (30, -25)
lx, ly = 40, 25
p2 = [pB, pC]
ΔN = 160
jslices = [2]*ΔN
shift(x) = [x[size(x,1)÷2+1:end, :]; x[1:size(x,1)÷2, :]]
yshift(x) = [x[:, size(x,2)÷2+1:end] x[:, 1:size(x,2)÷2]]
xhift(x) = yshift(shift(x))
##################################### 
var,scale = ∇bi, 1/M²₀;
x, y, z = nodes(var);
k = 202
cmap = Reverse(:grays)
rmin, rmax = 0, 1
hcolor, scolor = [:blue,:green], :black
fig = Figure(size = (640, 720))
gabc = fig[1, 1] = GridLayout()
axis_kwargs = (ylabel = "y (km)", aspect=1, limits = ((-50, 50), (-50, 50)))
ax_a = Axis(gabc[1,1]; titlealign = :left, title=L"\text{(a)}~\nabla{b}/M_0^2,~z=-2.8~\text{m}", axis_kwargs...)
ax_b = Axis(gabc[1,2]; titlealign = :left, title=L"\text{(b)}~\nabla{b}/M_0^2,~z=-25.3~\text{m}", aspect=1, limits = ((0, 100), (-50, 50)))
ax_c = Axis(gabc[3,1]; titlealign = :left, title=L"\text{(c)~Unstable front},~z=-25.3~\text{m}", ylabel = "y (km)", aspect=8/5, limits = ((pB[1], pB[1]+lx), (pB[2]-ly,pB[2]))) 
ax_d = Axis(gabc[3,2]; titlealign = :left, title=L"\text{(d)~Stable front},~z=-25.3~\text{m}", aspect=8/5, limits = ((pC[1], pC[1]+lx), (pC[2]-ly, pC[2]))) 
hm_a = heatmap!(ax_a, 1e-3x.-50, 1e-3y.-50, scale * xhift(interior(var,:,:,k+20)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
hm_b = heatmap!(ax_b, 1e-3x, 1e-3y.-50, scale * yshift(interior(var,:,:,k)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
hm_c = heatmap!(ax_c, 1e-3x.-50, 1e-3y.-50, scale * xhift(interior(var,:,:,k)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
hm_d = heatmap!(ax_d, 1e-3x, 1e-3y.-50, scale * yshift(interior(var,:,:,k)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
Colorbar(gabc[1,3], hm_b)
Colorbar(gabc[3,3], hm_d)
hidexdecorations!(ax_a, ticks = false)
hidexdecorations!(ax_b, ticks = false)
hidexdecorations!(ax_c, ticks = false)
hidexdecorations!(ax_d, ticks = false)
hlines!(ax_a, 0; color = :black, linestyle = :dash, linewidth = 0.8)
hlines!(ax_b, -50; color = :black, linestyle = :dash, linewidth = 0.8)
poly!(ax_a, Rect(pB[1], pB[2]-ly, lx, ly), color = (:white, 0.1), strokecolor = scolor, strokewidth = 0.5)
poly!(ax_b, Rect(pC[1], pC[2]-ly, lx, ly), color = (:white, 0.1), strokecolor = scolor, strokewidth = 0.5)
text!(ax_a, pB[1], pB[2], text = L"\text{Unstable}", color = :black, align = (:left, :top))
text!(ax_b, pC[1], pC[2], text = L"\text{Stable}", color = :black, align = (:left, :top))
hlines!(ax_c, -10; color = :black, linestyle = :dash, linewidth = 0.8)
hlines!(ax_d, -40; color = :black, linestyle = :dash, linewidth = 0.8)

zmin = -100
ja = 2ΔN
jc = 9*64
Δj = 64
kz = findlast(z .< zmin)
Nz = length(z)
axis_kwargs0 = (xlabel = "x (km)", ylabel = "z (m)", limits = ((-50, 50), (zmin, 0)))
axis_kwargs1 = NamedTuple{(:xlabel,:ylabel)}(axis_kwargs0)
ax_a = Axis(gabc[2,1]; titlealign = :left, title=L"y=0~\text{km}", axis_kwargs0...)
ax_b = Axis(gabc[2,2]; titlealign = :left, title=L"y=-50~\text{km}", xlabel = "x (km)", limits = ((0, 100), (zmin, 0)))
ax_c = Axis(gabc[4,1]; titlealign = :left, title=L"y=-10~\text{km}", limits = ((pB[1], pB[1]+lx), (zmin, 0)), axis_kwargs1...)
ax_d = Axis(gabc[4,2]; titlealign = :left, title=L"y=-40~\text{km}", xlabel = "x (km)", limits = ((pC[1], pC[1]+lx), (zmin, 0)))
vmin, vmax = rmin, rmax
hm_a = heatmap!(ax_a, 1e-3x.-50, z[kz:Nz], scale*shift(interior(var,:,2ja,kz:Nz)); rasterize = true, colormap = cmap, colorrange = (vmin, vmax))
hm_b = heatmap!(ax_b, 1e-3x, z[kz:Nz], scale*(interior(var,:,ja,kz:Nz)); rasterize = true, colormap = cmap, colorrange = (vmin, vmax))
hm_c = heatmap!(ax_c, 1e-3x.-50, z[kz:Nz], scale*shift(interior(var,:,jc,kz:Nz)); rasterize = true, colormap = cmap, colorrange = (vmin, vmax))
hm_d = heatmap!(ax_d, 1e-3x, z[kz:Nz], scale*(interior(var,:,ja+Δj,kz:Nz)); rasterize = true, colormap = cmap, colorrange = (vmin, vmax))
hideydecorations!(ax_b, ticks = false)
hideydecorations!(ax_d, ticks = false)
Colorbar(gabc[2,3], hm_b)
Colorbar(gabc[4,3], hm_d)
hlines!(ax_a, z[k]; color = :black, linestyle = :dash, linewidth = 0.8)
hlines!(ax_b, z[k]; color = :black, linestyle = :dash, linewidth = 0.8)
hlines!(ax_c, z[k]; color = :black, linestyle = :dash, linewidth = 0.8)
hlines!(ax_d, z[k]; color = :black, linestyle = :dash, linewidth = 0.8)
for (i,h) in enumerate([h0, h1])
    lines!(ax_a, 1e-3x.-50, -vec(shift(interior(h,:,2ja))); color = hcolor[i], linewidth = 0.8)
    lines!(ax_b, 1e-3x, -vec((interior(h,:,ja))); color = hcolor[i], linewidth = 0.8)
    lines!(ax_c, 1e-3x.-50, -vec(shift(interior(h,:,jc))); color = hcolor[i], linewidth = 0.8)
    lines!(ax_d, 1e-3x, -vec((interior(h,:,ja+Δj))); color = hcolor[i], linewidth = 0.8)
end
rowgap!(gabc, 3)
colgap!(gabc, 1, 15)
colgap!(gabc, 2, 3)
rowsize!(gabc, 2, Relative(0.13))
rowsize!(gabc, 3, Relative(0.28))
rowsize!(gabc, 4, Relative(0.14))
resize_to_layout!(fig)
save(filesave * "Dbhfields_" * fileparams * "_d6.pdf", fig; pt_per_unit = 1)

######################
T′ = compute!(Field(Ti - mean(Ti, dims = (1,2))));
u′ = compute!(Field(ui - mean(ui, dims = (1,2))));
v′ = compute!(Field(vi - mean(vi, dims = (1,2))));
w′ = compute!(Field(wi - mean(wi, dims = (1,2))));
x, y, z = nodes(Ti);
k = 202;
bmap, vmap = :diff,:delta
fig = Figure(size = (640, 750))
gabc = fig[1, 1] = GridLayout()
axis_kwargs = (ylabel = "y (km)", aspect=8/5, limits = ((pB[1], pB[1]+lx), (pB[2]-ly,pB[2])))
ax_a = Axis(gabc[1,1]; titlealign = :left, title=L"\text{(a)}~b\prime~\text{(10^{-3}m s^{-2})}", axis_kwargs...)
ax_b = Axis(gabc[1,2]; titlealign = :left, title=L"\text{(b)}~b\prime~\text{(10^{-3}m s^{-2})}", aspect=8/5, limits = ((pC[1], pC[1]+lx), (pC[2]-ly, pC[2])))
ax_c = Axis(gabc[2,1]; titlealign = :left, title=L"\text{(c)}~u\prime~\text{(m s^{-1})}", ylabel = "y (km)", aspect=8/5, limits = ((pB[1], pB[1]+lx), (pB[2]-ly,pB[2]))) 
ax_d = Axis(gabc[2,2]; titlealign = :left, title=L"\text{(d)}~u\prime~\text{(m s^{-1})}", aspect=8/5, limits = ((pC[1], pC[1]+lx), (pC[2]-ly, pC[2]))) 
ax_e = Axis(gabc[3,1]; titlealign = :left, title=L"\text{(e)}~v\prime~\text{(m s^{-1})}", ylabel = "y (km)", aspect=8/5, limits = ((pB[1], pB[1]+lx), (pB[2]-ly,pB[2]))) 
ax_f = Axis(gabc[3,2]; titlealign = :left, title=L"\text{(f)}~v\prime~\text{(m s^{-1})}", aspect=8/5, limits = ((pC[1], pC[1]+lx), (pC[2]-ly, pC[2]))) 
ax_g = Axis(gabc[4,1]; titlealign = :left, title=L"\text{(g)}~w\prime~\text{(m s^{-1})}", xlabel = "x (km)", ylabel = "y (km)", aspect=8/5, limits = ((pB[1], pB[1]+lx), (pB[2]-ly,pB[2]))) 
ax_h = Axis(gabc[4,2]; titlealign = :left, title=L"\text{(h)}~w\prime~\text{(m s^{-1})}", xlabel = "x (km)", aspect=8/5, limits = ((pC[1], pC[1]+lx), (pC[2]-ly, pC[2]))) 
rmin, rmax = -1, 1
hm_a = heatmap!(ax_a, 1e-3x.-50, 1e-3y.-50,1e3α * g * xhift(interior(T′,:,:,k)); rasterize = true, colormap = bmap, colorrange = (rmin, rmax))
hm_b = heatmap!(ax_b, 1e-3x, 1e-3y.-50,1e3α * g * yshift(interior(T′,:,:,k)); rasterize = true, colormap = bmap, colorrange = (rmin, rmax))
rmin, rmax = -0.2, 0.2
hm_c = heatmap!(ax_c, 1e-3x.-50, 1e-3y.-50, xhift(interior(u′,:,:,k)); rasterize = true, colormap = vmap, colorrange = (rmin, rmax))
hm_d = heatmap!(ax_d, 1e-3x, 1e-3y.-50, yshift(interior(u′,:,:,k)); rasterize = true, colormap = vmap, colorrange = (rmin, rmax))
hm_e = heatmap!(ax_e, 1e-3x.-50, 1e-3y.-50, xhift(interior(v′,:,:,k)); rasterize = true, colormap = vmap, colorrange = (rmin, rmax))
hm_f = heatmap!(ax_f, 1e-3x, 1e-3y.-50, yshift(interior(v′,:,:,k)); rasterize = true, colormap = vmap, colorrange = (rmin, rmax))
rmin, rmax = -0.01, 0.01
hm_g = heatmap!(ax_g, 1e-3x.-50, 1e-3y.-50, xhift(interior(w′,:,:,k)); rasterize = true, colormap = vmap, colorrange = (rmin, rmax))
hm_h = heatmap!(ax_h, 1e-3x, 1e-3y.-50, yshift(interior(w′,:,:,k)); rasterize = true, colormap = vmap, colorrange = (rmin, rmax))
Colorbar(gabc[1,3], hm_b)
Colorbar(gabc[2,3], hm_d)
Colorbar(gabc[3,3], hm_f)
Colorbar(gabc[4,3], hm_h)
hidexdecorations!(ax_a, ticks = false)
hidexdecorations!(ax_b, ticks = false)
hidexdecorations!(ax_c, ticks = false)
hidexdecorations!(ax_d, ticks = false)
hidexdecorations!(ax_e, ticks = false)
hidexdecorations!(ax_f, ticks = false)
Label(gabc[0, 1], L"\text{Unstable front},~z=-25.3~\text{m}", tellwidth = false)
Label(gabc[0, 2], L"\text{Stable front},~z=-25.3~\text{m}", tellwidth = false)
rowgap!(gabc, 3)
colgap!(gabc, 1, 10)
colgap!(gabc, 2, 1)
resize_to_layout!(fig)
save(filesave * "Tuvwpfields_" * fileparams * "_d6.pdf", fig; pt_per_unit = 1)

######################
x, y, z = nodes(Ti);
yrange1 = findfirst(p2[1][2]-ly .< 1e-3*y .- 100):findlast(1e-3*y .- 100 .<= p2[1][2]) 
yrange2 = findfirst(p2[2][2]-ly .< 1e-3*y .- 100):findlast(1e-3*y .- 100 .<= p2[2][2]) 
Tmap, vmap = :diff,:delta
zmin = -70
kz = findlast(z .< zmin)
Nz = length(z)
fig = Figure(size = (640, 750))
gabc = fig[1, 1] = GridLayout()
axis_kwargs = (aspect=8/5, ylabel = "z (m)", limits = ((pB[1], pB[1]+lx), (zmin,0)))
ax_a = Axis(gabc[1,1]; titlealign = :left, title=L"\text{(a)}~T-20~\text{({^\circ}C)}", axis_kwargs...)
ax_b = Axis(gabc[1,2]; titlealign = :left, title=L"\text{(b)}~T-20~\text{({^\circ}C)}", aspect=8/5, limits = ((pC[1], pC[1]+lx), (zmin,0))) 
ax_c = Axis(gabc[2,1]; titlealign = :left, title=L"\text{(c)}~u~\text{(m s^{-1})}", ylabel = "z (m)", aspect=8/5, limits = ((pB[1], pB[1]+lx), (zmin,0))) 
ax_d = Axis(gabc[2,2]; titlealign = :left, title=L"\text{(d)}~u~\text{(m s^{-1})}", aspect=8/5, limits = ((pC[1], pC[1]+lx), (zmin,0))) 
ax_e = Axis(gabc[3,1]; titlealign = :left, title=L"\text{(e)}~v~\text{(m s^{-1})}", ylabel = "z (m)", aspect=8/5, limits = ((pB[1], pB[1]+lx), (zmin,0))) 
ax_f = Axis(gabc[3,2]; titlealign = :left, title=L"\text{(f)}~v~\text{(m s^{-1})}", aspect=8/5, limits = ((pC[1], pC[1]+lx), (zmin,0))) 
ax_g = Axis(gabc[4,1]; titlealign = :left, title=L"\text{(g)}~w~\text{(10^{-3}m s^{-1})}", xlabel = "x (km)", ylabel = "z (m)", aspect=8/5, limits = ((pB[1], pB[1]+lx), (zmin,0))) 
ax_h = Axis(gabc[4,2]; titlealign = :left, title=L"\text{(h)}~w~\text{(10^{-3}m s^{-1})}", xlabel = "x (km)", aspect=8/5, limits = ((pC[1], pC[1]+lx), (zmin,0))) 
rmin, rmax = -0.5, 0.5
hm_a = heatmap!(ax_a, 1e-3x.-50, z[kz:Nz],shift(mean(interior(Ti,:,yrange1,kz:Nz),dims=2)[:,1,:]).-20; rasterize = true, colormap = Tmap, colorrange = (rmin, rmax))
hm_b = heatmap!(ax_b, 1e-3x, z[kz:Nz], (mean(interior(Ti,:,yrange2,kz:Nz),dims=2)[:,1,:]).-20; rasterize = true, colormap = Tmap, colorrange = (rmin, rmax))
rmin, rmax = -0.1, 0.1
hm_c = heatmap!(ax_c, 1e-3x.-50, z[kz:Nz], shift(mean(interior(ui,:,yrange1,kz:Nz),dims=2)[:,1,:]); rasterize = true, colormap = vmap, colorrange = (rmin, rmax))
hm_d = heatmap!(ax_d, 1e-3x, z[kz:Nz], (mean(interior(ui,:,yrange2,kz:Nz),dims=2)[:,1,:]); rasterize = true, colormap = vmap, colorrange = (rmin, rmax))
rmin, rmax = -0.2, 0.2
hm_e = heatmap!(ax_e, 1e-3x.-50, z[kz:Nz], shift(mean(interior(vi,:,yrange1,kz:Nz),dims=2)[:,1,:]); rasterize = true, colormap = vmap, colorrange = (rmin, rmax))
hm_f = heatmap!(ax_f, 1e-3x, z[kz:Nz], (mean(interior(vi,:,yrange2,kz:Nz),dims=2)[:,1,:]); rasterize = true, colormap = vmap, colorrange = (rmin, rmax))
rmin, rmax = -5, 5
hm_g = heatmap!(ax_g, 1e-3x.-50, z[kz:Nz], 1e3shift(mean(interior(wi,:,yrange1,kz:Nz),dims=2)[:,1,:]); rasterize = true, colormap = vmap, colorrange = (rmin, rmax))
hm_h = heatmap!(ax_h, 1e-3x, z[kz:Nz], 1e3(mean(interior(wi,:,yrange2,kz:Nz),dims=2)[:,1,:]); rasterize = true, colormap = vmap, colorrange = (rmin, rmax))

Colorbar(gabc[1,3], hm_b)
Colorbar(gabc[2,3], hm_d)
Colorbar(gabc[3,3], hm_f)
Colorbar(gabc[4,3], hm_h)
hidexdecorations!(ax_a, ticks = false)
hidexdecorations!(ax_b, ticks = false)
hidexdecorations!(ax_c, ticks = false)
hidexdecorations!(ax_d, ticks = false)
hidexdecorations!(ax_e, ticks = false)
hidexdecorations!(ax_f, ticks = false)
Label(gabc[0, 1], L"\text{Unstable front}", tellwidth = false)
Label(gabc[0, 2], L"\text{Stable front}", tellwidth = false)
rowgap!(gabc, 3)
colgap!(gabc, 1, 10)
colgap!(gabc, 2, 1)
resize_to_layout!(fig)
save(filesave * "Tuvwzfields_" * fileparams * "_d6.pdf", fig; pt_per_unit = 1)

######################
x, y, _ = nodes(h0);
cmap = :deep
fig = Figure(size = (640, 590))
gabc = fig[1, 1] = GridLayout()
axis_kwargs = (ylabel = "y (km)", aspect=8/5, limits = ((pB[1], pB[1]+lx), (pB[2]-ly,pB[2])))
ax_a = Axis(gabc[1,1]; titlealign = :left, title=L"\text{(a)~MLD}~h_D~\text{(m)}", axis_kwargs...)
ax_b = Axis(gabc[1,2]; titlealign = :left, title=L"\text{(b)~MLD}~h_D~\text{(m)}", aspect=8/5, limits = ((pC[1], pC[1]+lx), (pC[2]-ly, pC[2])))
ax_c = Axis(gabc[2,1]; titlealign = :left, title=L"\text{(c)~BLD}~h_K~\text{(m)}", ylabel = "y (km)", aspect=8/5, limits = ((pB[1], pB[1]+lx), (pB[2]-ly,pB[2]))) 
ax_d = Axis(gabc[2,2]; titlealign = :left, title=L"\text{(d)~BLD}~h_K~\text{(m)}", aspect=8/5, limits = ((pC[1], pC[1]+lx), (pC[2]-ly, pC[2]))) 
ax_e = Axis(gabc[3,1]; titlealign = :left, title=L"\text{(e)}~h_D-h_K~\text{(m)}", xlabel = "x (km)",ylabel = "y (km)", aspect=8/5, limits = ((pB[1], pB[1]+lx), (pB[2]-ly,pB[2]))) 
ax_f = Axis(gabc[3,2]; titlealign = :left, title=L"\text{(f)}~h_D-h_K~\text{(m)}", xlabel = "x (km)",aspect=8/5, limits = ((pC[1], pC[1]+lx), (pC[2]-ly, pC[2]))) 
rmin, rmax = 50, 70
hm_a = heatmap!(ax_a, 1e-3x.-50, 1e-3y.-50, xhift(interior(h0,:,:,1)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
hm_b = heatmap!(ax_b, 1e-3x, 1e-3y.-50, yshift(interior(h0,:,:,1)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
rmin, rmax = 20, 60
hm_c = heatmap!(ax_c, 1e-3x.-50, 1e-3y.-50, xhift(interior(h1,:,:,1)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
hm_d = heatmap!(ax_d, 1e-3x, 1e-3y.-50, yshift(interior(h1,:,:,1)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
rmin, rmax = -30, 30
cmap = :diff
hm_e = heatmap!(ax_e, 1e-3x.-50, 1e-3y.-50, xhift(interior(h0,:,:,1).-interior(h1,:,:,1)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
hm_f = heatmap!(ax_f, 1e-3x, 1e-3y.-50, yshift(interior(h0,:,:,1).-interior(h1,:,:,1)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
Colorbar(gabc[1,3], hm_b)
Colorbar(gabc[2,3], hm_d)
Colorbar(gabc[3,3], hm_f)
hidexdecorations!(ax_a, ticks = false)
hidexdecorations!(ax_b, ticks = false)
hidexdecorations!(ax_c, ticks = false)
hidexdecorations!(ax_d, ticks = false)
Label(gabc[0, 1], L"\text{Unstable front}", tellwidth = false)
Label(gabc[0, 2], L"\text{Stable front}", tellwidth = false)
rowgap!(gabc, 3)
colgap!(gabc, 1, 15)
colgap!(gabc, 2, 3)
resize_to_layout!(fig)
save(filesave * "MLDBLDfields_" * fileparams * "_d6.pdf", fig; pt_per_unit = 1)

########################
x, y, z = nodes(wi);
Nx, Ny, Nz = length(x), length(y), length(z)
Tfi = compute!(Field(@at (Center, Center, Face) snapshots[:T][i]));
fvz = compute!(Field(@at (Center, Center, Face) f * ∂z(vi)));
uzz = compute!(Field(@at (Center, Center, Face) ∂z(∂z(κi * ∂z(ui)))));

# coarse-graining
u̅i,v̅i = XFaceField(ui.grid),YFaceField(vi.grid);
T̅i = CenterField(Ti.grid);
cutoff = 300/2.4*2*π
coarse_graining!(ui, u̅i; kernel, cutoff)
coarse_graining!(vi, v̅i; kernel, cutoff)
coarse_graining!(Ti, T̅i; kernel, cutoff)
R̅i = compute!(Field(α * g * ∂z(T̅i)/(∂z(u̅i)^2 + ∂z(v̅i)^2)));
θ, τ, ρ₀=parameters.θ,parameters.τw,parameters.ρ₀
τ̅x,τ̅y=-τ*sind(θ),-τ*cosd(θ)
EBFτ̅ = compute!(Field(@at (Center, Center, Face) α * g * (τ̅y * ∂x(T̅i) - τ̅x * ∂y(T̅i))/f/ρ₀))

# Current feedback on stress
ρₐ = 1.225; # kg/m³
cd = 0.001; # Drag coefficient
Uₐ = sqrt(τ / (ρₐ * cd)); # Wind speed at the surface
Uₐx, Uₐy = Uₐ * cosd(θ), Uₐ * sind(θ) # Wind speed at the surface
Unorm = compute!(Field(sqrt((Uₐx + u̅i)^2 + (Uₐy + v̅i)^2))); # Wind speed at the surface
τx = compute!(Field(- ρₐ * cd * Unorm * (Uₐx + u̅i)));
τy = compute!(Field(- ρₐ * cd * Unorm * (Uₐy + v̅i)));
EBFτ′ = compute!(Field(-EBFτ̅ + α * g * (τy * ∂x(T̅i) - τx * ∂y(T̅i))/f/ρ₀))

shift(x) = [x[3size(x,1)÷4+1:end, :]; x[1:3size(x,1)÷4, :]]
xhift(x) = yshift(shift(x))
pA, pB, pC = (-45, 70), (-20, 0), (30, -25)
p3 = [pA, pB, pC]
alphabet = [letter for letter in 'a':'z'];
fig1 = Figure(size = (640, 600))
g1 = fig1[1, 1] = GridLayout()
fig2 = Figure(size = (640, 600))
g2 = fig2[1, 1] = GridLayout()
regiontitle = ["Unstable front", "Stable front"]
ufi = compute!(Field(@at (Center, Center, Face) ui))
for j = 2:3
    yrange = findfirst(p3[j][2]-ly .< 1e-3*y .- 50):findlast(1e-3*y .- 50 .<= p3[j][2]) 
    wij1,Tfij1 = xhift(interior(wi, :, :, 1))[:,yrange],xhift(interior(Tfi, :, :, 1))[:,yrange]
    hij = xhift(interior(hi, :, :, 1))[:,yrange]
    h̅ij = mean(hij, dims = 2);
    Nxij,Nyij = size(wij1)
    wij, Tfij = zeros(Nxij,Nyij,length(z)),zeros(Nxij,Nyij,length(z))
    vij = zeros(Nxij,Nyij,length(z)-1)
    ufij = zeros(Nxij,Nyij,length(z))
    fvzij = zeros(Nxij,Nyij,length(z))
    uzzij = zeros(Nxij,Nyij,length(z))
    wij[:, :, 1] = wij1;
    Tfij[:, :, 1] = Tfij1;
    ufij[:, :, 1] = xhift(interior(ufi, :, :, 1))[:,yrange];
    fvzij[:, :, 1] = xhift(interior(fvz, :, :, 1))[:,yrange];
    uzzij[:, :, 1] = xhift(interior(uzz, :, :, 1))[:,yrange];
    for k in 2:length(z)
        wij[:, :, k] = xhift(interior(wi, :, :, k))[:,yrange];
        Tfij[:, :, k] = xhift(interior(Tfi, :, :, k))[:,yrange];
        vij[:, :, k-1] = xhift(interior(vi, :, :, k-1))[:,yrange];
        ufij[:, :, k] = xhift(interior(ufi, :, :, k))[:,yrange];
        fvzij[:, :, k] = xhift(interior(fvz, :, :, k))[:,yrange];
        uzzij[:, :, k] = xhift(interior(uzz, :, :, k))[:,yrange];
    end
    w̅ij = mean(wij, dims = 2);
    T̅fij = mean(Tfij, dims = 2);
    u̅fij = mean(ufij, dims = 2);
    fvzij = mean(fvzij, dims = 2)[:,1,:];
    uzzij = mean(uzzij, dims = 2)[:,1,:];
    v̅ij = Field{Center, Face, Center}(wi.grid);
    set!(v̅ij, repeat(reshape(mean(vij, dims = 2), Nx, 1, Nz-1),1,Ny,1));
    fill_halo_regions!(v̅ij)
    b̅ij = Field{Center, Center, Face}(wi.grid);
    set!(b̅ij, α * g * repeat(reshape(T̅fij, Nx, 1, Nz),1,Ny,1))
    fill_halo_regions!(b̅ij)
    wbij = Field{Center, Nothing, Face}(wi.grid);
    set!(wbij, α * g * mean((wij .- w̅ij) .* (Tfij .- T̅fij), dims = 2))
    fill_halo_regions!(wbij)
    wuij = Field{Center, Center, Face}(wi.grid);
    set!(wuij, repeat(reshape(mean((wij .- w̅ij) .* (ufij .- u̅fij), dims = 2),Nx,1,Nz),1,Ny,1));
    fill_halo_regions!(wuij)
    wuzij = mean(compute!(Field(∂z(∂z(wuij)))), dims = 2);
    b̅ₓij = compute!(Field(@at (Center, Nothing, Face) mean(compute!(Field(∂x(b̅ij))), dims = 2)));
    b̅zij = compute!(Field(@at (Center, Nothing, Face) mean(compute!(Field(∂z(b̅ij))), dims = 2)));
    #Ψij = compute!(Field(wbij/b̅ₓij))
    # Nsq = mean(interior(b̅zij, xrange, 1, z .> -H))
    Ri = compute!(Field(f^2 * b̅zij / b̅ₓij^2))
    # MLI2 = 0.09H^2 * M⁴ / f * (1 .- (2z./H .+ 1).^2) ./ sqrt(1 + Ri)
    # lines!(ax3, 1e7MLI2, z)
    ωz = mean(compute!(Field(∂x(v̅ij))),dims=2);
    ωx = mean(compute!(Field(-∂z(v̅ij))),dims=2);
    PV = compute!(Field(b̅ₓij * ωx + (ωz + f) * b̅zij));

    xrange = findfirst(p3[j][1]+lx/4 .< 1e-3*x .- 25):findlast(1e-3*x .- 25 .<= p3[j][1]+3lx/4)
    ifront = zeros(Nxij,length(z))
    imin,imax = Nxij,1
    for k = 2:length(z)
        if z[k] .> -min(h̅ij...)
            b̅ₓks = interior(b̅ₓij, xrange, 1, k)*sign(j-2.5)
            b̅ₓmax = maximum(b̅ₓks)
            i1,i2 = findfirst(b̅ₓks .> 0.03*b̅ₓmax)+xrange[1]-1,findlast(b̅ₓks .> 0.03*b̅ₓmax)+xrange[1]-1
            ifront[i1:i2,k] .= 1
            imin, imax = min(imin, i1), max(imax, i2)
        end
    end

    title1 = "("*alphabet[3*(j-2)+1]*") " * regiontitle[j-1]
    title2 = "("*alphabet[3*(j-2)+2]*") " * regiontitle[j-1]
    title3 = "("*alphabet[3*(j-2)+3]*") " * regiontitle[j-1]
    axis_kwargs = (titlealign = :left, 
                   titlefont=texfont(), 
                   limits = ((p3[j][1], p3[j][1]+lx), (-100, 0)))
    # g1
    ax1 = Axis(g1[j-2,1]; title = title1, axis_kwargs...)
    ax2 = Axis(g1[j-2,2]; title = title2, axis_kwargs...)
    ax3 = Axis(g1[j-2,3]; title = title3, titlealign = :left, titlefont=texfont(),limits = (nothing, (-100, 0)))
    hideydecorations!(ax2, ticks = false)
    hideydecorations!(ax3, ticks = false)
    hm1 = heatmap!(ax1, (1e-3*x .- 25), z, 1e7interior(wbij, :, 1, :); rasterize = true, colormap = :balance, colorrange = (-1,1))
    hm2 = heatmap!(ax2, (1e-3*x .- 25), z, 1e7interior(b̅ₓij, :, 1, :); rasterize = true, colormap = :balance, colorrange = (-3,3))
    heatmap!(ax1, (1e-3*x .- 25), z, ifront; rasterize = true, colormap = (:grays,0.3), colorrange = (0,1))
    heatmap!(ax2, (1e-3*x .- 25), z, ifront; rasterize = true, colormap = (:grays,0.3), colorrange = (0,1))
    lines!(ax1, 1e-3x.-25, -vec(h̅ij); color = :black, linestyle = :dash, linewidth = 0.8)
    lines!(ax2, 1e-3x.-25, -vec(h̅ij); color = :black, linestyle = :dash, linewidth = 0.8)
    #vlines!(ax1, [10,30].+p3[j][1], -100, 0, color = :black, linewidth = 0.8)
    #vlines!(ax2, [10,30].+p3[j][1], -100, 0, color = :black, linewidth = 0.8)
 
    sfront = sum(ifront,dims=1)
    lines!(ax3, 1e7vec(sum(interior(wbij, :, 1, :).*ifront,dims=1)./sfront), z; label = L"\overline{w^\prime b^\prime}^{xy}")

    H = mean(h̅ij[imin:imax])
    μ = max.(0,(1 .- (2z./H .+ 1).^2) .* (1 .+ 5/21 * (2z./H .+ 1).^2))
    M⁴ = (sum(interior(b̅ₓij, :, 1, :).*ifront)/sum(ifront)).^2
    MLI1 = M⁴ * H^2 / f * μ
    band!(ax3, Point2f.(0.06*1e7MLI1, z), Point2f.(0.08*1e7MLI1, z), color = (:red, 0.5); label = "MLI")

    axislegend(ax3,labelsize=9, framevisible = false, position = :rb, padding = (0f0, 0f0, 0f0, 0f0), patchlabelgap = 3)
    ax1.ylabel = L"z~\text{(m)}"
    if j > 2
        Colorbar(g1[2,1], hm1, vertical = false, flipaxis = false, label = L"\overline{w^\prime b^\prime}~\text{(10^{-7}m^2 s^{-3})}")
        Colorbar(g1[2,2], hm2, vertical = false, flipaxis = false, label = L"\overline{b}_x~\text{(10^{-7} s^{-2})}")
        ax1.xlabel = L"x~\text{(km)}"
        ax2.xlabel = L"x~\text{(km)}"
        ax3.xlabel = L"\text{(10^{-7}m^2 s^{-3})}"
    end


    #g2
    T̅Wij = reshape(interior(b̅ₓij, :, 1, :),640,225) .- fvzij;
    T̅TWij = T̅Wij .- uzzij .+ reshape(interior(wuzij, :, 1, :),640,225);
    println(size(T̅Wij))
    ax1 = Axis(g2[j-2,1]; title = title1, axis_kwargs...)
    ax2 = Axis(g2[j-2,2]; title = title2, axis_kwargs...)
    ax3 = Axis(g2[j-2,3]; title = title3, titlealign = :left, xlabel = L"M_0^2",titlefont=texfont(),limits = ((-1,1), (-100, 0)))
    hideydecorations!(ax2, ticks = false)
    hideydecorations!(ax3, ticks = false)
    hm1 = heatmap!(ax1, (1e-3*x .- 25), z, T̅Wij./M²₀; rasterize = true, colormap = :balance, colorrange = (-1,1))
    hm2 = heatmap!(ax2, (1e-3*x .- 25), z, T̅TWij./M²₀; rasterize = true, colormap = :balance, colorrange = (-1,1))
    heatmap!(ax1, (1e-3*x .- 25), z, ifront; rasterize = true, colormap = (:grays,0.3), colorrange = (0,1))
    heatmap!(ax2, (1e-3*x .- 25), z, ifront; rasterize = true, colormap = (:grays,0.3), colorrange = (0,1))
    lines!(ax1, 1e-3x.-25, -vec(h̅ij); color = :black, linestyle = :dash, linewidth = 0.8)
    lines!(ax2, 1e-3x.-25, -vec(h̅ij); color = :black, linestyle = :dash, linewidth = 0.8)
    #vlines!(ax1, [10,30].+p3[j][1], -100, 0, color = :black, linewidth = 0.8)
    #vlines!(ax2, [10,30].+p3[j][1], -100, 0, color = :black, linewidth = 0.8)
    lines!(ax3, vec(sum(reshape(interior(b̅ₓij, :, 1, :),640,225).*ifront,dims=1)./sfront)./M²₀, z; label = L"\overline{b}_x^{xy}")
    lines!(ax3, vec(sum(fvzij.*ifront,dims=1)./sfront)./M²₀, z; label = L"f \overline{v}_z^{xy}")
    lines!(ax3, vec(sum(uzzij.*ifront,dims=1)./sfront)./M²₀, z; label = L"\overline{\kappa u_z}_{zz}^{xy}")
    lines!(ax3, vec(sum(reshape(interior(wuzij, :, 1, :),640,225).*ifront,dims=1)./sfront)./M²₀, z; label = L"\overline{w\prime u\prime}_{zz}^{xy}")

    axislegend(ax3,labelsize=9, framevisible = false, position = :rb, padding = (0f0, 0f0, 0f0, 0f0), patchlabelgap = 3)
    ax1.ylabel = L"z~\text{(m)}"
    if j > 2
        Colorbar(g2[2,1], hm1, vertical = false, flipaxis = false, label = L"\text{TW}/M_0^2")
        Colorbar(g2[2,2], hm2, vertical = false, flipaxis = false, label = L"\text{TTW}/M_0^2")
        ax1.xlabel = L"x~\text{(km)}"
        ax2.xlabel = L"x~\text{(km)}"
    end
end
rowgap!(g1, 5)
colsize!(g1, 3, Relative(0.2))
resize_to_layout!(fig1)
save(filesave * "MLIzfields_" * fileparams * "_d6.pdf", fig1; pt_per_unit = 1)
rowgap!(g2, 5)
colsize!(g2, 3, Relative(0.2))
resize_to_layout!(fig2)
save(filesave * "TTWzfields_" * fileparams * "_d6.pdf", fig2; pt_per_unit = 1)

#######################################
w̅i = ZFaceField(wi.grid); # compute!(Field(wi - w̅i));
T̅i = CenterField(Ti.grid);# compute!(Field(Ti - T̅i));
coarse_graining!(wi, w̅i; kernel, cutoff = ls)
coarse_graining!(Ti, T̅i; kernel, cutoff = ls)
∇b̅i = compute!(Field(α * g * (∂x(T̅i)^2 + ∂y(T̅i)^2)^0.5));
wˢ = compute!(Field(wi - w̅i));
Tˢ = compute!(Field(Ti - T̅i));
wˢbˢ = compute!(Field(α * g * Tˢ * wˢ));
w̅ˢ = ZFaceField(wi.grid); 
T̅ˢ = CenterField(Ti.grid);
coarse_graining!(wˢ, w̅ˢ; cutoff)
coarse_graining!(Tˢ, T̅ˢ; cutoff)
vbfi = compute!(Field(α * g * T̅ˢ * w̅ˢ));

######################################
# sound speed as a function of in situ temperature and depth
ρ₀ = parameters.ρ₀;
Ω = 7.2921e-5;
θ = asind(f/2/Ω)
_,_,z = nodes(snapshots[:T][1]);
p = gsw_p_from_z.(z, θ);
sa = 35;
i = 161;
Ti = snapshots[:T][i];
c = CenterField(snapshots[:T][1].grid);
set!(c,gsw_sound_speed_t_exact.(sa, interior(Ti, :, :, :), reshape(p, 1, 1, :)));
fill_halo_regions!(c)
c′ = compute!(Field(c - mean(c, dims = (1, 2))));

#####################
kwbs = [222, 202, 171];
T1 = snapshots[:T][3];
x, y, z = nodes(T1);
Δx, Δy = x[2]-x[1], y[2]-y[1];
Nx, Ny = length(x), length(y);
w1 = compute!(Field(@at (Center, Center, Center) snapshots[:w][3]));
b1 = compute!(Field(α * g * T1));
C1 = isotropic_powerspectrum(interior(b1, :, :, kwbs[1]), interior(w1, :, :, kwbs[1]), x, y);
freq = C1.freq
idx = freq .> 0
tsSt = zeros(length(times)÷2+1, length(kwbs), length(freq))
for i in 1:length(times)÷2+1
    println("Computing spectra at time $(times[2i+1]/24/3600) days...")
    for j in 1:3
        bij = α * g * interior(snapshots[:T][2i+1], :, :, kwbs[j])
        wij = interior(compute!(Field(@at (Center, Center, Center) snapshots[:w][2i+1])), :, :, kwbs[j])
        tsSt[i,j,:] = real.(isotropic_powerspectrum(bij, wij, x, y).spec)
    end
end

fig = Figure(size = (640, 750))
gl = fig[1, 1] = GridLayout()
titles = [L"\text{(a)}~k\hat{w}\hat{b}~\text{(m^2 s^{-3})},~z=-2.8~\text{m}", L"\text{(b)}~k\hat{w}\hat{b}~\text{(m^2 s^{-3})},~z=-25~\text{m}", L"\text{(c)}~k\hat{w}\hat{b}~\text{(m^2 s^{-3})},~z=-60~\text{m}"]
for i = 1:3
    ax = Axis(gl[i, 1]; ylabel = L"\text{Wavenumber}~k~\text{(m^{-1})}", yscale = log10,
                        titlealign = :left, title=titles[i],
                        limits = ((0.125,20),(6e-5, 3e-2)))
    Zi = tsSt[2:end,i,idx] .* reshape(freq[idx], (1,:)) 
    vbound = maximum(abs.(Zi))/2.
    hm = heatmap!(ax, times[3:2:end]/24/3600, freq[idx], Zi; rasterize = true, colorrange = (-vbound, vbound), colormap = :balance)
    Colorbar(gl[i, 2], hm)
    ylims!(ax, (6.5e-5, 2.5e-2))
    if i < 3
        hidexdecorations!(ax, ticks = false)
    else
        ax.xlabel = L"\text{Time}~t~\text{(days)}"
    end
end
colgap!(gl, 1, 3)
for i in 1:2
    rowgap!(gl, i, 3)
end
resize_to_layout!(fig)
save(filesave * "whbhts_" * fileparams * ".pdf", fig; pt_per_unit = 1)

###############################
# Compute spectral vertical boyancy flux 
Nz = length(z)
snapshot_number = 161
C1d = zeros(9, Nz, length(C1.spec));
for (n,i) in enumerate(-8:2:8)
    for k = 1:Nz
        if z[k] > -200
            bik = α * g * interior(snapshots[:T][snapshot_number+i], :, :, k)
            wik = interior(compute!(Field(@at (Center, Center, Center) snapshots[:w][snapshot_number+i])), :, :, k)
            C1d[n,k,:] = real.(isotropic_powerspectrum(bik, wik, x, y).spec)
            println("Spectral done at z = $(z[k]) m, snapshot number = $(snapshot_number+i)...")
        end
    end
end

Z1d = C1d .* reshape(C1.freq, (1,1,:));
fig = Figure(size = (640, 450))
g4 = fig[1, 1] = GridLayout()
axis_kwargs = (ylabel = L"z~\text{(m)}", xlabel = L"\text{Wavenumber (m^{-1})}",xscale = log10, ygridvisible = false, 
               limits = ((6e-5, 4e-2), (-150,0)),xticks = ([1e-4,1e-3,1e-2], [L"10^{-4}",L"10^{-3}",L"10^{-2}"]), xgridvisible = false,
               xminorticks = [4e-5:1e-5:9e-5; 2e-4:1e-4:9e-4; 2e-3:1e-3:9e-3; 2e-2:1e-2:5e-2],xminorticksvisible = true)
ax_a = Axis(g4[1,1]; titlealign = :left, title=L"\text{(a)}~k\hat{w}\hat{b}~\text{(m^2 s^{-3})}", axis_kwargs...)
ax_b = Axis(g4[1,2]; titlealign = :left, title=L"\text{(b) 1d average around day 10}", axis_kwargs...)
ax_c = Axis(g4[1,3]; titlealign = :left, title=L"\text{(c) submesoscale flux}", ylabel = L"z~\text{(m)}", xlabel = L"\langle w^s b^s \rangle \text{(10^{-9} m^2 s^{-3})}")
hideydecorations!(ax_b, ticks = false)
hideydecorations!(ax_c, ticks = false)
vbound = maximum(abs.(Z1d[5,:,:]))
hm_a = heatmap!(ax_a, C1.freq, z, Z1d[5,:,:]'; rasterize = true, colorrange = (-vbound, vbound), colormap = :balance)
hm_b = heatmap!(ax_b, C1.freq, z, mean(Z1d, dims = 1)[1,:,:]'; rasterize = true, colorrange = (-vbound, vbound), colormap = :balance)
vlines!(ax_a, 4e-4; color = :black, linewidth = 0.8)
vlines!(ax_b, 4e-4; color = :black, linewidth = 0.8)
vlines!(ax_a, 4e-3; color = :black, linewidth = 0.8)
vlines!(ax_b, 4e-3; color = :black, linewidth = 0.8)
Colorbar(g4[2, 1:2], hm_a, vertical = false)
for i = 1:2
    colgap!(g4, i, 5)
end
idx = 4e-4 .<= C1.freq .<= 4e-3
Δf = mean(diff(C1.freq[idx]))
lines!(ax_c, 1e9vec(sum(C1d[5,:,idx], dims = 2) .* Δf), z; linewidth = 1, label = "day 10")
lines!(ax_c, 1e9vec(mean(sum(C1d[:,:,idx], dims = 3), dims = 1) .* Δf), z; linewidth = 1, label = "1d average")
axislegend(ax_c,labelsize=9, framevisible = false, position = :rb, padding = (0f0, 0f0, 0f0, 0f0), patchlabelgap = 3)
resize_to_layout!(fig)
save(filesave * "whbhz_" * fileparams * "_d10.pdf", fig; pt_per_unit = 1)

################################

C[1, :] = real.(C1.spec)
# wc1 = (xhift(interior(wc, :, :, 1)))
# b1 = (xhift(interior(b, :, :, 1)))
# SCs = []
# for j = 1:3
#     xrange = findfirst(p3[j][1] .< 1e-3*xT .- 50):findlast(1e-3*xT .- 50 .<= p3[j][1]+l3[j]) 
#     yrange = findfirst(p3[j][2]-l3[j] .< 1e-3*yT .- 25):findlast(1e-3*yT .- 25 .<= p3[j][2]) 
#     SCj = isotropic_powerspectrum(b1[xrange,yrange], wc1[xrange,yrange], xT[xrange], yT[yrange];window=1)
#     push!(SCs, SCj)
# end
# Csub = zeros(3, Nz, length(SCs[1].spec))
# for j = 1:3
#     Csub[j, 1, :] = real.(SCs[j].spec)
# end
println("level 1 done.")
for k = 2:Nz
    Ck = isotropic_powerspectrum(interior(b, :, :, k), interior(wc, :, :, k), xT, yT)
    C[k,:] = real.(Ck.spec)
    # wck = (xhift(interior(wc, :, :, k)))
    # bk = (xhift(interior(b, :, :, k)))
    # for j = 1:3
    #     xrange = findfirst(p3[j][1] .< 1e-3*xT .- 50):findlast(1e-3*xT .- 50 .<= p3[j][1]+l3[j]) 
    #     yrange = findfirst(p3[j][2]-l3[j] .< 1e-3*yT .- 25):findlast(1e-3*yT .- 25 .<= p3[j][2]) 
    #     SCj = isotropic_powerspectrum(bk[xrange,yrange], wck[xrange,yrange], xT[xrange], yT[yrange];window=1)
    #     Csub[j, k, :] = real.(SCj.spec) 
    # end
    println("level $k done.")
end

fig = Figure(size = (640, 300))
g4 = fig[1, 1] = GridLayout()
axis_kwargs = (ylabel = L"z~\text{(m)}", xlabel = L"\text{Wavenumber (rad m^{-1})}",xscale = log10, ygridvisible = false, 
               limits = ((6e-5, 0.8e-2), (-150,0)),xticks = ([1e-4,1e-3], [L"10^{-4}",L"10^{-3}"]), xgridvisible = false,
               xminorticks = [4e-5:1e-5:9e-5; 2e-4:1e-4:9e-4; 2e-3:1e-3:9e-3],xminorticksvisible = true)
ax_a = Axis(g4[1,1]; titlealign = :left, title=L"\text{(a)}~\hat{w}\hat{b}~\text{(10^5 m^3 s^{-3})}", axis_kwargs...)
ax_b = Axis(g4[1,2]; titlealign = :left, title=L"\text{(b) Region A,}~5\times", axis_kwargs...)
ax_c = Axis(g4[1,3]; titlealign = :left, title=L"\text{(c) Region B,}~5\times", axis_kwargs...)
ax_d = Axis(g4[1,4]; titlealign = :left, title=L"\text{(d) Region C,}~5\times", axis_kwargs...)
hideydecorations!(ax_b, ticks = false)
hideydecorations!(ax_c, ticks = false)
hideydecorations!(ax_d, ticks = false)
hm_a = heatmap!(ax_a, C1.freq, zT, 1e-5*C'; rasterize = true, colormap = :balance, colorrange = (-2,2))
hm_b = heatmap!(ax_b, SCs[1].freq, zT, 5e-5*Csub[1,:,:]'; rasterize = true, colormap = :balance, colorrange = (-2,2))
hm_c = heatmap!(ax_c, SCs[2].freq, zT, 5e-5*Csub[2,:,:]'; rasterize = true, colormap = :balance, colorrange = (-2,2))
hm_d = heatmap!(ax_d, SCs[3].freq, zT, 5e-5*Csub[3,:,:]'; rasterize = true, colormap = :balance, colorrange = (-2,2))
Colorbar(g4[1, 5], hm_a)
for i = 1:3
    colgap!(g4, i, 5)
end
colgap!(g4, 4, 1)
resize_to_layout!(fig)
save(filesave * "spectral_vertical_boyancy_flux_d$(nday).pdf", fig; pt_per_unit = 1)

#################

vbf̅i  = CenterField(Ti.grid);
coarse_graining!(vbfi , vbf̅i; kernel, cutoff = lc)
vbf̅iz = MLaverage(snapshots,i,vbf̅i; kernel, scale=lc);
MLIi = compute!(Field(MLI(snapshots,i; kernel, scale=ls)));
#M̅i = CenterField(Ti.grid);
#coarse_graining!(MLi , M̅i ; cutoff = lc)

# Plot the fields
kw = 202
lA, lB, lC = 40, 40, 40
pA, pB, pC = (-45, 70), (-20, 15), (10, 70)
ΔN = 160
jslices = [2]*ΔN
Nsq = compute!(Field(N²(snapshots,i)));
#####################################
var,scale = vbfi,1e7;
x, y, z = nodes(var)
k = kw
cmap = :balance
rmin, rmax = -5, 5
hcolor, scolor = :gray, :black
h = compute!(MLD(snapshots,i; threshold = 0.1))
fig = Figure(size = (640, 785))
gabc = fig[1, 1] = GridLayout()
axis_kwargs = (ylabel = "y (km)", aspect=1, limits = ((-50, 50), (0, 100)))
ax_a = Axis(gabc[1,1]; titlealign = :left, title=L"\text{(a)}~\langle{w^s}\rangle \langle{b^s}\rangle~\text{(10^{-7} m^2 s^{-3})},~z=-25~\text{m}", axis_kwargs...)
ax_b = Axis(gabc[1,2]; titlealign = :left, title=L"\text{(b)~Region A}", aspect=1, limits = ((pA[1], pA[1]+lA), (pA[2]-lA, pA[2])))
ax_c = Axis(gabc[3,1]; titlealign = :left, title=L"\text{(c)~Region B}", ylabel = "y (km)", aspect=1, limits = ((pB[1], pB[1]+lB), (pB[2]-lB,pB[2]))) 
ax_d = Axis(gabc[3,2]; titlealign = :left, title=L"\text{(d)~Region C}", aspect=1, limits = ((pC[1], pC[1]+lC), (pC[2]-lC, pC[2]))) 
hm_a = heatmap!(ax_a, 1e-3x.-50, 1e-3y, scale * shift(interior(var,:,:,k)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
hm_b = heatmap!(ax_b, 1e-3x.-50, 1e-3y, scale * shift(interior(var,:,:,k)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
hm_c = heatmap!(ax_c, 1e-3x.-50, 1e-3y.-25, scale * xhift(interior(var,:,:,k)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
hm_d = heatmap!(ax_d, 1e-3x.-50, 1e-3y, scale * shift(interior(var,:,:,k)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
Colorbar(gabc[1,3], hm_b)
Colorbar(gabc[3,3], hm_d)
hidexdecorations!(ax_a, ticks = false)
hidexdecorations!(ax_b, ticks = false)
hidexdecorations!(ax_c, ticks = false)
hidexdecorations!(ax_d, ticks = false)
hlines!(ax_a, 1e-3*y[jslices]; color = :black, linestyle = :dash, linewidth = 0.8)
poly!(ax_a, Rect(pA[1], pA[2]-lA, lA, lA), color = (:white, 0.1), strokecolor = scolor, strokewidth = 0.5)
poly!(ax_a, Rect(pC[1], pC[2]-lC, lC, lC), color = (:white, 0.1), strokecolor = scolor, strokewidth = 0.5)
poly!(ax_a, Rect(pB[1], 75, lB, lB-pB[2]), color = (:white, 0.1), strokewidth = 0.)
poly!(ax_a, Rect(pB[1], 0, lB, pB[2]), color = (:white, 0.1), strokewidth = 0.)
vlines!(ax_a, [pB[1], pB[1]+lB]; ymin = 0.75, color = scolor, linewidth = 0.5)
vlines!(ax_a, [pB[1], pB[1]+lB]; ymax = 0.15, color = scolor, linewidth = 0.5)
hlines!(ax_a, [pB[2], 75]; xmin = 0.3, xmax = 0.7, color = scolor, linewidth = 0.5)
text!(ax_a, pA[1], pA[2], text = L"\text{A}", color = :black, align = (:left, :top))
text!(ax_a, pB[1], pB[2], text = L"\text{B}", color = :black, align = (:left, :top))
if cmap == :thermal
    text!(ax_a, pC[1], pC[2], text = L"\text{C}", color = :black, align = (:left, :top))
else
    text!(ax_a, pC[1], pC[2], text = L"\text{C}", color = :black, align = (:left, :top))
end

dΔN = 80
jslices = dΔN * [0] .+ 2ΔN
hlines!(ax_b, 1e-3*y[jslices].+6.25; color = :black, linestyle = :dash, linewidth = 0.8)
hlines!(ax_d, 1e-3*y[jslices].+6.25; color = :black, linestyle = :dash, linewidth = 0.8)
hlines!(ax_c, -10; color = :black, linestyle = :dash, linewidth = 0.8)

zmin = -100
ΔN = 160
ja = 2ΔN
jc = 9*64
Δj = 40
kz = findlast(z .< zmin)
Nz = length(z)
axis_kwargs0 = (xlabel = "x (km)", ylabel = "z (m)", limits = ((-50, 50), (zmin, 0)))
axis_kwargs1 = NamedTuple{(:xlabel,:ylabel)}(axis_kwargs0)
ax_a = Axis(gabc[2,1]; titlealign = :left, title=L"y=50~\text{km}", axis_kwargs0...)
ax_b = Axis(gabc[2,2]; titlealign = :left, title=L"y=56~\text{km}", xlabel = "x (km)", limits = ((pA[1], pA[1]+lA), (zmin, 0)))
ax_c = Axis(gabc[4,1]; titlealign = :left, title=L"y=-10~\text{km}", limits = ((pB[1], pB[1]+lB), (zmin, 0)), axis_kwargs1...)
ax_d = Axis(gabc[4,2]; titlealign = :left, title=L"y=56~\text{km}", xlabel = "x (km)", limits = ((pC[1], pC[1]+lC), (zmin, 0)))
if cmap == :thermal
    vmin, vmax = minimum(interior(var,:,:,kz:Nz)), maximum(interior(var,:,:,kz:Nz))
else
    vmin, vmax = rmin, rmax
end
hm_a = heatmap!(ax_a, 1e-3x.-50, z[kz:Nz], scale*shift(interior(var,:,ja,kz:Nz)); rasterize = true, colormap = cmap, colorrange = (vmin, vmax))
hm_b = heatmap!(ax_b, 1e-3x.-50, z[kz:Nz], scale*shift(interior(var,:,ja+Δj,kz:Nz)); rasterize = true, colormap = cmap, colorrange = (vmin, vmax))
hm_c = heatmap!(ax_c, 1e-3x.-50, z[kz:Nz], scale*shift(interior(var,:,jc,kz:Nz)); rasterize = true, colormap = cmap, colorrange = (vmin, vmax))
hm_d = heatmap!(ax_d, 1e-3x.-50, z[kz:Nz], scale*shift(interior(var,:,ja+Δj,kz:Nz)); rasterize = true, colormap = cmap, colorrange = (vmin, vmax))
hideydecorations!(ax_b, ticks = false)
hideydecorations!(ax_d, ticks = false)
Colorbar(gabc[2,3], hm_b)
Colorbar(gabc[4,3], hm_d)
hlines!(ax_a, z[k]; color = :black, linestyle = :dash, linewidth = 0.8)
hlines!(ax_b, z[k]; color = :black, linestyle = :dash, linewidth = 0.8)
hlines!(ax_c, z[k]; color = :black, linestyle = :dash, linewidth = 0.8)
hlines!(ax_d, z[k]; color = :black, linestyle = :dash, linewidth = 0.8)
lines!(ax_a, 1e-3x.-50, -vec(shift(interior(h,:,ja))); color = hcolor, linewidth = 0.8)
lines!(ax_b, 1e-3x.-50, -vec(shift(interior(h,:,ja+Δj))); color = hcolor, linewidth = 0.8)
lines!(ax_c, 1e-3x.-50, -vec(shift(interior(h,:,jc))); color = hcolor, linewidth = 0.8)
lines!(ax_d, 1e-3x.-50, -vec(shift(interior(h,:,ja+Δj))); color = hcolor, linewidth = 0.8)
rowgap!(gabc, 3)
colgap!(gabc, 1, 15)
colgap!(gabc, 2, 3)
for row = [2,4]
    rowsize!(gabc, row, Relative(0.1))
end
resize_to_layout!(fig)
save(filesave * "wasbasfields_" * fileparams * "_d10.pdf", fig; pt_per_unit = 1)

# Compute spectral vertical boyancy flux 
Nz = length(z)
C10 = zeros(Nz, length(C1.spec));
for k = 1:Nz
    if z[k] > -200
        bk = α * g * interior(T̅ˢ, :, :, k)
        wk = interior(compute!(Field(@at (Center, Center, Center) w̅ˢ)), :, :, k)
        C10[k,:] = real.(isotropic_powerspectrum(bk, wk, x, y).spec)
        println("Spectral done at z = $(z[k]) m...")
    end
end

Z10 = C10 .* reshape(C1.freq, (1,:));
fig = Figure(size = (640, 450))
g4 = fig[1, 1] = GridLayout()
axis_kwargs = (ylabel = L"z~\text{(m)}", xlabel = L"\text{Wavenumber (m^{-1})}",xscale = log10, ygridvisible = false, 
               limits = ((6e-5, 4e-2), (-150,0)),xticks = ([1e-4,1e-3,1e-2], [L"10^{-4}",L"10^{-3}",L"10^{-2}"]), xgridvisible = false,
               xminorticks = [4e-5:1e-5:9e-5; 2e-4:1e-4:9e-4; 2e-3:1e-3:9e-3; 2e-2:1e-2:5e-2],xminorticksvisible = true)
ax_a = Axis(g4[1,1]; titlealign = :left, title=L"\text{(a)}~k\hat{w}\hat{b}~\text{(m^2 s^{-3})}", axis_kwargs...)
ax_b = Axis(g4[1,2]; titlealign = :left, title=L"\text{(b) 1d average around day 10}", axis_kwargs...)
ax_c = Axis(g4[1,3]; titlealign = :left, title=L"\text{(c) submesoscale flux}", ylabel = L"z~\text{(m)}", xlabel = L"\langle w^s b^s \rangle \text{(10^{-9} m^2 s^{-3})}")
hideydecorations!(ax_b, ticks = false)
hideydecorations!(ax_c, ticks = false)
vbound = maximum(abs.(Z10[:,:]))
hm_a = heatmap!(ax_a, C1.freq, z, Z10[:,:]'; rasterize = true, colorrange = (-vbound, vbound), colormap = :balance)
#hm_b = heatmap!(ax_b, C1.freq, z, mean(Z1d, dims = 1)[1,:,:]'; rasterize = true, colorrange = (-vbound, vbound), colormap = :balance)
vlines!(ax_a, 4e-4; color = :black, linewidth = 0.8)
#vlines!(ax_b, 4e-4; color = :black, linewidth = 0.8)
vlines!(ax_a, 4e-3; color = :black, linewidth = 0.8)
#vlines!(ax_b, 4e-3; color = :black, linewidth = 0.8)
Colorbar(g4[2, 1:2], hm_a, vertical = false)
for i = 1:2
    colgap!(g4, i, 5)
end
idx = 4e-4 .<= C1.freq .<= 4e-3
Δf = mean(diff(C1.freq[idx]))
lines!(ax_c, 1e9vec(sum(C10[:,idx], dims = 2) .* Δf), z; linewidth = 1)
#lines!(ax_c, 1e9vec(mean(sum(C1d[:,:,idx], dims = 3), dims = 1) .* Δf), z; linewidth = 1, label = "1d average")
resize_to_layout!(fig)
save(filesave * "washbashz_" * fileparams * "_d10.pdf", fig; pt_per_unit = 1)

#############################
# Let's pick the last snapshot!
times = snapshots[:T].times
t0 = time()
idx = 33:2:length(times);
wbm = zeros(2,length(times));
MLms = zeros(2,length(times));
Cₑs = zeros(2,length(times));
for i = idx
    wi,Ti = snapshots[:w][i],snapshots[:T][i];
    for (j,l) in enumerate([8kilometer, 3.6kilometer])
        w̅i = ZFaceField(wi.grid); 
        T̅i = CenterField(Ti.grid);
        coarse_graining!(wi, w̅i; kernel, cutoff = l)
        coarse_graining!(Ti, T̅i; kernel, cutoff = l)
        wˢ = compute!(Field(wi - w̅i));
        Tˢ = compute!(Field(Ti - T̅i));
        w̅ˢ = ZFaceField(wi.grid); 
        T̅ˢ = CenterField(Ti.grid);
        coarse_graining!(wˢ, w̅ˢ; kernel, cutoff)
        coarse_graining!(Tˢ, T̅ˢ; kernel, cutoff)
        vbfi = compute!(Field(α * g * T̅ˢ * w̅ˢ));
        vbf̅i  = CenterField(Ti.grid);
        coarse_graining!(vbfi , vbf̅i; kernel, cutoff = l)
        vbf̅iz = MLaverage(snapshots,i,vbf̅i; kernel, scale=l);
        println("cg-ed wb down t = $((time()-t0)/60.0) minutes")

        MLIi = compute!(Field(MLI(snapshots,i; kernel, scale=l)));
        println("cg-ed MLI down t = $((time()-t0)/60.0) minutes")

        wbm[j,i] = median(interior(vbf̅iz,:,:,1))
        MLms[j,i] = median(interior(MLIi,:,:,1))
        Cₑs[j,i] = median(interior(vbf̅iz,:,:,1)./ interior(MLIi,:,:,1))
    end
    println("t = $(times[i]/3600/24) days, Cₑ1(t) = $(Cₑs[1,i]), Cₑ2(t) = $(Cₑs[2,i])")
end

fig = Figure(size = (640, 450))
gab = fig[1, 1] = GridLayout()
for (i,l) in enumerate([8kilometer, 3.2kilometer])
    ax1 = Axis(gab[i,1]; xlabel = "days",ylabel = L"(\text{m^2 s^{-3}})",title="MLI parameterization at $(l/1000)km")
    ax2 = Axis(gab[i,1], ylabel = L"C_e(t)", yticklabelcolor = :blue, yaxisposition = :right)
    hidespines!(ax2)
    hidexdecorations!(ax2)
    lines!(ax1, times[idx]/3600/24, wbm[i,idx], label = L"\langle \overline{{w^s b^s}}^z \rangle", color = :black)
    lines!(ax1, times[idx]/3600/24, Cₑs[i,idx] .* MLms[i,idx], label = L"C_e(t)\times \text{MLI}", color = :red)
    lines!(ax1, times[idx]/3600/24, mean(Cₑs[i,idx]) * MLms[i,idx], label = L"C_e \times \text{MLI}", color = :red, linestyle = :dash)
    lines!(ax2, times[idx]/3600/24, Cₑs[i,idx], label = L"C_e(t)", color = :blue)
    if i == 1
        hidexdecorations!(ax1)
        axislegend(ax1, labelsize=9, framevisible = false,
                padding = (0f0, 0f0, 0f0, 0f0), patchlabelgap = 3)
    end
end
resize_to_layout!(fig)
save(filesave * "ts_wb_MLI_Ce_" * fileparams * ".pdf", fig)

#############################
k = 1
Cₑi = median(interior(vbf̅iz,:,:,1)./ interior(MLIi,:,:,1));
vary = Cₑi*(xhift(interior(MLIi, :, :, k)));
varx = xhift(interior(vbf̅iz , :, :, k));

pos_indices = (varx .> 0);
fig = Figure(size = (560, 570))
gabc = fig[1, 1] = GridLayout()
limmin, limmax = 10^(-12), 10^(-6)
axis_kwargs = (xgridvisible = false, ygridvisible = false, xscale=log10, yscale=log10,limits = ((limmin, limmax), (limmin, limmax)))
ax_a = Axis(gabc[1,1]; titlealign = :left, title=L"\text{(a)~Full domain}", ylabel = L"C_e(t)\times\text{MLI}~\text{(m^2 s^{-3})}", axis_kwargs...)
ax_b = Axis(gabc[1,2]; titlealign = :left, title=L"\text{(b)~Region A}", axis_kwargs...) 
ax_c = Axis(gabc[2,1]; titlealign = :left, title=L"\text{(c)~Region B}", xlabel = L"| \overline{\langle w^s b^s \rangle}^z |~\text{(m^2 s^{-3})}", ylabel = L"C_e(t)\times\text{MLI}~\text{(m^2 s^{-3})}", axis_kwargs...) 
ax_d = Axis(gabc[2,2]; titlealign = :left, title=L"\text{(d)~Region C}", xlabel = L"| \overline{\langle w^s b^s \rangle}^z |~\text{( m^2 s^{-3})}", axis_kwargs...) 
scatter!(ax_a, vec(varx[pos_indices]), vec(vary[pos_indices]); markersize = 5, color = :red, rasterize = true, alpha = 0.5, label = L"\text{positive data}")
scatter!(ax_a, -vec(varx[.!pos_indices]), vec(vary[.!pos_indices]); markersize = 5, color = :blue, rasterize = true, alpha = 0.5, label = L"\text{negative data}")
lines!(ax_a, [limmin, limmax], [limmin, limmax], linestyle = :dash, color = :black, linewidth = 1, label = "1-to-1")
axislegend(ax_a, labelsize=9, framevisible = false, font = texfont(), position = :lt,patchsize = (15, 1), 
           padding = (0f0, 0f0, 0f0, 0f0), patchlabelgap = 3, rowgap = 1)
for (j,ax) in enumerate([ax_b,ax_c,ax_d])
    xrange = findfirst(p3[j][1] .< 1e-3*x .- 50):findlast(1e-3*x .- 50 .<= p3[j][1]+l3[j]) 
    yrange = findfirst(p3[j][2]-l3[j] .< 1e-3*y .- 25):findlast(1e-3*y .- 25 .<= p3[j][2]) 
    xj = vec((pos_indices.*varx)[xrange,yrange])
    yj = vec((pos_indices.*vary)[xrange,yrange])
    xj, yj = xj[xj.>0], yj[xj.>0]
    scatter!(ax, xj, yj; markersize = 5, color = :red, alpha = 0.5, rasterize = true, label = L"\text{data}")
    xj = vec((.!pos_indices.*varx)[xrange,yrange])
    yj = vec((.!pos_indices.*vary)[xrange,yrange])
    xj, yj = xj[xj.<0], yj[xj.<0]
    scatter!(ax, -xj, yj; markersize = 5, color = :blue, alpha = 0.5, rasterize = true, label = L"\text{data}")
    lines!(ax, [limmin, limmax], [limmin, limmax], linestyle = :dash, color = :black, linewidth = 1, label = L"|\overline{wb}^z|^2")
end
hidexdecorations!(ax_a, ticks = false)
hidexdecorations!(ax_b, ticks = false)
hideydecorations!(ax_b, ticks = false)
hideydecorations!(ax_d, ticks = false)
rowgap!(gabc, 3)
colgap!(gabc, 1, 15)
resize_to_layout!(fig)
save(filesave * "CeMLI_vbfz_" * fileparams * "_d10_cg3.6.pdf", fig; pt_per_unit = 1)

using StatsBase
x, y = vec(varx[pos_indices]), vec(vary[pos_indices])

# Create a 2D histogram
hist = fit(Histogram, (x, y), nbins=(30, 30))

# Plot the density
fig = Figure(size = (560, 570))
ax = Axis(fig[1, 1]; xlabel="X", ylabel="Y", title="2D Density Plot")
heatmap!(fig[1, 1], hist.edges[1], hist.edges[2], hist.weights, colormap=:viridis)
