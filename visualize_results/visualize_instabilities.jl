using LESStudySetup
using CairoMakie
using Printf, Dates
using Oceananigans: compute!
using Oceananigans.BoundaryConditions
using Oceananigans.Grids: xnodes, ynodes, znodes
using Oceananigans.Operators: div_xyᶜᶜᶜ
using Statistics: mean, std, median
using LESStudySetup.Diagnostics
using LESStudySetup.Oceananigans.Units
using LESStudySetup.Diagnostics: N², M², Bₕ, wb
using LESStudySetup.Diagnostics: load_snapshots, MLD, MLaverage
using LESStudySetup.Diagnostics: isotropic_powerspectrum, coarse_grained_fluxes, δ
using LESStudySetup.Diagnostics: MixedLayerN², MixedLayerDepth, MLI, spatial_filtering
using LESStudySetup.Diagnostics: subfilter_stress!, coarse_graining!, _horizontal_gauss_filter!
using MathTeXEngine,GibbsSeaWater
set_theme!(theme_latexfonts(), fontsize=12,figure_padding = 10)
shift(x) = [x[size(x,1)÷2+1:end, :]; x[1:size(x,1)÷2, :]]
yshift(x) = [x[:, 3size(x,2)÷4+1:end] x[:, 1:3size(x,2)÷4]]
xhift(x) = yshift(shift(x))

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

fileparams = "hydrostatic_twin_simulation"
filehead = "./"
filename = filehead * "hydrostatic_snapshots_" * fileparams * ".jld2"
metadata = filehead * "experiment_" * fileparams * "_metadata.jld2"
freename = filehead * "hydrostatic_free_surface_" * fileparams * ".jld2"
filesave = filehead * "results/"

# load all the data!!
println("Loading data from $filename...")
snapshots = load_snapshots(filename; metadata,variables = (:u, :v, :w, :T,:pHY′))
g = parameters.g;
α = parameters.α;
ls,lc = 20kilometer, 8kilometer;
kernel! = _horizontal_gauss_filter!;

# sound speed as a function of in situ temperature and depth
θ = asind(parameters.f/2/7.2921e-5)
_,_,z = nodes(snapshots[:T][1])
p = gsw_p_from_z.(z, θ);
sa = 35;
#c = gsw_sound_speed_t_exact.(sa, T, p)

# Let's pick the last snapshot!
times = snapshots[:T].times
t0 = now()
idx = 33:2:length(times);
MLsm = zeros(size(times));
MLm = zeros(size(times));
wbm = zeros(size(times));
Cₑs = zeros(size(times));
Cₑ = zeros(size(times));
for i = idx
    wi,Ti = snapshots[:w][i],snapshots[:T][i];
    w̅i = spatial_filtering(wi; smoothing_range = ls, kernel!);
    T̅i = spatial_filtering(Ti; smoothing_range = ls, kernel!);
    wˢ = compute!(Field(wi - w̅i));
    Tˢ = compute!(Field(Ti - T̅i));
    vbfi = compute!(Field(α * g * Tˢ * wˢ));
    vbfhi = MLaverage(snapshots,i,vbfi);
    vbfh̅i  = CenterField(Ti.grid);
    coarse_graining!(vbfhi , vbfh̅i ; cutoff = lc)
    println("cg-ed wb down t = $((now()-t0)/60) minutes")

    MLi = compute!(Field(MLI(snapshots,i)));
    M̅i = CenterField(Ti.grid);
    coarse_graining!(MLi , M̅i ; cutoff = lc)
    println("cg-ed MLI down t = $((now()-t0)/60) minutes")

    wbm[i] = median(interior(vbfh̅i,:,:,1))
    MLm[i] = median(interior(M̅i,:,:,1))
    Cₑ[i] = median(interior(vbfh̅i,:,:,1)./ interior(M̅i,:,:,1))
    println("t = $(times[i]/3600/24) days, Cₑ(t) = $(Cₑ[i]), Cₑs(t) = $(Cₑs[i])")
end

fig = Figure(size = (640, 300))
ga = fig[1, 1] = GridLayout()
ax1 = Axis(ga[1,1]; xlabel = "days",ylabel = L"\text(m^2 s^{-3})")
ax2 = Axis(ga[1, 1], ylabel = L"C_e(t)", yticklabelcolor = :blue, yaxisposition = :right)
hidespines!(ax2)
hidexdecorations!(ax2)
lines!(ax1, times[idx]/3600/24, wbm[idx], label = L"\langle \overline{{w^s b^s}}^z \rangle", color = :black)
lines!(ax1, times[idx]/3600/24, Cₑ[idx] .* MLm[idx], label = L"C_e(t)\times \text{MLI}", color = :red)
lines!(ax1, times[idx]/3600/24, mean(Cₑ[idx]) * MLm[idx], label = L"C_e \times \text{MLI}", color = :red, linestyle = :dash)
lines!(ax2, times[idx]/3600/24, Cₑ[idx], label = L"C_e(t)", color = :blue)
axislegend(ax1, labelsize=9, framevisible = false,
                padding = (0f0, 0f0, 0f0, 0f0), patchlabelgap = 3)
resize_to_layout!(fig)
save(filesave * "ts_wb_MLI_Ce_" * fileparams * ".pdf", fig)

i = 161;
cutoff = 300/2.4*2*π #785.40
wi,Ti = snapshots[:w][i],snapshots[:T][i];
w̅i = spatial_filtering(wi; smoothing_range = ls, kernel!);
T̅i = spatial_filtering(Ti; smoothing_range = ls, kernel!);
wˢ = ZFaceField(wi.grid);# compute!(Field(wi - w̅i));
Tˢ = CenterField(Ti.grid);# compute!(Field(Ti - T̅i));
coarse_graining!(compute!(Field(wi - w̅i)), wˢ; cutoff)
coarse_graining!(compute!(Field(Ti - T̅i)), Tˢ; cutoff)
vbfi = compute!(Field(α * g * Tˢ * wˢ));
vbfhi = MLaverage(snapshots,i,vbfi);
#vbfh̅i  = CenterField(Ti.grid);
#coarse_graining!(vbfhi , vbfh̅i ; cutoff = lc)
MLi = compute!(Field(MLI(snapshots,i)));
#M̅i = CenterField(Ti.grid);
#coarse_graining!(MLi , M̅i ; cutoff = lc)

k = 1
lA, lB, lC = 40, 40, 40
pA, pB, pC = (-45, 70), (-20, 15), (10, 70)
l3 = [lA, lB, lC]
p3 = [pA, pB, pC]
cmap = :amp 
rmin, rmax = 0, 5
scolor = :black
alphabet = [letter for letter in 'a':'z'];

Cₑi = median(interior(vbfh̅i,:,:,1)./ interior(M̅i,:,:,1));
var,scale = M̅i,1e8Cₑi;
x, y, z = nodes(var)
fig = Figure(size = (640, 600))
gabc = fig[1, 1] = GridLayout()
axis_kwargs = (ylabel = "y (km)", aspect=1, limits = ((-50, 50), (0, 100)))
ax_a = Axis(gabc[1,1]; titlealign = :left, title=L"\text{(a)}~C_e(t)\times \mathrm{MLI}~\text{(10^{-8} m^2 s^{-3})}", axis_kwargs...)
ax_b = Axis(gabc[1,2]; titlealign = :left, title=L"\text{(b)~Region A}", aspect=1, limits = ((pA[1], pA[1]+lA), (pA[2]-lA, pA[2])))
ax_c = Axis(gabc[2,1]; titlealign = :left, title=L"\text{(c)~Region B}", xlabel = L"x~\text{(km)}", ylabel = L"y~\text{(km)}", aspect=1, limits = ((pB[1], pB[1]+lB), (pB[2]-lB,pB[2]))) 
ax_d = Axis(gabc[2,2]; titlealign = :left, title=L"\text{(d)~Region C}", xlabel = L"x~\text{(km)}", aspect=1, limits = ((pC[1], pC[1]+lC), (pC[2]-lC, pC[2]))) 
hm_a = heatmap!(ax_a, 1e-3x.-50, 1e-3y, scale*shift(interior(var,:,:,k)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
hm_b = heatmap!(ax_b, 1e-3x.-50, 1e-3y, scale*shift(interior(var,:,:,k)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
hm_c = heatmap!(ax_c, 1e-3x.-50, 1e-3y.-25, scale*xhift(interior(var,:,:,k)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
hm_d = heatmap!(ax_d, 1e-3x.-50, 1e-3y, scale*shift(interior(var,:,:,k)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
Colorbar(gabc[1,3], hm_b)
Colorbar(gabc[2,3], hm_d)
poly!(ax_a, Rect(pA[1], pA[2]-lA, lA, lA), color = (:white, 0.1), strokecolor = scolor, strokewidth = 0.5)
poly!(ax_a, Rect(pC[1], pC[2]-lC, lC, lC), color = (:white, 0.1), strokecolor = scolor, strokewidth = 0.5)
poly!(ax_a, Rect(pB[1], 75, lB, lB-pB[2]), color = (:white, 0.1), strokewidth = 0.)
poly!(ax_a, Rect(pB[1], 0, lB, pB[2]), color = (:white, 0.1), strokewidth = 0.)
vlines!(ax_a, [pB[1], pB[1]+lB]; ymin = 0.75, color = scolor, linewidth = 0.5)
vlines!(ax_a, [pB[1], pB[1]+lB]; ymax = 0.15, color = scolor, linewidth = 0.5)
hlines!(ax_a, [pB[2], 75]; xmin = 0.3, xmax = 0.7, color = scolor, linewidth = 0.5)
text!(ax_a, pA[1], pA[2], text = L"\text{A}", color = :black, align = (:left, :top))
text!(ax_a, pB[1], pB[2], text = L"\text{B}", color = :black, align = (:left, :top))
text!(ax_a, pC[1], pC[2], text = L"\text{C}", color = :black, align = (:left, :top))
rowgap!(gabc, 3)
colgap!(gabc, 1, 15)
colgap!(gabc, 2, 5)
resize_to_layout!(fig)
save(filesave * "CeMLIfields_" * fileparams * "_d10.pdf", fig; pt_per_unit = 1)

#############################
k = 1
Cₑi = median(interior(vbfh̅i,:,:,1)./ interior(M̅i,:,:,1));
vary = Cₑi*(xhift(interior(M̅i, :, :, k)));
varx = xhift(interior(vbfh̅i , :, :, k));

pos_indices = (varx .> 0);
fig = Figure(size = (560, 570))
gabc = fig[1, 1] = GridLayout()
limmin, limmax = 1e-10, 1e-7
axis_kwargs = (xgridvisible = false, ygridvisible = false, xscale=log10, yscale=log10,limits = ((limmin, limmax), (limmin, limmax)))
ax_a = Axis(gabc[1,1]; titlealign = :left, title=L"\text{(a)~Full domain}", ylabel = L"C_e(t)\times\text{MLI}~\text{(m^2 s^{-3})}", axis_kwargs...)
ax_b = Axis(gabc[1,2]; titlealign = :left, title=L"\text{(b)~Region A}", axis_kwargs...) 
ax_c = Axis(gabc[2,1]; titlealign = :left, title=L"\text{(c)~Region B}", xlabel = L"|\langle \overline{w^s b^s}^z \rangle|~\text{(m^2 s^{-3})}", ylabel = L"C_e(t)\times\text{MLI}~\text{(m^2 s^{-3})}", axis_kwargs...) 
ax_d = Axis(gabc[2,2]; titlealign = :left, title=L"\text{(d)~Region C}", xlabel = L"|\langle \overline{w^s b^s}^z \rangle|~\text{( m^2 s^{-3})}", axis_kwargs...) 
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
save(filesave * "MLI_wbz_" * fileparams * "_d10.pdf", fig; pt_per_unit = 1)

