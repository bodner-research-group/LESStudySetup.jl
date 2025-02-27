using LESStudySetup
using CairoMakie
using Printf, Dates
using Oceananigans: compute!
using Oceananigans.BoundaryConditions
using Oceananigans.Grids: xnodes, ynodes, znodes
using Oceananigans.Operators: div_xyᶜᶜᶜ
using Statistics: mean, std
using LESStudySetup.Diagnostics
using LESStudySetup.Oceananigans.Units
using LESStudySetup.Diagnostics: N², M², Bₕ
using LESStudySetup.Diagnostics: load_snapshots, MLD
using LESStudySetup.Diagnostics: isotropic_powerspectrum, coarse_grained_fluxes, δ
using LESStudySetup.Diagnostics: MixedLayerN², MixedLayerDepth
using LESStudySetup.Diagnostics: subfilter_stress!, coarse_graining!
using MathTeXEngine
set_theme!(theme_latexfonts(), fontsize=12,figure_padding = 5)
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
filehead = "/orcd/data/abodner/001/shirui/LESStudySetup.jl/"
filename = filehead * "hydrostatic_snapshots_" * fileparams * ".jld2"
metadata = filehead * "experiment_" * fileparams * "_metadata.jld2"
#fileparams = "nonhydrostatic"
#filename = filehead * fileparams * "_snapshots_0.jld2"
#metadata = filehead * fileparams * "_experiment_metadata.jld2"
filesave = filehead * "results/"

# load all the data!!
println("Loading data from $filename...")
snapshots = load_snapshots(filename; metadata)

# Let's pick the last snapshot!
times = snapshots[:T].times
snapshot_number = 161
nday = @sprintf("%2.0f", (times[snapshot_number])/60^2/24)
println("Plotting snapshot $snapshot_number on day $(nday)...")
t0 = now()

α = parameters.α
g = parameters.g
T = snapshots[:T][snapshot_number]
u = snapshots[:u][snapshot_number]
v = snapshots[:v][snapshot_number]
w = snapshots[:w][snapshot_number]
h = compute!(MLD(snapshots,snapshot_number; threshold = 0.03))
println("Loading fields wall time: $((now() - t0).value/1e3) seconds.")
println("mean MLD: $(mean(interior(h)))m")

grid = T.grid

# Coordinate arrays
xu, yu, zu = nodes(u)
xv, yv, zv = nodes(v)
xw, yw, zw = nodes(w)
xT, yT, zT = nodes(T)

initfile = filehead * "hydrostatic_snapshots_hydrostatic_background.jld2"
initsnaps = load_snapshots(initfile)
Ub = initsnaps[:u][1]
Vb = initsnaps[:v][1]
iarrow = 1:32:640

##############################
T0 = snapshots[:T][1]
T5 = snapshots[:T][81]
v0 = snapshots[:v][1]
v5 = snapshots[:v][81]
h0 = compute!(MLD(snapshots,1; threshold = 0.03))
h5 = compute!(MLD(snapshots,81; threshold = 0.03))
fig = Figure(size = (640, 800))
gab = fig[1, 1] = GridLayout()
axis_kwargs1 = (titlealign = :left, xlabel = "x (km)", ylabel = "y (km)", aspect = 1, limits = ((-50, 50), (0, 100)))
axis_kwargs2 = NamedTuple{(:titlealign,:xlabel,:limits,:aspect)}(axis_kwargs1)
ax_a = Axis(gab[1,1]; title=L"\text{(a)}~T~\text{({^\circ}C)},~t=0~\text{d},~z=2.8~\text{m}", axis_kwargs1...)
ax_b = Axis(gab[1,2]; title=L"\text{(b)}~T~\text{({^\circ}C)},~t=5~\text{d},~z=2.8~\text{m}", axis_kwargs2...)
k = kT#length(zT)
#Tmin, Tmax = minimum(interior(T5,:,:,k)), maximum(interior(T0,:,:,k))
hm_a = heatmap!(ax_a, 1e-3xT.-50, 1e-3yT, shift(interior(T0,:,:,k)); rasterize = true, colormap = :thermal, colorrange = (Tmin, Tmax))
hm_b = heatmap!(ax_b, 1e-3xT.-50, 1e-3yT, shift(interior(T5,:,:,k)); rasterize = true, colormap = :thermal, colorrange = (Tmin, Tmax))
hideydecorations!(ax_b, ticks = false)
hidexdecorations!(ax_a, ticks = false)
hidexdecorations!(ax_b, ticks = false)
Colorbar(gab[1, 3], hm_b)#, label = "°C", labelfont = texfont(),labelrotation = 0)
arrows!(ax_a, [0], [50], [-15*sqrt(3)],[-15],linecolor = :black, arrowcolor =:black)
text!(ax_a, -15*sqrt(3), 34, text = L"\text{Wind}", color = :black, align = (:right, :top))

ax_a = Axis(gab[3,1]; title=L"\text{(c)}~v~\text{(m~s^{-1})},~t=0~\text{d},~z=2.8~\text{m}", axis_kwargs1...)
ax_b = Axis(gab[3,2]; title=L"\text{(d)}~v~\text{(m~s^{-1})},~t=5~\text{d},~z=2.8~\text{m}", axis_kwargs2...)
vbnd = maximum(abs.(interior(v0,:,:,k)))
hm_a = heatmap!(ax_a, 1e-3xT.-50, 1e-3yT, shift(interior(v0,:,:,k)); rasterize = true, colormap = :balance, colorrange = (-vbnd, vbnd))
hm_b = heatmap!(ax_b, 1e-3xT.-50, 1e-3yT, shift(interior(v5,:,:,k)); rasterize = true, colormap = :balance, colorrange = (-vbnd, vbnd))
hideydecorations!(ax_b, ticks = false)
hidexdecorations!(ax_a, ticks = false)
hidexdecorations!(ax_b, ticks = false)
Colorbar(gab[3, 3], hm_b)
arrows!(ax_a,xT[iarrow]/1e3.-50, yT[iarrow]/1e3, shift(interior(Ub,iarrow,iarrow,k)), shift(interior(Vb,iarrow,iarrow,k)), arrowsize = 3, lengthscale = 1e2,linecolor = :black, arrowcolor = :black, linewidth = 0.2)
text!(ax_a, 25, 25, text = L"\text{Warm eddy}", color = :red, align = (:center, :center))
text!(ax_a, 25, 75, text = L"\text{Cold eddy}", color = :blue, align = (:center, :center))
text!(ax_a, -25, 75, text = L"\text{Warm eddy}", color = :red, align = (:center, :center))
text!(ax_a, -25, 25, text = L"\text{Cold eddy}", color = :blue, align = (:center, :center))

zmin = -100
ΔN = 160
ja = 2ΔN
kz = findlast(zw .< zmin)
Nz = length(zT)
axis_kwargs0 = (xlabel = "x (km)", ylabel = "z (m)", limits = ((-50, 50), (zmin, 0)))
axis_kwargs1 = NamedTuple{(:xlabel,:limits)}(axis_kwargs0)
ax_a = Axis(gab[2,1]; titlealign = :left, title=L"y=50~\text{km}", axis_kwargs0...)
ax_b = Axis(gab[2,2]; titlealign = :left, title=L"y=50~\text{km}", axis_kwargs1...)
Tmin, Tmax = minimum(interior(T5,:,ja,kz:Nz)), maximum(interior(T0,:,ja,kz:Nz))
hm_a = heatmap!(ax_a, 1e-3xT.-50, zT[kz:Nz], shift(interior(T0,:,ja,kz:Nz)); rasterize = true, colormap = :thermal, colorrange = (Tmin, Tmax))
hm_b = heatmap!(ax_b, 1e-3xT.-50, zT[kz:Nz], shift(interior(T5,:,ja,kz:Nz)); rasterize = true, colormap = :thermal, colorrange = (Tmin, Tmax))
hideydecorations!(ax_b, ticks = false)
Colorbar(gab[2, 3], hm_b)
lines!(ax_a, 1e-3xT.-50, -vec(shift(interior(h0,:,ja))); color = :white, linewidth = 0.8)
lines!(ax_b, 1e-3xT.-50, -vec(shift(interior(h5,:,ja))); color = :white, linewidth = 0.8)

ax_a = Axis(gab[4,1]; titlealign = :left, title=L"y=50~\text{km}", axis_kwargs0...)
ax_b = Axis(gab[4,2]; titlealign = :left, title=L"y=50~\text{km}", axis_kwargs1...)
hm_a = heatmap!(ax_a, 1e-3xv.-50, zv[kz:Nz], shift(interior(v0,:,ja,kz:Nz)); rasterize = true, colormap = :balance, colorrange = (-vbnd, vbnd))
hm_b = heatmap!(ax_b, 1e-3xv.-50, zv[kz:Nz], shift(interior(v5,:,ja,kz:Nz)); rasterize = true, colormap = :balance, colorrange = (-vbnd, vbnd))
hideydecorations!(ax_b, ticks = false)
idxz = 1:5:Nz
#arrows!(ax_a,xv[iarrow]/1e3.-50, zv[idxz], shift(interior(Ub,iarrow,ja,idxz)), 0*shift(interior(Ub,iarrow,ja,idxz)), arrowsize = 3, lengthscale = 1e2,linecolor = :black, arrowcolor = :black, linewidth = 0.2)
Colorbar(gab[4, 3], hm_b)

rowgap!(gab, 3)
colgap!(gab, 1, 15)
colgap!(gab, 2, 3)
for row = [2,4]
    rowsize!(gab, row, Relative(0.1))
end
resize_to_layout!(fig)
save(filesave * "Tvfields_" * fileparams * "_d0d5.pdf", fig)

#################
# Load parameters and simulation results
T0 = snapshots[:T][1]
ρ₀ = parameters.ρ₀
threshold = 3e-4/g*ρ₀
surface,stratification=true,false
h    = compute!(MixedLayerDepth(grid, (; T=T0); ΔT = abs(threshold / ρ₀ / α), surface,stratification))
stratification = true
Nh² = compute!(MixedLayerN²(grid, (; T=T0); ΔT = abs(threshold / ρ₀ / α), surface,stratification))
us = sqrt(parameters.τw/ρ₀)
Qᵘ = parameters.Q/(ρ₀ * parameters.cp) * α * g
λₛ=compute!(Field(Qᵘ*h/us^3))
Rₕ=compute!(Field(Nh² * h^2/us^2))
Rₕw = compute!(Field(Rₕ/λₛ^(2/3)))
f = parameters.f
compute!(Field(f/Nh²^0.5))

fig = Figure(size = (640, 560))
gab = fig[1, 1] = GridLayout()
axis_kwargs1 = (xlabel = "x (km)", ylabel = "y (km)", aspect = 1,
            limits = ((0, 100), (0, 100)))
axis_kwargs2 = NamedTuple{(:xlabel,:limits,:aspect)}(axis_kwargs1)

ax_a = Axis(gab[1,1]; titlealign = :left, title=L"\text{(a)}~\lambda_s=\frac{-B_0 h}{u_*^3}", axis_kwargs1...)
ax_b = Axis(gab[1,3]; titlealign = :left, title=L"\text{(b)}~R_h=(\frac{N_h h}{u_*})^2", axis_kwargs2...)
ax_c = Axis(gab[2,1]; titlealign = :left, title=L"\text{(c)}~f/N_h", axis_kwargs1...)
ax_d = Axis(gab[2,3]; titlealign = :left, title=L"\text{(d)}~R^*_h=(\frac{N_h h}{w_*})^2", axis_kwargs2...)

hidexdecorations!(ax_a, ticks = false)
hidexdecorations!(ax_b, ticks = false)
hideydecorations!(ax_b, ticks = false)
hideydecorations!(ax_d, ticks = false)

k = 1
hm_a = heatmap!(ax_a, 1e-3xT, 1e-3yT, shift(interior(λₛ,:,:,k)); rasterize = true, colormap = :amp, colorscale = log10)
hm_b = heatmap!(ax_b, 1e-3xT, 1e-3yT, shift(interior(Rₕ,:,:,k)); rasterize = true, colormap = :balance, colorrange = (10^(2.27), 10^(3.73)), colorscale = log10)
hm_c = heatmap!(ax_c, 1e-3xT, 1e-3yT, parameters.f ./ sqrt.(shift(interior(Nh²,:,:,k))); rasterize = true, colormap = Reverse(:balance), colorrange = (10^(-2.11), 10^(-1.39)), colorscale = log10)
hm_d = heatmap!(ax_d, 1e-3xT, 1e-3yT, shift(interior(Rₕw,:,:,k)); rasterize = true, colormap = :balance, colorrange = (10^(2.37), 10^(3.63)),colorscale = log10)

Colorbar(gab[1, 2], hm_a)
Colorbar(gab[1, 4], hm_b)
Colorbar(gab[2, 2], hm_c)
Colorbar(gab[2, 4], hm_d)
colgap!(gab, 1, 3)
colgap!(gab, 3, 3)
colgap!(gab, 2, 5)
rowgap!(gab, 1, 3)

save(filesave * "fields_" * fileparams * "_mldass_d0.pdf", fig)

#################
k = 216
Bh = compute!(Field(Bₕ(snapshots, snapshot_number)))
Bv = compute!(Field(Bᵥ(snapshots, snapshot_number)))
fig = Figure(size = (640, 315))
gab = fig[1, 1] = GridLayout()
axis_kwargs1 = (xlabel = "x (km)", ylabel = "y (km)", aspect = 1, limits = ((0, 100), (0, 100)))
axis_kwargs2 = NamedTuple{(:xlabel,:limits,:aspect)}(axis_kwargs1)
ax_a = Axis(gab[1,1]; titlealign = :left, title=L"\text{(a)}~10^{16}\mathcal{B}_h~\text{({kg}^2 m^{-8}s^{-1})}", axis_kwargs1...)
ax_b = Axis(gab[1,2]; titlealign = :left, title=L"\text{(b)}~10^{16}\mathcal{B}_v~\text{({kg}^2 m^{-8}s^{-1})}", axis_kwargs2...)
hm_a = heatmap!(ax_a, 1e-3xT, 1e-3yT, 1e16*shift(interior(Bh,:,:,k)); rasterize = true, colormap = :amp, colorrange = (-0., 1))
hm_b = heatmap!(ax_b, 1e-3xT, 1e-3yT, 1e16*shift(interior(Bv,:,:,k)); rasterize = true, colormap = :amp, colorrange = (-0., 1))
hideydecorations!(ax_b, ticks = false)
Colorbar(gab[1, 3], hm_b)
colgap!(gab, 1, 5)
colgap!(gab, 2, 3)
save(filesave * "fields_" * fileparams * "_BhBv_day10.pdf", fig)

###################
# Plot the fields
ΔN = 160
kT, ks, kw = 222, 222, 202
jslices = [2]*ΔN
Tmin, Tmax = minimum(interior(T,:,:,kT)), maximum(interior(T,:,:,kT))
vbnd, wbnd = maximum(abs, interior(v)), 5
lA, lB, lC = 40, 40, 40
pA, pB, pC = (-45, 70), (-20, 15), (10, 70)

#####################################
var = T#compute!(Field(1e3*w))
x, y, z = nodes(var)
k = kT
cmap = :thermal
rmin, rmax = Tmin, Tmax
hcolor = :white
fig = Figure(size = (640, 785))
gabc = fig[1, 1] = GridLayout()
axis_kwargs = (ylabel = "y (km)", aspect=1, limits = ((-50, 50), (0, 100)))
ax_a = Axis(gabc[1,1]; titlealign = :left, title=L"\text{(a)}~T~\text{({^\circ}C)},~z=-2.8~\text{m}", axis_kwargs...)
ax_b = Axis(gabc[1,2]; titlealign = :left, title=L"\text{(b)~Region A}", aspect=1, limits = ((pA[1], pA[1]+lA), (pA[2]-lA, pA[2])))
ax_c = Axis(gabc[3,1]; titlealign = :left, title=L"\text{(c)~Region B}", ylabel = "y (km)", aspect=1, limits = ((pB[1], pB[1]+lB), (pB[2]-lB,pB[2]))) 
ax_d = Axis(gabc[3,2]; titlealign = :left, title=L"\text{(d)~Region C}", aspect=1, limits = ((pC[1], pC[1]+lC), (pC[2]-lC, pC[2]))) 
hm_a = heatmap!(ax_a, 1e-3x.-50, 1e-3y, shift(interior(var,:,:,k)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
hm_b = heatmap!(ax_b, 1e-3x.-50, 1e-3y, shift(interior(var,:,:,k)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
hm_c = heatmap!(ax_c, 1e-3x.-50, 1e-3y.-25, xhift(interior(var,:,:,k)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
hm_d = heatmap!(ax_d, 1e-3x.-50, 1e-3y, shift(interior(var,:,:,k)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
Colorbar(gabc[1,3], hm_b)
Colorbar(gabc[3,3], hm_d)
hidexdecorations!(ax_a, ticks = false)
hidexdecorations!(ax_b, ticks = false)
hidexdecorations!(ax_c, ticks = false)
hidexdecorations!(ax_d, ticks = false)
hlines!(ax_a, 1e-3*yT[jslices]; color = :black, linestyle = :dash, linewidth = 0.8)
poly!(ax_a, Rect(pA[1], pA[2]-lA, lA, lA), color = (:white, 0.1), strokecolor = :black, strokewidth = 0.5)
poly!(ax_a, Rect(pC[1], pC[2]-lC, lC, lC), color = (:white, 0.1), strokecolor = :white, strokewidth = 0.5)
poly!(ax_a, Rect(pB[1], 75, lB, lB-pB[2]), color = (:white, 0.1), strokewidth = 0.)
poly!(ax_a, Rect(pB[1], 0, lB, pB[2]), color = (:white, 0.1), strokewidth = 0.)
vlines!(ax_a, [pB[1], pB[1]+lB]; ymin = 0.75, color = :white, linewidth = 0.5)
vlines!(ax_a, [pB[1], pB[1]+lB]; ymax = 0.15, color = :white, linewidth = 0.5)
hlines!(ax_a, [pB[2], 75]; xmin = 0.3, xmax = 0.7, color = :white, linewidth = 0.5)
text!(ax_a, pA[1], pA[2], text = L"\text{A}", color = :black, align = (:left, :top))
text!(ax_a, pB[1], pB[2], text = L"\text{B}", color = :black, align = (:left, :top))
if cmap == :thermal
    text!(ax_a, pC[1], pC[2], text = L"\text{C}", color = :white, align = (:left, :top))
else
    text!(ax_a, pC[1], pC[2], text = L"\text{C}", color = :black, align = (:left, :top))
end

dΔN = 80
jslices = dΔN * [0] .+ 2ΔN
hlines!(ax_b, 1e-3*yv[jslices].+6.25; color = :black, linestyle = :dash, linewidth = 0.8)
hlines!(ax_d, 1e-3*yv[jslices].+6.25; color = :black, linestyle = :dash, linewidth = 0.8)
hlines!(ax_c, -10; color = :black, linestyle = :dash, linewidth = 0.8)

zmin = -100
ΔN = 160
ja = 2ΔN
jc = 9*64
Δj = 40
kz = findlast(zw .< zmin)
Nz = length(zT)
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
hm_a = heatmap!(ax_a, 1e-3x.-50, z[kz:Nz], shift(interior(var,:,ja,kz:Nz)); rasterize = true, colormap = cmap, colorrange = (vmin, vmax))
hm_b = heatmap!(ax_b, 1e-3x.-50, z[kz:Nz], shift(interior(var,:,ja+Δj,kz:Nz)); rasterize = true, colormap = cmap, colorrange = (vmin, vmax))
hm_c = heatmap!(ax_c, 1e-3x.-50, z[kz:Nz], shift(interior(var,:,jc,kz:Nz)); rasterize = true, colormap = cmap, colorrange = (vmin, vmax))
hm_d = heatmap!(ax_d, 1e-3x.-50, z[kz:Nz], shift(interior(var,:,ja+Δj,kz:Nz)); rasterize = true, colormap = cmap, colorrange = (vmin, vmax))
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
save(filesave * "Tfields_" * fileparams * "_d$(nday).pdf", fig; pt_per_unit = 1)
println("Finished plotting fields, wall time: $((now() - t0).value/1e3) seconds.")

#####################
kTs = [222, 202, 171]
var = interior(snapshots[:T][1], :, :, kTs[1])
St0 = isotropic_powerspectrum(var, var, xT, yT)
freq = St0.freq
idx = freq .> 0
tsSt = zeros(length(times)÷2+1, length(kTs), length(freq))
tsvT = zeros(length(times)÷2+1, length(kTs))
for i in 1:length(times)÷2+1
    println("Computing spectra at time $(times[2i-1]/24/3600) days...")
    for j in 1:length(kTs)
        vari = interior(snapshots[:T][2i-1], :, :, kTs[j])
        tsSt[i,j,:] = Real.(isotropic_powerspectrum(vari, vari, xT, yT).spec./St0.spec[idx][1])
        tsvT[i,j] = mean((vari .- mean(vari)).^2)
    end
end
fig = Figure(size = (640, 750))
gl = fig[1, 1] = GridLayout()
ax = Axis(gl[1, 1]; titlealign = :left, title=L"\text{(a)}~\langle T^{\prime 2} \rangle(t,z)/\langle T^{\prime 2} \rangle(0,z=-2.8~\text{m})", limits = ((0.125,20),nothing))
lines!(ax, times[3:2:end]/24/3600, tsvT[2:end,1]/tsvT[1,1]; label = L"z=-2.8~\text{m}", linewidth = 1)
lines!(ax, times[3:2:end]/24/3600, tsvT[2:end,2]/tsvT[1,1]; label = L"z=-25~\text{m}", linewidth = 1)
lines!(ax, times[3:2:end]/24/3600, tsvT[2:end,3]/tsvT[1,1]; label = L"z=-60~\text{m}", linewidth = 1)
axislegend(ax, labelsize=9, framevisible = false, orientation = :horizontal,
           padding = (0f0, 0f0, 0f0, 0f0), patchlabelgap = 3)
hidexdecorations!(ax, ticks = false)
titles = [L"\text{(b)}~E_T(k,t,z=-2.8~\text{m})/E_T(k_{\min},0,-2.8~\text{m})", L"\text{(c)}~E_T(k,t,z=-25~\text{m})/E_T(k_{\min},0,-2.8~\text{m})", L"\text{(d)}~E_T(k,t,z=-60~\text{m})/E_T(k_{\min},0,-2.8~\text{m})"]
for i = 1:3
    ax = Axis(gl[1+i, 1]; ylabel = L"\text{Wavenumber}~k~\text{(rad~m^{-1})}", yscale = log10,
                        titlealign = :left, title=titles[i],
                        limits = ((0.125,20),(6e-5, 3e-2)))
    hm = heatmap!(ax, times[3:2:end]/24/3600, freq[idx], tsSt[2:end,i,idx]; colorrange = (1e-10,1), colormap = :amp, colorscale = log10)
    Colorbar(gl[1+i, 2], hm)
    ylims!(ax, (6.5e-5, 2.5e-2))
    if i < 3
        hidexdecorations!(ax, ticks = false)
    else
        ax.xlabel = L"\text{Time}~t~\text{(days)}"
    end
end
colgap!(gl, 1, 3)
for i in 1:3
    rowgap!(gl, i, 3)
end
rowsize!(gl, 1, Relative(0.2))
resize_to_layout!(fig)
save(filesave * "ETts_" * fileparams * ".pdf", fig; pt_per_unit = 1)

#####################
# Compute the horizontal spectrum of T, u, v, w 
axis_kwargs1 = (xlabel = L"\text{Wavenumber (rad m^{-1})}", xgridvisible = false,
                ylabel = L"E_i(k)/E_{T,v}(k_{\min},z=-2.8~\text{m})", ygridvisible = false,
                xscale = log10, yscale = log10,
                limits = ((4e-5, 6e-2), (1e-9,1e2)),
                xminorticks = [4e-5:1e-5:9e-5; 2e-4:1e-4:9e-4; 2e-3:1e-3:9e-3; 2e-2:1e-2:6e-2],xminorticksvisible = true)
axis_kwargs2 = NamedTuple{(:xlabel,:xscale,:yscale,:limits,:xgridvisible,:ygridvisible,:xminorticks,:xminorticksvisible)}(axis_kwargs1)
axis_kwargs3 = NamedTuple{(:ylabel,:xscale,:yscale,:limits,:xgridvisible,:ygridvisible,:xminorticks,:xminorticksvisible)}(axis_kwargs1)
axis_kwargs4 = NamedTuple{(:xscale,:yscale,:limits,:xgridvisible,:ygridvisible,:xminorticks,:xminorticksvisible)}(axis_kwargs1)
Nr = length(xv) ÷ 2

fig = Figure(size = (640, 800))
g33 = fig[1, 1] = GridLayout()
S0s = []
alphabet = [letter for letter in 'a':'z'];
l3 = [lA, lB, lC]
p3 = [pA, pB, pC]
for (i,klev) in enumerate([222, 202, 171])
    println("Plotting spectra at z = $(zT[klev])m...")
    if i == 3
        ax = Axis(g33[1, i]; titlealign = :left, title=L"\text{(i)}~z=-60~\text{m}", axis_kwargs4...)
        hideydecorations!(ax, ticks = false)
    else
        if i == 2
            ax = Axis(g33[1, i]; titlealign = :left, title=L"\text{(e)}~z=-25~\text{m}", axis_kwargs4...)
            hideydecorations!(ax, ticks = false)
        else
            ax = Axis(g33[1, 1]; titlealign = :left, title=L"\text{(a)}~z=-2.8~\text{m}", axis_kwargs3...)
            vlines!(ax, 1/300; color = :black, linewidth = 0.8)
        end
    end

    hidexdecorations!(ax, ticks = false)
    Su = isotropic_powerspectrum(interior(u, :, :, klev), interior(u, :, :, klev), xu, yu)
    Sv = isotropic_powerspectrum(interior(v, :, :, klev), interior(v, :, :, klev), xv, yv)
    wk = (interior(w, :, :, klev)+interior(w, :, :, klev+1))/2
    Sw = isotropic_powerspectrum(wk, wk, xw, yw)
    St = isotropic_powerspectrum(interior(T, :, :, klev), interior(T, :, :, klev), xT, yT)

    if i == 1
        global Sv0,St0 = Sv,St
    end

    idx = Su.freq .> 0
    lines!(ax, Su.freq[idx], (10^-6.5)Su.freq[idx].^-2, linestyle = :dash, color = :black)
    text!(ax, 10^-3, 0.7; text = L"k^{-2}")
    idxf = idx#0 .< Su.freq .<= 1e-3
    lines!(ax, Su.freq[idxf], 1e-13Su.freq[idxf].^-3, linestyle = :dash, color = :gray)
    text!(ax, 10^-4, 1e-3; text = L"k^{-3}")
    #idxf = Su.freq .> 1e-3
    #lines!(ax, Su.freq[idxf], 1e-9Su.freq[idxf].^(-5/3), linestyle = :dash, color = :gray)
    #text!(ax, 10^-3, 10^-5.5; text = L"k^{-5/3}")
    lines!(ax, St.freq[idx], Real.(St.spec[idx]./St0.spec[idx][1]), color = :red, label = L"E_T")
    lines!(ax, Su.freq[idx], Real.(Su.spec[idx]./Sv0.spec[idx][1]), color = :blue, label = L"E_u")
    lines!(ax, Sv.freq[idx], Real.(Sv.spec[idx]./Sv0.spec[idx][1]), color = :green, label = L"E_v")
    lines!(ax, Sw.freq[idx], Real.(Sw.spec[idx]./Sv0.spec[idx][1]), color = :black, label = L"E_w")
    xlims!(ax, (4e-5, 3.5e-2))
    ylims!(ax, (1e-9, 1e1))
    if i == 1
        axislegend(ax, labelsize=9, patchsize = (15, 1), framevisible = false,
                   padding = (0f0, 0f0, 0f0, 0f0), patchlabelgap = 3, rowgap = 1)
    end

    uk = (xhift(interior(u, :, :, klev)))
    vk = (xhift(interior(v, :, :, klev)))
    wk = ((xhift(interior(w, :, :, klev))+xhift(interior(w, :, :, klev+1)))/2)
    Tk = (xhift(interior(T, :, :, klev)))
    for j = 1:3
        println([i,j])
        window = 1
        nalpha = 1+j+(i-1)*4
        title = "("*alphabet[nalpha]*") Region " * uppercase(alphabet[j])

        ax = Axis(g33[j+1, i]; titlealign = :left, title=title, titlefont=texfont(),axis_kwargs4...)
        if j < 3
            hidexdecorations!(ax, ticks = false)
        else
            ax.xlabel = L"\text{Wavenumber (rad m^{-1})}"
        end

        if i > 1
            hideydecorations!(ax, ticks = false)
        else
            ax.ylabel = L"E_i(k)/E_{T,v}(k_{\min},z=-2.8~\text{m})"
        end
        jr = 1-j
        xrange = findfirst(p3[j][1] .< 1e-3*xT .- 50):findlast(1e-3*xT .- 50 .<= p3[j][1]+l3[j]) 
        yrange = findfirst(p3[j][2]-l3[j] .< 1e-3*yT .- 25):findlast(1e-3*yT .- 25 .<= p3[j][2]) 
        #xrange = 1+(j-1)*Nr÷2:(j-1)*Nr÷2+Nr
        #yrange = j != 2 ? (1+Nr:2Nr) : (1:Nr)

        Su = isotropic_powerspectrum(uk[xrange,yrange], uk[xrange,yrange], xu[xrange], yu[yrange];window)
        Sv = isotropic_powerspectrum(vk[xrange,yrange], vk[xrange,yrange], xv[xrange], yv[yrange];window)
        Sw = isotropic_powerspectrum(wk[xrange,yrange], wk[xrange,yrange], xw[xrange], yw[yrange];window)
        St = isotropic_powerspectrum(Tk[xrange,yrange], Tk[xrange,yrange], xT[xrange], yT[yrange];window)

        idx = Su.freq .> 0
        if i == 1
            push!(S0s, [St.spec[idx][1], Sv.spec[idx][1]])
        end
        lines!(ax, Su.freq[idx], (10^-6.5)Su.freq[idx].^-2, linestyle = :dash, color = :black)
        idxf = idx#0 .< Su.freq .<= 1e-3
        lines!(ax, Su.freq[idxf], 1e-13Su.freq[idxf].^-3, linestyle = :dash, color = :gray)
        #idxf = Su.freq .> 1e-3
        #lines!(ax, Su.freq[idxf], 1e-9Su.freq[idxf].^(-5/3), linestyle = :dash, color = :gray)
        lines!(ax, St.freq[idx], Real.(St.spec[idx]./(S0s[j][1])), color = :red, label = L"E_T")
        lines!(ax, Su.freq[idx], Real.(Su.spec[idx]./(S0s[j][2])), color = :blue, label = L"E_u")
        lines!(ax, Sv.freq[idx], Real.(Sv.spec[idx]./(S0s[j][2])), color = :green, label = L"E_v")
        lines!(ax, Sw.freq[idx], Real.(Sw.spec[idx]./(S0s[j][2])), color = :black, label = L"E_w")
        xlims!(ax, (4e-5, 3.5e-2))
        ylims!(ax, (1e-9, 1e1))
    end

end
rowgap!(g33, 3)
resize_to_layout!(fig)
save(filesave * "spectra_" * fileparams * "_d$(nday).pdf", fig; pt_per_unit = 1)
println("Finished plotting spectra, wall time: $((now() - t0).value/1e3) seconds.")

###################################
# Plot x-y slices of ζ/f and δ/f Bh
Bh = compute!(Field(Bₕ(snapshots, snapshot_number)))
#title=L"\text{(a)}~10^{16}\mathcal{B}_h~\text{({kg}^2 m^{-8}s^{-1})}", axis_kwargs1...)
f = parameters.f
u̅  = XFaceField(u.grid)
v̅  = YFaceField(v.grid)
cutoff = 300/2.4*2*π #785.40
coarse_graining!(u , u̅ ; cutoff)
coarse_graining!(v , v̅ ; cutoff)
fill_halo_regions!(u̅)
fill_halo_regions!(v̅)
R̅o = compute!(Field((∂x(v̅) - ∂y(u̅))/f))
D̅ = compute!(Field(KernelFunctionOperation{Center, Center, Center}(div_xyᶜᶜᶜ, grid, u̅, v̅)/f))
alphabet = [letter for letter in 'a':'z'];
∇b = compute!(Field(α * g * (∂x(T)^2 + ∂y(T)^2)^0.5))
k = 222
rbnd,dbnd = 2,1
M²₀ = parameters.M²₀
b̅ = CenterField(T.grid)
coarse_graining!(compute!(Field(α * g * T)), b̅; cutoff)
fill_halo_regions!(b̅)
w̅  = ZFaceField(w.grid)
coarse_graining!(w , w̅ ; cutoff)
fill_halo_regions!(w̅)
ωz = ∂x(v̅) - ∂y(u̅) + f
ωx = ∂y(w̅) - ∂z(v̅)
ωy = ∂z(u̅) - ∂x(w̅)
PV = compute!(Field(ωx * ∂x(b̅) + ωy * ∂y(b̅) + ωz * ∂z(b̅)))
θ, τ, ρ₀=parameters.θ,parameters.τ,parameters.ρ₀
τx,τy=-τ*sind(θ),τ*cosd(30.)
EBF = compute!(Field((τx * ∂x(b̅) + τy * ∂y(b̅))/f/ρ₀))
B̅h = compute!(Field(-(∂x(b̅)^2 * ∂x(u̅) + ∂y(b̅)^2 * ∂y(v̅))-∂x(b̅)*∂y(b̅)*(∂x(v̅) + ∂y(u̅))))
∇b̅ = compute!(Field((∂x(b̅)^2 + ∂y(b̅)^2)^0.5))

var = PV#B̅h
x, y, z = nodes(var)
cmap = :balance #Reverse(:grays)
rmin, rmax = -1, 1
fig = Figure(size = (640, 610))
gabc = fig[1, 1] = GridLayout()
axis_kwargs = (ylabel = "y (km)", aspect=1, limits = ((-50, 50), (0, 100)))
ax_a = Axis(gabc[1,1]; titlealign = :left, title=L"\text{(a)}~\tilde{q}~\text{(10^{-8} s^{-3})},~z=-2.8~\text{m}", axis_kwargs...)
ax_b = Axis(gabc[1,2]; titlealign = :left, title=L"\text{(b)~Region A}", aspect=1, limits = ((-50, 0), (25, 75)))
ax_c = Axis(gabc[2,1]; titlealign = :left, title=L"\text{(c)~Region B}", xlabel = "x (km)", ylabel = "y (km)", aspect=1, limits = ((-25, 25), (-25,25))) 
ax_d = Axis(gabc[2,2]; titlealign = :left, title=L"\text{(d)~Region C}", xlabel = "x (km)", aspect=1, limits = ((0, 50), (25, 75))) 
hm_a = heatmap!(ax_a, 1e-3x.-50, 1e-3y, 1e8shift(interior(var,:,:,k)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
hm_b = heatmap!(ax_b, 1e-3x.-50, 1e-3y, 1e8shift(interior(var,:,:,k)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
hm_c = heatmap!(ax_c, 1e-3x.-50, 1e-3y.-75, 1e8xhift(interior(var,:,:,k)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
hm_d = heatmap!(ax_d, 1e-3x.-50, 1e-3y, 1e8shift(interior(var,:,:,k)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
Colorbar(gabc[1,3], hm_b)
Colorbar(gabc[2,3], hm_d)
poly!(ax_a, [Rect(-50i, 25, 50, 50) for i in 0:1], color = (:white, 0.1), strokecolor = :black, strokewidth = 0.5)
poly!(ax_a, [Rect(-25, 75i, 50, 25) for i in 0:1], color = (:white, 0.1), strokewidth = 0.)
vlines!(ax_a, [-25, 25]; ymin = 0.75, color = :black, linewidth = 0.5)
vlines!(ax_a, [-25, 25]; ymax = 0.25, color = :black, linewidth = 0.5)
text!(ax_a, -50, 75, text = L"\text{A}", color = :black, align = (:left, :top))
text!(ax_a, -25, 25, text = L"\text{B}", color = :black, align = (:left, :top))
text!(ax_a, 0, 75, text = L"\text{C}", color = :black, align = (:left, :top))
rowgap!(gabc, 3)
colgap!(gabc, 1, 15)
colgap!(gabc, 2, 5)
resize_to_layout!(fig)
save(filesave * "PVfields_" * fileparams * "_d$(nday).pdf", fig; pt_per_unit = 1)

ro = compute!(Field(ζ(snapshots, snapshot_number)/f));
rd = compute!(Field(δ(snapshots, snapshot_number)/f));

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

Nr = length(xv) ÷ 2
fig = Figure(size = (560, 570))
gabc = fig[1, 1] = GridLayout()
axis_kwargs = (xgridvisible = false, ygridvisible = false,limits = ((-7, 20), (-7, 20)))
ax_a = Axis(gabc[1,1]; titlealign = :left, title=L"\text{(a)}~z=-2.8~\text{m}", ylabel = L"\tilde{\mathcal{B}}_h/(M_0^4 f)", axis_kwargs...)
ax_b = Axis(gabc[1,2]; titlealign = :left, title=L"\text{(b)~Region A}", axis_kwargs...) 
ax_c = Axis(gabc[2,1]; titlealign = :left, title=L"\text{(c)~Region B}", xlabel = L"-| \nabla_h \overline{b}|^2 \delta/(M_0^4 f)", ylabel = L"\tilde{\mathcal{B}}_h/(M_0^4 f)", axis_kwargs...) 
ax_d = Axis(gabc[2,2]; titlealign = :left, title=L"\text{(d)~Region C}", xlabel = L"-| \nabla_h \overline{b}|^2 \delta/(M_0^4 f)", axis_kwargs...) 
vary = (xhift(interior(B̅h, :, :, k)))/(f*M²₀^2)
varx = -xhift(interior(D̅ , :, :, k) .* interior(∇b̅, :, :, k).^2)/(M²₀^2)
scatter!(ax_a, vec(varx), vec(vary); markersize = 5, color = :blue, rasterize = true, label = L"\text{data}")
slope = rma_slope(vec(varx), vec(vary))# vec(varx)' * vec(vary)/(vec(varx)' * vec(varx))
lines!(ax_a, [-5, 25], slope*[-5, 25]; color = :red, linewidth = 1, label = "y = "*string(round(slope, sigdigits=2))*"x")
axislegend(ax_a, labelsize=9, framevisible = false, font = texfont(), position = :lt,patchsize = (15, 1), 
padding = (0f0, 0f0, 0f0, 0f0), patchlabelgap = 3, rowgap = 1)
for (j,ax) in enumerate([ax_b,ax_c,ax_d])
    xrange = 1+(j-1)*Nr÷2:(j-1)*Nr÷2+Nr
    yrange = j != 2 ? (1:Nr) : (1+Nr:2Nr)
    scatter!(ax, vec(varx[xrange,yrange]), vec(vary[xrange,yrange]); markersize = 5, color = :blue, rasterize = true, label = L"\text{data}")
    slope = rma_slope(vec(varx[xrange,yrange]), vec(vary[xrange,yrange]))#vec(varx[xrange,yrange])' * vec(vary[xrange,yrange])/(vec(varx[xrange,yrange])' * vec(varx[xrange,yrange]))
    lines!(ax, [-5, 25], slope*[-5, 25]; color = :red, linewidth = 1, label = "y = "*string(round(slope, sigdigits=2))*"x")
    axislegend(ax, labelsize=9, framevisible = false, font = texfont(), position = :lt,patchsize = (15, 1), 
    padding = (0f0, 0f0, 0f0, 0f0), patchlabelgap = 3, rowgap = 1)
end
hidexdecorations!(ax_a, ticks = false)
hidexdecorations!(ax_b, ticks = false)
hideydecorations!(ax_b, ticks = false)
hideydecorations!(ax_d, ticks = false)
rowgap!(gabc, 3)
colgap!(gabc, 1, 15)
resize_to_layout!(fig)
save(filesave * "Bh_Dbd_" * fileparams * "_d$(nday).pdf", fig; pt_per_unit = 1)

#####################
# Compute the horizontal spectrum of ro, rd, Db
axis_kwargs1 = (xlabel = L"\text{Wavenumber (rad m^{-1})}", xgridvisible = false,
                ylabel = L"E_i(k)/E_{i}(k_{\min},z=-2.8~\text{m})", ygridvisible = false,
                xscale = log10, yscale = log10,
                limits = ((4e-5, 3.5e-2), (1e-3,0.5e4)),
                xminorticks = [4e-5:1e-5:9e-5; 2e-4:1e-4:9e-4; 2e-3:1e-3:9e-3; 2e-2:1e-2:6e-2],xminorticksvisible = true)
axis_kwargs2 = NamedTuple{(:xlabel,:xscale,:yscale,:limits,:xgridvisible,:ygridvisible,:xminorticks,:xminorticksvisible)}(axis_kwargs1)
axis_kwargs3 = NamedTuple{(:ylabel,:xscale,:yscale,:limits,:xgridvisible,:ygridvisible,:xminorticks,:xminorticksvisible)}(axis_kwargs1)
axis_kwargs4 = NamedTuple{(:xscale,:yscale,:limits,:xgridvisible,:ygridvisible,:xminorticks,:xminorticksvisible)}(axis_kwargs1)
Nr = length(xv) ÷ 2

fig = Figure(size = (640, 800))
g33 = fig[1, 1] = GridLayout()
S0s = []
alphabet = [letter for letter in 'a':'z'];
for (i,klev) in enumerate([222, 202, 171])
    println("Plotting spectra at z = $(zT[klev])m...")
    if i == 3
        ax = Axis(g33[1, i]; titlealign = :left, title=L"\text{(i)}~z=-60~\text{m}", axis_kwargs4...)
        hideydecorations!(ax, ticks = false)
    else
        if i == 2
            ax = Axis(g33[1, i]; titlealign = :left, title=L"\text{(e)}~z=-25~\text{m}", axis_kwargs4...)
            hideydecorations!(ax, ticks = false)
        else
            ax = Axis(g33[1, 1]; titlealign = :left, title=L"\text{(a)}~z=-2.8~\text{m}", axis_kwargs3...)
        end
    end

    hidexdecorations!(ax, ticks = false)
    Sro = isotropic_powerspectrum(interior(ro, :, :, klev), interior(ro, :, :, klev), x, y)
    Srd = isotropic_powerspectrum(interior(rd, :, :, klev), interior(rd, :, :, klev), x, y)
    Sdb = isotropic_powerspectrum(interior(∇b, :, :, klev), interior(∇b, :, :, klev), x, y)
    SR̅o = isotropic_powerspectrum(interior(R̅o, :, :, klev), interior(R̅o, :, :, klev), x, y)
    SD̅d = isotropic_powerspectrum(interior(D̅, :, :, klev), interior(D̅, :, :, klev), x, y)
    Sd̅b = isotropic_powerspectrum(interior(∇b̅, :, :, klev), interior(∇b̅, :, :, klev), x, y)

    if i == 1
        global So0,Sd0,Sb0 = Sro,Srd,Sdb
    end

    idx = Sro.freq .> 0
    lines!(ax, Sro.freq[idx], 1.5Sro.freq[idx].^0, linestyle = :dash, color = :black)
    text!(ax, 5e-4, 1; text = L"k^0",align = (:left, :top))
    idxf = Sro.freq .> 0
    lines!(ax, Sro.freq[idxf], 2e-2Sro.freq[idxf].^-1, linestyle = :dash, color = :gray)
    text!(ax, 5e-4, 1e2; text = L"k^{-1}",align = (:left, :bottom))
    lines!(ax, Sdb.freq[idx], Real.(Sd̅b.spec[idx]./Sb0.spec[idx][1]), color = :red, label = L"E_{|\nabla_h \overline{b}|}")
    lines!(ax, Sro.freq[idx], Real.(SR̅o.spec[idx]./So0.spec[idx][1]), color = :blue, label = L"E_{\tilde{\zeta}}")
    lines!(ax, Srd.freq[idx], Real.(SD̅d.spec[idx]./Sd0.spec[idx][1]), color = :green, label = L"E_{\tilde{\delta}}")
    lines!(ax, Sdb.freq[idx], Real.(Sdb.spec[idx]./Sb0.spec[idx][1]), color = :red,  linestyle = :dash, label = L"E_{|\nabla_h b|}")
    lines!(ax, Sro.freq[idx], Real.(Sro.spec[idx]./So0.spec[idx][1]), color = :blue, linestyle = :dash, label = L"E_{\zeta}")
    lines!(ax, Srd.freq[idx], Real.(Srd.spec[idx]./Sd0.spec[idx][1]), color = :green, linestyle = :dash, label = L"E_{\delta}")
    if i == 1
        axislegend(ax, labelsize=9, patchsize = (15, 1), framevisible = false, position = :lt,
                   padding = (0f0, 0f0, 0f0, 0f0), patchlabelgap = 3, rowgap = 1)
    end

    dbk = (xhift(interior(∇b, :, :, klev)))
    rok = (xhift(interior(ro, :, :, klev)))
    rdk = (xhift(interior(rd, :, :, klev)))
    R̅ok = (xhift(interior(R̅o, :, :, klev)))
    D̅dk = (xhift(interior(D̅, :, :, klev)))
    db̅k = (xhift(interior(∇b̅, :, :, klev)))
    for j = 1:3
        println([i,j])
        window = 1
        nalpha = 1+j+(i-1)*4
        title = "("*alphabet[nalpha]*") Region " * uppercase(alphabet[j])

        ax = Axis(g33[j+1, i]; titlealign = :left, title=title, titlefont=texfont(),axis_kwargs4...)
        if j < 3
            hidexdecorations!(ax, ticks = false)
        else
            ax.xlabel = L"\text{Wavenumber (rad m^{-1})}"
        end

        if i > 1
            hideydecorations!(ax, ticks = false)
        else
            ax.ylabel = L"E_i(k)/E_{T,v}(k_{\min},z=-2.8~\text{m})"
        end

        xrange = 1+(j-1)*Nr÷2:(j-1)*Nr÷2+Nr
        yrange = j != 2 ? (1:Nr) : (1+Nr:2Nr)
        Sdb = isotropic_powerspectrum(dbk[xrange,yrange], dbk[xrange,yrange], x[xrange], y[yrange];window)
        Sro = isotropic_powerspectrum(rok[xrange,yrange], rok[xrange,yrange], x[xrange], y[yrange];window)
        Srd = isotropic_powerspectrum(rdk[xrange,yrange], rdk[xrange,yrange], x[xrange], y[yrange];window)
        SR̅o = isotropic_powerspectrum(R̅ok[xrange,yrange], R̅ok[xrange,yrange], x[xrange], y[yrange];window)
        SD̅d = isotropic_powerspectrum(D̅dk[xrange,yrange], D̅dk[xrange,yrange], x[xrange], y[yrange];window)
        Sd̅b = isotropic_powerspectrum(db̅k[xrange,yrange], db̅k[xrange,yrange], x[xrange], y[yrange];window)

        idx = Sro.freq .> 0
        if i == 1
            push!(S0s, [Sdb.spec[idx][1], Sro.spec[idx][1], Srd.spec[idx][1]])
        end
        lines!(ax, Sro.freq[idx], 1.5Sro.freq[idx].^0, linestyle = :dash, color = :black)
        text!(ax, 5e-4, 1; text = L"k^0",align = (:left, :top))
        idxf = Sro.freq .> 0
        lines!(ax, Sro.freq[idxf], 2e-2Sro.freq[idxf].^-1, linestyle = :dash, color = :gray)
        text!(ax, 5e-4, 1e2; text = L"k^{-1}",align = (:left, :bottom))
        lines!(ax, Sdb.freq[idx], Real.(Sd̅b.spec[idx]./S0s[j][1]), color = :red, label = L"E_{|\nabla_h b|}")
        lines!(ax, Sro.freq[idx], Real.(SR̅o.spec[idx]./S0s[j][2]), color = :blue, label = L"E_{\tilde{\zeta}}")
        lines!(ax, Srd.freq[idx], Real.(SD̅d.spec[idx]./S0s[j][3]), color = :green, label = L"E_{\tilde{\delta}}")
        lines!(ax, Sdb.freq[idx], Real.(Sdb.spec[idx]./S0s[j][1]), color = :red, linestyle = :dash)
        lines!(ax, Sro.freq[idx], Real.(Sro.spec[idx]./S0s[j][2]), color = :blue, linestyle = :dash)
        lines!(ax, Srd.freq[idx], Real.(Srd.spec[idx]./S0s[j][3]), color = :green, linestyle = :dash)
    end
end
rowgap!(g33, 3)
resize_to_layout!(fig)
save(filesave * "dspectra_" * fileparams * "_d$(nday).pdf", fig; pt_per_unit = 1)
println("Finished plotting spectra, wall time: $((now() - t0).value/1e3) seconds.")

###############################
# Compute spectral vertical boyancy flux 
Nz = length(zT)
wc = compute!(Field(@at (Center, Center, Center) snapshots[:w][snapshot_number]))
b = compute!(Field(α * g * T))
C1 = isotropic_powerspectrum(interior(b, :, :, 1), interior(wc, :, :, 1), xT, yT)
C = zeros(Nz, length(C1.spec))
C[1, :] = real.(C1.spec)
wc1 = (xhift(interior(wc, :, :, 1)))
b1 = (xhift(interior(b, :, :, 1)))
SCs = []
for j = 1:3
    xrange = 1+(j-1)*Nr÷2:(j-1)*Nr÷2+Nr
    yrange = j != 2 ? (1:Nr) : (1+Nr:2Nr)
    SCj = isotropic_powerspectrum(b1[xrange,yrange], wc1[xrange,yrange], xT[xrange], yT[yrange];window=1)
    push!(SCs, SCj)
end
Csub = zeros(3, Nz, length(SCs[1].spec))
for j = 1:3
    Csub[j, 1, :] = real.(SCs[j].spec)
end
println("level 1 done.")
for k = 2:Nz
    Ck = isotropic_powerspectrum(interior(b, :, :, k), interior(wc, :, :, k), xT, yT)
    C[k,:] = real.(Ck.spec)
    wck = (xhift(interior(wc, :, :, k)))
    bk = (xhift(interior(b, :, :, k)))
    for j = 1:3
        xrange = 1+(j-1)*Nr÷2:(j-1)*Nr÷2+Nr
        yrange = j != 2 ? (1:Nr) : (1+Nr:2Nr)
        SCj = isotropic_powerspectrum(bk[xrange,yrange], wck[xrange,yrange], xT[xrange], yT[yrange];window=1)
        Csub[j, k, :] = real.(SCj.spec) 
    end
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

####################################
# Plot x-y slices of KE and Π

# Compute the coarse-grained cross-scale fluxes 
l0 = 4
u̅l, v̅l, w̅l, τ̅uul, τ̅vvl, τ̅wwl, Πhl, Πδl, Πvl, Πgl = coarse_grained_fluxes(snapshots, snapshot_number; cutoff = l0*1kilometer)
#Πhl = compute!(Field(Πhl))
#Πvl = compute!(Field(Πvl))

ke = compute!(Field(0.5 * (u^2 + v^2)))
kτ = compute!(Field(0.5 * (τ̅uul + τ̅vvl)))
#k̅e = compute!(Field(0.5 * (u̅l^2 + v̅l^2 + w̅l^2)))
# Coordinate arrays
x, y, z = nodes(ke)

# Fh = compute!(Field(-((u+u0) * ∂x(ke) + (v+v0) * ∂y(ke))))
# Fv = compute!(Field(-(w * ∂z(ke))))
# Fb = compute!(Field(w*(α * g * T)))

k = 222
var = Πhl
x, y, z = nodes(var)
cmap = :balance
rmin, rmax = -2, 2
fig = Figure(size = (640, 600))
gabc = fig[1, 1] = GridLayout()
axis_kwargs = (ylabel = "y (km)", aspect=1, limits = ((-50, 50), (0, 100)))
ax_a = Axis(gabc[1,1]; titlealign = :left, title=L"\text{(a)}~10^9 \Pi_h^4~\text{(m^2 s^{-3})},~z=-2.8~\text{m}", axis_kwargs...)
ax_b = Axis(gabc[1,2]; titlealign = :left, title=L"\text{(b)~Region A}", aspect=1, limits = ((-50, 0), (25, 75)))
ax_c = Axis(gabc[2,1]; titlealign = :left, title=L"\text{(c)~Region B}", xlabel = "x (km)", ylabel = "y (km)", aspect=1, limits = ((-25, 25), (-25,25))) 
ax_d = Axis(gabc[2,2]; titlealign = :left, title=L"\text{(d)~Region C}", xlabel = "x (km)", aspect=1, limits = ((0, 50), (25, 75))) 
hm_a = heatmap!(ax_a, 1e-3x.-50, 1e-3y, 1e9shift(interior(var,:,:,k)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
hm_b = heatmap!(ax_b, 1e-3x.-50, 1e-3y, 1e9shift(interior(var,:,:,k)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
hm_c = heatmap!(ax_c, 1e-3x.-50, 1e-3y.-75, 1e9xhift(interior(var,:,:,k)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
hm_d = heatmap!(ax_d, 1e-3x.-50, 1e-3y, 1e9shift(interior(var,:,:,k)); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
Colorbar(gabc[1,3], hm_b)
Colorbar(gabc[2,3], hm_d)
poly!(ax_a, [Rect(-50i, 25, 50, 50) for i in 0:1], color = (:white, 0.1), strokecolor = :black, strokewidth = 0.5)
poly!(ax_a, [Rect(-25, 75i, 50, 25) for i in 0:1], color = (:white, 0.1), strokewidth = 0.)
vlines!(ax_a, [-25, 25]; ymin = 0.75, color = :black, linewidth = 0.5)
vlines!(ax_a, [-25, 25]; ymax = 0.25, color = :black, linewidth = 0.5)
text!(ax_a, -50, 75, text = L"\text{A}", color = :black, align = (:left, :top))
text!(ax_a, -25, 25, text = L"\text{B}", color = :black, align = (:left, :top))
if cmap == :thermal
    text!(ax_a, 0, 75, text = L"\text{C}", color = :white, align = (:left, :top))
else
    text!(ax_a, 0, 75, text = L"\text{C}", color = :black, align = (:left, :top))
end
rowgap!(gabc, 3)
colgap!(gabc, 1, 15)
colgap!(gabc, 2, 5)
resize_to_layout!(fig)
save(filesave * "Pifields_" * fileparams * "_d$(nday).pdf", fig; pt_per_unit = 1)

using JLD2 
cgfilename = filesave * "cgfluxes_" * fileparams * "_d$(nday).jld2"
Δh = parameters.Δh/1e3
Nr = length(xu) ÷ 2
ls = [Δh/4:Δh/4:5Δh; 4:10; 12:2:30; 40:10:60]
Πhls = zeros(4,length(ls), length(zu));
Πδls = zeros(4,length(ls), length(zu));
Πvls = zeros(4,length(ls), length(zu));
Πgls = zeros(4,length(ls), length(zu));
for (i, l) in enumerate(ls)
    if l>1.5Δh
        u̅l, v̅l, w̅l, τ̅uul, τ̅vvl, τ̅wwl, Πhl, Πδl, Πvl, Πgl = coarse_grained_fluxes(snapshots, snapshot_number; cutoff = l*1kilometer)
        Πhls[1,i, :] = mean(interior(Πhl, :, :, :), dims = (1,2))
        Πδls[1,i, :] = mean(interior(Πδl, :, :, :), dims = (1,2))
        Πvls[1,i, :] = mean(interior(Πvl, :, :, :), dims = (1,2))
        Πgls[1,i, :] = mean(interior(Πgl, :, :, :), dims = (1,2))
        for k = 1:length(zu)
            Πhk = xhift(interior(Πhl, :, :, k))
            Πδk = xhift(interior(Πδl, :, :, k))
            Πvk = xhift(interior(Πvl, :, :, k))
            Πgk = xhift(interior(Πgl, :, :, k))
            for j = 1:3
                xrange = 1+(j-1)*Nr÷2:(j-1)*Nr÷2+Nr
                yrange = j != 2 ? (1:Nr) : (1+Nr:2Nr)
                Πhls[j+1,i,k] = mean(Πhk[xrange,yrange])
                Πδls[j+1,i,k] = mean(Πδk[xrange,yrange])
                Πvls[j+1,i,k] = mean(Πvk[xrange,yrange])
                Πgls[j+1,i,k] = mean(Πgk[xrange,yrange])
            end
        end
        println("Coarse-graining to $l km...")
        jldsave(cgfilename; ls, z=zu, Πhls, Πδls, Πvls, Πgls)
    end
end


ls = jldopen(cgfilename)["ls"]
z = jldopen(cgfilename)["z"]
Πhls = jldopen(cgfilename)["Πhls"];
Πδls = jldopen(cgfilename)["Πδls"];
Πvls = jldopen(cgfilename)["Πvls"];
Πgls = jldopen(cgfilename)["Πgls"];

axis_kwargs = (limits = ((1/60, 1/0.15), (-100, 0)), xgridvisible = false, ygridvisible = false, 
               ylabel = L"z~\text{(m)}", titlealign = :left,                
               xscale = log10, xticks = ([1/50,1/10,1], ["50⁻¹","10⁻¹","1"]),
               xminorticks = [3e-2:1e-2:9e-2; 0.2:0.1:0.9; 2:9],xminorticksvisible = true)
axis_kwargs1 = NamedTuple{(:xscale,:xticks,:xminorticks,:xminorticksvisible, :ygridvisible,:titlealign)}(axis_kwargs)
fig = Figure(size = (640, 350))
gab = fig[1, 1] = GridLayout()
axlimits = ((1/60, 1/0.15), (-0.3,0.4))
ax = Axis(gab[1,1]; limits = axlimits, title=L"\text{(a)}", ylabel = L"\text{(10^{11} m^2 s^{-2})}", axis_kwargs1...)
ax_b = Axis(gab[2,1]; title=L"\Pi_H~\text{(10^{11} m^2 s^{-2})}",  xlabel = L"l^{-1}~\text{(km^{-1})}", axis_kwargs...)
#ax_c = Axis(gab[3,1]; title=L"10^9\Pi_V", xlabel = L"l^{-1}~\text{(km^{-1})}", axis_kwargs...)
#ax_c = Axis(gab[1,3]; titlealign = :left, limits = ((1/55, 2), (-100, 0)), title=L"\text{(c)}~10^8\Pi_vg^l", axis_kwargs2...)

hm_b = heatmap!(ax_b, ls.^(-1), z, 1e11Πhls[1,:,:]; rasterize = true, colormap = :balance, colorrange = (-3, 3))
hidexdecorations!(ax, ticks = false)
#Colorbar(gab[2, 2], hm_b)
#Colorbar(gab[3, 2], hm_c)

#lines!(ax, ls.^(-1), 1e9vec(mean((Πhs+Πvs)[1,:,z.>=-60],dims=2)), label = L"\Pi", color = :black)
lines!(ax, ls.^(-1), 1e11vec(mean(Πhls[1,:,z.>=-60],dims=2)), label = L"\Pi_H", color = :blue)
lines!(ax, ls.^(-1), 1e11vec(mean((Πhls-Πδls)[1,:,z.>=-60],dims=2)), label = L"\Pi_\alpha", color = :orange)
lines!(ax, ls.^(-1), 1e11vec(mean(Πδls[1,:,z.>=-60],dims=2)), label = L"\Pi_\delta", color = :green)
hlines!(ax, 0, color = :black, linestyle = :dash)
axislegend(ax, position = :lt, patchsize = (15, 1), framevisible = false, fontsize = 9,
padding = (0f0, 0f0, 0f0, 0f0), patchlabelgap = 3, rowgap = 1)

alphabet = [letter for letter in 'a':'z']
for i = 1:3
    nalpha = i+1
    title = "("*alphabet[nalpha]*") Region " * uppercase(alphabet[i])
    ax = Axis(gab[1,i+1]; limits = axlimits, title=title, titlefont=texfont(), axis_kwargs1...)
    ax_b = Axis(gab[2,i+1]; title=L"\Pi_H~\text{(10^{11} m^2 s^{-2})}", xlabel = L"l^{-1}~\text{(km^{-1})}", axis_kwargs...)
    #ax_c = Axis(gab[3,i+1]; title=L"10^8\Pi_V", xlabel = L"l^{-1}~\text{(km^{-1})}", axis_kwargs...)
    #ax_c = Axis(gab[1,3]; titlealign = :left, limits = ((1/55, 2), (-100, 0)), title=L"\text{(c)}~10^8\Pi_vg^l", axis_kwargs2...)

    hm_b = heatmap!(ax_b, ls.^(-1), z, 1e11Πhls[i+1,:,:]; rasterize = true, colormap = :balance, colorrange = (-3, 3))
    #hm_c = heatmap!(ax_c, ls.^(-1), z, 1e9Πvs[i,:,:]; rasterize = true, colormap = :balance, colorrange = (-3, 3))
    #hm_c = heatmap!(ax_c, ls.^(-1), zu, 1e8Πgls; rasterize = true, colormap = :balance, colorrange = (-3, 3))
    #hidexdecorations!(ax_b, ticks = false)
    hidexdecorations!(ax, ticks = false)
    hideydecorations!(ax, ticks = false)
    #hideydecorations!(ax_c, ticks = false)
    hideydecorations!(ax_b, ticks = false)

    #lines!(ax, ls.^(-1), 1e8vec(mean((Πhs+Πvs)[1+i,:,z.>=-50],dims=2)), label = L"\Pi", color = :black)
    lines!(ax, ls.^(-1), 1e11vec(mean(Πhls[1+i,:,z.>=-60],dims=2)), label = L"\Pi_H", color = :blue)
    lines!(ax, ls.^(-1), 1e11vec(mean((Πhls-Πδls)[1+i,:,z.>=-50],dims=2)), label = L"\Pi_\alpha", color = :orange)
    lines!(ax, ls.^(-1), 1e11vec(mean(Πδls[1+i,:,z.>=-60],dims=2)), label = L"\Pi_\delta", color = :green)
    #lines!(ax, ls.^(-1), 1e8vec(mean(Πvs[1+i,:,z.>=-50],dims=2)), label = L"\Pi_V", color = :red)
    #lines!(ax, ls.^(-1), 1e8vec(mean(Πgls[:,zu.>=-50],dims=2)), color = :green)
    hlines!(ax, 0, color = :black, linestyle = :dash)
    colgap!(gab, i, 5)
end
Colorbar(gab[2, 5], hm_b)
#Colorbar(gab[3, 5], hm_c)
colgap!(gab, 4, 3)
rowgap!(gab, 1, 3)
#rowgap!(gab, 2, 3)
# using LESStudySetup.Diagnostics: Ah
# nfactor = 1000
# Ahi,zh,fh = Ah(snapshots, snapshot_number; u0,v0,nfactor)
# println("Computed spectral transfer Ah...")
# Δz = zh[2] - zh[1]
# h = 50
# A_h = sum(Ahi[zh .> -h,:], dims = 1) * Δz / h
# Πhf = cumtrapz(fh, vec(A_h))
# Πhf .-= Πhf[end]
# lines!(ax, 1e3fh/3, 1e6Πhf, color = :green)
rowsize!(gab, 1, Relative(0.3))
save(filesave * "Pihvfields_" * fileparams * "_d$(nday).pdf", fig)

havg(field) = mean(interior(field, :, :, :), dims = (1,2))
Π1 = compute!(Field(-(τuu * ∂x(u̅))));
Π2 = compute!(Field(-(τvv * ∂y(v̅))));
Π4 = compute!(Field(-(τvu * ∂x(v̅))));
Π3 = compute!(Field(-(τuv * ∂y(u̅))));
Π̅1 = mean(interior(Π1, :, :, :), dims = (1,2));
Π̅2 = mean(interior(Π2, :, :, :), dims = (1,2));
Π̅3 = mean(interior(Π3, :, :, :), dims = (1,2));
Π̅4 = mean(interior(Π4, :, :, :), dims = (1,2));
fig = Figure()
ax = Axis(fig[1,1][1,1];aspect=1,xlabel=L"x~\text{(km)}",ylabel=L"y~\text{(km)}")
hm = heatmap!(ax, xT/1e3,yT/1e3,1e9shift(interior(Π3,:,:,217));rasterize=true,colormap=:balance,colorrange=(-5,5))
Colorbar(fig[1,1][1, 2], hm)
ax = Axis(fig[1,2];xlabel=L"\text{(10^9 m^2 s^{-3})}",ylabel = L"z~\text{(m)}")
lines!(ax,1e9 * vec(Π̅1), zu; label = L"-\tau_{uu} \overline{u}_x")
lines!(ax,1e9 * vec(Π̅2), zu; label = L"-\tau_{vv} \overline{v}_y")
lines!(ax,1e9 * vec(Π̅3), zu; label = L"-\tau_{uv} \overline{u}_y")
lines!(ax,1e9 * vec(Π̅4), zu; label = L"-\tau_{vu} \overline{v}_x")
lines!(ax,1e9 * vec(Π̅4+Π̅3+Π̅2+Π̅1), zu; color=:black,label = L"\Pi_h")
hlines!(ax,zT[217];linestyle=:dash)
axislegend(ax,position = :rb)
save(filesave * "Pivtest600_" * fileparams * "_d$(nday).pdf", fig)
# Let's pick the last snapshot!
times = snapshots[:T].times
k = 127
ts, es = [], []
for snapshot_number = 1:4:length(times)
    nday = @sprintf("%2.0f", (times[snapshot_number])/60^2/24)
    println("Plotting snapshot $snapshot_number on day $(nday)...")

    u = snapshots[:u][snapshot_number]
    v = snapshots[:v][snapshot_number]
    ke = compute!(Field(0.5 * (u^2 + v^2)))

    push!(ts, (times[snapshot_number])/60^2/24)
    push!(es, mean(interior(ke, :, :, k)))
end

fig = Figure(size = (640, 560))
gab = fig[1, 1] = GridLayout()
axis_kwargs1 = (xlabel = "x (km)", ylabel = "y (km)", aspect = 1,
            limits = ((0, 100), (0, 100)))
axis_kwargs2 = NamedTuple{(:xlabel,:limits,:aspect)}(axis_kwargs1)

ax_a = Axis(gab[1,1]; titlealign = :left, title=L"\text{(a) Initial mixed layer depth (m)}", axis_kwargs1...)
ax_b = Axis(gab[1,3]; titlealign = :left, title=L"\text{(b) Day10 mixed layer depth (m)}", axis_kwargs2...)
ax_c = Axis(gab[2,1]; titlealign = :left, title=L"\text{(c) Day20 mixed layer depth (m)}", axis_kwargs1...)
ax_d = Axis(gab[2,3]; titlealign = :left, title=L"\text{(d) Day30 mixed layer depth (m)}", axis_kwargs2...)

hidexdecorations!(ax_a, ticks = false)
hidexdecorations!(ax_b, ticks = false)
hideydecorations!(ax_b, ticks = false)
hideydecorations!(ax_d, ticks = false)

hm_a = heatmap!(ax_a, 1e-3xT, 1e-3yT, h[1,:,:]; rasterize = true, colormap = :deep)
hm_b = heatmap!(ax_b, 1e-3xT, 1e-3yT, h[findfirst(t.==10),:,:]; rasterize = true, colormap = :deep)
hm_c = heatmap!(ax_c, 1e-3xT, 1e-3yT, h[findfirst(t.==20),:,:]; rasterize = true, colormap = :deep)
hm_d = heatmap!(ax_d, 1e-3xT, 1e-3yT, h[findfirst(t.==30),:,:]; rasterize = true, colormap = :deep, colorrange = (0, 250))

Colorbar(gab[1, 2], hm_a)
Colorbar(gab[1, 4], hm_b)
Colorbar(gab[2, 2], hm_c)
Colorbar(gab[2, 4], hm_d)
colgap!(gab, 1, 3)
colgap!(gab, 3, 3)
colgap!(gab, 2, 5)
rowgap!(gab, 1, 3)

save(filesave * "insfields_" * fileparams * "_d$(nday).pdf", fig)

##############################

N̅² = mean(compute!(Field(N²(snapshots, snapshot_number))), dims=(1,2))
Γx, Γy = M²(snapshots, snapshot_number)
Γx, Γy = mean(compute!(Field(Γx)), dims=(1,2)), mean(compute!(Field(Γy)), dims=(1,2))
Γxz = compute!(Field(@at (Nothing, Nothing, Center) ∂z(Γx)))
Γyz = compute!(Field(@at (Nothing, Nothing, Center) ∂z(Γy)))

α = parameters.α
g = parameters.g
T = snapshots[:T][snapshot_number]
u = compute!(Field(snapshots[:u][snapshot_number] + u0))
v = compute!(Field(snapshots[:v][snapshot_number] + v0))
w = snapshots[:w][snapshot_number]
xT, yT, zT = nodes(T)
_, _, zw = nodes(w)
B = compute!(Field(α * g * T))
xT, yT = xT .- mean(xT), yT' .- mean(yT)

bz = deepcopy(mean(B, dims=(1,2)))
bx = deepcopy(mean(B, dims=(2)))
by = deepcopy(mean(B, dims=(1)))
bxz = deepcopy(mean(B, dims=(2)))
byz = deepcopy(mean(B, dims=(1)))
d = interior(bz)
d[1,1,:] = cumtrapz(zw, vec(interior(N̅²)))[2:end]
set!(bz, d)
fill_halo_regions!(bz)
d = interior(bx)
set!(bx, xT .* interior(Γx))
fill_halo_regions!(bx)
set!(by, reshape(yT,(1,:,1)) .* interior(Γy))
fill_halo_regions!(by)
set!(bxz, xT .* interior(Γxz))
fill_halo_regions!(bxz)
set!(byz, reshape(yT,(1,:,1)) .* interior(Γyz))
fill_halo_regions!(byz)
b̃ = compute!(Field(B - bz - bx - by))

b̃ᵖ = compute!(Field(b̃ - mean(b̃, dims = (1,2))))
uᵖ = compute!(Field(u - mean(u, dims = (1,2))))
vᵖ = compute!(Field(v - mean(v, dims = (1,2))))

println("compute buoyacy variance term 1...")
var1 = mean(compute!(Field(b̃ᵖ * (u * ∂x(b̃ᵖ) + v * ∂y(b̃ᵖ) + w * ∂z(b̃ᵖ)))), dims = (1,2))
println("compute buoyacy variance term 2...")
var2 = compute!(Field(mean(uᵖ * b̃ᵖ, dims = (1,2)) * Γx + mean(vᵖ * b̃ᵖ, dims = (1,2)) * Γy))
w = compute!(Field(@at (Center, Center, Center) snapshots[:w][snapshot_number]))
wᵖ = compute!(Field(w - mean(w, dims = (1,2))))
wbp = mean(wᵖ * b̃ᵖ, dims = (1,2))
println("compute buoyacy variance term 3...")
var3 = compute!(Field((N̅² + mean(∂z(b̃ᵖ), dims = (1,2))) * wbp))
var4 = mean(compute!(Field(b̃ᵖ * w * (bxz + byz))), dims = (1,2))

# Plot the buoyancy variance terms
fig = Figure(size = (500, 500))
gab = fig[1, 1] = GridLayout()
axis_kwargs = (ylabel = "z (m)", limits = (nothing, (-250, 0)))
ax = Axis(gab[1, 1]; titlealign = :left, title=L"\text{(a)}", xlabel = L"\langle w^\prime \tilde{b}^\prime \rangle~\text{(m^2/s^3)}", axis_kwargs...)
lines!(ax, vec(wbp), zT)
ax = Axis(gab[1, 2]; titlealign = :left, title=L"\text{(b)}", xlabel = "buoyancy variance terms (m²/s)", axis_kwargs...)
lines!(ax, vec(var1), zT; label = L"\langle \mathbf{u} \cdot \nabla \frac{1}{2}\tilde{b}^{\prime 2} \rangle")
lines!(ax, vec(var2), zT; label = L"\langle \mathbf{u}^\prime \tilde{b}^\prime \rangle \cdot \mathbf{\Gamma}")
lines!(ax, vec(var1 + var2), zT; label = L"\langle \mathbf{u} \cdot \nabla \frac{1}{2}\tilde{b}^{\prime 2} \rangle +\langle \mathbf{u}' \tilde{b}' \rangle \cdot \mathbf{\Gamma}")
lines!(ax, vec(var3), zw; label = L"( N^2 +\langle \partial_z \tilde{b}^\prime \rangle) \langle w^\prime \tilde{b}^\prime \rangle")
lines!(ax, vec(var4), zT; label = L"\langle \tilde{b}^{\prime} w \partial_z \mathbf{\Gamma} \cdot \mathbf{x} \rangle")
hideydecorations!(ax, ticks = false, grid = true)
axislegend(ax, position = :rb, framevisible = false)
colgap!(gab, 1, 3)
save(filesave * "bvaraince_" * fileparams * "_d$(nday).pdf", fig)