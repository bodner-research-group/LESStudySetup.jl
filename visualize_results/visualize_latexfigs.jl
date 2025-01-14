using LESStudySetup
using CairoMakie
using SixelTerm
using Printf, Dates
using Oceananigans: compute!
using Oceananigans.BoundaryConditions
using Oceananigans.Grids: xnodes, ynodes, znodes
using Statistics: mean, std
using LESStudySetup.Diagnostics
using LESStudySetup.Oceananigans.Units
using LESStudySetup.Diagnostics: N², M²
using LESStudySetup.Diagnostics: load_snapshots, MLD
using LESStudySetup.Diagnostics: isotropic_powerspectrum, coarse_grained_fluxes, δ
using LESStudySetup.Diagnostics: MixedLayerN², MixedLayerDepth
set_theme!(theme_latexfonts(), fontsize=12,figure_padding = 5)
shift(x) = [x[size(x,1)÷2+1:end, :]; x[1:size(x,1)÷2, :]]

# Load parameters and simulation results
α  = parameters.α
ρ₀ = parameters.ρ₀
grid = T.grid
xT, yT, zT = nodes(T)
g = parameters.g
threshold = 3e-4/g*ρ₀
surface,stratification=false,false
h    = compute!(MixedLayerDepth(grid, (; T); ΔT = abs(threshold / ρ₀ / α), surface,stratification))
stratification = true
Nh² = compute!(MixedLayerN²(grid, (; T); ΔT = abs(threshold / ρ₀ / α), surface,stratification))
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
ax_d = Axis(gab[2,3]; titlealign = :left, title=L"\text{(c)}~R^*_h=(\frac{N_h h}{w_*})^2", axis_kwargs2...)

hidexdecorations!(ax_a, ticks = false)
hidexdecorations!(ax_b, ticks = false)
hideydecorations!(ax_b, ticks = false)
hideydecorations!(ax_d, ticks = false)

k = 1
hm_a = heatmap!(ax_a, 1e-3xT, 1e-3yT, interior(λₛ,:,:,k); rasterize = true, colormap = :amp, colorscale = log10)
hm_b = heatmap!(ax_b, 1e-3xT, 1e-3yT, interior(Rₕ,:,:,k); rasterize = true, colormap = :balance, colorrange = (10^(2.27), 10^(3.73)), colorscale = log10)
hm_c = heatmap!(ax_c, 1e-3xT, 1e-3yT, parameters.f ./ sqrt.(interior(Nh²,:,:,k)); rasterize = true, colormap = Reverse(:balance), colorrange = (10^(-2.11), 10^(-1.39)), colorscale = log10)
hm_d = heatmap!(ax_d, 1e-3xT, 1e-3yT, interior(Rₕw,:,:,k); rasterize = true, colormap = :balance, colorrange = (10^(2.37), 10^(3.63)),colorscale = log10)

Colorbar(gab[1, 2], hm_a)
Colorbar(gab[1, 4], hm_b)
Colorbar(gab[2, 2], hm_c)
Colorbar(gab[2, 4], hm_d)
colgap!(gab, 1, 3)
colgap!(gab, 3, 3)
colgap!(gab, 2, 5)
rowgap!(gab, 1, 3)

save("mldassessfields.pdf", fig)

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
filehead = "./"
filename = filehead * "hydrostatic_snapshots_" * fileparams * ".jld2"
metadata = filehead * "experiment_" * fileparams * "_metadata.jld2"
filesave = "./results/"

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

T0 = snapshots[:T][1]
fig = Figure(size = (640, 315))
gab = fig[1, 1] = GridLayout()
axis_kwargs1 = (xlabel = "x (km)", ylabel = "y (km)", aspect = 1, limits = ((0, 100), (0, 100)))
axis_kwargs2 = NamedTuple{(:xlabel,:limits,:aspect)}(axis_kwargs1)
ax_a = Axis(gab[1,1]; titlealign = :left, title=L"\text{(a)}", axis_kwargs1...)
ax_b = Axis(gab[1,2]; titlealign = :left, title=L"\text{(b)}", axis_kwargs2...)
k = length(zT)
Tmin, Tmax = minimum(interior(T,:,:,k)), maximum(interior(T0,:,:,k))
hm_a = heatmap!(ax_a, 1e-3xT, 1e-3yT, shift(interior(T0,:,:,k)); rasterize = true, colormap = :thermal, colorrange = (Tmin, Tmax))
hm_b = heatmap!(ax_b, 1e-3xT, 1e-3yT, shift(interior(T,:,:,k)); rasterize = true, colormap = :thermal, colorrange = (Tmin, Tmax))
hideydecorations!(ax_b, ticks = false)
Colorbar(gab[1, 3], hm_b)
colgap!(gab, 1, 5)
colgap!(gab, 2, 3)
arrows!(ax_a, [50], [50], [-15*sqrt(3)],[-15],arrowcolor =:black)
save(filesave * "fields_" * fileparams * "_T0T10.pdf", fig)

k = 216
b = compute!(Field(α * g * T))
Bₕ = compute!(Field(-(∂x(b)^2 * ∂x(u) + ∂y(b)^2 * ∂y(v))-∂x(b)*∂y(b)*(∂x(v) + ∂y(u))))
Bᵥ = compute!(Field(-(∂z(b) * (∂x(w) * ∂x(b) + ∂y(w) * ∂y(b)))))
fig = Figure(size = (640, 315))
gab = fig[1, 1] = GridLayout()
axis_kwargs1 = (xlabel = "x (km)", ylabel = "y (km)", aspect = 1, limits = ((0, 100), (0, 100)))
axis_kwargs2 = NamedTuple{(:xlabel,:limits,:aspect)}(axis_kwargs1)
ax_a = Axis(gab[1,1]; titlealign = :left, title=L"\text{(a)}~10^{16}\mathcal{B}_h~\text{({kg}^2 m^{-8}s^{-1})}", axis_kwargs1...)
ax_b = Axis(gab[1,2]; titlealign = :left, title=L"\text{(b)}~10^{16}\mathcal{B}_v~\text{({kg}^2 m^{-8}s^{-1})}", axis_kwargs2...)
hm_a = heatmap!(ax_a, 1e-3xT, 1e-3yT, 1e16*shift(interior(Bₕ,:,:,k)); rasterize = true, colormap = :amp, colorrange = (-0., 1))
hm_b = heatmap!(ax_b, 1e-3xT, 1e-3yT, 1e16*shift(interior(Bᵥ,:,:,k)); rasterize = true, colormap = :amp, colorrange = (-0., 1))
hideydecorations!(ax_b, ticks = false)
Colorbar(gab[1, 3], hm_b)
colgap!(gab, 1, 5)
colgap!(gab, 2, 3)
save(filesave * "fields_" * fileparams * "_BhBv_day10.pdf", fig)

# Plot the fields
ΔN = 160
kT, ks, kw = 222, 222, 202
jslices = [1,2,3]*ΔN
Tmin, Tmax = minimum(interior(T,:,:,kT)), maximum(interior(T,:,:,kT))
vbnd, wbnd = maximum(abs, interior(v)), 0.005

#####################################
var = w
x, y, z = nodes(var)
k = kw
cmap = :balance
rmin, rmax = -wbnd,wbnd#Tmin, Tmax
hcolor = :gray
fig = Figure(size = (640, 410))
gabc = fig[1, 1] = GridLayout()
axis_kwargs = (ylabel = "y (km)", aspect=1, limits = ((0, 100), (0, 100)))
ax_a = Axis(gabc[1,1]; titlealign = :left, title=L"\text{(a)}~w,~z=-25.9~\text{m}", axis_kwargs...)
ax_b = Axis(gabc[1,2]; titlealign = :left, title=L"\text{(b) Region A}", aspect=1, limits = ((50, 100), (25, 75)))
ax_c = Axis(gabc[1,3]; titlealign = :left, title=L"\text{(c) Region B}", aspect=1, limits = ((0, 50), (0, 50))) 
hm_a = heatmap!(ax_a, 1e-3x, 1e-3y, interior(var,:,:,k); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
hm_b = heatmap!(ax_b, 1e-3x, 1e-3y, interior(var,:,:,k); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
hm_c = heatmap!(ax_c, 1e-3x, 1e-3y, interior(var,:,:,k); rasterize = true, colormap = cmap, colorrange = (rmin, rmax))
Colorbar(gabc[1,4], hm_c)
hidexdecorations!(ax_a, ticks = false)
hidexdecorations!(ax_b, ticks = false)
hidexdecorations!(ax_c, ticks = false)
hlines!(ax_a, 1e-3*yT[jslices]; color = :black, linewidth = 0.5)
poly!(ax_a, [Rect(50i, 25i, 50, 50) for i in 0:1], color = (:white, 0.1), strokecolor = :white, strokewidth = 0.5)
text!(ax_a, 50, 75, text = L"\text{A}", color = :black, align = (:left, :top))
if cmap == :thermal
    text!(ax_a, 0, 50, text = L"\text{B}", color = :white, align = (:left, :top))
else
    text!(ax_a, 0, 50, text = L"\text{B}", color = :black, align = (:left, :top))
end

dΔN = 80
jslices = dΔN * [-1,0,1] .+ 2ΔN
hlines!(ax_b, 1e-3*yv[jslices]; color = :black, linewidth = 0.5)
jslices = dΔN * [-1,0,1] .+ ΔN
hlines!(ax_c, 1e-3*yw[jslices]; color = :black, linewidth = 0.5)

zmin = -100
kz = findlast(zw .< zmin)
Nz = length(z)
axis_kwargs0 = (xlabel = "x (km)", ylabel = "z (m)", limits = ((0, 100), (zmin, 0)))
axis_kwargs2 = NamedTuple{(:ylabel,:limits)}(axis_kwargs0)
for j in [1,2,3]
    ja, jb, jc = ΔN*(4-j), dΔN*(2-j)+2ΔN, dΔN*(2-j)+ΔN
    if j == 3
        local ax_a = Axis(gabc[j+1,1]; titlealign = :left, title=L"y=25~\text{km}", axis_kwargs0...)
        local ax_b = Axis(gabc[j+1,2]; titlealign = :left, title=L"y=37.5~\text{km}", xlabel = "x (km)", limits = ((50, 100), (zmin, 0)))
        local ax_c = Axis(gabc[j+1,3]; titlealign = :left, title=L"y=12.5~\text{km}", xlabel = "x (km)", limits = ((0, 50), (zmin, 0)))
    elseif j == 2
        local ax_a = Axis(gabc[j+1,1]; titlealign = :left, title=L"y=50~\text{km}", axis_kwargs2...)
        local ax_b = Axis(gabc[j+1,2]; titlealign = :left, title=L"y=50~\text{km}", limits = ((50, 100), (zmin, 0)))
        local ax_c = Axis(gabc[j+1,3]; titlealign = :left, title=L"y=25~\text{km}", limits = ((0, 50), (zmin, 0)))
    else
        local ax_a = Axis(gabc[j+1,1]; titlealign = :left, title=L"y=75~\text{km}", axis_kwargs2...)
        local ax_b = Axis(gabc[j+1,2]; titlealign = :left, title=L"y=62.5~\text{km}", limits = ((50, 100), (zmin, 0)))
        local ax_c = Axis(gabc[j+1,3]; titlealign = :left, title=L"y=37.5~\text{km}", limits = ((0, 50), (zmin, 0)))
    end
    
    hm_a = heatmap!(ax_a, 1e-3x, z[kz:Nz], interior(var,:,ja,kz:Nz); rasterize = true, colormap = cmap)
    hm_b = heatmap!(ax_b, 1e-3x, z[kz:Nz], interior(var,:,jb,kz:Nz); rasterize = true, colormap = cmap)
    hm_c = heatmap!(ax_c, 1e-3x, z[kz:Nz], interior(var,:,jc,kz:Nz); rasterize = true, colormap = cmap)
    if j < 3
        hidexdecorations!(ax_a, ticks = false)
        hidexdecorations!(ax_b, ticks = false)
        hidexdecorations!(ax_c, ticks = false)
    end
    hideydecorations!(ax_b, ticks = false)
    hideydecorations!(ax_c, ticks = false)
    if cmap == :balance
        hm_a.colorrange = (rmin,rmax)
        hm_b.colorrange = (rmin,rmax)
        hm_c.colorrange = (rmin,rmax)
        if Nz == 224
            Colorbar(gabc[j+1, 4], hm_c, ticks = -0.2:0.2:0.2)
        else
            Colorbar(gabc[j+1, 4], hm_c)
        end
    else
        Colorbar(gabc[j+1, 4], hm_c)
    end
    hlines!(ax_a, [z[k]]; color = :black, linewidth = 0.5)
    hlines!(ax_b, [z[k]]; color = :black, linewidth = 0.5)
    hlines!(ax_c, [z[k]]; color = :black, linewidth = 0.5)
    lines!(ax_a, x/1e3, -interior(h,:,ja); color = hcolor, linewidth = 0.8)
    lines!(ax_b, x/1e3, -interior(h,:,jb); color = hcolor, linewidth = 0.8)
    lines!(ax_c, x/1e3, -interior(h,:,jc); color = hcolor, linewidth = 0.8)
    println(mean(-interior(h,:,ja)))
end
rowgap!(gabc, 3)
colgap!(gabc, 3, 3)
colgap!(gabc, 1, 10)
colgap!(gabc, 2, 10)
for row = 2:4
    rowsize!(gabc, row, Relative(0.15))
end
resize_to_layout!(fig)
save(filesave * "wfields_" * fileparams * "_d$(nday).pdf", fig; pt_per_unit = 1)
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
ax = Axis(gl[1, 1]; titlealign = :left, title=L"(a)~\langle T^{\prime 2} \rangle(t)/\langle T^{\prime 2} \rangle(0,z=-2.8~\text{m})", limits = ((0.125,20),nothing))
lines!(ax, times[3:2:end]/24/3600, tsvT[2:end,1]/tsvT[1,1]; label = L"z=-2.8~\text{m}", linewidth = 1)
lines!(ax, times[3:2:end]/24/3600, tsvT[2:end,2]/tsvT[1,1]; label = L"z=-25~\text{m}", linewidth = 1)
lines!(ax, times[3:2:end]/24/3600, tsvT[2:end,3]/tsvT[1,1]; label = L"z=-60~\text{m}", linewidth = 1)
axislegend(ax, labelsize=9, framevisible = false, orientation = :horizontal,
           padding = (0f0, 0f0, 0f0, 0f0), patchlabelgap = 3)
hidexdecorations!(ax, ticks = false)
titles = [L"(b)~E_T(k,t,z=-2.8~\text{m})/E_T(k_{\min},0,-2.8~\text{m})", L"(c)~E_T(k,t,z=-25~\text{m})/E_T(k_{\min},0,-2.8~\text{m})", L"(d)~E_T(k,t,z=-60~\text{m})/E_T(k_{\min},0,-2.8~\text{m})"]
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
axis_kwargs1 = (xlabel = "Wavenumber (rad⋅m⁻¹)", xgridvisible = false,
                ylabel = L"E_i(k)/E_{T,v}(k_{\min},z=-2.8~\text{m})", ygridvisible = false,
                xscale = log10, yscale = log10,
                limits = ((4e-5, 6e-2), (1e-9,1e2)),
                xminorticks = [4e-5:1e-5:9e-5; 2e-4:1e-4:9e-4; 2e-3:1e-3:9e-3; 2e-2:1e-2:6e-2],xminorticksvisible = true)
axis_kwargs2 = NamedTuple{(:xlabel,:xscale,:yscale,:limits,:xgridvisible,:ygridvisible,:xminorticks,:xminorticksvisible)}(axis_kwargs1)
axis_kwargs3 = NamedTuple{(:ylabel,:xscale,:yscale,:limits,:xgridvisible,:ygridvisible,:xminorticks,:xminorticksvisible)}(axis_kwargs1)
axis_kwargs4 = NamedTuple{(:xscale,:yscale,:limits,:xgridvisible,:ygridvisible,:xminorticks,:xminorticksvisible)}(axis_kwargs1)
Nr = length(xv) ÷ 2
using MathTeXEngine
fig = Figure(size = (640, 660))
g33 = fig[1, 1] = GridLayout()
S0s = []
alphabet = [letter for letter in 'a':'z']
for (i,klev) in enumerate([222, 202, 171])
    println("Plotting spectra at z = $(zT[klev])m...")
    if i == 3
        ax = Axis(g33[i, 1]; titlealign = :left, title=L"\text{(g)}~z=-60~\text{m}", axis_kwargs1...)
    else
        if i == 2
            ax = Axis(g33[i, 1]; titlealign = :left, title=L"\text{(d)}~z=-25~\text{m}", axis_kwargs3...)
        else
            ax = Axis(g33[i, 1]; titlealign = :left, title=L"\text{(a)}~z=-2.8~\text{m}", axis_kwargs3...)
        end
        hidexdecorations!(ax, ticks = false)
    end

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
    idxf = 0 .< Su.freq .<= 1e-3
    lines!(ax, Su.freq[idxf], 1e-13Su.freq[idxf].^-3, linestyle = :dash, color = :gray)
    text!(ax, 10^-4, 1e-3; text = L"k^{-3}")
    idxf = Su.freq .> 1e-3
    lines!(ax, Su.freq[idxf], 1e-9Su.freq[idxf].^(-5/3), linestyle = :dash, color = :gray)
    text!(ax, 10^-3, 10^-5.5; text = L"k^{-5/3}")
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

    for j = 0:1
        println([i,j])
        window = 1
        nalpha = j+2+(i-1)*3
        title = "("*alphabet[nalpha]*") Region " * uppercase(alphabet[j+1])
        if i == 3
            ax = Axis(g33[i, j+2]; titlealign = :left, title=title, titlefont=texfont(),axis_kwargs2...)
        else
            ax = Axis(g33[i, j+2]; titlealign = :left, title=title, titlefont=texfont(),axis_kwargs4...)
            hidexdecorations!(ax, ticks = false)
        end
        hideydecorations!(ax, ticks = false)
        jr = 1-j
        Su = isotropic_powerspectrum(interior(u, 1+jr*Nr:(jr+1)*Nr, 1+jr*Nr÷2:(jr+2)*Nr÷2, klev), interior(u, 1+jr*Nr:(jr+1)*Nr, 1+jr*Nr÷2:(jr+2)*Nr÷2, klev), xu[1+jr*Nr:(jr+1)*Nr], yu[1+jr*Nr÷2:(jr+2)*Nr÷2];window)
        Sv = isotropic_powerspectrum(interior(v, 1+jr*Nr:(jr+1)*Nr, 1+jr*Nr÷2:(jr+2)*Nr÷2, klev), interior(v, 1+jr*Nr:(jr+1)*Nr, 1+jr*Nr÷2:(jr+2)*Nr÷2, klev), xv[1+jr*Nr:(jr+1)*Nr], yv[1+jr*Nr÷2:(jr+2)*Nr÷2];window)
        wk = (interior(w, 1+jr*Nr:(jr+1)*Nr, 1+jr*Nr÷2:(jr+2)*Nr÷2, klev)+interior(w, 1+jr*Nr:(jr+1)*Nr, 1+jr*Nr÷2:(jr+2)*Nr÷2, klev+1))/2
        Sw = isotropic_powerspectrum(wk, wk, xw[1+jr*Nr:(jr+1)*Nr], yw[1+jr*Nr÷2:(jr+2)*Nr÷2];window)
        St = isotropic_powerspectrum(interior(T, 1+jr*Nr:(jr+1)*Nr, 1+jr*Nr÷2:(jr+2)*Nr÷2, klev), interior(T, 1+jr*Nr:(jr+1)*Nr, 1+jr*Nr÷2:(jr+2)*Nr÷2, klev), xT[1+jr*Nr:(jr+1)*Nr], yT[1+jr*Nr÷2:(jr+2)*Nr÷2];window)

        idx = Su.freq .> 0
        if i == 1
            push!(S0s, [St.spec[idx][1], Sv.spec[idx][1]])
        end
        lines!(ax, Su.freq[idx], (10^-6.5)Su.freq[idx].^-2, linestyle = :dash, color = :black)
        idxf = 0 .< Su.freq .<= 1e-3
        lines!(ax, Su.freq[idxf], 1e-13Su.freq[idxf].^-3, linestyle = :dash, color = :gray)
        idxf = Su.freq .> 1e-3
        lines!(ax, Su.freq[idxf], 1e-9Su.freq[idxf].^(-5/3), linestyle = :dash, color = :gray)
        lines!(ax, St.freq[idx], Real.(St.spec[idx]./(S0s[j+1][1])), color = :red, label = L"E_T")
        lines!(ax, Su.freq[idx], Real.(Su.spec[idx]./(S0s[j+1][2])), color = :blue, label = L"E_u")
        lines!(ax, Sv.freq[idx], Real.(Sv.spec[idx]./(S0s[j+1][2])), color = :green, label = L"E_v")
        lines!(ax, Sw.freq[idx], Real.(Sw.spec[idx]./(S0s[j+1][2])), color = :black, label = L"E_w")
        xlims!(ax, (4e-5, 3.5e-2))
        ylims!(ax, (1e-9, 1e1))
    end
end
rowgap!(g33, 3)
resize_to_layout!(fig)
save(filesave * "spectra_" * fileparams * "_d$(nday).pdf", fig; pt_per_unit = 1)
println("Finished plotting spectra, wall time: $((now() - t0).value/1e3) seconds.")

###################################
# Plot x-y slices of ζ/f and δ/f
alphabet = [letter for letter in 'a':'z']
f = parameters.f 
ro = compute!(Field(ζ(snapshots, snapshot_number)/f))
rd = compute!(Field(δ(snapshots, snapshot_number)/f))
∇b = compute!(Field(α * g * (∂x(T)^2 + ∂y(T)^2)^0.5))
# Coordinate arrays
x, y, z = nodes(ro)
k = 222
Nr = length(x) ÷ 2
using MathTeXEngine
fig = Figure(size = (640, 595))
gabc = fig[1, 1] = GridLayout()
axis_kwargs1 = (xlabel = "x (km)", ylabel = "y (km)", aspect = 1,
            limits = ((0, 100), (0, 100)))
axis_kwargs3 = NamedTuple{(:ylabel,:limits,:aspect)}(axis_kwargs1)
axis_kwargs2 = [(xlabel = "x (km)",aspect=1, limits = ((50, 100), (25, 75))),(xlabel = "x (km)",aspect=1, limits = ((0, 50), (0, 50)))]
axis_kwargs4 = [(aspect=1, limits = ((50, 100), (25, 75))),(aspect=1, limits = ((0, 50), (0, 50)))]

ax_db = Axis(gabc[1,1]; titlealign = :left, title=L"\text{(a)}~10^7 |\nabla_h b|, ~z=-2.8~\text{m}", axis_kwargs3...)
ax_ro = Axis(gabc[2,1]; titlealign = :left, title=L"\text{(d)}~\zeta/f, ~z=-2.8~\text{m}", axis_kwargs3...)
ax_rd = Axis(gabc[3,1]; titlealign = :left, title=L"\text{(g)}~\delta/f,~z=-2.8~\text{m}", axis_kwargs1...)

axs = [ax_db, ax_ro, ax_rd]
cmaps = [:binary, :balance, :balance]
ctext = [:black, :white, :white]
cranges = [(0, 6), (-2, 2), (-1, 1)]
for (i, dvar) in enumerate([∇b, ro, rd])
    a = i==1 ? 1e7 : 1
    hm = heatmap!(axs[i], 1e-3x, 1e-3y, a*interior(dvar,:,:,k); rasterize = true, colormap = cmaps[i], colorrange = cranges[i])
    poly!(axs[i], [Rect(50i, 25i, 50, 50) for i in 0:1], color = (:white, 0.1), strokecolor = ctext[i], strokewidth = 0.5)
    text!(axs[i], 50, 75, text = L"\text{A}", color = ctext[i], align = (:left, :top))
    text!(axs[i], 0, 50,  text = L"\text{B}", color = ctext[i], align = (:left, :top))

    for j = 0:1
        nalpha = j+2+(i-1)*3
        title = "("*alphabet[nalpha]*") Region " * uppercase(alphabet[j+1])
        if i < 3
            ax = Axis(gabc[i, 2+j]; titlealign = :left, title=title, titlefont=texfont(), axis_kwargs4[j+1]...)
            hidexdecorations!(ax, ticks = false)
            hidexdecorations!(axs[i], ticks = false)
        else
            ax = Axis(gabc[i, 2+j]; titlealign = :left, title=title, titlefont=texfont(), axis_kwargs2[j+1]...)
        end

        jr = 1-j
        hm = heatmap!(ax, 1e-3x[1+jr*Nr:(jr+1)*Nr], 1e-3y[1+jr*Nr÷2:(jr+2)*Nr÷2], a*interior(dvar,1+jr*Nr:(jr+1)*Nr,1+jr*Nr÷2:(jr+2)*Nr÷2,k); rasterize = true, colormap = cmaps[i], colorrange = cranges[i])
    end
    Colorbar(gabc[i, 4], hm)
end

for i = 1:2
    colgap!(gabc, i,10)
end
colgap!(gabc, 3, 3)
rowgap!(gabc, 3)
resize_to_layout!(fig)
save(filesave * "dfields_" * fileparams * "_d$(nday).pdf", fig)

##############################################
# Compute the horizontal spectrum of |∇b| and ζ/f, δ/f
axis_kwargs1 = (xlabel = "Wavenumber (rad⋅m⁻¹)", xgridvisible = false,
                ylabel = L"E_i(k)/E_{i}(k_{\min},z=-2.8~\text{m})", ygridvisible = false,
                xscale = log10, yscale = log10,
                limits = ((4e-5, 3.5e-2), (1e-3,1e4)),
                xminorticks = [4e-5:1e-5:9e-5; 2e-4:1e-4:9e-4; 2e-3:1e-3:9e-3; 2e-2:1e-2:6e-2],xminorticksvisible = true)
axis_kwargs2 = NamedTuple{(:xlabel,:xscale,:yscale,:limits,:xgridvisible,:ygridvisible,:xminorticks,:xminorticksvisible)}(axis_kwargs1)
axis_kwargs3 = NamedTuple{(:ylabel,:xscale,:yscale,:limits,:xgridvisible,:ygridvisible,:xminorticks,:xminorticksvisible)}(axis_kwargs1)
axis_kwargs4 = NamedTuple{(:xscale,:yscale,:limits,:xgridvisible,:ygridvisible,:xminorticks,:xminorticksvisible)}(axis_kwargs1)
Nr = length(xv) ÷ 2
fig = Figure(size = (640, 660))
g33 = fig[1, 1] = GridLayout()
S0s = []
alphabet = [letter for letter in 'a':'z']
for (i,klev) in enumerate([222, 202, 171])
    println("Plotting spectra at z = $(zT[klev])m...")
    if i == 3
        ax = Axis(g33[i, 1]; titlealign = :left, title=L"\text{(g)}~z=-60~\text{m}", axis_kwargs1...)
    else
        if i == 2
            ax = Axis(g33[i, 1]; titlealign = :left, title=L"\text{(d)}~z=-25~\text{m}", axis_kwargs3...)
        else
            ax = Axis(g33[i, 1]; titlealign = :left, title=L"\text{(a)}~z=-2.8~\text{m}", axis_kwargs3...)
        end
        hidexdecorations!(ax, ticks = false)
    end

    Sro = isotropic_powerspectrum(interior(ro, :, :, klev), interior(ro, :, :, klev), x, y)
    Srd = isotropic_powerspectrum(interior(rd, :, :, klev), interior(rd, :, :, klev), x, y)
    Sdb = isotropic_powerspectrum(interior(∇b, :, :, klev), interior(∇b, :, :, klev), x, y)

    if i == 1
        global So0,Sd0,Sb0 = Sro,Srd,Sdb
    end

    idx = Sro.freq .> 0
    lines!(ax, Sro.freq[idx], 1.5Sro.freq[idx].^0, linestyle = :dash, color = :black)
    text!(ax, 2e-3, 1; text = L"k^0",align = (:left, :top))
    idxf = Sro.freq .<= 1e-3
    lines!(ax, Sro.freq[idxf], 2e-2Sro.freq[idxf].^-1, linestyle = :dash, color = :gray)
    text!(ax, 10^-4, 1e2; text = L"k^{-1}",align = (:left, :top))
    idxf = Sro.freq .> 1e-3
    lines!(ax, Sro.freq[idxf], 4e2Sro.freq[idxf].^(1/3), linestyle = :dash, color = :gray)
    text!(ax, 10^-2.4, 0.8e2; text = L"k^{1/3}")
    lines!(ax, Sdb.freq[idx], Real.(Sdb.spec[idx]./Sb0.spec[idx][1]), color = :red, label = L"E_{|\nabla_h b|}")
    lines!(ax, Sro.freq[idx], Real.(Sro.spec[idx]./So0.spec[idx][1]), color = :blue, label = L"E_{\zeta}")
    lines!(ax, Srd.freq[idx], Real.(Srd.spec[idx]./Sd0.spec[idx][1]), color = :green, label = L"E_{\delta}")
    #xlims!(ax, (8e-5, 4e-2))
    if i == 1
        axislegend(ax, labelsize=9, patchsize = (15, 1), framevisible = false, position = :lt,
                   padding = (0f0, 0f0, 0f0, 0f0), patchlabelgap = 3, rowgap = 1)
    end

    for j = 0:1
        println([i,j])
        window = 1
        nalpha = j+2+(i-1)*3
        title = "("*alphabet[nalpha]*") Region " * uppercase(alphabet[j+1])
        if i == 3
            ax = Axis(g33[i, j+2]; titlealign = :left, title=title, titlefont=texfont(),axis_kwargs2...)
        else
            ax = Axis(g33[i, j+2]; titlealign = :left, title=title, titlefont=texfont(),axis_kwargs4...)
            hidexdecorations!(ax, ticks = false)
        end
        hideydecorations!(ax, ticks = false)
        jr = 1-j
        Sro = isotropic_powerspectrum(interior(ro, 1+jr*Nr:(jr+1)*Nr, 1+jr*Nr÷2:(jr+2)*Nr÷2, klev), interior(ro, 1+jr*Nr:(jr+1)*Nr, 1+jr*Nr÷2:(jr+2)*Nr÷2, klev), x[1+jr*Nr:(jr+1)*Nr], y[1+jr*Nr÷2:(jr+2)*Nr÷2];window)
        Srd = isotropic_powerspectrum(interior(rd, 1+jr*Nr:(jr+1)*Nr, 1+jr*Nr÷2:(jr+2)*Nr÷2, klev), interior(rd, 1+jr*Nr:(jr+1)*Nr, 1+jr*Nr÷2:(jr+2)*Nr÷2, klev), x[1+jr*Nr:(jr+1)*Nr], y[1+jr*Nr÷2:(jr+2)*Nr÷2];window)
        Sdb = isotropic_powerspectrum(interior(∇b, 1+jr*Nr:(jr+1)*Nr, 1+jr*Nr÷2:(jr+2)*Nr÷2, klev), interior(∇b, 1+jr*Nr:(jr+1)*Nr, 1+jr*Nr÷2:(jr+2)*Nr÷2, klev), x[1+jr*Nr:(jr+1)*Nr], y[1+jr*Nr÷2:(jr+2)*Nr÷2];window)

        idx = Sro.freq .> 0
        if i == 1
            push!(S0s, [Sro.spec[idx][1], Srd.spec[idx][1], Sdb.spec[idx][1]])
        end
        idx = Sro.freq .> 0
        lines!(ax, Sro.freq[idx], 1.5Sro.freq[idx].^0, linestyle = :dash, color = :black)
        idxf = Sro.freq .<= 1e-3
        lines!(ax, Sro.freq[idxf], 2e-2Sro.freq[idxf].^-1, linestyle = :dash, color = :gray)
        idxf = Sro.freq .> 1e-3
        lines!(ax, Sro.freq[idxf], 4e2Sro.freq[idxf].^(1/3), linestyle = :dash, color = :gray)

        lines!(ax, Sdb.freq[idx], Real.(Sdb.spec[idx]./(S0s[j+1][3])), color = :red, label = L"E_{|\nabla_h b|}")
        lines!(ax, Sro.freq[idx], Real.(Sro.spec[idx]./(S0s[j+1][1])), color = :blue, label = L"E_{\zeta}")
        lines!(ax, Srd.freq[idx], Real.(Srd.spec[idx]./(S0s[j+1][2])), color = :green, label = L"E_{\delta}")
        #xlims!(ax, (4e-5, 3.5e-2))
        #ylims!(ax, (1e-9, 1e1))
    end
end
rowgap!(g33, 3)
resize_to_layout!(fig)
save(filesave * "dspectra_" * fileparams * "_d$(nday).pdf", fig; pt_per_unit = 1)
println("Finished plotting dspectra, wall time: $((now() - t0).value/1e3) seconds.")

####################################
# Plot x-y slices of KE and Π
initfile = filehead * "hydrostatic_snapshots_hydrostatic_background_dh_100_dz_1.jld2"
initsnaps = load_snapshots(initfile)
u0 = initsnaps[:u][1]
v0 = initsnaps[:v][1]

# Compute the coarse-grained cross-scale fluxes 
l0 = 4
u̅l, v̅l, w̅l, τ̅uul, τ̅vvl, τ̅wwl, Πhl, Πvl = coarse_grained_fluxes(snapshots, snapshot_number; u0, v0, cutoff = l0*1kilometer)
Πhl = compute!(Field(Πhl))
Πvl = compute!(Field(Πvl))

ke = compute!(Field(0.5 * (u^2 + v^2)))
kτ = compute!(Field(0.5 * (τ̅uul + τ̅vvl)))
#k̅e = compute!(Field(0.5 * (u̅l^2 + v̅l^2 + w̅l^2)))
# Coordinate arrays
x, y, z = nodes(ke)

# Fh = compute!(Field(-((u+u0) * ∂x(ke) + (v+v0) * ∂y(ke))))
# Fv = compute!(Field(-(w * ∂z(ke))))
# Fb = compute!(Field(w*(α * g * T)))

vh, vv, vb = 1, 1, 1e-5 
fig = Figure(size = (640, 560))
gab = fig[1, 1] = GridLayout()
axis_kwargs1 = (xlabel = "x (km)", ylabel = "y (km)", aspect = 1,
            limits = ((0, 100), (0, 100)))
axis_kwargs2 = NamedTuple{(:xlabel,:limits,:aspect)}(axis_kwargs1)

ax_a = Axis(gab[1,1]; titlealign = :left, title=L"\text{(a)}~\mathcal{E}, ~z=-3~\text{m}", axis_kwargs1...)
ax_b = Axis(gab[1,3]; titlealign = :left, title=L"\text{(b)}~\mathcal{E}^{\prime 4}, ~z=-3~\text{m}", axis_kwargs2...)
ax_c = Axis(gab[2,1]; titlealign = :left, title=L"\text{(c)}~10^7 \Pi_h^4, ~z=-3~\text{m}", axis_kwargs1...)
ax_d = Axis(gab[2,3]; titlealign = :left, title=L"\text{(b)}~10^7 \Pi_v^4, ~z=-3~\text{m}", axis_kwargs2...)

hidexdecorations!(ax_a, ticks = false)
hidexdecorations!(ax_b, ticks = false)
hideydecorations!(ax_b, ticks = false)
hideydecorations!(ax_d, ticks = false)

k = 198
hm_a = heatmap!(ax_a, 1e-3x, 1e-3y, interior(ke,:,:,k); rasterize = true, colormap = :oslo, colorrange = (1e-3, 1e-1), colorscale = log10)
hm_b = heatmap!(ax_b, 1e-3x, 1e-3y, interior(kτ,:,:,k); rasterize = true, colormap = :oslo, colorrange = (1e-3, 1e-1), colorscale = log10)
hm_c = heatmap!(ax_c, 1e-3x, 1e-3y, 1e7interior(Πhl,:,:,k); rasterize = true, colormap = :balance, colorrange = (-vh, vh))
hm_d = heatmap!(ax_d, 1e-3x, 1e-3y, 1e7interior(Πvl,:,:,k); rasterize = true, colormap = :balance, colorrange = (-vv, vv))

Colorbar(gab[1, 2], hm_a)
Colorbar(gab[1, 4], hm_b)
Colorbar(gab[2, 2], hm_c)
Colorbar(gab[2, 4], hm_d)
colgap!(gab, 1, 3)
colgap!(gab, 3, 3)
colgap!(gab, 2, 5)
rowgap!(gab, 1, 3)

save(filesave * "kefields_" * fileparams * "_d$(nday).pdf", fig)

ls = [0.3; 0.5; 1:60]
Πhls = zeros(length(ls), length(zu))
Πδls = zeros(length(ls), length(zu))
Πvls = zeros(length(ls), length(zu))
Πgls = zeros(length(ls), length(zu))
for (i, l) in enumerate(ls)
    u̅l, v̅l, w̅l, τ̅uul, τ̅vvl, τ̅wwl, Πhl, Πδl, Πvl, Πgl = coarse_grained_fluxes(snapshots, snapshot_number; u0, v0, cutoff = l*1kilometer)
    Πhl = compute!(Field(Πhl))
    Πδl = compute!(Field(Πδl))
    Πvl = compute!(Field(Πvl))
    Πgl = compute!(Field(Πgl))
    Πhls[i, :] = mean(interior(Πhl, :, :, :), dims = (1,2))
    Πδls[i, :] = mean(interior(Πδl, :, :, :), dims = (1,2))
    Πvls[i, :] = mean(interior(Πvl, :, :, :), dims = (1,2))
    Πgls[i, :] = mean(interior(Πgl, :, :, :), dims = (1,2))
    println("Coarse-graining to $l km...")
end

using JLD2 
jldsave("cgfluxes_d12.jld2"; ls, z=zu, Πhls, Πδls, Πvls, Πgls)

ls = jldopen("cgfluxes_d12.jld2")["ls"]
z = jldopen("cgfluxes_d12.jld2")["z"]
Πhls = jldopen("cgfluxes_d12.jld2")["Πhls"]
Πδls = jldopen("cgfluxes_d12.jld2")["Πδls"]
Πvls = jldopen("cgfluxes_d12.jld2")["Πvls"]

fig = Figure(size = (400, 500))
gab = fig[1, 1] = GridLayout()
axis_kwargs = (limits = ((1/55, 4), (-100, 0)), xgridvisible = false, ygridvisible = false, 
               ylabel = L"z~\text{(m)}", titlealign = :left,                
               xscale = log10, xticks = ([1/50,1/10,1], ["50⁻¹","10⁻¹","1"]),
               xminorticks = [3e-2:1e-2:9e-2; 0.2:0.1:0.9],xminorticksvisible = true)
axis_kwargs1 = NamedTuple{(:xscale,:xticks,:xminorticks,:xminorticksvisible,:titlealign)}(axis_kwargs)
ax = Axis(gab[1,1]; limits = ((1/55, 4), nothing), title=L"\text{(a)}", ylabel = L"\text{(10^8 m^2 s^{-2})}", axis_kwargs1...)
ax_b = Axis(gab[2,1]; title=L"\text{(b)}~10^8\Pi_H", axis_kwargs...)
ax_c = Axis(gab[3,1]; title=L"\text{(c)}~10^8\Pi_V", xlabel = L"l^{-1}~\text{(km^{-1})}", axis_kwargs...)
#ax_c = Axis(gab[1,3]; titlealign = :left, limits = ((1/55, 2), (-100, 0)), title=L"\text{(c)}~10^8\Pi_vg^l", axis_kwargs2...)

hm_b = heatmap!(ax_b, ls.^(-1), z, 1e8Πhls; rasterize = true, colormap = :balance, colorrange = (-3, 3))
hm_c = heatmap!(ax_c, ls.^(-1), z, 1e8Πvls; rasterize = true, colormap = :balance, colorrange = (-3, 3))
#hm_c = heatmap!(ax_c, ls.^(-1), zu, 1e8Πgls; rasterize = true, colormap = :balance, colorrange = (-3, 3))
hidexdecorations!(ax_b, ticks = false)
hidexdecorations!(ax, ticks = false)
Colorbar(gab[2, 2], hm_b)
Colorbar(gab[3, 2], hm_c)
colgap!(gab, 1, 3)
rowgap!(gab, 1, 3)
rowgap!(gab, 2, 3)

lines!(ax, ls.^(-1), 1e8vec(mean((Πhls+Πvls)[:,z.>=-50],dims=2)), label = L"\Pi", color = :black)
lines!(ax, ls.^(-1), 1e8vec(mean(Πhls[:,z.>=-50],dims=2)), label = L"\Pi_H", color = :blue)
lines!(ax, ls.^(-1), 1e8vec(mean((Πhls-Πδls)[:,z.>=-50],dims=2)), label = L"\Pi_\alpha", color = :orange)
lines!(ax, ls.^(-1), 1e8vec(mean(Πδls[:,z.>=-50],dims=2)), label = L"\Pi_\delta", color = :green)
lines!(ax, ls.^(-1), 1e8vec(mean(Πvls[:,z.>=-50],dims=2)), label = L"\Pi_V", color = :red)
#lines!(ax, ls.^(-1), 1e8vec(mean(Πgls[:,zu.>=-50],dims=2)), color = :green)
hlines!(ax, 0, color = :black, linestyle = :dash)
axislegend(ax, position = :rb, framevisible = false, orientation = :horizontal)

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

save(filesave * "Pifields_" * fileparams * "_d$(nday).pdf", fig)

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