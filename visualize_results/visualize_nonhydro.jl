using LESStudySetup
using CairoMakie, Makie
using Printf, Dates
using Statistics: mean, std, quantile
using LESStudySetup.Diagnostics
using LESStudySetup.Diagnostics: load_distributed_checkpoint,load_subdomain_snapshot
using LESStudySetup.Diagnostics: isotropic_powerspectrum, coarse_grained_fluxes
using LESStudySetup.Diagnostics: coarse_graining!, TKE, MLD, along_front_averages, MixedLayerDepth
using Oceananigans.Fields: interpolate
using MathTeXEngine
set_theme!(theme_latexfonts(), fontsize=12, figure_padding = 10)
using JLD2, CUDA
set_value!(; Δh = 4.8828125)
filehead = "/orcd/data/abodner/002/shared_datasets/nhyles_output/" 
filesave = "results/"
Q, h₀, ρ₀, cₚ, α, g = 40, 60, parameters.ρ₀, parameters.cp, parameters.α, parameters.g
wₛ = (α * g * Q * h₀ / (ρ₀ * cₚ))^(1/3)
# --- Helper Functions ---
function get_iterations_regex(filehead, fileparam, directory="."; subdirparam="subdomains", rank = 1008)
    subdir = joinpath(directory, filehead * subdirparam)
    
    if !isdir(subdir)
        return Int[]
    end
    
    files = readdir(subdir)
    if subdirparam == "subdomains"
        pattern = Regex("^" * escape_string(fileparam) * "_snapshot_iter(\\d+)\\.jld2\$")
    elseif subdirparam[1:4] == "iter"
        pattern = Regex("^" * escape_string(fileparam) * "_$(rank)_iteration(\\d+)\\.jld2\$")
    end
    
    iterations = Int[]
    
    for file in files
        m = match(pattern, file)
        if m !== nothing
            push!(iterations, parse(Int, m.captures[1]))
        end
    end
    
    return sort(iterations)
end

##########################
f = parameters.f;
α = parameters.α;
g = parameters.g;
# shift(x) = [x[size(x,1)÷2+1:end, :]; x[1:size(x,1)÷2, :]]
# filename0 = "./hydrostatic_snapshots_init.jld2"
# snapshots = load_snapshots(filename0);
# v0 = snapshots[:v][1];
# u0 = snapshots[:u][1];
# T0 = snapshots[:T][1];
# xu, yu, zu = nodes(u0);
# xv, yv, zv = nodes(v0);
# xT, yT, zT = nodes(T0);
# _ , _ , zw = nodes(snapshots[:w][1])
# h0 = MLD(snapshots,1; threshold = 0.09)
# hfront = shift(interior(h0,:,:,1))[320-5:326,:]
# @info "hfront: mean = $(mean(hfront)/60), 10th = $(quantile(vec(hfront), 0.1)/60), 90th = $(quantile(vec(hfront), 0.9)/60)"
# ζ₀ = compute!(Field(∂x(v0) - ∂y(u0)));
# x, y, z = nodes(ζ₀);
# Nz = length(z)
# initfile = "./hydrostatic_snapshots_free.jld2"
# initsnaps = load_snapshots(initfile)
# Ub = initsnaps[:u][1];
# Vb = compute!(Field(initsnaps[:v][1] - snapshots[:v][1]));
# iUb, iVb = interior(Ub), interior(Vb);
# iUb, iVb = [iUb[size(xu,1)÷2+1:end, :, :]; iUb[1:size(xu,1)÷2, :, :]],[iVb[size(xv,1)÷2+1:end, :, :]; iVb[1:size(xv,1)÷2, :, :]];
# ipU, ipV = (xu .- 5e4,vec(yu)), (xv .- 5e4,vec(yv))
# σn = compute!(Field(-(∂x(Ub)-∂y(Vb))/2));
# σs = compute!(Field((∂y(Ub)+∂x(Vb))/2));
# ζb = compute!(Field(∂x(Vb)-∂y(Ub)));
# @info "σn extrema: $(extrema(interior(σn)))"
# @info "σs extrema: $(extrema(interior(σs)))"
# iarrow = 2:32:640
# k = length(zT)

meshgrid(x::AbstractVector, y::AbstractVector) =
        repeat(x, 1, length(y)), repeat(y', length(x), 1)

iterations = Int[]
fileparam = "nonhydrostatic_checkpoint";

# ---------- load the data ----------
rank = 992
Hx,Hy,Hz = 6,6,6
Δh = parameters.Δh
nx = 640;ny = 640;Nz = 224;Lz = 252;
downsample = 8
xgrid = collect(range(-5e4+Δh/2, 5e4-Δh/2; length = 32*nx÷downsample))/1e3;
ygrid = (collect(range(Δh/2, 1e5-Δh/2; length = 32*ny÷downsample)))/1e3;
zgrid = collect(range(-Lz, 0; length = Nz+1));
zgridc = 0.5*(zgrid[1:end-1] + zgrid[2:end]);

output_filename = filehead * "subdomains/sublevels_snapshot_iter164410.jld2";
snapshot = load_subdomain_snapshot(output_filename; variables = ("T","v"), level=224);
top_T = interior(snapshot[:T], 1:downsample:32*nx, 1:downsample:32*ny, 1)'; 
output_filename = filehead * "subdomains/slice_y0_iter164410.jld2";
snapshot = load_subdomain_snapshot(output_filename; variables = ("T","v"));
front_T = interior(snapshot[:T], 1:downsample:32*nx, 1, :);
output_filename = filehead * "subdomains/slice_x5e4_iter164410.jld2";
snapshot = load_subdomain_snapshot(output_filename; variables = ("T","v"));
right_T = interior(snapshot[:T], nx, 1:downsample:32*ny, :);

# ----------------------------- figure & axis ---------------------------------
fig = Figure(size = (640, 320))

# common keyword bundle -------------------------------------------------------
cm_temp = :thermal

clims = (minimum(front_T)+1.0,maximum(top_T))
kw_top  = (colorrange = clims, colormap = cm_temp,
            rasterize = true, shading=NoShading)
kw_front = (colormap = cm_temp, rasterize = true, colorrange = clims, shading=NoShading)
kw_right = (colormap  = cm_temp, rasterize = true, colorrange = clims, shading=NoShading)

# --------------------------- build the slice grids ---------------------------
# 1) constant-z slices (top)
y_bt, x_bt     = meshgrid(ygrid, xgrid);            # Ny × Nx 
z0  = zeros(size(x_bt));  

# 2) constant-y slices (y = Ly/2)
x_xt, z_xt     = meshgrid(xgrid, zgrid[zgrid .> -100]);           # Nx × Nz
y_front        = fill(ygrid[1], size(x_xt));

# 3) constant-x slices (x = 0)
y_yt, z_yt     = meshgrid(ygrid, zgridc[zgridc .> -100]);           # Ny × Nz
x_left         = fill(xgrid[end], size(y_yt));

ax  = Axis3(fig[2, 1]; 
        ylabel = L"x~\text{(km)}", xlabel = L"y~\text{(km)}", zlabel = L"z~\text{(m)}",
        aspect = (1.1, 1, 0.1),
        limits = ((minimum(ygrid), maximum(ygrid)),
                    (minimum(xgrid), maximum(xgrid)),
                    (-100, maximum(zgrid))),
        elevation = 0.15π, azimuth = 0.24π,
        xspinesvisible = false, yspinesvisible = false, zspinesvisible = false,
        xgridvisible = false, ygridvisible = false, zgridvisible = false,
        perspectiveness = 0.7,protrusions = (40,10,0,0))
ax.xreversed = true

# ----------------------------- draw the surfaces -----------------------------  
sft = surface!(ax, y_bt, x_bt, z0;         color = top_T,    kw_top...);   # top
sfw = surface!(ax, y_front,  x_xt,  z_xt;     color = front_T,      kw_front...);   # front 
sfv = surface!(ax, y_yt,    x_left,  z_yt;     color = right_T,      kw_right...);   # right
Colorbar(fig[1, 1], sft, vertical = false, label = L"T~\text{({^\circ}C)}", width = Relative(0.3), tellheight=false)
lines!(ax, [0, 0], [0, 0], [-100, 0]; color = :black, linewidth = 0.5, linestyle = :dash)
lines!(ax, [0, 0], [-3.125, -3.125], [-100, 0]; color = :black, linewidth = 0.5, linestyle = :dash)
lines!(ax, [0, 3.125], [0, 0], [0, 0]; color = :black, linewidth = 0.5, linestyle = :dash)
lines!(ax, [3.125, 3.125], [0, -3.125], [0, 0]; color = :black, linewidth = 0.5, linestyle = :dash)
lines!(ax, [0, 3.125], [-3.125, -3.125], [0, 0]; color = :black, linewidth = 0.5, linestyle = :dash)

xgrid = collect(range(-Δh*nx+Δh/2, -Δh/2; length = nx))/1e3;
ygrid = (collect(range(Δh/2, Δh*ny+Δh/2; length = ny)))/1e3;
filename = filehead * "iteration16x/nonhydrostatic_checkpoint_";
file = jldopen(filename * "$(rank)_iteration164410.jld2");
top2_T = file["NonhydrostaticModel/T/data"][Hx+1:end-Hx, Hy+1:end-Hy, end-Hz]'; 
front2_T = file["NonhydrostaticModel/T/data"][Hx+1:end-Hx, Hy+1, Hz+1:end-Hz][:,zgridc .> -100];
right2_T = file["NonhydrostaticModel/T/data"][end-Hx, Hy+1:end-Hy, Hz+1:end-Hz][:,zgridc .> -100];
close(file)
# --------------------------- build the slice grids ---------------------------
# 1) constant-z slices (top)
y_bt, x_bt     = meshgrid(ygrid, xgrid);            # Ny × Nx 
z0  = zeros(size(x_bt));  

# 2) constant-y slices (y = Ly/2)
x_xt, z_xt     = meshgrid(xgrid, zgrid[zgrid .> -100]);           # Nx × Nz
y_front        = fill(ygrid[1], size(x_xt));

# 3) constant-x slices (x = 0)
y_yt, z_yt     = meshgrid(ygrid, zgridc[zgridc .> -100]);           # Ny × Nz
x_left         = fill(xgrid[end], size(y_yt));

ax  = Axis3(fig[2, 2]; 
        xlabelvisible = false, ylabelvisible = false, zlabelvisible = false,
        aspect = (1.1, 1, 0.75),
        limits = ((minimum(ygrid), maximum(ygrid)),
                    (minimum(xgrid), maximum(xgrid)),
                    (-100, maximum(zgrid))),
        elevation = 0.1π, azimuth = 0.24π,
        xspinesvisible = false, yspinesvisible = false, zspinesvisible = false,
        xgridvisible = false, ygridvisible = false, zgridvisible = false,
        perspectiveness = 0.7,protrusions = (40,10,0,0))
ax.xreversed = true

# ----------------------------- draw the surfaces -----------------------------  
sft = surface!(ax, y_bt, x_bt, z0;         color = top2_T,    kw_top...)   # top
sfw = surface!(ax, y_front,  x_xt,  z_xt;     color = front2_T,      kw_front...)   # w 
sfv = surface!(ax, y_yt,    x_left,  z_yt;     color = right2_T,      kw_right...)   # v
rowgap!(fig.layout, 1, Relative(-0.3))
colgap!(fig.layout, 1, Relative(-0.1))
rowsize!(fig.layout, 1, Relative(0.1))
colsize!(fig.layout, 2, Relative(0.2))
resize_to_layout!(fig)
save(filesave * "T_rank$(rank)_iter164410.pdf", fig; pt_per_unit = 1)

########################################
# fig = Figure(size = (640, 750))
# gab = fig[1, 1] = GridLayout()
# axis_kwargs1 = (titlealign = :left, xlabel = "x (km)", ylabel = "y (km)", aspect = 1, limits = ((-50, 50), (0, 100)))
# axis_kwargs2 = NamedTuple{(:titlealign,:xlabel,:limits,:aspect)}(axis_kwargs1)
# ax_σ = Axis(gab[1,1]; title=L"\text{(a)}~\sigma_n/f,~z=-0.56~\text{m}", axis_kwargs1...)
# ax_ζ = Axis(gab[1,3]; title=L"\text{(b)}~\zeta_e/f,~z=-0.56~\text{m}", axis_kwargs1...)
# ax_T = Axis(gab[3,1]; title=L"\text{(c)}~T_i~\text{({^\circ}C)},~z=-0.56~\text{m}", axis_kwargs1...)
# ax_v = Axis(gab[3,3]; title=L"\text{(d)}~v_0~\text{(m~s^{-1})},~z=-0.56~\text{m}", axis_kwargs1...)
# σmin, σmax = extrema(interior(σn,:,:,k))
# hm_σ = heatmap!(ax_σ, 1e-3xT.-50, 1e-3yT, shift(interior(σn,:,:,k))/f; rasterize = true, colormap = :diff, colorrange = (1.1σmin/f, 1.1σmax/f))
# Colorbar(gab[1, 2], hm_σ)
# hidexdecorations!(ax_σ, ticks = false)
# arrows!(ax_σ,xT[iarrow]/1e3.-50, yT[iarrow]/1e3, shift(interior(Ub,iarrow,iarrow,k)), shift(interior(Vb,iarrow,iarrow,k)), arrowsize = 3, lengthscale = 1e2,linecolor = :white, arrowcolor = :white, linewidth = 0.2)
# hm_ζ = heatmap!(ax_ζ, 1e-3xT.-50, 1e-3yT, shift(interior(ζb,:,:,k))/f; rasterize = true, colormap = :curl)
# Colorbar(gab[1, 4], hm_ζ)
# hidexdecorations!(ax_ζ, ticks = false)
# hideydecorations!(ax_ζ, ticks = false)
# Tmin, Tmax = minimum(interior(T0,:,:,k))-0.5, maximum(interior(T0,:,:,k))
# hm_T = heatmap!(ax_T, 1e-3xT.-50, 1e-3yT, shift(interior(T0,:,:,k)); rasterize = true, colormap = :thermal, colorrange = (Tmin, Tmax))
# hidexdecorations!(ax_T, ticks = false)
# Colorbar(gab[3, 2], hm_T)
# arrows!(ax_T, [0], [50], [-15*sqrt(3)],[-15],linecolor = :black, arrowcolor =:black)
# text!(ax_T, -10*sqrt(3), 50, text = L"\text{Wind}", color = :black, align = (:right, :top))
# vbnd = maximum(abs.(interior(v0,:,:,k)))
# hm_v = heatmap!(ax_v, 1e-3xv.-50, 1e-3yv, shift(interior(v0,:,:,k)); rasterize = true, colormap = :balance, colorrange = (-vbnd, vbnd))
# hideydecorations!(ax_v, ticks = false)
# hidexdecorations!(ax_v, ticks = false)
# Colorbar(gab[3, 4], hm_v)
# arrows!(ax_v,xT[iarrow]/1e3.-50, yT[iarrow]/1e3, shift(interior(Ub,iarrow,iarrow,k)), shift(interior(Vb,iarrow,iarrow,k)), arrowsize = 3, lengthscale = 1e2,linecolor = :black, arrowcolor = :black, linewidth = 0.2)
# text!(ax_T, 25, 25, text = L"\text{Warm eddy}", color = :red, align = (:center, :center))
# text!(ax_T, 25, 75, text = L"\text{Cold eddy}", color = :blue, align = (:center, :center))
# text!(ax_T, -25, 75, text = L"\text{Warm eddy}", color = :red, align = (:center, :center))
# text!(ax_T, -25, 25, text = L"\text{Cold eddy}", color = :blue, align = (:center, :center))

# zmin = -100
# ΔN = 160
# ja = 2ΔN
# kz = findlast(zw .< zmin)
# Nz = length(zT)
# axis_kwargs0 = (xlabel = "x (km)", ylabel = "z (m)", limits = ((-50, 50), (zmin, 0)))
# axis_kwargs1 = NamedTuple{(:xlabel,:limits)}(axis_kwargs0)
# ax_a = Axis(gab[2,1]; titlealign = :left, title=L"y=0~\text{km}", axis_kwargs0...)
# ax_b = Axis(gab[2,3]; titlealign = :left, title=L"y=50~\text{km}", axis_kwargs0...)
# ax_c = Axis(gab[4,1]; titlealign = :left, title=L"y=50~\text{km}", axis_kwargs0...)
# ax_d = Axis(gab[4,3]; titlealign = :left, title=L"y=50~\text{km}", axis_kwargs0...)
# hm_a = heatmap!(ax_a, 1e-3xT.-50, zT[kz:Nz], shift(interior(σn, :, 1, kz:Nz))/f; rasterize = true, colormap = :diff, colorrange = (1.1σmin/f, 1.1σmax/f))
# Colorbar(gab[2, 2], hm_a)
# hm_b = heatmap!(ax_b, 1e-3xT.-50, zT[kz:Nz], shift(interior(ζb, :, ΔN, kz:Nz))/f; rasterize = true, colormap = :curl)
# hideydecorations!(ax_b, ticks = false)
# Colorbar(gab[2, 4], hm_b)
# Tmin, Tmax = minimum(interior(T0,:,ja,kz:Nz)), maximum(interior(T0,:,ja,kz:Nz))
# hm_c = heatmap!(ax_c, 1e-3xT.-50, zT[kz:Nz], shift(interior(T0,:,ja,kz:Nz)); rasterize = true, colormap = :thermal, colorrange = (Tmin, Tmax))
# Colorbar(gab[4, 2], hm_c)
# lines!(ax_c, 1e-3xT.-50, -vec(shift(interior(h0,:,ja))); color = :white, linewidth = 0.8)
# hm_d = heatmap!(ax_d, 1e-3xv.-50, zv[kz:Nz], shift(interior(v0,:,ja,kz:Nz)); rasterize = true, colormap = :balance, colorrange = (-vbnd, vbnd))
# hideydecorations!(ax_d, ticks = false)
# idxz = 1:10:Nz
# arrows!(ax_a,xu[iarrow]/1e3.-50, zu[idxz], shift(interior(Ub,iarrow,1,idxz)), 0*shift(interior(Ub,iarrow,1,idxz)), arrowsize = 3, lengthscale = 1e2,linecolor = :white, arrowcolor = :white, linewidth = 0.2)
# arrows!(ax_d,xu[iarrow]/1e3.-50, zu[idxz], shift(interior(Ub,iarrow,ja,idxz)), 0*shift(interior(Ub,iarrow,ja,idxz)), arrowsize = 3, lengthscale = 1e2,linecolor = :black, arrowcolor = :black, linewidth = 0.2)
# Colorbar(gab[4, 4], hm_d,ticks = -0.2:0.2:0.2)

# rowgap!(gab, 1, 0)
# rowgap!(gab, 3, 0)
# rowgap!(gab, 2, 5)
# colgap!(gab, 1, 1)
# colgap!(gab, 3, 1)
# colgap!(gab, 2, 10)
# for row = [2,4]
#     rowsize!(gab, row, Relative(0.1))
# end
# resize_to_layout!(fig)
# save(filesave * "SnZeT0v0fields.pdf", fig)

#################
# # Load parameters and simulation results
# T0 = snapshots[:T][1]
# ρ₀ = parameters.ρ₀
# threshold = 3e-4/g*ρ₀
# surface,stratification=true,false
# h    = MixedLayerDepth(T0.grid, (; T=T0); ΔT = abs(threshold / ρ₀ / α), surface,stratification)
# stratification = true
# Nh² = MixedLayerDepth(T0.grid, (; T=T0); ΔT = abs(threshold / ρ₀ / α), surface,stratification)
# us = sqrt(parameters.τw/ρ₀)
# Qᵘ = parameters.Q/(ρ₀ * parameters.cp) * α * g
# λₛ=compute!(Field(Qᵘ*h/us^3))
# Rₕ=compute!(Field(Nh² * h^2/us^2))
# Rₕw = compute!(Field(Rₕ/λₛ^(2/3)))

# fig = Figure(size = (640, 590))
# gab = fig[1, 1] = GridLayout()
# axis_kwargs1 = (xlabel = "x (km)", ylabel = "y (km)", aspect = 1,
#                 limits = ((-50, 50), (0, 100)))
# axis_kwargs2 = NamedTuple{(:xlabel,:limits,:aspect)}(axis_kwargs1)

# ax_a = Axis(gab[1,1]; titlealign = :left, title=L"\text{(a)}~\lambda_s=\frac{-B_0 h}{u_*^3}", axis_kwargs1...)
# ax_b = Axis(gab[1,3]; titlealign = :left, title=L"\text{(b)}~R_h=(\frac{N_h h}{u_*})^2", axis_kwargs2...)
# ax_c = Axis(gab[2,1]; titlealign = :left, title=L"\text{(c)}~f/N_h", axis_kwargs1...)
# ax_d = Axis(gab[2,3]; titlealign = :left, title=L"\text{(d)}~R^*_h=(\frac{N_h h}{w_*})^2", axis_kwargs2...)

# hidexdecorations!(ax_a, ticks = false)
# hidexdecorations!(ax_b, ticks = false)
# hideydecorations!(ax_b, ticks = false)
# hideydecorations!(ax_d, ticks = false)

# hm_a = heatmap!(ax_a, 1e-3xT.-50, 1e-3yT, shift(interior(λₛ,:,:,1)); rasterize = true, colormap = :amp, colorscale = log10)
# hm_b = heatmap!(ax_b, 1e-3xT.-50, 1e-3yT, shift(interior(Rₕ,:,:,1)); rasterize = true, colormap = :balance, colorrange = (10^(2.27), 10^(3.73)), colorscale = log10)
# hm_c = heatmap!(ax_c, 1e-3xT.-50, 1e-3yT, f ./ sqrt.(shift(interior(Nh²,:,:,1))); rasterize = true, colormap = Reverse(:balance), colorrange = (10^(-2.11), 10^(-1.39)), colorscale = log10)
# hm_d = heatmap!(ax_d, 1e-3xT.-50, 1e-3yT, shift(interior(Rₕw,:,:,1)); rasterize = true, colormap = :balance, colorrange = (10^(2.76), 10^(4.2)),colorscale = log10)

# Colorbar(gab[1, 2], hm_a)
# Colorbar(gab[1, 4], hm_b)
# Colorbar(gab[2, 2], hm_c)
# Colorbar(gab[2, 4], hm_d)
# colgap!(gab, 1, 0)
# colgap!(gab, 3, 0)
# colgap!(gab, 2, 5)
# rowgap!(gab, 1, 3)

# save(filesave * "fields_init_mldass_d0.pdf", fig)

# use_gpu = true
# A=rand(Float32,20480,20480,1);
# if use_gpu && CUDA.functional()
#     p = CUDA.CUFFT.plan_rfft(CuArray(A), (1,2));
#     pk = CUDA.CUFFT.plan_rfft(CuArray(A[:,:,1:1]), (1,2));
#     ip = CUDA.CUFFT.plan_irfft(p * CuArray(A), 20480, (1,2));
# else
#     @error "CUDA not available"
# end
# fileparam = "sublevels"
# output_filename = filehead * "subdomains/" * fileparam * "_snapshot_iter37003.jld2"
# snapshot = load_subdomain_snapshot(output_filename;variables = ("v","u"),level=224);
# grid = snapshot[:u].grid
# u̅  = XFaceField(grid);
# v̅  = YFaceField(grid);
# _, _, zi = nodes(snapshot[:v]);
# Nzi = length(zi)
# coarse_graining!(snapshot[:v] , v̅; kernel=:gaussian, cutoff=300, border = :circular, method = :spectral, use_gpu,plans=(p,ip));
# coarse_graining!(snapshot[:u] , u̅; kernel=:gaussian, cutoff=300, border = :circular, method = :spectral, use_gpu,plans=(p,ip));
# ζ̄i = compute!(Field((∂x(v̅-mean(v̅,dims=2))-∂y(u̅))));

# filename = "./hydrostatic_snapshots_40_Q.jld2"
# snapshots = load_snapshots(filename);
# times = snapshots[:T].times
# v30 = snapshots[:v][length(times)];
# u30 = snapshots[:u][length(times)];
# ζ̄i = compute!(Field((∂x(v30-mean(v30,dims=2))-∂y(u30))));

# fig = Figure(size = (640, 600))
# g4 = fig[1, 1] = GridLayout()
# aspect = 1
# crange = (-3, 3)
# axis_kwargs = (xlabel = L"x~\text{(km)}", ylabel = L"y~\text{(km)}", limits = ((-50,50),(0,100)), aspect=aspect)
# ax = Axis(g4[1,1]; title=L"{\zeta}^s/f,~t=30~\text{h},~z=0.56~\text{m}", axis_kwargs...) 
# xi, yi, zi = nodes(ζ̄i);
# Nzi = length(zi)
# hm = heatmap!(ax, 1e-3xi .- 50, 1e-3yi, shift(interior(ζ̄i,:,:,Nzi))./parameters.f; rasterize = true, colormap = :curl, colorrange = crange)
# arrows!(ax,xT[iarrow]/1e3 .- 50, yT[iarrow]/1e3, shift(interior(Ub,iarrow,iarrow,k)), shift(interior(Vb,iarrow,iarrow,k)), arrowsize = 3, lengthscale = 1e2,linecolor = :black, arrowcolor = :black, linewidth = 0.2,alpha=0.5)
# Colorbar(g4[1,2], hm)
# colgap!(g4, 1, 1)
# resize_to_layout!(fig)
# save(filesave * "curls_surface_hydrostatic_t30h.pdf", fig; pt_per_unit = 1)

# use_gpu = true
# A=rand(Float32,10240,20480,1);
# if use_gpu && CUDA.functional()
#     p = CUDA.CUFFT.plan_rfft(CuArray(A), (1,2));
#     pk = CUDA.CUFFT.plan_rfft(CuArray(A[:,:,1:1]), (1,2));
#     ip = CUDA.CUFFT.plan_irfft(p * CuArray(A), 10240, (1,2));
# else
#     @error "CUDA not available"
# end
# fileparam = "xband1sublevels"
# fig = Figure(size = (640, 580))
# g4 = fig[1, 1] = GridLayout()
# aspect = 0.25
# crange = (-1, 1)
# axis_kwargs = (xlabel = L"x~\text{(km)}", limits = ((-12.5,12.5),(0,100)), aspect=aspect)
# ax_a = Axis(g4[1,1]; titlealign = :left, title=L"\text{(a)}~\delta_0/f", ylabel = L"y~\text{(km)}", axis_kwargs...)
# ax_b = Axis(g4[1,2]; titlealign = :left, title=L"\text{(b)}~\overline{\delta}/f,~t=16~\text{h}", axis_kwargs...)
# ax_c = Axis(g4[1,3]; titlealign = :left, title=L"\text{(c)}~\overline{\delta}/f,~t=30~\text{h}", axis_kwargs...) 
# ax_d = Axis(g4[1,4]; titlealign = :left, title=L"\text{(d)}~\overline{\delta}/f,~t=44~\text{h}", axis_kwargs...) 
# hideydecorations!(ax_b, ticks = false)
# hideydecorations!(ax_c, ticks = false)
# hideydecorations!(ax_d, ticks = false)
# ζ₀ = compute!(Field(∂y(v0) + ∂x(u0)));
# hm_a = heatmap!(ax_a, 1e-3x .- 50, 1e-3y, shift(interior(ζ₀,:,:,Nz))./f; rasterize = true, colormap = :diff, colorrange = crange)
# arrows!(ax_a,xT[iarrow]/1e3.-50, yT[iarrow]/1e3, shift(interior(Ub,iarrow,iarrow,k)), shift(interior(Vb,iarrow,iarrow,k)), arrowsize = 3, lengthscale = 1e2,linecolor = :black, arrowcolor = :black, linewidth = 0.2,alpha=0.5)
# axs = [ax_b, ax_c, ax_d]
# iterations = [26066,37003,49086]
# for (i,iteration) in enumerate(iterations)
#     output_filename = filehead * "subdomains/" * fileparam * "_snapshot_iter$(iteration).jld2"
#     snapshot = load_subdomain_snapshot(output_filename;variables = ("v","u"),level=224);
#     grid = snapshot[:u].grid
#     u̅  = XFaceField(grid);
#     v̅  = YFaceField(grid);
#     _, _, zi = nodes(snapshot[:v]);
#     Nzi = length(zi)
#     coarse_graining!(snapshot[:v] , v̅; kernel=:gaussian, cutoff=300, border = :ycircular, method = :spectral, use_gpu,plans=(p,ip));
#     coarse_graining!(snapshot[:u] , u̅; kernel=:gaussian, cutoff=300, border = :ycircular, method = :spectral, use_gpu,plans=(p,ip));
#     # ζ̄i = compute!(Field((∂x(v̅)-∂y(u̅))));
#     ζ̄i = compute!(Field((∂y(v̅)+∂x(u̅))));
#     xi, yi, _ = nodes(ζ̄i);
#     hm_i = heatmap!(axs[i], 1e-3xi, 1e-3yi, (interior(ζ̄i,:,:,Nzi))./f; rasterize = true, colormap = :diff, colorrange = crange)
#     arrows!(axs[i],xT[iarrow]/1e3.-50, yT[iarrow]/1e3, shift(interior(Ub,iarrow,iarrow,k)), shift(interior(Vb,iarrow,iarrow,k)), arrowsize = 3, lengthscale = 1e2,linecolor = :black, arrowcolor = :black, linewidth = 0.2,alpha=0.5)
# end
# Colorbar(g4[1,5], hm_a)
# colgap!(g4, 1, 5)
# colgap!(g4, 2, 5)
# colgap!(g4, 3, 5)
# colgap!(g4, 4, 1)
# resize_to_layout!(fig)
# save(filesave * "div_" * fileparam * "_0_16_30_44h_cg3hm.pdf", fig; pt_per_unit = 1)

# meshgrid(x::AbstractVector, y::AbstractVector) =
#         repeat(x, 1, length(y)), repeat(y', length(x), 1)

# iterations = Int[]
# fileparam = "nonhydrostatic_checkpoint";
# for i = 3:7
#     subdirparam = "iteration$(i)x";
#     iterations = [iterations; get_iterations_regex(filehead, fileparam;subdirparam)]
# end

# # ----------------------------- figure & axis ---------------------------------
# fig = Figure(size = (640, 320))

# # common keyword bundle -------------------------------------------------------
# cm_temp = :thermal
# cm_w    = :delta
# cm_uv   = :balance

# clims = (19.7,20.2)
# kw_temp  = (colorrange = clims, colormap = cm_temp,
#             rasterize = true, shading=NoShading)
# kw_w   = (colormap = cm_w, rasterize = true, colorrange = (-0.02,0.02), shading=NoShading)
# kw_uv   = (colormap   = cm_uv, rasterize = true, colorrange = (-0.1,0.1), shading=NoShading)
# Δh = parameters.Δh
# nx = 640;ny = 640;Nz = 224;Lz = 252;
# xgrid = collect(range(-Δh*nx+Δh/2, -Δh/2; length = nx))/1e3;
# ygrid = (collect(range(Δh/2, Δh*ny+Δh/2; length = ny)) .+ parameters.Ly/2)/1e3;
# zgrid = collect(range(-Lz, 0; length = Nz+1));
# zgridc = 0.5*(zgrid[1:end-1] + zgrid[2:end]);
# Hx,Hy,Hz = 6,6,6

# # --------------------------- build the slice grids ---------------------------
# # 1) constant-z slices (top)
# y_bt, x_bt     = meshgrid(ygrid, xgrid);            # Ny × Nx 
# z0  = zeros(size(x_bt));  

# # 2) constant-y slices (y = Ly/2)
# x_xt, z_xt     = meshgrid(xgrid, zgrid[zgrid .> -100]);           # Nx × Nz
# y_front        = fill(ygrid[1], size(x_xt));

# # 3) constant-x slices (x = 0)
# y_yt, z_yt     = meshgrid(ygrid, zgridc[zgridc .> -100]);           # Ny × Nz
# x_left         = fill(xgrid[end], size(y_yt));

# n = Observable(1)
# ax  = Axis3(fig[2, 2]; 
#         ylabel = L"\text{Cross-front}~x~\text{(km)}", xlabel = L"\text{Along-front}~y~\text{(km)}", zlabel = L"\text{Vertical}~z~\text{(m)}",
#         aspect = (1.1, 1, 0.3),
#         limits = ((minimum(ygrid), maximum(ygrid)),
#                     (minimum(xgrid), maximum(xgrid)),
#                     (-100, maximum(zgrid))),
#         elevation = 0.3, azimuth = 0.24π,
#         xspinesvisible = false, yspinesvisible = false, zspinesvisible = false,
#         xgridvisible = false, ygridvisible = false, zgridvisible = false,
#         perspectiveness = 0.7,protrusions = (40,10,0,0))
# ax.xreversed = true

# # ---------- load the data ----------
# filename = @lift filehead * "iteration" * "$(iterations[$n])"[1] * "x/nonhydrostatic_checkpoint_";
# file = @lift jldopen($filename * "1008_iteration$(iterations[$n]).jld2");

# top_T = @lift $file["NonhydrostaticModel/T/data"][Hx+1:end-Hx, Hy+1:end-Hy, end-Hz]'; 
# front_w = @lift $file["NonhydrostaticModel/w/data"][Hx+1:end-Hx, Hy+1, Hz+1:end-Hz][:,zgrid .> -100];
# right_v = @lift $file["NonhydrostaticModel/v/data"][end-Hx, Hy+1:end-Hy, Hz+1:end-Hz][:,zgridc .> -100];

# # ----------------------------- draw the surfaces -----------------------------  
# sft = surface!(ax, y_bt, x_bt, z0;         color = top_T,    kw_temp...)   # top
# sfw = surface!(ax, y_front,  x_xt,  z_xt;     color = front_w,      kw_w...)   # w 
# sfv = surface!(ax, y_yt,    x_left,  z_yt;     color = right_v,      kw_uv...)   # v
# Colorbar(fig[1, 2], sft, vertical = false, label = @lift("t = $($n+23) h"))
# Colorbar(fig[2, 1], sfw, height = Relative(0.3), tellheight=false, ticks = -0.02:0.01:0.02)
# Colorbar(fig[2, 3], sfv, height = Relative(0.3), tellheight=false)
# #title = L"\text{Buoyancy evolution}~b^0(x,z,t)"
# #fig[1, 1:2] = Label(fig, title; tellwidth = false, padding = (0, 0, -120, 0))
# rowgap!(fig.layout, 1, Relative(-0.3))
# #colgap!(fig.layout, 1, Relative(0))    
# resize_to_layout!(fig)
# frames = 1:length(iterations)
# t0 = now()
# @info "Making a neat animation of rank 1008..."
# record(fig, filesave * "rank1008_" * fileparam * ".mp4", frames, framerate=4) do i
#     println("Loading fields $i wall time: $((now() - t0).value/1e3) seconds.")
#     n[] = i
# end

# save(filesave * "Twv_rank$(rank)_.pdf", fig; pt_per_unit = 1)
##############################
# # --- Simulation and Subdomain Parameters ---
# fileparam = "xband1surf9"
iteration = 37003

# # 1. Define the filename of the saved snapshot
# output_filename = filehead * "subdomains/" * fileparam * "_snapshot_iter$(iteration).jld2"

# # 2. Load the snapshot using the new function
# snapshot = load_subdomain_snapshot(output_filename; variables = ("u", "v", "w"));

# u̅, v̅, w̅, τuu, τvv, τww = TKE(snapshot; cutoff=300, border=:ycircular, Lx = snapshot[:grid].Lx, Ly = snapshot[:grid].Ly,
#                                        method=:spectral,use_gpu,plans=(p,pk,ip));
# @info "τuu extrema: $(extrema(interior(τuu)))"
# @info "τvv extrema: $(extrema(interior(τvv)))"
# @info "τww extrema: $(extrema(interior(τww)))"

# # 3. Plot the TKE fields
# Q, h₀, ρ₀, cₚ, α, g = 40, 60, parameters.ρ₀, parameters.cp, parameters.α, parameters.g
# wₛ = (α * g * Q * h₀ / (ρ₀ * cₚ))^(1/3)
# TKEₛ = compute!(Field((τvv + τuu + τww)/2/wₛ^2));
# SKEₛ = compute!(Field(((u̅ - mean(u̅;dims=(1,2)))^2 + (v̅ - mean(v̅;dims=(1,2)))^2 + (w̅ - mean(w̅;dims=(1,2)))^2)/2/wₛ^2));
# using Makie
# fig = Figure(size = (640, 640))
# gab = fig[1, 1] = GridLayout()
# idxEm = argmax(interior(TKEₛ))
# var = interior(TKEₛ, :, :, idxEm[3]);
# x, y, z = nodes(TKEₛ);
# axis_kwargs = (titlealign = :left,xlabel=L"x~\text{(km)}",limits=((-12.5,12.5),(0,100)))
# ax_a = Axis(gab[1,1]; title=L"\text{(a) TKE}/w_*^2", ylabel=L"y~\text{(km)}", axis_kwargs...)
# ax_b = Axis(gab[1,3]; title=L"\text{(b) SKE}/w_*^2", axis_kwargs...)
# hm_a = heatmap!(ax_a, 1e-3x, 1e-3y, var; rasterize = true, colormap = :amp, colorrange = (0, max(var...)))
# Colorbar(gab[1,2], hm_a)
# iEmax, jEmax = idxEm[1], idxEm[2]
# scatter!(ax_a, 1e-3x[iEmax], 1e-3y[jEmax]; marker = :star4, markersize = 10, color = :black)
# hideydecorations!(ax_b, ticks = false)
# varb = interior(SKEₛ, :, :, idxEm[3]);
# x, y, z = nodes(SKEₛ);
# hm_b = heatmap!(ax_b, 1e-3x, 1e-3y, varb; rasterize = true, colormap = :amp, colorrange = (0, max(varb...)))
# Colorbar(gab[1,4], hm_b)
# idxSKEm = argmax(varb)
# scatter!(ax_b, 1e-3x[idxSKEm[1]], 1e-3y[idxSKEm[2]]; marker = :star4, markersize = 10, color = :black)
# colgap!(gab, 1, 1)
# colgap!(gab, 3, 1)
# resize_to_layout!(fig)
# save(filesave * "TKESKE_" * fileparam * "_44h_iter$(iteration)_cg3hm.pdf", fig; pt_per_unit = 1)
# println("Finished plotting TKE SKE fields")
# @info "Memory usage: " * (Sys.free_memory() |> Base.format_bytes) * " GB free"

# fileparam = "sublevels"
# iterations = get_iterations_regex(filehead, fileparam)
# mMLDs = zeros(length(iterations), 3)
# h10ps = zeros(length(iterations), 3)
# h90ps = zeros(length(iterations), 3)
# # MEw3s = zeros(length(iterations), 3)
# for (i,iteration) in enumerate(iterations)
#     # --- Load the snapshot ---
#     output_filename = filehead * "subdomains/" * fileparam * "_snapshot_iter$(iteration).jld2"
#     snapshot = load_subdomain_snapshot(output_filename; variables = ("MLD", "MLD2", "MLD3"));
#     mMLDs[i, 1] = mean(snapshot[:MLD])
#     h10ps[i, 1] = quantile(vec(snapshot[:MLD]), 0.1) # 10th percentile mixed layer depth
#     h90ps[i, 1] = quantile(vec(snapshot[:MLD]), 0.9) # 90th percentile mixed layer depth
#     mMLDs[i, 2] = mean(snapshot[:MLD2])
#     h10ps[i, 2] = quantile(vec(snapshot[:MLD2]), 0.1) # 10th percentile mixed layer depth
#     h90ps[i, 2] = quantile(vec(snapshot[:MLD2]), 0.9) # 90th percentile mixed layer depth
#     mMLDs[i, 3] = mean(snapshot[:MLD3])
#     h10ps[i, 3] = quantile(vec(snapshot[:MLD3]), 0.1) # 10th percentile mixed layer depth
#     h90ps[i, 3] = quantile(vec(snapshot[:MLD3]), 0.9) # 90th percentile mixed layer depth
#     # MEW3s[i, 1] = maximum(snapshot[:Ew])
#     # MEW3s[i, 2] = maximum(snapshot[:Ew2])
#     # MEW3s[i, 3] = maximum(snapshot[:Ew3])
#     @info "Iteration $iteration: mMLDs = $(mMLDs[i, :])"
# end

# fig = Figure(size = (640, 360))
# t = Array(24:2:length(iterations)*2+22)
# t[iterations .> 110000] .+= 24
# ax1 = Axis(fig[1, 1], xlabel = "Time (d)", ylabel = L"h/h_0")
# # ax2 = Axis(fig[2, 1], ylabel = "mMLD2 (m)")
# # ax3 = Axis(fig[2, 1], xlabel = "Time (h)", ylabel = L"\overline{(w^2/2)}^h_{\max}/w_*^2")
# lines!(ax1, t/24, mMLDs[:, 1]/h₀, label = "MLD")
# lines!(ax1, t/24, mMLDs[:, 2]/h₀, label = "MLD2")
# lines!(ax1, t/24, mMLDs[:, 3]/h₀, label = "MLD3")
# fill_between!(ax1, t/24, h10ps[:, 1]/h₀, h90ps[:,1]/h₀; alpha = 0.3)
# fill_between!(ax1, t/24, h10ps[:, 2]/h₀, h90ps[:,2]/h₀; alpha = 0.3)
# fill_between!(ax1, t/24, h10ps[:, 3]/h₀, h90ps[:,3]/h₀; alpha = 0.3)
# # lines!(ax3, t, MEW3s[:, 1]/wₛ^2, label = "Ew1")
# # lines!(ax3, t, MEW3s[:, 2]/wₛ^2, label = "Ew2")
# # lines!(ax3, t, MEW3s[:, 3]/wₛ^2, label = "Ew3")
# # hidexdecorations!(ax1, grid = false)
# # hidexdecorations!(ax2, grid = false)
# resize_to_layout!(fig)
# save(filesave * "mp10p90MLD_" * fileparam * "_1d_nd.pdf", fig; pt_per_unit = 1)

# cutoff = 100
# for (i,iteration) in enumerate(iterations)
#     # --- Load the snapshot ---
#     output_filename = filehead * "subdomains/" * fileparam * "_snapshot_iter$(iteration).jld2"
#     snapshot = load_subdomain_snapshot(output_filename; variables = ("u", "v"), level=224);
#     u̅  = XFaceField(snapshot[:grid], Float32);
#     v̅  = YFaceField(snapshot[:grid], Float32);

#     @time coarse_graining!(snapshot[:u] , u̅ ; kernel=:gaussian, cutoff)
#     @time coarse_graining!(snapshot[:v] , v̅ ; kernel=:gaussian, cutoff, method=:spectral)

#     @time ζ̄i = compute!(Field((∂x(v̅)-∂y(u̅))));
#     @info "Iteration $iteration: max ζ̄i/f = $(max(abs.(interior(ζ̄i))...)/parameters.f)"
# end

# #######################
# # 1. Define the filename of the saved snapshot
# output_filename = filehead * "subdomains/" * fileparam * "_snapshot_iter$(iteration).jld2"

# # 2. Load the snapshot using the new function
# snapshot = load_subdomain_snapshot(output_filename; variables = ("u", "v", "w"),level=221);

# ###################################
# window = nothing #:xhann
# axis_kwargs1 = (xlabel = L"\text{Wavenumber (m^{-1})}", xgridvisible = false,
#                 ylabel = L"E_i(k)/E_{v}(k_{\min})", ygridvisible = false,
#                 xscale = log10, yscale = log10,
#                 limits = ((3e-5, 7e-1), (1e-8,1e0)),
#                 xminorticks = [2e-5:1e-5:9e-5; 2e-4:1e-4:9e-4; 2e-3:1e-3:9e-3; 2e-2:1e-2:9e-2; 2e-1:1e-1:9e-1],xminorticksvisible = true,
#                 yminorticks = [2e-8:1e-8:9e-8;2e-7:1e-7:9e-7;2e-6:1e-6:9e-6; 2e-5:1e-5:9e-5; 2e-4:1e-4:9e-4; 2e-3:1e-3:9e-3; 2e-2:1e-2:9e-2; 2e-1:1e-1:9e-1; 2e0:1e0:9e0; 2e1:1e1:9e1],yminorticksvisible = true);

# fig = Figure(size = (640, 270))
# g3 = fig[1, 1] = GridLayout()
# alphabet = [letter for letter in 'a':'z'];
# klev = 1
# L_filter = 300
# ax = Axis(g3[1, 1]; title=L"\text{Nonhydrostatic LES spectra at 44~h}~z=-4~\text{m}", axis_kwargs1...)
# vlines!(ax, 1/200; color = :black, linewidth = 0.8)
# #vlines!(ax, 2*3.14/300; color = :red, linewidth = 0.8)
# vlines!(ax, 2*3.14/400; color = :red, linewidth = 0.8)

# xidx = 1:20480#5120:15360
# @time Su = isotropic_powerspectrum(interior(snapshot[:u], xidx, :, klev), interior(snapshot[:u], xidx, :, klev);window);
# @time Sv = isotropic_powerspectrum(interior(snapshot[:v], xidx, :, klev), interior(snapshot[:v], xidx, :, klev);window,L_filter)
# wk = (interior(snapshot[:w], xidx, :, klev));#+interior(snapshot[:w], :, :, klev+1))/2
# Sw = isotropic_powerspectrum(wk, wk;window)
# #St = isotropic_powerspectrum(interior(snapshot[:T], :, :, klev), interior(snapshot[:T], :, :, klev), xT, yT;window)
# idx = Su.freq .> 0
# lines!(ax, Su.freq[idx], 0.1*(Su.freq[idx]./Su.freq[idx][1]).^-2, linestyle = :dash, color = :black)
# #text!(ax, 10^-3, 0.7; text = L"k^{-2}")
# lines!(ax, Su.freq[idx], 5e2*(Su.freq[idx]./Su.freq[idx][1]).^(-5/3), linestyle = :dash, color = :gray)
# #text!(ax, 10^-3, 10^-5.5; text = L"k^{-5/3}")
# #lines!(ax, St.freq[idx], Real.(St.spec[idx]./St.spec[idx][1]), color = :red, label = L"E_T")
# lines!(ax, Su.freq[idx], Real.(Su.spec[idx]./Sv.spec[idx][1]), color = :blue, label = L"E_u")
# lines!(ax, Sv.freq[idx], Real.(Sv.spec[idx]./Sv.spec[idx][1]), color = :green, label = L"E_v")
# lines!(ax, Sv.freq[idx], Real.(Sv.specf[idx]./Sv.spec[idx][1]), color = :green, linestyle = :dash, label = L"E_{\overline{v}}")
# lines!(ax, Sw.freq[idx], Real.(Sw.spec[idx]./Sv.spec[idx][1]), color = :black, label = L"E_w")
# axislegend(ax, labelsize=9, patchsize = (15, 1), framevisible = false,
#             padding = (0f0, 0f0, 0f0, 0f0), patchlabelgap = 3, rowgap = 1)
# resize_to_layout!(fig)
# save(filesave * "spectra_" * fileparam * "_44h4m_iter$(iteration).pdf", fig; pt_per_unit = 1)

#############################
# uᵃ, vᵃ, wᵃ, Bᵃ = 0.0, 0.0, 0.0, 0.0
# nsubdomains = 16
# for i = 1:nsubdomains
#     fileparam = "subdomain" * string(i)
#     @info "Computing and saving data for $fileparam"
#     # 1. Define the filename of the saved snapshot
#     output_filename = filehead * "subdomains/" * fileparam * "_snapshot_iter$(iteration).jld2"

#     # 2. Load the snapshot using the new function
#     snapshot = load_subdomain_snapshot(output_filename; variables = ("u", "v", "w", "T"));

#     x0, y0, z0 = nodes(snapshot[:T])
#     yc = 0.5 * (y0[640] + y0[641])
#     to_grid = RectilinearGrid(snapshot[:grid].architecture,Float32;
#                               size = (2048*2, 512*2, length(z0)),
#                               x = (-10000,10000),
#                               y = (yc-parameters.Δh*512,yc+parameters.Δh*512),
#                               z = (-81,0),
#                               topology = (Bounded, Bounded, Bounded))

#     uᵃi, vᵃi, wᵃi, Bᵃi = along_front_averages(snapshot; to_grid, cutoff=300, border=:reflect, Lx = snapshot[:grid].Lx, Ly = snapshot[:grid].Ly);
#     global uᵃ = uᵃ .+ interior(uᵃi)
#     global vᵃ = vᵃ .+ interior(vᵃi)
#     global wᵃ = wᵃ .+ interior(wᵃi)
#     global Bᵃ = Bᵃ .+ interior(Bᵃi)
# end
# uᵃ ./= nsubdomains
# vᵃ ./= nsubdomains
# wᵃ ./= nsubdomains
# Bᵃ ./= nsubdomains

# jldopen(filehead * "subdomains/Vavgs_iter$(iteration)_afront.jld2", "w") do file
#    file["fields/ba"] = Bᵃ
#    file["fields/ua"] = uᵃ
#    file["fields/va"] = vᵃ
#    file["fields/wa"] = wᵃ
# end

file = jldopen(filehead * "subdomains/Vavgs_iter$(iteration)_afront.jld2", "r");
Bᵃ = file["fields/ba"]
uᵃ = file["fields/ua"]
vᵃ = file["fields/va"]
wᵃ = file["fields/wa"]
close(file)

# fig = Figure(size = (640, 400))
# gab = fig[1, 1] = GridLayout()
# fileparam = "subdomain1"
# slices = jldopen(filehead * "subdomains/Vslices_"*fileparam*"_iter$(iteration)_afront.jld2", "r")
# xi, zi = slices["metadata/x"], slices["metadata/z"]
# close(slices)

# ax1 = Axis(gab[1,1]; titlealign = :left, title="(a)", titlefont=texfont(), ylabel=L"z~\text{(m)}",limits=(nothing,(-80,0)))
# ax2 = Axis(gab[1,3]; titlealign = :left, title="(b)", titlefont=texfont(), limits=(nothing,(-80,0)))
# ax3 = Axis(gab[2,1]; titlealign = :left, title="(c)", titlefont=texfont(), ylabel=L"z~\text{(m)}",limits=(nothing,(-80,0)))
# ax4 = Axis(gab[2,3]; titlealign = :left, title="(d)", titlefont=texfont(), limits=(nothing,(-80,0)))
# hm1 = heatmap!(ax1, 1e-3xi, zi, Bᵃ[:,1,:]/α/g; rasterize = true, colormap = :thermal)
# hm2 = heatmap!(ax2, 1e-3xi, zi, vᵃ[:,1,:]; rasterize = true, colorrange = (-0.12,0.12), colormap = :balance)
# hm3 = heatmap!(ax3, 1e-3xi, zi, uᵃ[:,1,:]; rasterize = true, colorrange = (-0.08,0.08), colormap = :balance)
# hm4 = heatmap!(ax4, 1e-3xi, zi, wᵃ[:,1,:]; rasterize = true, colorrange = (-4e-4,4e-4), colormap = :delta)

# levels = Makie.get_tickvalues(Makie.LinearTicks(20), extrema(Bᵃ[:,1,:])...)
# contour!(ax1, 1e-3xi, zi, Bᵃ[:,1,:]; color=:black)
# contour!(ax2, 1e-3xi, zi, Bᵃ[:,1,:]; color=:black) 
# contour!(ax3, 1e-3xi, zi, Bᵃ[:,1,:]; color=:black) 
# contour!(ax4, 1e-3xi, zi, Bᵃ[:,1,:]; color=:black) 
# hideydecorations!(ax2, ticks = false)
# hideydecorations!(ax4, ticks = false)
# hidexdecorations!(ax1, ticks = false)
# hidexdecorations!(ax2, ticks = false)
# ax3.xlabel = L"x~\text{(km)}"
# ax4.xlabel = L"x~\text{(km)}"
# Colorbar(gab[1,2], hm1)
# Colorbar(gab[1,4], hm2)
# Colorbar(gab[2,2], hm3)
# Colorbar(gab[2,4], hm4)
# Label(gab[1, 1, Top()], L"T^a~({^\circ}C)", valign = :bottom,font = texfont(),padding = (0, 0, 5, 0))
# Label(gab[1, 3, Top()], L"v^a~(\text{m s^{-1}})", valign = :bottom,font = texfont(),padding = (0, 0, 5, 0))
# Label(gab[2, 1, Top()], L"u^a~(\text{m s^{-1}})", valign = :bottom,font = texfont(),padding = (0, 0, 5, 0))
# Label(gab[2, 3, Top()], L"w^a~(\text{m s^{-1}})", valign = :bottom,font = texfont(),padding = (0, 0, 5, 0))
# rowgap!(gab, 3)
# colgap!(gab, 1, 0)
# colgap!(gab, 2, 5)
# colgap!(gab, 3, 0)
# resize_to_layout!(fig)
# save(filesave * "afrontfields_vslices_30h_iter$(iteration).pdf", fig; pt_per_unit = 1)

# varᵃ = (uᵃ, vᵃ, wᵃ, Bᵃ)
# for i = 4:4
#     fileparam = "subdomain" * string(i)
#     @info "Computing and saving data for $fileparam"
#     # 1. Define the filename of the saved snapshot
#     output_filename = filehead * "subdomains/" * fileparam * "_snapshot_iter$(iteration).jld2"

#     # 2. Load the snapshot using the new function
#     snapshot = load_subdomain_snapshot(output_filename; variables = ("u", "v", "w", "T", "MLD3"));

#     x0, y0, z0 = nodes(snapshot[:T])
#     yc = 0.5 * (y0[640] + y0[641])
#     to_grid = RectilinearGrid(snapshot[:grid].architecture,Float32;
#                               size = (2048*2, 512*2, length(z0)),
#                               x = (-10000,10000),
#                               y = (yc-parameters.Δh*512,yc+parameters.Δh*512),
#                               z = (-81,0),
#                               topology = (Bounded, Bounded, Bounded))

#     u̅, v̅, w̅, B̅, uˢ, vˢ, wˢ, Bˢ, τuu, τvv, τww, τwb, Πₕ, Πᵥ, Πᵥg, Pᵃ, Pᵃᵥg, Pˢ, Pˢᵥg, Pᵀ, wˢbˢ = coarse_grained_fluxes(snapshot, Ub, Vb, varᵃ; to_grid, cutoff=300, border=:reflect, Lx = snapshot[:grid].Lx, Ly = snapshot[:grid].Ly);

#     xi, _, _ = nodes(CenterField(to_grid,Float32))
#     itp2nodes(field) = Array([interpolate((x, yc, z), field) for x in xi, z in z0'])
#     jldopen(filehead * "subdomains/Vslices_"*fileparam*"_iter$(iteration)_afront.jld2", "w") do file
#         file["fields/bcg"] = B̅[1]
#         file["fields/ucg"] = u̅[1]
#         file["fields/vcg"] = v̅[1]
#         file["fields/wcg"] = w̅[1]
#         file["fields/bcga"] = interior(B̅[2])
#         file["fields/ucga"] = interior(u̅[2])
#         file["fields/vcga"] = interior(v̅[2])
#         file["fields/wcga"] = interior(w̅[2])
#         file["fields/us"]  = itp2nodes(uˢ)
#         file["fields/vs"]  = itp2nodes(vˢ)
#         file["fields/ws"]  = itp2nodes(wˢ)
#         file["fields/bs"]  = itp2nodes(Bˢ)
#         println(mean(uˢ, dims=2))
#         file["fields/usa"]  = interior(mean(uˢ, dims=2))
#         file["fields/vsa"]  = interior(mean(vˢ, dims=2))
#         file["fields/wsa"]  = interior(mean(wˢ, dims=2))
#         file["fields/bsa"]  = interior(mean(Bˢ, dims=2))
#         file["fields/τwb"] = τwb[1]
#         file["fields/τwba"] = interior(τwb[2])
#         file["fields/wbs"] = itp2nodes(wˢbˢ)
#         file["fields/wbsa"] = interior(mean(wˢbˢ, dims=2))
#         file["fields/τuu"] = τuu[1]
#         file["fields/τww"] = τww[1]
#         file["fields/τvv"] = τvv[1]
#         file["fields/τuua"] = interior(τuu[2])
#         file["fields/τwwa"] = interior(τww[2])
#         file["fields/τvva"] = interior(τvv[2])
#         file["fields/Ph"] = itp2nodes(Πₕ)#, itp2nodes(mean(Πₕ, dims=2))] #./ (parameters.f * wₛ^2)
#         file["fields/Pv"] = itp2nodes(Πᵥ)#, itp2nodes(mean(Πᵥ, dims=2))] #./ (parameters.f * wₛ^2)
#         file["fields/Pvg"] = itp2nodes(Πᵥg)
#         file["fields/Pas"] = itp2nodes(Pᵃ)#, itp2nodes(mean(Pᵃ, dims=2))] #./ (parameters.f * wₛ^2)
#         file["fields/Pasg"] = itp2nodes(Pᵃᵥg)
#         file["fields/Pss"] = itp2nodes(Pˢ)#, itp2nodes(mean(Pˢ, dims=2))] #./ (parameters.f * wₛ^2)
#         file["fields/Pssg"] = itp2nodes(Pˢᵥg)
#         file["fields/PTs"] = itp2nodes(Pᵀ)#, itp2nodes(mean(Pᵀ, dims=2))] #./ (parameters.f * wₛ^2)
#         file["fields/Pha"] = interior(mean(Πₕ, dims=2))
#         file["fields/Pva"] = interior(mean(Πᵥ, dims=2))
#         file["fields/Pasa"] = interior(mean(Pᵃ, dims=2))
#         file["fields/Pssa"] = interior(mean(Pˢ, dims=2))
#         file["fields/PTsa"] = interior(mean(Pᵀ, dims=2))
        
#         file["metadata/iteration"] = iteration
#         file["metadata/x"] = xi
#         file["metadata/y"] = yc
#         file["metadata/z"] = z0
#     end
# end
# println("Finished plotting TKE, SKE, and P fields")

# 3. Plot the TKE fields
# alphabet = [letter for letter in 'a':'z'];
# cmap = :balance
# ylabels = [L"y=18.75~\text{km}",L"y=43.75~\text{km}",L"y=68.75~\text{km}",L"y=93.75~\text{km}"]
# iarrow = [320-32, 320+32]
# idxz = findfirst(zu .> -80) : 10 : length(zu)
# skwargs = (; markersize = 15, strokewidth = 1, color = :transparent, strokecolor = :black)
# skwargs1 = (; markersize=15, marker=:xcross, color=:black)
# skwargs3 = (; markersize=3, color=:black)
# xmarkers = vec([-5*ones(4) 5*ones(4)])
# zmarkers = vec(repeat(-70.:20.:-10., 2, 1))
# fig = Figure(size = (640, 680))
# gab = fig[1, 1] = GridLayout()
# for i = 1:4
#     fileparam = "subdomain$(5-i)"
#     @info "Ploting vslices for $fileparam"
#     # 1. Define the filename of the saved snapshot
#     output_filename = filehead * "subdomains/" * fileparam * "_snapshot_iter$(iteration).jld2"

#     # 2. Load the snapshot using the new function
#     snapshot = load_subdomain_snapshot(output_filename; variables = ("MLD", "MLD3"));
#     slices = jldopen(filehead * "subdomains/Vslices_"*fileparam*"_iter$(iteration)_afront.jld2", "r")
#     xi, yi, zi = slices["metadata/x"], slices["metadata/y"], slices["metadata/z"]
#     # TKEₛi = (slices["fields/τvv"] + slices["fields/τuu"] + slices["fields/τww"])/2/wₛ^2;
#     # SKEₛi = (slices["fields/us"].^2 + slices["fields/vs"].^2 + slices["fields/ws"].^2)/2/wₛ^2;
#     var1i = slices["fields/Ph"]/(parameters.f * wₛ^2)
#     var2i = slices["fields/Pv"]/(parameters.f * wₛ^2)
#     bcgi = slices["fields/bcg"][:,1,:]
#     close(slices)

#     title1 = "("*alphabet[2i-1]*")" 
#     title2 = "("*alphabet[2i]*")" 
#     ax1 = Axis(gab[i,1]; titlealign = :left, title=title1, titlefont=texfont(), ylabel=L"z~\text{(m)}",limits=(nothing,(-80,0)))
#     ax2 = Axis(gab[i,3]; titlealign = :left, title=title2, titlefont=texfont(), limits=(nothing,(-80,0)))
#     lim1 = maximum(abs.(var1i[1025:3096,1,:]))
#     lim2 = maximum(abs.(var2i[1025:3096,1,:]))
#     hm1 = heatmap!(ax1, 1e-3xi, zi, var1i[:,1,:]; rasterize = true, colormap = cmap, colorrange = (-lim1, lim1))
#     hm2 = heatmap!(ax2, 1e-3xi, zi, var2i[:,1,:]; rasterize = true, colormap = cmap, colorrange = (-lim2, lim2))

#     if (5-i) % 2 == 0
#         arrows!(ax1,xu[iarrow]/1e3.-50, zu[idxz], shift(interior(Ub,iarrow,160(5-i)-40,idxz)), 0*shift(interior(Ub,iarrow,160i,idxz)), arrowsize = 6, lengthscale = 1e2,linecolor = :black, arrowcolor = :black, linewidth = 0.8)
#         arrows!(ax2,xu[iarrow]/1e3.-50, zu[idxz], shift(interior(Ub,iarrow,160(5-i)-40,idxz)), 0*shift(interior(Ub,iarrow,160i,idxz)), arrowsize = 6, lengthscale = 1e2,linecolor = :black, arrowcolor = :black, linewidth = 0.8) 
#     else
#         scatter!(ax1,xmarkers, zmarkers; skwargs...)
#         scatter!(ax2,xmarkers, zmarkers; skwargs...)
#     end
#     if (5-i) == 1
#         scatter!(ax1,xmarkers, zmarkers; skwargs1...)
#         scatter!(ax2,xmarkers, zmarkers; skwargs1...)
#     elseif (5-i) == 3
#         scatter!(ax1,xmarkers, zmarkers; skwargs3...)
#         scatter!(ax2,xmarkers, zmarkers; skwargs3...)
#     end
#     # hm2 = heatmap!(ax2, 1e-3x, z, var2i[:,1,:]; rasterize = true, colormap = cmap, colorrange = (-50,50))
#     # hm1 = heatmap!(ax1, 1e-3x, z, var1i[:,1,:]; rasterize = true, colormap = cmap, colorrange = (-50,50))
#     lines!(ax1, 1e-3xi, -interior(snapshot[:MLD3],513:4608,640,1), color = :blue, linewidth = 1)
#     lines!(ax2, 1e-3xi, -interior(snapshot[:MLD3],513:4608,640,1), color = :blue, linewidth = 1)
#     # We can use Makie's tick finders to get some nice looking contour levels:
#     levels = Makie.get_tickvalues(Makie.LinearTicks(20), extrema(bcgi)...)
#     contour!(ax1, 1e-3xi, zi, bcgi; color = :black)
#     contour!(ax2, 1e-3xi, zi, bcgi; color = :black) 
#     contour!(ax1, 1e-3xi, zi, Bᵃ[:,1,:]; color=:black, linestyle=:dash)
#     contour!(ax2, 1e-3xi, zi, Bᵃ[:,1,:]; color=:black, linestyle=:dash) 
#     hideydecorations!(ax2, ticks = false)
#     if i<4
#         hidexdecorations!(ax1, ticks = false)
#         hidexdecorations!(ax2, ticks = false)
#     else
#         ax1.xlabel = L"x~\text{(km)}"
#         ax2.xlabel = L"x~\text{(km)}"
#     end
#     Colorbar(gab[i,2], hm1)
#     Colorbar(gab[i,4], hm2)
#     Label(gab[i, 1, Top()], ylabels[5-i], valign = :bottom,font = texfont(),padding = (0, 0, 5, 0))
#     Label(gab[i, 3, Top()], ylabels[5-i], valign = :bottom,font = texfont(),padding = (0, 0, 5, 0))
#     Label(gab[i, 2, Top()], L"P_H/(fw_*^2)", valign = :bottom, font = texfont(), padding = (0, 0, 5, 0))
#     Label(gab[i, 4, Top()], L"P_V/(fw_*^2)", valign = :bottom, font = texfont(), padding = (0, 0, 5, 0))
# end
# rowgap!(gab, 3)
# colgap!(gab, 1, 0)
# colgap!(gab, 2, 10)
# colgap!(gab, 3, 0)
# resize_to_layout!(fig)
# save(filesave * "PhPv_vslices_30h_iter$(iteration)_afront.pdf", fig; pt_per_unit = 1)

# cmap = :balance
# fig = Figure(size = (640, 540))
# gab = fig[1, 1] = GridLayout()
# for i = 1:4
#     fileparam = "subdomain$(5-i)"
#     @info "Ploting vslices for $fileparam"
#     # 1. Define the filename of the saved snapshot
#     output_filename = filehead * "subdomains/" * fileparam * "_snapshot_iter$(iteration).jld2"

#     # 2. Load the snapshot using the new function
#     snapshot = load_subdomain_snapshot(output_filename; variables = ("MLD", "MLD3"));
#     slices = jldopen(filehead * "subdomains/Vslices_"*fileparam*"_iter$(iteration)_afront.jld2", "r")
#     xi, yi, zi = slices["metadata/x"], slices["metadata/y"], slices["metadata/z"]
#     @info "yi = $(1e-3y)km"
#     # TKEₛi = (slices["fields/τuv"] + slices["fields/τuu"] + slices["fields/τww"])/2/wₛ^2;
#     # SKEₛi = (slices["fields/us"].^2 + slices["fields/vs"].^2 + slices["fields/ws"].^2)/2/wₛ^2;
#     var1i = slices["fields/Ph"]/(parameters.f * wₛ^2)
#     var2i = slices["fields/Pv"]/(parameters.f * wₛ^2)
#     var3i = slices["fields/Pvg"]/(parameters.f * wₛ^2)
#     bcgi = slices["fields/bcg"][:,1,:]
#     close(slices)

#     title1 = "("*alphabet[3i-2]*")" 
#     title2 = "("*alphabet[3i-1]*")" 
#     title3 = "("*alphabet[3i]*")" 
#     ax1 = Axis(gab[i,1]; titlealign = :left, title=title1, titlefont=texfont(), ylabel=L"z~\text{(m)}",limits=(nothing,(-80,0)))
#     ax2 = Axis(gab[i,3]; titlealign = :left, title=title2, titlefont=texfont(), limits=(nothing,(-80,0)))
#     ax3 = Axis(gab[i,5]; titlealign = :left, title=title3, titlefont=texfont(), limits=(nothing,(-80,0)))
#     lim1 = maximum(abs.(var1i[1025:3096,1,:]))
#     lim2 = maximum(abs.(var2i[1025:3096,1,:]))
#     lim3 = maximum(abs.(var3i[1025:3096,1,:]))
#     hm1 = heatmap!(ax1, 1e-3xi, zi, var1i[:,1,:]; rasterize = true, colormap = cmap, colorrange = (-lim1,lim1))
#     hm2 = heatmap!(ax2, 1e-3xi, zi, var2i[:,1,:]; rasterize = true, colormap = cmap, colorrange = (-lim2,lim2))
#     hm3 = heatmap!(ax3, 1e-3xi, zi, var3i[:,1,:]; rasterize = true, colormap = cmap, colorrange = (-lim3,lim3))
#     lines!(ax1, 1e-3xi, -interior(snapshot[:MLD3],513:4608,640,1), color = :blue, linewidth = 1)
#     lines!(ax2, 1e-3xi, -interior(snapshot[:MLD3],513:4608,640,1), color = :blue, linewidth = 1)
#     lines!(ax3, 1e-3xi, -interior(snapshot[:MLD3],513:4608,640,1), color = :blue, linewidth = 1)
#     if (5-i) % 2 == 0
#         arrows!(ax1,xu[iarrow]/1e3.-50, zu[idxz], shift(interior(Ub,iarrow,160(5-i)-40,idxz)), 0*shift(interior(Ub,iarrow,160i,idxz)), arrowsize = 6, lengthscale = 1e2,linecolor = :black, arrowcolor = :black, linewidth = 0.8)
#         arrows!(ax2,xu[iarrow]/1e3.-50, zu[idxz], shift(interior(Ub,iarrow,160(5-i)-40,idxz)), 0*shift(interior(Ub,iarrow,160i,idxz)), arrowsize = 6, lengthscale = 1e2,linecolor = :black, arrowcolor = :black, linewidth = 0.8) 
#         arrows!(ax3,xu[iarrow]/1e3.-50, zu[idxz], shift(interior(Ub,iarrow,160(5-i)-40,idxz)), 0*shift(interior(Ub,iarrow,160i,idxz)), arrowsize = 6, lengthscale = 1e2,linecolor = :black, arrowcolor = :black, linewidth = 0.8) 
#     else
#         scatter!(ax1,xmarkers, zmarkers; skwargs...)
#         scatter!(ax2,xmarkers, zmarkers; skwargs...)
#         scatter!(ax3,xmarkers, zmarkers; skwargs...)
#     end
#     if (5-i) == 1
#         scatter!(ax1,xmarkers, zmarkers; skwargs1...)
#         scatter!(ax2,xmarkers, zmarkers; skwargs1...)
#         scatter!(ax3,xmarkers, zmarkers; skwargs1...)
#     elseif (5-i) == 3
#         scatter!(ax1,xmarkers, zmarkers; skwargs3...)
#         scatter!(ax2,xmarkers, zmarkers; skwargs3...)
#         scatter!(ax3,xmarkers, zmarkers; skwargs3...)
#     end
#     # We can use Makie's tick finders to get some nice looking contour levels:
#     levels = Makie.get_tickvalues(Makie.LinearTicks(20), extrema(bcgi)...)
#     contour!(ax1, 1e-3xi, zi, bcgi; color = :black)
#     contour!(ax2, 1e-3xi, zi, bcgi; color = :black) 
#     contour!(ax3, 1e-3xi, zi, bcgi; color = :black) 
#     contour!(ax1, 1e-3xi, zi, Bᵃ[:,1,:]; color=:black, linestyle=:dash)
#     contour!(ax2, 1e-3xi, zi, Bᵃ[:,1,:]; color=:black, linestyle=:dash) 
#     contour!(ax3, 1e-3xi, zi, Bᵃ[:,1,:]; color=:black, linestyle=:dash)
#     hideydecorations!(ax2, ticks = false)
#     hideydecorations!(ax3, ticks = false)
#     if i<4
#         hidexdecorations!(ax1, ticks = false)
#         hidexdecorations!(ax2, ticks = false)
#         hidexdecorations!(ax3, ticks = false)
#     else
#         ax1.xlabel = L"x~\text{(km)}"
#         ax2.xlabel = L"x~\text{(km)}"
#         ax3.xlabel = L"x~\text{(km)}"
#     end
#     Colorbar(gab[i,2], hm1)
#     Colorbar(gab[i,4], hm2)
#     Colorbar(gab[i,6], hm3)
#     Label(gab[i, 1, Top()], ylabels[5-i], valign = :bottom,font = texfont(),padding = (0, 0, 5, 0))
#     Label(gab[i, 3, Top()], ylabels[5-i], valign = :bottom,font = texfont(),padding = (0, 0, 5, 0))
#     Label(gab[i, 5, Top()], ylabels[5-i], valign = :bottom,font = texfont(),padding = (0, 0, 5, 0))
#     Label(gab[i, 2, Top()], L"P_H/(fw_*^2)", valign = :bottom, font = texfont(), padding = (0, 0, 5, 0))
#     Label(gab[i, 4, Top()], L"P_V/(fw_*^2)", valign = :bottom, font = texfont(), padding = (0, 0, 5, 0))
#     Label(gab[i, 6, Top()], L"P_{Vg}/(fw_*^2)", valign = :bottom, font = texfont(), padding = (0, 0, 5, 0))
# end
# rowgap!(gab, 3)
# colgap!(gab, 1, 0)
# colgap!(gab, 2, 8)
# colgap!(gab, 3, 0)
# colgap!(gab, 4, 8)
# colgap!(gab, 5, 0)
# resize_to_layout!(fig)
# save(filesave * "PhPvPvg_vslices_30h_iter$(iteration)_afront.pdf", fig; pt_per_unit = 1)

    # TKEₛj = file["fields/TKEs"];
    # Πₕj = file["fields/PHs"];
    # Πᵥj = file["fields/PVs"];
    # τwbj = file["fields/Bs"];
    # x = file["metadata/x"];
    # z = file["metadata/z"];
    # nxj,nzj = length(x), length(z)
    # titles = [L"\text{(a) TKE}/w_*^2",L"\text{(b) }P_H/(fw_*^2)",L"\text{(c) }P_V/(fw_*^2)",L"\text{(d) }B/(fw_*^2)"]
    # for (i, var0) in enumerate([TKEₛj, Πₕj, Πᵥj, τwbj])
    #     cmap = i==1 ? :amp : :balance
    #     #x, y, z = nodes(var0);
    #     var = var0#interior(var0, :, 640,:) * (i==1 ? 1 : 1e4/wₛ^2);
    #     nx, nz = size(var)
    #     crange = i==1 ? (0, max(var...)) : (-50,50)
    #     ax_a = Axis(gab[i,1]; titlealign = :left, title=titles[i], xlabel=L"x~\text{(km)}", ylabel=L"z~\text{(m)}",limits=(nothing,(-80,0)))
    #     ax_b = Axis(gab[i,3]; titlealign = :left, limits=(nothing,(-80,0)))
    #     hm_a = heatmap!(ax_a, 1e-3x, z, var[1:nxj, nz-nzj+1:nz]; rasterize = true, colormap = cmap, colorrange = crange)
    #     Colorbar(gab[i,2], hm_a)
    #     idxEm = argmax(var[641:end-640,:])
    #     iEmax, jEmax = idxEm[1]+640, idxEm[2]
    #     @info "max at $(x[iEmax]/1e3), $(z[jEmax]): $(var[iEmax,jEmax])"
    #     scatter!(ax_a, 1e-3x[iEmax], z[jEmax]; marker = :star4, markersize = 10, color = :black)
    #     idxs = [argmin(vec(mean(var;dims=2))[641:end-640])+640,argmax(vec(mean(var;dims=2))[641:end-640])+640,iEmax]
    #     vlines!(ax_a, 1e-3x[idxs], color = Makie.wong_colors()[2:4], linewidth = 0.8)
    #     lines!(ax_a, 1e-3x, -interior(snapshot[:MLD3],:,640,1), color = :blue, linewidth = 1, label="MLD")
    #     hideydecorations!(ax_b, ticks = false)
    #     lines!(ax_b, vec(mean(var, dims =(1))), z; linewidth = 2, label="mean")
    #     lines!(ax_b, var[idxs[1],nz-nzj+1:nz], z; linewidth = 1)
    #     lines!(ax_b, var[idxs[2],nz-nzj+1:nz], z; linewidth = 1)
    #     lines!(ax_b, var[idxs[3],nz-nzj+1:nz], z; linewidth = 1)
    #     # hlines!(ax_b, -interior(snapshot[:BLD], :, 1497, 1)[[2004,1004,3004]], linestyle = :dash, color = [:orange, :green, :purple], linewidth = 0.8)
    #     colsize!(gab, 3, Relative(0.3))
    #     colgap!(gab, 1, 1)
    #     colgap!(gab, 2, 5)
    #     resize_to_layout!(fig)
    #     if i==1
    #         axislegend(ax_a, labelsize=10,position = :lb,patchsize = (15, 5), patchlabelgap = 3, rowgap = 1)
    #         axislegend(ax_b,  labelsize=10, position = :rb,patchsize = (15, 3), patchlabelgap = 3)
    #     end
    #     if i<4
    #         hidexdecorations!(ax_a, ticks = false)
    #     else
    #         ax_b.xlabel = xlabel=L"\text{mean vs. local}"
    #     end
    # end
    # rowgap!(gab, 1)
    # save(filesave * "TKE3Ps_" * fileparam * "_44h_iter$(iteration)_cg3hm.pdf", fig; pt_per_unit = 1)

# file = jldopen(filehead * "TKEphysical_"*fileparam*"_iter$(iteration).jld2")
# TKEw = file["fields/τww"]/2
# zmin, zmax = file["metadata/zlims"]
# Nx, Nz = size(TKEw)
# Δz = (zmax-zmin)/Nz
# x = range(-12500, 12500, length = 640*8)
# z = zmin+Δz/2:Δz:zmax-Δz/2
# BLDw = zeros(Nx);
# r = 0.3
# for i = 1:Nx
#     imax = argmax(TKEw[i,:])
#     maxi = TKEw[i,imax]
#     iBLD = findlast(TKEw[i,1:imax-1] .< r*maxi)
#     BLDw[i] = -(z[iBLD] + (z[iBLD+1]-z[iBLD])*(r*maxi-TKEw[i,iBLD])/(TKEw[i,iBLD+1]-TKEw[i,iBLD]))
#     @info "x = $(1e-3x[i])km, TKEw = $maxi, BLD = $(BLDw[i])m"
# end

# TKEs = (interior(τuu, 1:5120, 640, :) .+ interior(τuu,2:5121, 640, :))/4;
# TKEs .+= (interior(τvv, :, 640, :) .+ interior(τvv, :, 641, :))/4;
# TKEs .+= (interior(τww, :, 640, 1:224) .+ interior(τww, :, 640, 2:225))/4;

#################################
#iterations = [32207,52543,72635]#72635 #52543 #32207 
# limits=(nothing,(-80,0))
# f = 1e-4
# fig = Figure(size = (640, 320))
# gab = fig[1, 1] = GridLayout()
# ax_a = Axis(gab[1,1]; titlealign = :left, title=L"\text{(a)}", xlabel=L"\langle u \rangle~\text{(m~s^{-1})}", ylabel=L"z~\text{(m)}",limits=limits)
# ax_b = Axis(gab[1,2]; titlealign = :left, title=L"\text{(d)}", xlabel=L"\langle v \rangle~\text{(m~s^{-1})}", limits=limits)
# ax_c = Axis(gab[1,3]; titlealign = :left, title=L"\text{(c)}", xlabel=L"10^6\langle\text{KE}_w\rangle~\text{(m^2~s^{-2})}",limits=limits)
# ax_d = Axis(gab[1,4]; titlealign = :left, title=L"\text{(d)}", xlabel=L"\text{max}_{(x,y)} \zeta/f", limits=limits)
# colgap!(gab, 3)
# hideydecorations!(ax_b, ticks = false)
# hideydecorations!(ax_c, ticks = false)
# hideydecorations!(ax_d, ticks = false)
# for iteration in iterations
#     # 1. Define the filename of the saved snapshot
#     output_filename = filehead * "subdomains/" * fileparam * "_snapshot_iter$(iteration).jld2"

#     # 2. Load the snapshot using the new function
#     snapshot = load_subdomain_snapshot(output_filename)
#     _, _, zu = nodes(snapshot[:u]);
#     _, _, zw = nodes(snapshot[:w]);
#     lines!(ax_a, vec(mean(compute!(Field(snapshot[:u])), dims =(1,2))), zu; linewidth = 1)
#     lines!(ax_b, vec(mean(compute!(Field(snapshot[:v])), dims =(1,2))), zu; linewidth = 1)
#     lines!(ax_c, 1e6*vec(mean(compute!(Field(snapshot[:w]^2/2)), dims =(1,2))), zw; linewidth = 1)
#     Ro = compute!(Field((∂x(snapshot[:v])-∂y(snapshot[:u]))/f));
#     x, y, z = nodes(Ro);
#     nx, ny = length(x), length(y)
#     lines!(ax_d, vec(maximum(interior(Ro,2:nx-1,2:ny-1,:), dims =(1,2))), z; linewidth = 1)
#     resize_to_layout!(fig)
#     save(filesave * "uvEwRo_" * fileparam * ".pdf", fig; pt_per_unit = 1)
# end
# println("Finished plotting uvEw fields")

#################################
# iteration = 72635
# # 1. Define the filename of the saved snapshot
# output_filename = filehead * "subdomains/" * fileparam * "_snapshot_iter$(iteration).jld2"

# # 2. Load the snapshot using the new function
# snapshot = load_subdomain_snapshot(output_filename)
# fig = Figure(size = (640, 320))
# gab = fig[1, 1] = GridLayout()
# x, y, z = nodes(snapshot[:w]);
# ax_a = Axis(gab[1,1]; titlealign = :left, title=L"\text{(a)}~w~\text{(mm s^{-1})}", xlabel=L"x~\text{(km)}", ylabel=L"z~\text{(m)}",limits=(nothing,(-90,0)))
# ax_b = Axis(gab[1,3]; titlealign = :left, title=L"\text{(b)}", xlabel=L"10^6\langle\text{KE}_w\rangle~\text{(m^2~s^{-2})}",limits=(nothing,(-90,0)))
# hm_a = heatmap!(ax_a, 1e-3x, z, 1e3interior(snapshot[:w], :, 1497, :); rasterize = true, colormap = :delta, colorrange = (-20, 20))
# Colorbar(gab[1,2], hm_a)
# lines!(ax_a, 1e-3x, -interior(snapshot[:BLD], :, 1497, 1), color = :black, linewidth = 0.5, alpha=0.8)
# #lines!(ax_a, 1e-3x, -interior(snapshot[:MLD], :, 1497, 1), color = :red, linewidth = 1)
# #lines!(ax_a, 1e-3x, -interior(snapshot[:MLD2], :, 1497, 1), color = :blue, linewidth = 1)
# lines!(ax_a, 1e-3x, -interior(snapshot[:MLD3], :, 1497, 1), color = :green, linewidth = 1)
# vlines!(ax_a, 1e-3x[[2004,1004,3004]], color = [:orange, :green, :purple], linewidth = 0.8)
# hideydecorations!(ax_b, ticks = false)
# lines!(ax_b, 1e6*vec(mean(interior(snapshot[:w], :, 1497, :).^2/2, dims =(1))), z; linewidth = 1)
# lines!(ax_b, 1e6*interior(snapshot[:w], 2004, 1497, :).^2/2, z; linewidth = 1)
# lines!(ax_b, 1e6*interior(snapshot[:w], 1004, 1497, :).^2/2, z; linewidth = 1)
# lines!(ax_b, 1e6*interior(snapshot[:w], 3004, 1497, :).^2/2, z; linewidth = 1)
# hlines!(ax_b, -interior(snapshot[:BLD], :, 1497, 1)[[2004,1004,3004]], linestyle = :dash, color = [:orange, :green, :purple], linewidth = 0.8)
# colsize!(gab, 3, Relative(0.3))
# colgap!(gab, 1, 1)
# resize_to_layout!(fig)
# save(filesave * "w_" * fileparam * "_2d_iter$(iteration).pdf", fig; pt_per_unit = 1)
# println("Finished plotting w fields")

##############################
# iteration = 72635
# # 1. Define the filename of the saved snapshot
# output_filename = filehead * "subdomains/" * fileparam * "_snapshot_iter$(iteration).jld2"

# # 2. Load the snapshot using the new function
# snapshot = load_subdomain_snapshot(output_filename)
# x, y, z = nodes(snapshot[:T]);
# _, _, zw = nodes(snapshot[:w]);
# k = 77
# Tmap, wmap, vmap = :thermal,:delta,:balance
# wmax,umax,vmax=0.01,0.15,0.2
# #fig = Figure(size = (640, 450))
# fig = Figure(size = (640, 750))
# gabc = fig[1, 1] = GridLayout()
# aspect = 1
# axis_kwargs = (ylabel = L"y~\text{(km)}", aspect=aspect)
# ax_a = Axis(gabc[1,1]; titlealign = :left, title=L"\text{(a)}~T~\text{({^\circ}C)}", axis_kwargs...)
# ax_b = Axis(gabc[1,3]; titlealign = :left, title=L"\text{(b)}~w~\text{(m s^{-1})}", aspect=aspect)
# ax_c = Axis(gabc[3,1]; titlealign = :left, title=L"\text{(c)}~u~\text{(m s^{-1})}", axis_kwargs...) 
# ax_d = Axis(gabc[3,3]; titlealign = :left, title=L"\text{(d)}~v~\text{(m s^{-1})}", aspect=aspect) 
# hm_a = heatmap!(ax_a, 1e-3x, 1e-3y, (interior(snapshot[:T],:,:,k)); rasterize = true, colormap = Tmap)
# hm_b = heatmap!(ax_b, 1e-3x, 1e-3y, (interior(snapshot[:w],:,:,k)); rasterize = true, colormap = wmap, colorrange = (-wmax, wmax))
# hm_c = heatmap!(ax_c, 1e-3x, 1e-3y, (interior(snapshot[:u],:,:,k)); rasterize = true, colormap = vmap, colorrange = (-umax, umax))
# hm_d = heatmap!(ax_d, 1e-3x, 1e-3y, (interior(snapshot[:v],:,:,k)); rasterize = true, colormap = vmap, colorrange = (-vmax, vmax))
# Colorbar(gabc[1,2], hm_a)
# Colorbar(gabc[1,4], hm_b)
# Colorbar(gabc[3,2], hm_c)
# Colorbar(gabc[3,4], hm_d)
# hidexdecorations!(ax_a, ticks = false)
# hidexdecorations!(ax_b, ticks = false)
# hidexdecorations!(ax_c, ticks = false)
# hidexdecorations!(ax_d, ticks = false)
# hideydecorations!(ax_b, ticks = false)
# hideydecorations!(ax_d, ticks = false)

# # T̄ = (mean(snapshot[:T], dims = 2));
# # w̄ = (mean(snapshot[:w], dims = 2));
# # ū = (mean(snapshot[:u], dims = 2));
# # v̄ = (mean(snapshot[:v], dims = 2));
# zmin = -75
# kz = findfirst(z .≥ zmin)
# Nz = length(z)
# axis_kwargs0 = (xlabel = L"x~\text{(km)}", ylabel = L"z~\text{(m)}", limits = (nothing, (zmin, 0)))
# axis_kwargs1 = NamedTuple{(:xlabel,:ylabel)}(axis_kwargs0)
# ax_a = Axis(gabc[2,1]; titlealign = :left, axis_kwargs0...)
# ax_b = Axis(gabc[2,3]; titlealign = :left, xlabel = L"x~\text{(km)}", limits = (nothing, (zmin, 0)))
# ax_c = Axis(gabc[4,1]; titlealign = :left,  limits = (nothing, (zmin, 0)), axis_kwargs1...)
# ax_d = Axis(gabc[4,3]; titlealign = :left, xlabel = L"x~\text{(km)}", limits = (nothing, (zmin, 0)))
# wmax,umax,vmax=0.01,0.1,0.2
# hm_a = heatmap!(ax_a, 1e-3x, z[kz:Nz], (interior(snapshot[:T],:,3*640,kz:Nz)); rasterize = true, colormap = Tmap)#, colorrange = (vmin, vmax))
# hm_b = heatmap!(ax_b, 1e-3x, zw[kz:Nz], (interior(snapshot[:w],:,3*640,kz:Nz)); rasterize = true, colormap = wmap, colorrange = (-wmax, wmax))
# hm_c = heatmap!(ax_c, 1e-3x, z[kz:Nz], (interior(snapshot[:u],:,3*640,kz:Nz)); rasterize = true, colormap = vmap, colorrange = (-umax, umax))
# hm_d = heatmap!(ax_d, 1e-3x, z[kz:Nz], (interior(snapshot[:v],:,3*640,kz:Nz)); rasterize = true, colormap = vmap, colorrange = (-vmax, vmax))
# hideydecorations!(ax_b, ticks = false)
# hideydecorations!(ax_d, ticks = false)
# Colorbar(gabc[2,2], hm_a)
# Colorbar(gabc[2,4], hm_b)
# Colorbar(gabc[4,2], hm_c)
# Colorbar(gabc[4,4], hm_d)
# rowgap!(gabc, 4)
# colgap!(gabc, 1, 1)
# colgap!(gabc, 3, 1)
# colgap!(gabc, 2, 5)
# for row = [2,4]
#     rowsize!(gabc, row, Relative(0.14))
# end
# resize_to_layout!(fig)
# save(filesave * "Twuv_" * fileparam * "_3d_iter$(iteration).pdf", fig; pt_per_unit = 1)
# println("Finished plotting Twuv fields")

########################
# Compute the horizontal spectrum of T, u, v, w 
# window = :hann
# xT, yT, zT = nodes(snapshot[:T]);
# xu, yu, _ = nodes(snapshot[:u]);
# xv, yv, _ = nodes(snapshot[:v]);
# grid = snapshot[:u].grid
# u̅  = XFaceField(grid);
# v̅  = YFaceField(grid);
# cutoff = 300
# set_value!(; Δh = 4.8828125)
# t0 = now()
# coarse_graining!(snapshot[:u] , u̅ ; kernel=:lanczos, cutoff, levels = [77,58,27])
# coarse_graining!(snapshot[:v] , v̅ ; kernel=:lanczos, cutoff, levels = [77,58,27])
# @info "Coarse graining took $((now() - t0).value/1e3) seconds."
# axis_kwargs1 = (xlabel = L"\text{Wavenumber (m^{-1})}", xgridvisible = false,
#                 ylabel = L"E_i(k)/E_{T,v}(k_{\min},z=-3.9~\text{m})", ygridvisible = false,
#                 xscale = log10, yscale = log10,
#                 limits = ((1e-4, 6e-1), (1e-7,1e1)),
#                 xminorticks = [2e-4:1e-4:9e-4; 2e-3:1e-3:9e-3; 2e-2:1e-2:9e-2; 2e-1:1e-1:9e-1],xminorticksvisible = true)
# axis_kwargs2 = NamedTuple{(:xlabel,:xscale,:yscale,:limits,:xgridvisible,:ygridvisible,:xminorticks,:xminorticksvisible)}(axis_kwargs1)

# fig = Figure(size = (640, 250))
# g3 = fig[1, 1] = GridLayout()
# S0s = []
# alphabet = [letter for letter in 'a':'z'];
# for (i,klev) in enumerate([77, 58, 27])
#     println("Plotting spectra at z = $(zT[klev])m...")
#     if i == 3
#         ax = Axis(g3[1, i]; titlealign = :left, title=L"\text{(c)}~z=-60~\text{m}", axis_kwargs2...)
#         hideydecorations!(ax, ticks = false)
#     else
#         if i == 2
#             ax = Axis(g3[1, i]; titlealign = :left, title=L"\text{(b)}~z=-25~\text{m}", axis_kwargs2...)
#             hideydecorations!(ax, ticks = false)
#         else
#             ax = Axis(g3[1, 1]; titlealign = :left, title=L"\text{(a)}~z=-3.9~\text{m}", axis_kwargs1...)
#         end
#     end
#     vlines!(ax, 2π/300; color = :black, linewidth = 0.8)
#     Su = isotropic_powerspectrum(interior(snapshot[:u], :, :, klev), interior(snapshot[:u], :, :, klev), xu, yu;window)
#     Sv = isotropic_powerspectrum(interior(snapshot[:v], :, :, klev), interior(snapshot[:v], :, :, klev), xv, yv;window)
#     wk = (interior(snapshot[:w], :, :, klev)+interior(snapshot[:w], :, :, klev+1))/2
#     Sw = isotropic_powerspectrum(wk, wk, xT, yT;window)
#     St = isotropic_powerspectrum(interior(snapshot[:T], :, :, klev), interior(snapshot[:T], :, :, klev), xT, yT;window)
#     Su̅ = isotropic_powerspectrum(interior(u̅, :, :, klev), interior(u̅, :, :, klev), xu, yu;window)
#     Sv̅ = isotropic_powerspectrum(interior(v̅, :, :, klev), interior(v̅, :, :, klev), xv, yv;window)

#     if i == 1
#         global Sv0,St0 = Sv,St
#     end

#     idx = Su.freq .> 0
#     lines!(ax, Su.freq[idx], 0.5*(Su.freq[idx]./Su.freq[idx][1]).^-2, linestyle = :dash, color = :black)
#     #text!(ax, 10^-3, 0.7; text = L"k^{-2}")
#     #idxf = idx#0 .< Su.freq .<= 1e-3
#     #lines!(ax, Su.freq[idxf], 1e-13Su.freq[idxf].^-3, linestyle = :dash, color = :gray)
#     #text!(ax, 10^-4, 1e-3; text = L"k^{-3}")
#     lines!(ax, Su.freq[idx], 50*(Su.freq[idx]./Su.freq[idx][1]).^(-5/3), linestyle = :dash, color = :gray)
#     #text!(ax, 10^-3, 10^-5.5; text = L"k^{-5/3}")
#     lines!(ax, St.freq[idx], Real.(St.spec[idx]./St0.spec[idx][1]), color = :red, label = L"E_T")
#     lines!(ax, Su.freq[idx], Real.(Su.spec[idx]./Sv0.spec[idx][1]), color = :blue, label = L"E_u")
#     lines!(ax, Sv.freq[idx], Real.(Sv.spec[idx]./Sv0.spec[idx][1]), color = :green, label = L"E_v")
#     lines!(ax, Sw.freq[idx], Real.(Sw.spec[idx]./Sv0.spec[idx][1]), color = :black, label = L"E_w")
#     lines!(ax, Su̅.freq[idx], Real.(Su̅.spec[idx]./Sv0.spec[idx][1]), color = :blue, linestyle = :dash)
#     lines!(ax, Sv̅.freq[idx], Real.(Sv̅.spec[idx]./Sv0.spec[idx][1]), color = :green, linestyle = :dash)
#     if i == 1
#         axislegend(ax, labelsize=9, patchsize = (15, 1), framevisible = false,
#                    padding = (0f0, 0f0, 0f0, 0f0), patchlabelgap = 3, rowgap = 1)
#     end
#     resize_to_layout!(fig)
#     save(filesave * "spectra_" * fileparam * "_3d_iter$(iteration).pdf", fig; pt_per_unit = 1)
# end
# colgap!(g3, 3)
# resize_to_layout!(fig)
# save(filesave * "spectra_" * fileparam * "_3d_iter$(iteration).pdf", fig; pt_per_unit = 1)
# println("Finished plotting spectra.")

###########################
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