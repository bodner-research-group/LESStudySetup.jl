using LESStudySetup
using CairoMakie
using Printf, Dates
using Statistics: mean, std
using LESStudySetup.Diagnostics
using LESStudySetup.Diagnostics: load_distributed_checkpoint,load_subdomain_snapshot
using LESStudySetup.Diagnostics: isotropic_powerspectrum
using LESStudySetup.Diagnostics: coarse_graining!, TKE, MLD
set_theme!(theme_latexfonts(), fontsize=12, figure_padding = 10)
using JLD2, CUDA
set_value!(; Δh = 4.8828125)
filehead = "/orcd/data/abodner/002/shared_datasets/nhyles_output/" 
filesave = "results/"
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
shift(x) = [x[size(x,1)÷2+1:end, :]; x[1:size(x,1)÷2, :]]
filename0 = "./hydrostatic_snapshots_init.jld2"
snapshots = load_snapshots(filename0);
v0 = snapshots[:v][1];
u0 = snapshots[:u][1];
T0 = snapshots[:T][1];
xu, yu, zu = nodes(u0);
xv, yv, zv = nodes(v0);
xT, yT, zT = nodes(T0);
_ , _ , zw = nodes(snapshots[:w][1])
h0 = MLD(snapshots,1; threshold = 0.09)
f = parameters.f;
ζ₀ = compute!(Field(∂x(v0) - ∂y(u0)));
x, y, z = nodes(ζ₀);
Nz = length(z)
initfile = "./hydrostatic_snapshots_free.jld2"
initsnaps = load_snapshots(initfile)
Ub = initsnaps[:u][1];
Vb = compute!(Field(initsnaps[:v][1] - snapshots[:v][1]));
σn = compute!(Field(-(∂x(Ub)-∂y(Vb))/2));
σs = compute!(Field((∂y(Ub)+∂x(Vb))/2));
@info "σn extrema: $(extrema(interior(σn)))"
@info "σs extrema: $(extrema(interior(σs)))"
iarrow = 2:32:640
k = length(zT)

# fig = Figure(size = (640, 270))
# gab = fig[1, 1] = GridLayout()
# axis_kwargs1 = (titlealign = :left, xlabel = "x (km)", ylabel = "y (km)", aspect = 1, limits = ((-50, 50), (0, 100)))
# axis_kwargs2 = NamedTuple{(:titlealign,:xlabel,:limits,:aspect)}(axis_kwargs1)
# ax_σ = Axis(gab[1,1]; title=L"\text{(a)}~S_n~\text{(10^{-6}s^{-1})},~z=-0.56~\text{m}", axis_kwargs1...)
# ax_T = Axis(gab[1,3]; title=L"\text{(b)}~T_i~\text{({^\circ}C)},~z=-0.56~\text{m}", axis_kwargs1...)
# ax_v = Axis(gab[1,5]; title=L"\text{(c)}~v_0~\text{(m~s^{-1})},~z=-0.56~\text{m}", axis_kwargs1...)
# σmin, σmax = 1.1e6*minimum(interior(σn,:,:,k)), 1.1e6*maximum(interior(σn,:,:,k))
# hm_σ = heatmap!(ax_σ, 1e-3xT.-50, 1e-3yT, 1e6*shift(interior(σn,:,:,k)); rasterize = true, colormap = :diff, colorrange = (σmin, σmax))
# Colorbar(gab[1, 2], hm_σ)
# hidexdecorations!(ax_σ, ticks = false)
# arrows!(ax_σ,xT[iarrow]/1e3.-50, yT[iarrow]/1e3, shift(interior(Ub,iarrow,iarrow,k)), shift(interior(Vb,iarrow,iarrow,k)), arrowsize = 3, lengthscale = 1e2,linecolor = :white, arrowcolor = :white, linewidth = 0.2)
# Tmin, Tmax = minimum(interior(T0,:,:,k))-0.5, maximum(interior(T0,:,:,k))
# hm_T = heatmap!(ax_T, 1e-3xT.-50, 1e-3yT, shift(interior(T0,:,:,k)); rasterize = true, colormap = :thermal, colorrange = (Tmin, Tmax))
# hidexdecorations!(ax_T, ticks = false)
# hideydecorations!(ax_T, ticks = false)
# Colorbar(gab[1, 4], hm_T)
# arrows!(ax_T, [0], [50], [-15*sqrt(3)],[-15],linecolor = :black, arrowcolor =:black)
# text!(ax_T, -10*sqrt(3), 50, text = L"\text{Wind}", color = :black, align = (:right, :top))
# vbnd = maximum(abs.(interior(v0,:,:,k)))
# hm_v = heatmap!(ax_v, 1e-3xv.-50, 1e-3yv, shift(interior(v0,:,:,k)); rasterize = true, colormap = :balance, colorrange = (-vbnd, vbnd))
# hideydecorations!(ax_v, ticks = false)
# hidexdecorations!(ax_v, ticks = false)
# Colorbar(gab[1, 6], hm_v)
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
# ax_c = Axis(gab[2,5]; titlealign = :left, title=L"y=50~\text{km}", axis_kwargs0...)
# hm_a = heatmap!(ax_a, 1e-3xT.-50, zT[kz:Nz], 1e6*shift(interior(σn, :, 1, kz:Nz)); rasterize = true, colormap = :diff, colorrange = (σmin, σmax))
# Colorbar(gab[2, 2], hm_a)
# Tmin, Tmax = minimum(interior(T0,:,ja,kz:Nz)), maximum(interior(T0,:,ja,kz:Nz))
# hm_b = heatmap!(ax_b, 1e-3xT.-50, zT[kz:Nz], shift(interior(T0,:,ja,kz:Nz)); rasterize = true, colormap = :thermal, colorrange = (Tmin, Tmax))
# hideydecorations!(ax_b, ticks = false)
# Colorbar(gab[2, 4], hm_b)
# lines!(ax_b, 1e-3xT.-50, -vec(shift(interior(h0,:,ja))); color = :white, linewidth = 0.8)
# hm_c = heatmap!(ax_c, 1e-3xv.-50, zv[kz:Nz], shift(interior(v0,:,ja,kz:Nz)); rasterize = true, colormap = :balance, colorrange = (-vbnd, vbnd))
# hideydecorations!(ax_c, ticks = false)
# idxz = 1:10:Nz
# arrows!(ax_a,xu[iarrow]/1e3.-50, zu[idxz], shift(interior(Ub,iarrow,1,idxz)), 0*shift(interior(Ub,iarrow,1,idxz)), arrowsize = 3, lengthscale = 1e2,linecolor = :white, arrowcolor = :white, linewidth = 0.2)
# arrows!(ax_c,xu[iarrow]/1e3.-50, zu[idxz], shift(interior(Ub,iarrow,ja,idxz)), 0*shift(interior(Ub,iarrow,ja,idxz)), arrowsize = 3, lengthscale = 1e2,linecolor = :black, arrowcolor = :black, linewidth = 0.2)
# Colorbar(gab[2, 6], hm_c,ticks = -0.2:0.2:0.2)

# rowgap!(gab, 0)
# colgap!(gab, 1, 1)
# colgap!(gab, 3, 1)
# colgap!(gab, 5, 1)
# colgap!(gab, 2, 10)
# colgap!(gab, 4, 10)
# rowsize!(gab, 2, Relative(0.22))
# resize_to_layout!(fig)
# save(filesave * "sMT0v0fields.pdf", fig)

use_gpu = true
A=rand(Float32,10240,20480,1);
if use_gpu && CUDA.functional()
    p = CUDA.CUFFT.plan_rfft(CuArray(A), (1,2));
    ip = CUDA.CUFFT.plan_irfft(p * CuArray(A), 10240, (1,2));
else
    @error "CUDA not available"
end
fileparam = "xband1sublevels"
fig = Figure(size = (640, 580))
g4 = fig[1, 1] = GridLayout()
aspect = 0.25
crange = (-10, 10)
axis_kwargs = (xlabel = L"x~\text{(km)}", limits = ((-12.5,12.5),(0,100)), aspect=aspect)
ax_a = Axis(g4[1,1]; titlealign = :left, title=L"\text{(a)}~\zeta_0/f", ylabel = L"y~\text{(km)}", axis_kwargs...)
ax_b = Axis(g4[1,2]; titlealign = :left, title=L"\text{(b)}~\overline{\zeta}/f,~t=16~\text{h}", axis_kwargs...)
ax_c = Axis(g4[1,3]; titlealign = :left, title=L"\text{(c)}~\overline{\zeta}/f,~t=30~\text{h}", axis_kwargs...) 
ax_d = Axis(g4[1,4]; titlealign = :left, title=L"\text{(d)}~\overline{\zeta}/f,~t=44~\text{h}", axis_kwargs...) 
hideydecorations!(ax_b, ticks = false)
hideydecorations!(ax_c, ticks = false)
hideydecorations!(ax_d, ticks = false)
hm_a = heatmap!(ax_a, 1e-3x .- 50, 1e-3y, shift(interior(ζ₀,:,:,Nz))./f; rasterize = true, colormap = :curl, colorrange = crange)
arrows!(ax_a,xT[iarrow]/1e3.-50, yT[iarrow]/1e3, shift(interior(Ub,iarrow,iarrow,k)), shift(interior(Vb,iarrow,iarrow,k)), arrowsize = 3, lengthscale = 1e2,linecolor = :black, arrowcolor = :black, linewidth = 0.2,alpha=0.5)
axs = [ax_b, ax_c, ax_d]
iterations = [26066,37003,49086]
for (i,iteration) in enumerate(iterations)
    output_filename = filehead * "subdomains/" * fileparam * "_snapshot_iter$(iteration).jld2"
    snapshot = load_subdomain_snapshot(output_filename;variables = ("v","u"),level=224);
    grid = snapshot[:u].grid
    u̅  = XFaceField(grid);
    v̅  = YFaceField(grid);
    _, _, zi = nodes(snapshot[:v]);
    Nzi = length(zi)
    coarse_graining!(snapshot[:v] , v̅; kernel=:gaussian, cutoff=100, border = :ycircular, method = :spectral, use_gpu,plans=(p,ip));
    coarse_graining!(snapshot[:u] , u̅; kernel=:gaussian, cutoff=100, border = :ycircular, method = :spectral, use_gpu,plans=(p,ip));
    ζ̄i = compute!(Field((∂x(v̅)-∂y(u̅))));
    xi, yi, _ = nodes(ζ̄i);
    hm_i = heatmap!(axs[i], 1e-3xi, 1e-3yi, (interior(ζ̄i,:,:,Nzi))./f; rasterize = true, colormap = :curl, colorrange = crange)
    arrows!(axs[i],xT[iarrow]/1e3.-50, yT[iarrow]/1e3, shift(interior(Ub,iarrow,iarrow,k)), shift(interior(Vb,iarrow,iarrow,k)), arrowsize = 3, lengthscale = 1e2,linecolor = :black, arrowcolor = :black, linewidth = 0.2,alpha=0.5)
end
Colorbar(g4[1,5], hm_a)
colgap!(g4, 1, 5)
colgap!(g4, 2, 5)
colgap!(g4, 3, 5)
colgap!(g4, 4, 1)
resize_to_layout!(fig)
save(filesave * "curl_" * fileparam * "_0_16_30_44h.pdf", fig; pt_per_unit = 1)

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
# fileparam = "sublevels"

# iterations = get_iterations_regex(filehead, fileparam)
# cutoff = 100

# mMLDs = zeros(length(iterations), 3)
# for (i,iteration) in enumerate(iterations)
#     # --- Load the snapshot ---
#     output_filename = filehead * "subdomains/" * fileparam * "_snapshot_iter$(iteration).jld2"
#     snapshot = load_subdomain_snapshot(output_filename; variables = ("MLD", "MLD2", "MLD3"));
#     mMLDs[i, 1] = mean(snapshot[:MLD])
#     mMLDs[i, 2] = mean(snapshot[:MLD2])
#     mMLDs[i, 3] = mean(snapshot[:MLD3])
#     @info "Iteration $iteration: mMLDs = $(mMLDs[i, :])"
# end

# fig = Figure(size = (640, 320))
# t = 24:2:length(iterations)*2+22
# ax1 = Axis(fig[1, 1], ylabel = "mMLD (m)")
# ax2 = Axis(fig[2, 1], ylabel = "mMLD2 (m)")
# ax3 = Axis(fig[3, 1], xlabel = "Time (h)", ylabel = "mMLD3 (m)")
# lines!(ax1, t, mMLDs[:, 1], label = "MLD")
# lines!(ax2, t, mMLDs[:, 2], label = "MLD2")
# lines!(ax3, t, mMLDs[:, 3], label = "MLD3")
# hidexdecorations!(ax1, grid = false)
# hidexdecorations!(ax2, grid = false)
# resize_to_layout!(fig)
# save(filesave * "mMLD_" * fileparam * "_1d_nd.pdf", fig; pt_per_unit = 1)

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
# snapshot = load_subdomain_snapshot(output_filename; level = 198, variables = ("u", "v", "w"));

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
# L_filter = 100
# ax = Axis(g3[1, 1]; title=L"\text{Nonhydrostatic LES spectra at 48~h}~z=-4~\text{m}", axis_kwargs1...)
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
# lines!(ax, Sv.freq[idx], Real.(Sv.specf[idx]./Sv.spec[idx][1]), color = :green, linestyle = :dash)
# lines!(ax, Sw.freq[idx], Real.(Sw.spec[idx]./Sv.spec[idx][1]), color = :black, label = L"E_w")
# axislegend(ax, labelsize=9, patchsize = (15, 1), framevisible = false,
#             padding = (0f0, 0f0, 0f0, 0f0), patchlabelgap = 3, rowgap = 1)
# resize_to_layout!(fig)
# save(filesave * "spectra_" * fileparam * "_2d4m_iter$(iteration).pdf", fig; pt_per_unit = 1)

# τuu, τvv, τww, _, _, _ = TKE(snapshot; border=:reflect, Lx = snapshot[:grid].Lx, Ly = snapshot[:grid].Ly);

# jldopen(filehead * "TKEphysical_"*fileparam*"_iter$(iteration).jld2", "w") do file
#     file["fields/τuu"] = (interior(τuu, 1:5120, 640, :) .+ interior(τuu,2:5121, 640, :))/2
#     file["fields/τvv"] = (interior(τvv, :, 640, :) .+ interior(τvv, :, 641, :))/2
#     file["fields/τww"] = (interior(τww, :, 640, 1:224) .+ interior(τww, :, 640, 2:225))/2
    
#     file["metadata/iteration"] = iteration
#     file["metadata/xlims"] = (-12500,12500)
#     file["metadata/ylims"] = (5e4, 5e4)
#     file["metadata/zlims"] = (-252,0)
#     file["metadata/Nx"] = 640*8
#     file["metadata/Nz"] = 224
# end

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

# 3. Plot the TKE fields
# using Makie
# fig = Figure(size = (640, 320))
# gab = fig[1, 1] = GridLayout()
# #x, y, z = nodes(τww);
# var = TKEw;#interior(τww, :, 640,:);
# ax_a = Axis(gab[1,1]; titlealign = :left, title=L"\text{(a) TKE}_w~\text{(mm^2 s^{-2})}", xlabel=L"x~\text{(km)}", ylabel=L"z~\text{(m)}",limits=(nothing,(-90,0)))
# ax_b = Axis(gab[1,3]; titlealign = :left, title=L"\text{(b)}", xlabel=L"\langle\text{TKE}_w\rangle~\text{(mm^2 s^{-2})}",limits=(nothing,(-90,0)))
# hm_a = heatmap!(ax_a, 1e-3x, z, 1e6*var; rasterize = true, colormap = :amp, colorrange = (0, 70))
# Colorbar(gab[1,2], hm_a)
# # lines!(ax_a, 1e-3x, -interior(snapshot[:BLD], :, 1497, 1), color = :black, linewidth = 0.5, alpha=0.8)
# # #lines!(ax_a, 1e-3x, -interior(snapshot[:MLD], :, 1497, 1), color = :red, linewidth = 1)
# # #lines!(ax_a, 1e-3x, -interior(snapshot[:MLD2], :, 1497, 1), color = :blue, linewidth = 1)
# # lines!(ax_a, 1e-3x, -interior(snapshot[:MLD3], :, 1497, 1), color = :green, linewidth = 1)
# idxs = [1400,2834,3240]
# vlines!(ax_a, 1e-3x[idxs], color = Makie.wong_colors()[2:4], linewidth = 0.8)
# lines!(ax_a, 1e-3x, -BLDw, color = :black, linewidth = 1, label="BLD")
# lines!(ax_a, 1e-3x, -interior(snapshot[:MLD3],:,640,1), color = :blue, linewidth = 1, label="MLD3")
# axislegend(ax_a, labelsize=10,position = :lb,patchsize = (15, 5), patchlabelgap = 3, rowgap = 1)
# hideydecorations!(ax_b, ticks = false)
# lines!(ax_b, 1e6*vec(mean(var, dims =(1))), z; linewidth = 2, label="mean")
# lines!(ax_b, 1e6*var[idxs[1],:], z; linewidth = 1)
# lines!(ax_b, 1e6*var[idxs[2],:], z; linewidth = 1)
# lines!(ax_b, 1e6*var[idxs[3],:], z; linewidth = 1)
# axislegend(ax_b,  labelsize=10, position = :rb,patchsize = (15, 3), patchlabelgap = 3)
# # hlines!(ax_b, -interior(snapshot[:BLD], :, 1497, 1)[[2004,1004,3004]], linestyle = :dash, color = [:orange, :green, :purple], linewidth = 0.8)
# colsize!(gab, 3, Relative(0.3))
# colgap!(gab, 1, 1)
# resize_to_layout!(fig)
# save(filesave * "TKEw_" * fileparam * "_2d_iter$(iteration).pdf", fig; pt_per_unit = 1)
# println("Finished plotting TKE fields")

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