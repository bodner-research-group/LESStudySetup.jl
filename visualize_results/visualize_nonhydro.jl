using LESStudySetup
using CairoMakie
using Printf, Dates
using Statistics: mean, std
using LESStudySetup.Diagnostics
using LESStudySetup.Diagnostics: load_distributed_checkpoint,load_subdomain_snapshot
using LESStudySetup.Diagnostics: isotropic_powerspectrum
using LESStudySetup.Diagnostics: coarse_graining!
set_theme!(theme_latexfonts(), fontsize=12, figure_padding = 10)

# --- Simulation and Subdomain Parameters ---
filehead = "/orcd/data/abodner/002/nhyles_output/" 
fileparam = "subdomain3"
filesave = "results/"
iterations = [32207,52543,72635]#72635 #52543 #32207 
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

# fig = Figure(size = (640, 320))
# gab = fig[1, 1] = GridLayout()
# x, y, z = nodes(snapshot[:w]);
# ax_a = Axis(gab[1,1]; titlealign = :left, title=L"\text{(a)}~w~\text{(mm~s^{-1})}", xlabel=L"x~\text{(km)}", ylabel=L"z~\text{(m)}",limits=(nothing,(-250,0)))
# ax_b = Axis(gab[1,3]; titlealign = :left, title=L"\text{(b)}", xlabel=L"10^6\langle\text{KE}_w\rangle~\text{(m^2~s^{-2})}",limits=(nothing,(-250,0)))
# hm_a = heatmap!(ax_a, 1e-3x, z, 1e3*interior(snapshot[:w], :, 1, :); rasterize = true, colormap = :delta, colorrange = (-10, 10))
# Colorbar(gab[1,2], hm_a)
# hideydecorations!(ax_b, ticks = false)
# lines!(ax_b, 1e6*vec(mean(compute!(Field(snapshot[:w]^2/2)), dims =(1,2))), z; linewidth = 1)
# colsize!(gab, 3, Relative(0.3))
# colgap!(gab, 1, 1)
# resize_to_layout!(fig)
# save(filesave * "w_" * fileparam * "_36h_iter$(iteration).pdf", fig; pt_per_unit = 1)
# println("Finished plotting w fields")

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

###################################
# window = :xhann
# xT, yT, _ = nodes(snapshot[:T]);
# xu, yu, _ = nodes(snapshot[:u]);
# xv, yv, _ = nodes(snapshot[:v]);
# axis_kwargs1 = (xlabel = L"\text{Wavenumber (m^{-1})}", xgridvisible = false,
#                 ylabel = L"E_i(k)/E_{T,v}(k_{\min})", ygridvisible = false,
#                 xscale = log10, yscale = log10,
#                 limits = ((4e-5, 3e-1), (1e-6,1e1)),
#                 xminorticks = [4e-5:1e-5:9e-5; 2e-4:1e-4:9e-4; 2e-3:1e-3:9e-3; 2e-2:1e-2:9e-2; 2e-1:1e-1:3e-1],xminorticksvisible = true,
#                 yminorticks = [2e-6:1e-6:9e-6; 2e-5:1e-5:9e-5; 2e-4:1e-4:9e-4; 2e-3:1e-3:9e-3; 2e-2:1e-2:9e-2; 2e-1:1e-1:9e-1; 2e0:1e0:9e0],yminorticksvisible = true)

# fig = Figure(size = (640, 270))
# g3 = fig[1, 1] = GridLayout()
# alphabet = [letter for letter in 'a':'z'];
# klev = 5
# ax = Axis(g3[1, 1]; titlealign = :left, title=L"\text{(a)}~z=-3.9~\text{m}", axis_kwargs1...)
# vlines!(ax, 2π/300; color = :black, linewidth = 0.8)

# Su = isotropic_powerspectrum(interior(snapshot[:u], :, :, klev), interior(snapshot[:u], :, :, klev), xu, yu;window)
# Sv = isotropic_powerspectrum(interior(snapshot[:v], :, :, klev), interior(snapshot[:v], :, :, klev), xv, yv;window)
# wk = (interior(snapshot[:w], :, :, klev)+interior(snapshot[:w], :, :, klev+1))/2
# Sw = isotropic_powerspectrum(wk, wk, xT, yT;window)
# St = isotropic_powerspectrum(interior(snapshot[:T], :, :, klev), interior(snapshot[:T], :, :, klev), xT, yT;window)
# idx = Su.freq .> 0
# lines!(ax, Su.freq[idx], (Su.freq[idx]./Su.freq[idx][1]).^-2, linestyle = :dash, color = :black)
# text!(ax, 10^-3, 0.7; text = L"k^{-2}")
# lines!(ax, Su.freq[idx], 10*(Su.freq[idx]./Su.freq[idx][1]).^(-5/3), linestyle = :dash, color = :gray)
# text!(ax, 10^-3, 10^-5.5; text = L"k^{-5/3}")
# lines!(ax, St.freq[idx], Real.(St.spec[idx]./St.spec[idx][1]), color = :red, label = L"E_T")
# lines!(ax, Su.freq[idx], Real.(Su.spec[idx]./Sv.spec[idx][1]), color = :blue, label = L"E_u")
# lines!(ax, Sv.freq[idx], Real.(Sv.spec[idx]./Sv.spec[idx][1]), color = :green, label = L"E_v")
# lines!(ax, Sw.freq[idx], Real.(Sw.spec[idx]./Sv.spec[idx][1]), color = :black, label = L"E_w")
# axislegend(ax, labelsize=9, patchsize = (15, 1), framevisible = false,
#             padding = (0f0, 0f0, 0f0, 0f0), patchlabelgap = 3, rowgap = 1)
# resize_to_layout!(fig)
# save(filesave * "spectra_" * fileparam * "_3d4m_iter$(iteration).pdf", fig; pt_per_unit = 1)

##########################
shift(x) = [x[size(x,1)÷2+1:end, :]; x[1:size(x,1)÷2, :]]
filename0 = "./hydrostatic_snapshots_init.jld2"
snapshots = load_snapshots(filename0);
v = snapshots[:v][1];
u = snapshots[:u][1];
f = parameters.f;
ζ₀ = compute!(Field(∂x(v) - ∂y(u)));
x, y, z = nodes(ζ₀);
Nz = length(z)
fig = Figure(size = (640, 580))
g4 = fig[1, 1] = GridLayout()
aspect = 0.25
crange = (-10, 10)
axis_kwargs = (xlabel = L"x~\text{(km)}", limits = ((-12.5,12.5),(0,100)), aspect=aspect)
ax_a = Axis(g4[1,1]; titlealign = :left, title=L"\text{(a)}~\zeta_0/f", ylabel = L"y~\text{(km)}", axis_kwargs...)
ax_b = Axis(g4[1,2]; titlealign = :left, title=L"\text{(b)}~\overline{\zeta_1}/f", axis_kwargs...)
ax_c = Axis(g4[1,3]; titlealign = :left, title=L"\text{(c)}~\overline{\zeta_2}/f", axis_kwargs...) 
ax_d = Axis(g4[1,4]; titlealign = :left, title=L"\text{(d)}~\overline{\zeta_3}/f", axis_kwargs...) 
hideydecorations!(ax_b, ticks = false)
hideydecorations!(ax_c, ticks = false)
hideydecorations!(ax_d, ticks = false)
hm_a = heatmap!(ax_a, 1e-3x .- 50, 1e-3y, shift(interior(ζ₀,:,:,Nz))./f; rasterize = true, colormap = :curl, colorrange = crange)
axs = [ax_b, ax_c, ax_d]
cutoff = 300
set_value!(; Δh = 4.8828125)
for (i,iteration) in enumerate(iterations)
    output_filename = filehead * "subdomains/" * fileparam * "_snapshot_iter$(iteration).jld2"
    snapshot = load_subdomain_snapshot(output_filename)
    grid = snapshot[:u].grid
    u̅  = XFaceField(grid);
    v̅  = YFaceField(grid);
    xi, yi, zi = nodes(snapshot[:T]);
    Nzi = length(zi)
    coarse_graining!(snapshot[:u] , u̅ ; kernel=:lanczos, cutoff, levels = [Nzi], window=:xhann)
    coarse_graining!(snapshot[:v] , v̅ ; kernel=:lanczos, cutoff, levels = [Nzi], window=:xhann)
    ζ̄i = compute!(Field((∂x(v̅)-∂y(u̅))));
    hm_i = heatmap!(axs[i], 1e-3xi, 1e-3yi, (interior(ζ̄i,:,:,Nzi))./f; rasterize = true, colormap = :curl, colorrange = crange)
end
Colorbar(g4[1,5], hm_a)
colgap!(g4, 1, 5)
colgap!(g4, 2, 5)
colgap!(g4, 3, 5)
colgap!(g4, 4, 1)
resize_to_layout!(fig)
save(filesave * "curl_" * fileparam * "_0123d_4m.pdf", fig; pt_per_unit = 1)
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