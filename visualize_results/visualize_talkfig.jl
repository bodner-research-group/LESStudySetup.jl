using LESStudySetup
using CairoMakie
using SixelTerm
using Printf, Dates
using Oceananigans: compute!
using Oceananigans.BoundaryConditions
using Oceananigans.Grids: xnodes, ynodes, znodes
using Statistics: mean
using LESStudySetup.Diagnostics
using LESStudySetup.Oceananigans.Units
using LESStudySetup.Diagnostics: N², M²
using LESStudySetup.Diagnostics: load_snapshots, MLD
set_theme!(theme_latexfonts(), fontsize=15,figure_padding = 10)

fig = Figure(size = (800, 360))
ax = Axis3(fig[2, 1],
           aspect=(1, 1, 0.5),
           xlabel = "x (km)",
           ylabel = "y (km)",
           zlabel = "z (m)",
           xlabeloffset = 50,
           ylabeloffset = 50,
           zlabeloffset = 50,
           limits = ((x[1]/1e3, x[end]/1e3), (y[1]/1e3, y[end]/1e3), (z[1], z[end])),
           elevation = 0.45,
           #azimuth = 6.8,
           xspinesvisible = false,
           zgridvisible = false,
           protrusions = 40,
           perspectiveness = 0.7)
           
kwargs = (colorrange=(Tmin,Tmax), colormap=:thermal, shading=NoShading,rasterize = true)
surface!(ax, x_yz_west/1e3, y_yz/1e3, z_yz;  color = T_slices.west, kwargs...)
surface!(ax, x_xz/1e3, y_xz_south/1e3, z_xz; color = T_slices.south, kwargs...)
sf = surface!(ax, x_xy/1e3, y_xy/1e3, z_xy_top;  color = T_slices.top, kwargs...)
Colorbar(fig[2, 2], sf, height = Relative(0.3), tellheight=false)
title = "Initial temperature (°C)"
fig[1, 1:2] = Label(fig, title; fontsize = 20, tellwidth = false, padding = (0, 0, -120, 0))

ax = Axis3(fig[2, 3],
                  aspect=(1, 1, 0.5),
                  xlabel = "x (km)",
                  ylabel = "y (km)",
                  zlabel = "z (m)",
                  xlabeloffset = 50,
                  ylabeloffset = 50,
                  zlabeloffset = 50,
                  limits = ((x[1]/1e3, x[end]/1e3), (y[1]/1e3, y[end]/1e3), (z[1], z[end])),
                  elevation = 0.45,
                  #azimuth = 6.8,
                  xspinesvisible = false,
                  zgridvisible = false,
                  protrusions = 40,
                  perspectiveness = 0.7)

kwargs = (colorrange=(-clim,clim), colormap=:balance, shading=NoShading,rasterize = true)
surface!(ax, x_yz_west/1e3, y_yz/1e3, z_yz;  color = v_slices.west, kwargs...)
surface!(ax, x_yz_west2/1e3, y_yz/1e3, z_yz;  color = v_slices.west2, kwargs...)
sf = surface!(ax, x_xy/1e3, y_xy/1e3, z_xy_top;  color = v_slices.top, kwargs...)
Colorbar(fig[2, 4], sf, height = Relative(0.3), tellheight=false)
title = "Initial frontal velocity (m s⁻¹)"
hidezdecorations!(ax, ticks = false)
fig[1, 3:4] = Label(fig, title; fontsize = 20, tellwidth = false, padding = (0, 0, -120, 0))
arrows!(ax,x[1:10:200]/1e3, y[1:10:200]/1e3, interior(ub,1:10:200,1:10:200,ks), interior(vb,1:10:200,1:10:200,ks), arrowsize = 3, lengthscale = 1e2,linecolor = :gray, arrowcolor = :black, linewidth = 0.2)

rowgap!(fig.layout, 1, Relative(-0.3))
colgap!(fig.layout, 1, Relative(-0.05))
colgap!(fig.layout, 2, Relative(0))
colgap!(fig.layout, 3, Relative(-0.05))
save("initial_Tv.pdf", fig)