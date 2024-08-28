using LESStudySetup
using CairoMakie
using SixelTerm
using Printf, Dates
using Oceananigans: compute!
using Oceananigans.Grids: nodes
using Statistics: mean
using LESStudySetup.Diagnostics
using LESStudySetup.Diagnostics: load_snapshots, BLD1D
set_theme!(Theme(fontsize = 12))

filename = "turbulence_generator_output.jld2"
metadata = "turbulence_generator_metadata.jld2"

# load all the data!!
println("Loading data from $filename...")
snapshots = load_snapshots(filename; metadata)

# Let's pick the last snapshot!
times = snapshots[:T].times
snapshot_number = length(times)

fig = Figure(size = (800, 300))
hs = zeros(snapshot_number)
for i = 2:snapshot_number
    hs[i] = mean(compute!(BLD1D(snapshots,i)))
end
ax = Axis(fig[1,1];xlabel="Time (hours)", ylabel="Boundary layer depth (m)")
lines!(ax, times[2:end]/3600, hs[2:end])
save("results/initial_turbulence_hb.pdf", fig)

u = snapshots[:u][snapshot_number]
v = snapshots[:v][snapshot_number]
w = snapshots[:w][snapshot_number]

grid = T.grid

# Coordinate arrays
xu, yu, zu = nodes(u)
xv, yv, zv = nodes(v)
x, y, z = nodes(w)

Nx, Ny, Nz = size(grid)
Nz += 1

# Slices
x_xz = repeat(x, 1, Nz)
y_xz_north = y[end] * ones(Nx, Nz)
z_xz = repeat(reshape(z, 1, Nz), Nx, 1)

x_yz_east = x[end] * ones(Ny, Nz)
y_yz = repeat(y, 1, Nz)
z_yz = repeat(reshape(z, 1, Nz), grid.Ny, 1)

x_xy = x
y_xy = y
z_xy_top = z[end] * ones(grid.Nx, grid.Ny)

fig = Figure(size = (800, 700))

ax = Axis3(fig[2, 1],
           aspect=(1, 1, 1),
           xlabel = "x (m)",
           ylabel = "y (m)",
           zlabel = "z (m)",
           xlabeloffset = 50,
           ylabeloffset = 50,
           zlabeloffset = 50,
           limits = ((x[1], x[end]), (y[1], y[end]), (z[1], z[end])),
           elevation = 0.45,
           azimuth = 6.8,
           xspinesvisible = false,
           zgridvisible = false,
           protrusions = 40,
           perspectiveness = 0.7)

w_slices = (east   = interior(w, 1, :, :),
            north  = interior(w, :, 1, :),
            top    = interior(w, :, :, Nz-2))

u = snapshots[:u][7];
u = compute!(Field(u - mean(u)))
u_slices = (east   = interior(u, 1, :, :),
            north  = interior(u, :, 1, :),
            top    = interior(u, :, :, Nz))

clim = 0.5 * maximum(abs.(w_slices.top))

kwargs = (colorrange=(-clim,clim), colormap=:balance, shading=NoShading,rasterize = true)

surface!(ax, x_yz_east, y_yz, z_yz;  color = w_slices.east, kwargs...)
surface!(ax, x_xz, y_xz_north, z_xz; color = w_slices.north, kwargs...)
sf = surface!(ax, x_xy, y_xy, z_xy_top;  color = w_slices.top, kwargs...)

Colorbar(fig[2, 2], sf, label = "m s⁻¹", height = Relative(0.4), tellheight=false)

title = "Vertical velocity at t = 4 hours"
fig[1, 1:2] = Label(fig, title; fontsize = 24, tellwidth = false, padding = (0, 0, -120, 0))

rowgap!(fig.layout, 1, Relative(-0.2))
colgap!(fig.layout, 1, Relative(-0.1))

save("results/initial_turbulence_w3d4h.pdf", fig)