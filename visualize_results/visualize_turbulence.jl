using LESStudySetup
using CairoMakie
using SixelTerm
using Printf, Dates
using Oceananigans: compute!
using Oceananigans.Grids: nodes
using Statistics: mean
using LESStudySetup.Diagnostics
using LESStudySetup.Diagnostics: load_snapshots, BLD1D, N², KE
using LinearAlgebra: \
set_theme!(Theme(fontsize = 12))

filename = "turbulence_generator_output.jld2"
metadata = "turbulence_generator_metadata.jld2"

# load all the data!!
println("Loading data from $filename...")
snapshots = load_snapshots(filename; metadata)

# Let's pick the last snapshot!
times = snapshots[:T].times
snapshot_number = length(times)

x, y, z = nodes(snapshots[:w][snapshot_number])
xu, yu, zu = nodes(snapshots[:u][snapshot_number])
A = [ones(length(zu)) zu]
A2 = A'*A
hbs = zeros(snapshot_number)
hms = zeros(snapshot_number)
hes = zeros(snapshot_number)
t0 = now()
for i = 2:snapshot_number
    println("Loading fields $i wall time: $((now() - t0).value/1e3) seconds.")
    hbs[i] = mean(compute!(BLD1D(snapshots,i)))
    hms[i] = -z[argmax(vec(mean(N²(snapshots, i),dims = (1,2))))]
    mKE = vec(mean((KE(snapshots, i)), dims = (1,2)))
    kKE = findlast(mKE .< 1e-20)
    hKE = kKE == nothing ? -z[1] : -z[kKE]
    hes[i] = min(hKE, 2/(A2\(A' * log.(mKE)))[2])
end
fig = Figure(size = (800, 300))
ax = Axis(fig[1,1];xticks = 0:10,yticks = 0:25:75, xlabel="Time (hours)", ylabel="Depth (m)",limits=((0,10), (0, 75)))
lines!(ax, times[2:end]/3600, hbs[2:end]; label = "boundary layer")
lines!(ax, times[2:end]/3600, hms[2:end]; label = "max N² layer")
lines!(ax, times[2:end]/3600, hes[2:end]; label = "Ekman layer")
axislegend(ax; position = :ct, orientation = :horizontal)
save("results/initial_turbulence_h.pdf", fig)

α  = parameters.α
g  = parameters.g
fig = Figure(size = (500, 400))
axis_kwargs1 = (xlabel = "buoyancy", ylabel = "z (m)",
                limits = ((0, 100), (0, 100))) 
axis_kwargs2 = (xlabel = "velocity²", 
                limits = ((0, 100), (0, 100))) 
axis_kwargs3 = (xlabel = "Ri", 
                limits = ((0, 100), (0, 100))) 
ax_B = Axis(fig[2, 1]; title = "Buoyancy", axis_kwargs...)
ax_u = Axis(fig[2, 2]; title = "Velocity", axis_kwargs...)
ax_ri = Axis(fig[2, 3]; title = "Velocity", axis_kwargs...)
n = Observable(1)

B  = @lift vec(mean(compute!(Field(α * g * snapshots[:T][$n])), dims = (1, 2)))
uc = @lift snapshots[:u][$n]
vc = @lift snapshots[:v][$n]
wc = @lift compute!(Field(@at (Center, Center, Center) snapshots[:w][$n]))
Ut² = @lift (var($uc,dims=(1, 2))+var($vc,dims=(1, 2))+var($wc,dims=(1, 2)))/3
uc = @lift vec(mean($uc, dims = (1, 2)))
vc = @lift vec(mean($vc, dims = (1, 2)))
wc = @lift vec(mean($wc, dims = (1, 2)))
Δb = @lift $B[end] .- $B
Δu² = @lift (($uc[end] .- $uc).^2 .+ ($vc[end] .- $vc).^2 .+ ($wc[end] .- $wc).^2)/3
Ri = @lift -zT.* $Δb ./ ($Δu² .+ $Ut²)

lines!(ax_B, Δb, zT)
lines!(ax_u, Δu², zT)
lines!(ax_u, Ut², zT)
lines!(ax_ri, Ri, zT)
hideydecorations!(ax_u, ticks = false)
hideydecorations!(ax_ri, ticks = false)

frames = 1:2:length(times)
@info "Making a neat animation of Ri"
record(fig, "results/initial_turbulence_ri.mp4", frames, framerate=4) do i
    println("Loading fields $i wall time: $((now() - t0).value/1e3) seconds.")
    n[] = i
end

u = snapshots[:u][snapshot_number]
v = snapshots[:v][snapshot_number]
w = snapshots[:w][43]

grid = u.grid

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

u = snapshots[:u][7];
u = compute!(Field(u - mean(u)))
u_slices = (east   = interior(u, 1, :, :),
            north  = interior(u, :, 1, :),
            top    = interior(u, :, :, Nz-1))

w_slices = (east   = interior(w, 1, :, :),
            north  = interior(w, :, 1, :),
            top    = interior(w, :, :, Nz-2))

clim = 0.5 * maximum(abs.(w_slices.top))

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

kwargs = (colorrange=(-clim,clim), colormap=:balance, shading=NoShading,rasterize = true)

surface!(ax, x_yz_east, y_yz, z_yz;  color = w_slices.east, kwargs...)
surface!(ax, x_xz, y_xz_north, z_xz; color = w_slices.north, kwargs...)
sf = surface!(ax, x_xy, y_xy, z_xy_top;  color = w_slices.top, kwargs...)

Colorbar(fig[2, 2], sf, label = "m s⁻¹", height = Relative(0.4), tellheight=false)

title = "Vertical velocity at t = 7 hours"
fig[1, 1:2] = Label(fig, title; fontsize = 24, tellwidth = false, padding = (0, 0, -120, 0))

rowgap!(fig.layout, 1, Relative(-0.2))
colgap!(fig.layout, 1, Relative(-0.1))

save("results/initial_turbulence_w3d7h.pdf", fig)