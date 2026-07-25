#######################################################################
# windmixing_hydro_test.jl
#
# End‑to‑end demo:
#   0) Run a short NH “wind mixing + convection” simulation
#   1) Coarse‑grain the NH fields
#   2) Restart a HydrostaticFreeSurfaceModel with the coarse fields
#######################################################################

using Oceananigans, Oceananigans.OutputWriters
using JLD2, Statistics, Printf, Random
using Oceananigans.Units: minute, minutes, hour

# ────────────────────────────────────────────────────────────────────
# Stage 0:  WIND‑MIXING  NON‑HYDROSTATIC  TEST CASE
# Adapted from the Oceananigans example docs :contentReference[oaicite:0]{index=0}
# ────────────────────────────────────────────────────────────────────

@info "Stage 0 ▶ Running a short wind‑mixing NH simulation"

### the grid
Nx = Ny = 32     # number of points in each of horizontal directions
Nz = 24          # number of points in the vertical direction

Lx = Ly = 64     # (m) domain horizontal extents
Lz = 32          # (m) domain depth

refinement = 1.2 # controls spacing near surface (higher means finer spaced)
stretching = 12  # controls rate of stretching at bottom

# Normalized height ranging from 0 to 1
h(k) = (k - 1) / Nz

# Linear near-surface generator
ζ₀(k) = 1 + (h(k) - 1) / refinement

# Bottom-intensified stretching function
Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))

# Generating function
z_faces(k) = Lz * (ζ₀(k) * Σ(k) - 1)

grid = RectilinearGrid(CPU(), size = (Nx, Nx, Nz),
                          x = (0, Lx),
                          y = (0, Ly),
                          z = z_faces)

### set buoyancy
buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(thermal_expansion = 2e-4,
                                                                    haline_contraction = 8e-4))
                                                                    
### set boundary conditions
# temperature
Qʰ = 200.0  # W m⁻², surface _heat_ flux
ρₒ = 1026.0 # kg m⁻³, average density at the surface of the world ocean
cᴾ = 3991.0 # J K⁻¹ kg⁻¹, typical heat capacity for seawater
Qᵀ = Qʰ / (ρₒ * cᴾ) # K m s⁻¹, surface _temperature_ flux
dTdz = 0.01 # K m⁻¹
T_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵀ),
                                bottom = GradientBoundaryCondition(dTdz))

# velocity                   
u₁₀ = 10    # m s⁻¹, average wind velocity 10 meters above the ocean
cᴰ = 2.5e-3 # dimensionless drag coefficient
ρₐ = 1.225  # kg m⁻³, average density of air at sea-level
Qᵘ = - ρₐ / ρₒ * cᴰ * u₁₀ * abs(u₁₀) # m² s⁻²
u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵘ))

# salinity
@inline Qˢ(x, y, t, S, evaporation_rate) = - evaporation_rate * S # [salinity unit] m s⁻¹
evaporation_rate = 1e-3 / hour # m s⁻¹
evaporation_bc = FluxBoundaryCondition(Qˢ, field_dependencies=:S, parameters=evaporation_rate)
S_bcs = FieldBoundaryConditions(top=evaporation_bc)

### initialize model
model = NonhydrostaticModel(; grid, buoyancy,
                            advection = UpwindBiased(order=5),
                            tracers = (:T, :S),
                            coriolis = FPlane(f=1e-4),
                            closure = AnisotropicMinimumDissipation(),
                            boundary_conditions = (u=u_bcs, T=T_bcs, S=S_bcs))

### initial conditions
# Random noise damped at top and bottom
Ξ(z) = randn() * z / model.grid.Lz * (1 + z / model.grid.Lz) # noise

# Temperature initial condition: a stable density gradient with random noise superposed.
Tᵢ(x, y, z) = 20 + dTdz * z + dTdz * model.grid.Lz * 1e-6 * Ξ(z)

# Velocity initial condition: random noise scaled by the friction velocity.
uᵢ(x, y, z) = sqrt(abs(Qᵘ)) * 1e-3 * Ξ(z)

# `set!` the `model` fields using functions or constants:
set!(model, u=uᵢ, w=uᵢ, T=Tᵢ, S=35)

# set-up a simulation
simulation = Simulation(model, Δt=10.0, stop_time=40minutes)
wizard = TimeStepWizard(cfl=1.0, max_change=1.1, max_Δt=1minute)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

# Print a progress message
progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, max(|w|) = %.1e ms⁻¹, wall time: %s\n",
                                iteration(sim), prettytime(sim), prettytime(sim.Δt),
                                maximum(abs, sim.model.velocities.w), prettytime(sim.run_wall_time))

add_callback!(simulation, progress_message, IterationInterval(20))

### output
# Create a NamedTuple with eddy viscosity
# filename = "windmixing_hydro_test"

# simulation.output_writers[:slices] =
#     JLD2OutputWriter(model, merge(model.velocities, model.tracers),
#                      filename = filename * ".jld2",
#                      indices = (:, grid.Ny/2, :),
#                      schedule = TimeInterval(15minute),
#                      overwrite_existing = true)

run!(simulation)

# Snapshot and save fields
u_nh = Array(model.velocities.u)
v_nh = Array(model.velocities.v)
T_nh = Array(model.tracers.T)

# ────────────────────────────────────────────────────────────────────
# Stage 1:  COARSE‑GRAINING  BY BOX FILTER
# ────────────────────────────────────────────────────────────────────

@info "Stage 1 ▶ Coarse‑graining NH fields"

function box_filter(A::Array{T,3}, r::Int) where T
    B = similar(A)
    M,N,O = size(A)
    CI = CartesianIndices(A)
    @inbounds for I in CI
        i,j,k = Tuple(I)
        i0,i1 = max(i-r,1), min(i+r,M)
        j0,j1 = max(j-r,1), min(j+r,N)
        k0,k1 = max(k-r,1), min(k+r,O)
        B[i,j,k] = mean(@view A[i0:i1, j0:j1, k0:k1])
    end
    return B
end

r = 2                   # half‑width = 2 cells → 5³ stencil
u_cg = box_filter(u_nh, r)
v_cg = box_filter(v_nh, r)
T_cg = box_filter(T_nh, r)

# ────────────────────────────────────────────────────────────────────
# Stage 2:  HYDROSTATIC  RESTART
# ────────────────────────────────────────────────────────────────────

@info "Stage 2 ▶ Starting HydrostaticFreeSurfaceModel with coarse fields"

model_H = HydrostaticFreeSurfaceModel(; grid, buoyancy,
                                      tracers = (:T, :S),
                                      coriolis = FPlane(f=1e-4),
                                      boundary_conditions = (u=u_bcs, T=T_bcs, S=S_bcs))

set!(model_H, u = u_cg, v = v_cg, T = T_cg, S = 35)

sim_H = Simulation(model_H, Δt=10.0, stop_time=40minutes)
wizard = TimeStepWizard(cfl=0.25, max_change=1.1, max_Δt=1minute)
sim_H.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

# Print a progress message
progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, max(|w|) = %.1e ms⁻¹, wall time: %s\n",
                                iteration(sim), prettytime(sim), prettytime(sim.Δt),
                                maximum(abs, sim.model.velocities.w), prettytime(sim.run_wall_time))

add_callback!(sim_H, progress_message, IterationInterval(20))
run!(sim_H)

@info @sprintf "Finished!  NH→Hydro pipeline executed in %.2f s wall‑clock." sim_H.clock.timer.elapsed

#######################################################################
# End of script
#######################################################################
