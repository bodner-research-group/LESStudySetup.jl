# AGENTS.md

Agent-specific guidance for LESStudySetup.jl — a Julia package orchestrating large-eddy simulations of submesoscale ocean fronts on top of [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl).

---

## Quick Reference

| Task | Command |
|------|---------|
| Instantiate deps | `julia --project -e 'using Pkg; Pkg.instantiate()'` |
| Run tests | `julia --project -e 'using Pkg; Pkg.test()'` |
| Run single test | `julia --project test/runtests.jl` |
| Format (if needed) | `julia --project -e 'using JuliaFormatter; format(".")'` |

**CI**: Julia 1.9 and 1.10 on `ubuntu-latest`. Tests must be deterministic and GPU-free.

---

## Repository Layout

```
LESStudySetup.jl/
├── src/
│   ├── LESStudySetup.jl      # Main module, re-exports Oceananigans
│   ├── parameters.jl          # Mutable ProblemConstants singleton
│   ├── initial_conditions.jl  # Initial velocity/temperature fields (uᵢ, vᵢ, Tᵢ)
│   ├── background_field_forcing.jl
│   ├── model_setup.jl         # Model builders, kernel dispatchers
│   ├── idealized_setup.jl     # Simulation factories (idealized_setup, turbulence_generator_setup)
│   └── Diagnostics/
│       ├── Diagnostics.jl     # Submodule entry, exports loaders/writers
│       ├── load_distributed_snapshot.jl  # MPI rank-stitching utilities
│       ├── mixed_layer.jl     # MLD calculations
│       ├── spectra.jl         # Spectral analysis
│       ├── filtering.jl       # Spatial filtering
│       ├── boundary_layer.jl  # BL diagnostics
│       └── pointwise_diagnostics.jl
├── experiments/               # Standalone simulation drivers
├── test/                      # Integration/unit tests
├── visualize_results/         # Makie plotting scripts
├── job.sh, prof.sh            # HPC batch scripts (Perlmutter)
└── setup_perlmutter.sh        # HPC environment setup
```

---

## Core Patterns

### 1. Module Organization

The main module uses explicit re-exports and includes:

```julia
module LESStudySetup

export parameters
export idealized_setup, turbulence_generator_setup
export set_value!, set!

using Reexport
@reexport using Oceananigans  # Re-exports entire Oceananigans API
using Oceananigans.Units

include("parameters.jl")
include("initial_conditions.jl")
# ... other includes
include("Diagnostics/Diagnostics.jl")

using .Diagnostics  # Uses the submodule

end
```

**When adding new files**: Include them in `src/LESStudySetup.jl` and add exports there.

### 2. Parameters (Mutable Singleton Pattern)

Global simulation constants live in `parameters::ProblemConstants`:

```julia
# Reading parameters
Δh = parameters.Δh
Lx = parameters.Lx

# Setting parameters (two equivalent methods)
set_value!(; Δh = 1000.0, Δz = 10.0)
set!(parameters; Δh = 1000.0, Δz = 10.0)

# Setting single value
set_value!(:Δh, 1000.0)
```

**Key fields**: `Δh`, `Δz`, `Lx`, `Ly`, `Lz`, `f`, `τw`, `Q`, `α`, `N²s`, `N²T`, `M²₀`, `T₀`, `m₀`

**In tests**: Always configure parameters at the start using `set_value!` before creating simulations.

### 3. Initial Condition Functions

Initial conditions use `@inline` functions that read from `parameters`:

```julia
@inline function Tᵢ(x, y, z)
    Le = parameters.Le
    R = 25e3
    # ... compute temperature
end
```

**Convention**: 
- `uᵢ`, `vᵢ`, `wᵢ` — velocity components
- `Tᵢ` — temperature
- `vᶠ`, `Tᶠ` — pure frontal profiles
- `vᵢᶠ` — combined eddy + frontal

### 4. Model Setup Dispatch

Model type selection uses `Val` dispatch:

```julia
model_type(::Val{true})  = HydrostaticFreeSurfaceModel
model_type(::Val{false}) = NonhydrostaticModel

# Usage
ModelType = model_type(Val(hydrostatic_approximation))
model = ModelType(; grid, coriolis, buoyancy, boundary_conditions, settings...)
```

### 5. Kernel Patterns (KernelAbstractions)

GPU-compatible kernels use `@kernel` with `launch!`:

```julia
using KernelAbstractions: @index, @kernel
using Oceananigans.Utils: launch!

@kernel function _compute_v_from_continuity!(v, grid, u)
    i, k = @index(Global, NTuple)
    @inbounds v[i, 1, k] = v[i, 0, k] - Δyᶜᶜᶜ(i, 0, k, grid) * ∂xᶜᶜᶜ(i, 0, k, grid, u)
    for j in 2:size(grid, 2)
        @inbounds v[i, j, k] = v[i, j-1, k] - Δyᶜᶜᶜ(i, j-1, k, grid) * ∂xᶜᶜᶜ(i, j-1, k, grid, u)
    end
end

function compute_v_from_continuity!(v_background, arch, grid, u_background)
    launch!(architecture(grid), grid, :xz, _compute_v_from_continuity!, v_background, grid, u_background)
    fill_halo_regions!(v_background)
    return nothing
end
```

**Rules**:
- Use `@inbounds` inside tight loops
- Always call `fill_halo_regions!` after mutating fields
- Use `architecture(grid)` for CPU/GPU dispatch

### 6. Simulation Factory Pattern

`idealized_setup` is the main entry point:

```julia
simulation = idealized_setup(arch;
    stop_time = 10days,
    hydrostatic_approximation = false,
    background_forcing = true)
```

Returns a configured `Simulation` with:
- Model (Hydrostatic or Nonhydrostatic)
- `TimeStepWizard` callback for CFL-based Δt
- Progress callback (every 100 iterations)

### 7. Field Operations

```julia
# Creating fields
u = XFaceField(grid)
v = YFaceField(grid)
T = CenterField(grid)

# Setting initial conditions
set!(model, u = uᵢ, v = vᵢ, T = Tᵢ)

# Accessing interior data (no halos)
ui = interior(u)
Ti = interior(T)

# After mutation, fill halos
fill_halo_regions!(u)
```

---

## Experiment Script Pattern

All experiments follow this structure:

```julia
using LESStudySetup
using LESStudySetup.Oceananigans.Units

# 1. Configure architecture
arch = Distributed(GPU(), partition = Partition(32, 32))

# 2. Set parameters
set_value!(; Δh = 5.0, Δz = 1.125, Q = 40.0, τw = 0.1)

# 3. Create simulation
simulation = idealized_setup(arch; stop_time = 10days, background_forcing = true)

# 4. Attach outputs
model = simulation.model
output_fields = merge(model.velocities, model.tracers)
simulation.output_writers[:snapshots] = JLD2OutputWriter(model, output_fields;
    schedule = TimeInterval(30minutes),
    filename = "snapshots_$(arch.local_rank)")

# 5. Run
run!(simulation)
```

---

## Test Patterns

Tests use clear setup/verify/cleanup phases:

```julia
@testset "LESStudySetup.jl" begin
    @testset "Specific feature" begin
        # 1. Configure parameters
        set_value!(; Δh = 1000.0, Δz = 10.0, Lx = 10000.0)

        # 2. Create simulation/model
        simulation = idealized_setup(CPU(); stop_time = 1)
        model = simulation.model

        # 3. Build test data
        full_snapshot = Dict{Symbol, Any}()
        full_snapshot[:u] = model.velocities.u
        # ...

        # 4. Test operations
        subdomain = extract_subdomain(full_snapshot; xlims = ...)
        @test subdomain[:grid].Nx == 5

        # 5. Verify data matches
        @test interior(subdomain[:T]) ≈ interior(model.tracers.T)[3:7, 3:7, 3:7]

        # 6. Save/load round-trip (if applicable)
        test_file = tempname() * ".jld2"
        save_subdomain_snapshot(test_file, subdomain; iteration = 0)
        loaded = load_subdomain_snapshot(test_file)
        @test interior(loaded[:T]) ≈ interior(subdomain[:T])

        # 7. Cleanup
        rm(test_file, force = true)
    end
end
```

**Adding new tests**: Wire them through `test/runtests.jl`.

---

## Diagnostics Submodule

### Loading Data

```julia
using LESStudySetup.Diagnostics

# Single-file time series
snapshots = load_snapshots("snapshots.jld2"; architecture = CPU())

# Distributed (MPI rank-partitioned) data
snapshot = load_distributed_checkpoint(filename, iteration; architecture = CPU())
snapshot = load_distributed_snapshot(filename, iteration; level = 10)

# Subdomain extraction
subdomain = load_distributed_checkpoint_subdomain(prefix, iteration;
    xlims = (2000.0, 7000.0),
    ylims = (2000.0, 7000.0),
    zlims = (-80.0, -30.0),
    getMLD = 3,
    getEw = true)
```

### Saving/Loading Subdomains

```julia
# Save
save_subdomain_snapshot("output.jld2", subdomain;
    iteration = 100,
    xlims = (x1, x2),
    ylims = (y1, y2),
    zlims = (z1, z2))

# Load
loaded = load_subdomain_snapshot("output.jld2")
```

### Computed Diagnostics

```julia
# Available functions (all operate on snapshot dictionaries)
ζ(snapshots)      # Relative vorticity
ub(snapshots)     # Buoyancy flux (u*b)
vb(snapshots)     # Buoyancy flux (v*b)
wb(snapshots)     # Vertical buoyancy flux
uw(snapshots)     # Momentum flux
vw(snapshots)     # Momentum flux
KE(snapshots)     # Kinetic energy
MLD(snapshots)    # Mixed layer depth
BLD1D(snapshots)  # Boundary layer depth (1D)
PV(snapshots)     # Potential vorticity
```

---

## Naming Conventions

| Entity | Convention | Examples |
|--------|------------|----------|
| Functions | `snake_case` | `idealized_setup`, `compute_v_from_continuity!` |
| Types/Structs | `PascalCase` | `ProblemConstants` |
| Mutating functions | `!` suffix | `set_value!`, `set!`, `fill_halo_regions!` |
| Grid dimensions | `N` prefix | `Nx`, `Ny`, `Nz` |
| Domain lengths | `L` prefix | `Lx`, `Ly`, `Lz` |
| Spacing | `Δ` prefix | `Δh`, `Δz` |
| Initial conditions | `ᵢ` subscript | `uᵢ`, `vᵢ`, `Tᵢ` |
| Physical constants | Unicode | `ρ₀`, `α`, `σ²`, `τw` |

---

## Critical Do's and Don'ts

### DO

- Use `set_value!` to configure parameters before creating simulations
- Call `fill_halo_regions!` after mutating fields
- Use `interior(field)` for array access in tests
- Use `@inbounds` inside performance-critical kernel loops
- Keep tests deterministic and small-grid
- Follow existing patterns in `src/` for new features
- Use `architecture(grid)` for CPU/GPU dispatch

### DON'T

- Don't create new global constants — use the `parameters` singleton
- Don't access field data without `interior()` in tests
- Don't forget halo filling after field mutations
- Don't add GPU-dependent tests to CI
- Don't introduce new dependencies without explicit request
- Don't run long simulations in tests (use small grids, short times)
- Don't delete user data or result directories
- Don't refactor unrelated code during bug fixes

---

## Boundary Conditions

```julia
# Surface flux boundary conditions
u_top = FluxBoundaryCondition(τw * cosd(θ) / ρ₀)
v_top = FluxBoundaryCondition(τw * sind(θ) / ρ₀)
T_top = FluxBoundaryCondition(Q / ρ₀ / cₚ)  # Positive = cooling

u_bcs = FieldBoundaryConditions(top = u_top)
v_bcs = FieldBoundaryConditions(top = v_top)
T_bcs = FieldBoundaryConditions(top = T_top)

boundary_conditions = (u = u_bcs, v = v_bcs, T = T_bcs)
```

---

## Dependencies

Key dependencies from `Project.toml`:

| Package | Purpose |
|---------|---------|
| `Oceananigans` (0.95.7) | Core ocean modeling |
| `KernelAbstractions` | GPU/CPU kernel dispatch |
| `JLD2` | Data I/O |
| `MPI` | Distributed computing |
| `CairoMakie` | Visualization |
| `CUDA` | GPU support |
| `FFTW` | Spectral analysis |

---

## HPC (Perlmutter)

```bash
# Environment setup
source setup_perlmutter.sh

# Production run
sbatch job.sh    # Runs experiments/nonhydrostatic_experiment.jl

# Benchmark/profiling
sbatch prof.sh   # Runs experiments/benchmark_nonhydrostatic.jl
```

---

## Common Gotchas

1. **Halo regions**: Always `fill_halo_regions!` after `set!` or direct interior manipulation
2. **Positive flux = negative direction**: In Oceananigans, positive top flux = cooling/downward
3. **MPI barriers**: Use `MPI.Barrier(MPI.COMM_WORLD)` when synchronizing distributed kernels
4. **Index conventions**: Oceananigans uses 1-based indexing; halos extend into negative indices
5. **Unicode in Julia**: Unicode variable names (`Δh`, `ρ₀`) are fully supported and preferred
6. **TimeStepWizard**: Hydrostatic models need smaller CFL (~0.25) than nonhydrostatic (~0.75)
