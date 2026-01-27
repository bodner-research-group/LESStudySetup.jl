# LESStudySetup.jl

[![CI](https://github.com/bodner-research-group/LESStudySetup.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/bodner-research-group/LESStudySetup.jl/actions/workflows/CI.yml)
[![Julia](https://img.shields.io/badge/Julia-1.9%2B-blue.svg)](https://julialang.org/)
[![Oceananigans](https://img.shields.io/badge/Oceananigans-0.95.7-purple.svg)](https://github.com/CliMA/Oceananigans.jl)

LESStudySetup.jl orchestrates large-eddy simulations (LES) of submesoscale ocean fronts built on [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl). If you're familiar with Oceananigans, think of this package as a **pre-configured simulation factory** that bundles physical parameters, initial conditions, boundary conditions, and diagnostics for studying submesoscale-boundary layer interactions.

## Scientific Context

This package supports research on multiscale upper ocean dynamics—specifically the interplay between:
- **Mesoscale eddies** (~100 km, weeks–months)
- **Submesoscale fronts** (~1–10 km, hours–days)
- **Boundary layer turbulence** (~10–1000 m, minutes–hours)

The simulations use a 100 km domain with meter-scale resolution to simultaneously resolve boundary layer turbulence energized by surface forcing (wind stress, convective cooling) and submesoscale frontal dynamics modulated by a prescribed mesoscale eddy field.

**Reference**: For the scientific methodology and results, see:

> S. Peng, S. Silvestri, A. Bodner. *Submesoscale and boundary layer turbulence under mesoscale forcing in the upper ocean*. arXiv:2601.10441, 2026.
> [https://arxiv.org/abs/2601.10441](https://arxiv.org/abs/2601.10441)

## For Oceananigans Users

If you already know Oceananigans, here's how LESStudySetup maps to familiar concepts:

| Oceananigans Concept | LESStudySetup Equivalent |
|---------------------|--------------------------|
| `RectilinearGrid(...)` | Auto-constructed from `parameters.Δh`, `parameters.Lx`, etc. |
| `NonhydrostaticModel(...)` | `idealized_setup(arch; hydrostatic_approximation=false)` |
| `HydrostaticFreeSurfaceModel(...)` | `idealized_setup(arch; hydrostatic_approximation=true)` |
| `set!(model, u=..., T=...)` | Pre-configured initial conditions (`uᵢ`, `vᵢ`, `Tᵢ`) |
| `FluxBoundaryCondition(...)` | Auto-configured from `parameters.τw`, `parameters.Q` |
| `Simulation(model; Δt, ...)` | Returned by `idealized_setup()` with `TimeStepWizard` |
| `FieldTimeSeries(...)` | `Diagnostics.load_snapshots(...)` |

**Key difference**: Instead of building everything from scratch, you configure the `parameters` singleton and call `idealized_setup()`. The package handles grid construction, boundary conditions, initial conditions, and time-stepping configuration.

```julia
# Oceananigans way (verbose)
grid = RectilinearGrid(GPU(), size=(100,100,64), x=(0,1e5), y=(0,1e5), z=(-256,0))
model = NonhydrostaticModel(; grid, coriolis=FPlane(f=1e-4), ...)
set!(model, u=my_u_init, T=my_T_init)
simulation = Simulation(model; Δt=10, stop_time=86400)

# LESStudySetup way (streamlined)
set_value!(; Δh=1000.0, Δz=4.0, Lx=100e3, Ly=100e3, Lz=256.0)
simulation = idealized_setup(GPU(); stop_time=1days)
```

## Requirements

- **Julia**: 1.9 or later
- **Oceananigans.jl**: 0.95.7
- **Optional**: CUDA-capable GPU for accelerated simulations, MPI for distributed runs

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/bodner-research-group/LESStudySetup.jl")
```

Or clone and develop locally:

```bash
git clone https://github.com/bodner-research-group/LESStudySetup.jl.git
cd LESStudySetup.jl
julia --project -e 'using Pkg; Pkg.instantiate()'
```

## Quick Start

```julia
using LESStudySetup
using LESStudySetup.Oceananigans.Units

# Configure physical parameters
set_value!(;
    Δh = 100.0,       # Horizontal grid spacing (m)
    Δz = 5.0,         # Vertical grid spacing (m)
    Lx = 10kilometers,
    Ly = 10kilometers,
    Lz = 200meters,
    Q  = 40.0,        # Surface heat flux (W/m², positive = cooling)
    τw = 0.1          # Wind stress (N/m²)
)

# Create and run simulation
simulation = idealized_setup(CPU();
    stop_time = 1days,
    hydrostatic_approximation = false,
    background_forcing = true)

run!(simulation)
```

## Repository Layout

```
LESStudySetup.jl/
├── src/
│   ├── LESStudySetup.jl          # Main module (re-exports Oceananigans)
│   ├── parameters.jl              # Mutable ProblemConstants singleton
│   ├── initial_conditions.jl      # Initial velocity/temperature fields (uᵢ, vᵢ, Tᵢ)
│   ├── background_field_forcing.jl
│   ├── model_setup.jl             # Model builders and kernel dispatchers
│   ├── idealized_setup.jl         # Simulation factories
│   └── Diagnostics/               # Submodule for analysis tools
├── experiments/                   # Production simulation scripts
├── test/                          # Integration/unit tests
├── visualize_results/             # Makie plotting scripts
└── job.sh, setup_perlmutter.sh    # HPC batch scripts (Perlmutter)
```

## API Reference

### Main Exports

| Function | Description |
|----------|-------------|
| `idealized_setup(arch; ...)` | Create a configured `Simulation` for idealized LES |
| `turbulence_generator_setup(arch; ...)` | Create a small-domain simulation for generating initial turbulence |
| `set_value!(; kwargs...)` | Configure global simulation parameters |
| `parameters` | Global `ProblemConstants` singleton |

### idealized_setup

```julia
simulation = idealized_setup(arch;
    stop_time = 100days,
    stop_iteration = Inf,
    hydrostatic_approximation = false,
    background_forcing = true)
```

**Arguments**:
- `arch`: Architecture—`CPU()`, `GPU()`, or `Distributed(GPU(), partition=Partition(nx, ny))`
- `hydrostatic_approximation`: `true` → `HydrostaticFreeSurfaceModel`, `false` → `NonhydrostaticModel`
- `background_forcing`: Include mesoscale eddies as background forcing

**Returns**: A `Simulation` with `TimeStepWizard` and progress callbacks pre-configured.

## Parameter Reference

Configure via `set_value!(; param=value, ...)` before calling `idealized_setup()`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `Δh` | 1000 m | Horizontal grid spacing |
| `Δz` | 4 m | Vertical grid spacing |
| `Lx`, `Ly` | 100 km | Horizontal domain dimensions |
| `Lz` | 256 m | Domain depth |
| `f` | 1e-4 s⁻¹ | Coriolis parameter |
| `τw` | 0.1 N/m² | Surface wind stress magnitude |
| `θ` | 30° | Wind stress angle (0° = zonal) |
| `Q` | 10 W/m² | Surface heat flux (positive = cooling) |
| `α` | 2e-4 K⁻¹ | Thermal expansion coefficient |
| `N²s` | 5e-7 s⁻² | Surface stratification |
| `N²T` | 1e-4 s⁻² | Pycnocline stratification |
| `M²₀` | 5e-7 s⁻² | Frontal horizontal buoyancy gradient |
| `m₀` | 50 m | Initial mixed layer depth |
| `T₀` | 5 °C | Surface reference temperature |

**Tip**: Use `@info parameters` to print current values.

## Diagnostics

The `Diagnostics` submodule provides tools for loading, analyzing, and saving simulation output.

### Loading Data

```julia
using LESStudySetup.Diagnostics

# Single-file time series (standard Oceananigans output)
snapshots = load_snapshots("snapshots.jld2"; architecture=CPU())

# Distributed checkpoint (MPI rank-partitioned files)
snapshot = load_distributed_checkpoint("checkpoint_prefix", iteration)

# Extract subdomain from distributed data
subdomain = load_distributed_checkpoint_subdomain("prefix", iteration;
    xlims = (20e3, 80e3),
    ylims = (20e3, 80e3),
    zlims = (-100.0, 0.0))
```

### Saving/Loading Subdomains

```julia
save_subdomain_snapshot("output.jld2", subdomain; iteration=100)
loaded = load_subdomain_snapshot("output.jld2")
```

### Computed Diagnostics

| Function | Returns |
|----------|---------|
| `ζ(snapshots)` | Relative vorticity ζ = ∂v/∂x - ∂u/∂y |
| `ub`, `vb`, `wb` | Buoyancy fluxes |
| `uw`, `vw` | Vertical momentum fluxes |
| `KE` | Kinetic energy ½(u² + v² + w²) |
| `MLD` | Mixed layer depth |
| `BLD1D` | 1D boundary layer depth |
| `PV` | Ertel potential vorticity |

## Running Experiments

### Local (CPU or single GPU)

```julia
using LESStudySetup
using LESStudySetup.Oceananigans.Units

set_value!(; Δh=100.0, Δz=5.0, Q=40.0, τw=0.1)
simulation = idealized_setup(GPU(); stop_time=1days)
run!(simulation)
```

### Distributed (MPI + Multi-GPU)

```julia
using MPI
MPI.Init()
using LESStudySetup
using LESStudySetup.Oceananigans.Units

arch = Distributed(GPU(), partition=Partition(4, 4))  # 16 GPUs
set_value!(; Δh=5.0, Δz=1.125)
simulation = idealized_setup(arch; stop_time=10days)

# Per-rank output
model = simulation.model
simulation.output_writers[:snapshots] = JLD2OutputWriter(model,
    merge(model.velocities, model.tracers);
    schedule = TimeInterval(30minutes),
    filename = "snapshots_$(arch.local_rank)")

run!(simulation)
```

### HPC (Perlmutter)

```bash
source setup_perlmutter.sh
sbatch job.sh
```

## Testing

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

## Working with AI Agents

This repository includes an `AGENTS.md` file that provides AI coding assistants (like Claude, GPT, or Cursor) with project-specific context. Here's how to collaborate effectively:

### Getting Started with AI Agents

1. **Use a capable tool**: AI agents work best in environments like [Cursor](https://cursor.sh), [GitHub Copilot Chat](https://github.com/features/copilot), or API-based assistants that can read files and execute commands.

2. **Point the agent to AGENTS.md**: When starting a session, tell the agent:
   > "Read AGENTS.md for project conventions before making changes."

3. **Be specific about what you want**: Good prompts include context:
   > "Add a new diagnostic function to compute vertical buoyancy flux variance. Follow the pattern in `src/Diagnostics/pointwise_diagnostics.jl`."

### Effective Prompts for This Project

| Task | Example Prompt |
|------|----------------|
| **Run simulation** | "Set up a 10km domain with 50m resolution and run for 1 day on CPU" |
| **Add diagnostic** | "Add a function to compute Rossby number Ro = ζ/f to the Diagnostics module" |
| **Load data** | "Load the distributed checkpoint at iteration 1000 and extract a 20km×20km subdomain centered at (50km, 50km)" |
| **Debug** | "The simulation crashes with a CFL error. Check the TimeStepWizard configuration in idealized_setup.jl" |
| **Modify physics** | "Change the surface boundary condition to use a diurnally-varying heat flux" |

### What Agents Do Well Here

- **Finding patterns**: "Show me how initial conditions are defined" → Agent searches `initial_conditions.jl`
- **Explaining code**: "What does `background_forcing=true` actually do?" → Agent traces through `idealized_setup.jl` and `background_field_forcing.jl`
- **Writing boilerplate**: "Create an experiment script like `nonhydrostatic_experiment.jl` but for a smaller domain"
- **Running tests**: "Run the test suite and explain any failures"

### What to Watch For

- **Always verify physics**: Agents can write syntactically correct code that's physically wrong. Check units, signs (positive flux = cooling!), and boundary condition orientations.
- **Test changes**: Ask the agent to run `julia --project -e 'using Pkg; Pkg.test()'` after modifications.
- **Halo regions**: After any field mutation, ensure `fill_halo_regions!()` is called—agents sometimes forget this.
- **Parameter singleton**: The `parameters` object is mutable and global. Changes persist across function calls within a session.

### Example Session

```
You: I want to add mixed layer restratification by submesoscale eddies as a
     parameterization. Where should this go?

Agent: Based on the codebase structure, parameterizations that modify the model
       equations belong in `src/model_setup.jl`. Looking at the existing code,
       you'd add it as a `forcing` term in the model constructor. Here's the
       pattern used for background forcing...

You: Implement it following Fox-Kemper et al. (2008).

Agent: I'll add the MLE parameterization. First, let me check how the existing
       closures are configured... [reads files, writes code, runs tests]
```

## Citation

If you use this package in your research, please cite:

```bibtex
@article{peng2026submesoscale,
  title={Submesoscale and boundary layer turbulence under mesoscale forcing in the upper ocean},
  author={Peng, S. and Silvestri, S. and Bodner, A.},
  journal={arXiv preprint arXiv:2601.10441},
  year={2026}
}
```

## Authors

- Simone Silvestri (MIT, Politecnico di Torino)
- Shirui Peng (MIT)
- Abigail Bodner (MIT)

## License

See [LICENSE](LICENSE) for details.
