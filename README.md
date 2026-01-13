# LESStudySetup.jl

[![CI](https://github.com/YOUR_ORG/LESStudySetup.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/YOUR_ORG/LESStudySetup.jl/actions/workflows/CI.yml)
[![Julia](https://img.shields.io/badge/Julia-1.9%2B-blue.svg)](https://julialang.org/)
[![Oceananigans](https://img.shields.io/badge/Oceananigans-0.95.7-purple.svg)](https://github.com/CliMA/Oceananigans.jl)

LESStudySetup.jl is a Julia package that orchestrates large-eddy simulations (LES) of submesoscale ocean fronts built on top of [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl). The package provides parameter management, model setup utilities, diagnostics, experiment drivers, and visualization recipes for running idealized hydrostatic or nonhydrostatic studies.

## Requirements

- **Julia**: 1.9 or later
- **Oceananigans.jl**: 0.95.7
- **Optional**: CUDA-capable GPU for accelerated simulations, MPI for distributed runs

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/YOUR_ORG/LESStudySetup.jl")
```

Or clone and develop locally:

```bash
git clone https://github.com/YOUR_ORG/LESStudySetup.jl.git
cd LESStudySetup.jl
julia --project -e 'using Pkg; Pkg.instantiate()'
```

## Quick Start

```julia
using LESStudySetup
using LESStudySetup.Oceananigans.Units

# Configure physical parameters
set_value!(;
    Δh = 100.0,      # Horizontal grid spacing (m)
    Δz = 5.0,        # Vertical grid spacing (m)
    Lx = 10kilometers,
    Ly = 10kilometers,
    Lz = 200meters,
    Q  = 40.0,       # Surface heat flux (W/m²)
    τw = 0.1         # Wind stress (N/m²)
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
│   ├── initial_conditions.jl      # Initial velocity/temperature fields
│   ├── background_field_forcing.jl
│   ├── model_setup.jl             # Model builders and kernel dispatchers
│   ├── idealized_setup.jl         # Simulation factories
│   └── Diagnostics/               # Submodule for analysis tools
├── experiments/                   # Production simulation scripts
├── test/                          # Integration/unit tests
├── visualize_results/             # Makie plotting scripts
└── job.sh, setup_perlmutter.sh    # HPC batch scripts
```

## API Reference

### Main Exports

| Function | Description |
|----------|-------------|
| `idealized_setup(arch; ...)` | Create a configured `Simulation` for idealized LES |
| `turbulence_generator_setup(arch; ...)` | Create a small-domain simulation for generating initial turbulence |
| `set_value!(; kwargs...)` | Configure global simulation parameters |
| `set!(parameters; kwargs...)` | Alternative parameter setter |
| `parameters` | Global `ProblemConstants` singleton |

### idealized_setup

```julia
simulation = idealized_setup(arch;
    stop_time = 100days,
    stop_iteration = Inf,
    hydrostatic_approximation = false,
    background_forcing = true)
```

- `arch`: Architecture (`CPU()`, `GPU()`, or `Distributed(...)`)
- `hydrostatic_approximation`: Use `HydrostaticFreeSurfaceModel` (true) or `NonhydrostaticModel` (false)
- `background_forcing`: Include eddies as background forcing

## Parameter Reference

Key parameters accessible via `parameters` singleton:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `Δh` | 1000 m | Horizontal grid spacing |
| `Δz` | 4 m | Vertical grid spacing |
| `Lx`, `Ly` | 100 km | Domain dimensions |
| `Lz` | 256 m | Domain depth |
| `f` | 1e-4 s⁻¹ | Coriolis parameter |
| `τw` | 0.1 N/m² | Surface wind stress |
| `θ` | 30° | Wind stress angle |
| `Q` | 10 W/m² | Surface heat flux (positive = cooling) |
| `α` | 2e-4 K⁻¹ | Thermal expansion coefficient |
| `N²s` | 5e-7 s⁻² | Surface stratification |
| `N²T` | 1e-4 s⁻² | Pycnocline stratification |
| `M²₀` | 5e-7 s⁻² | Frontal density gradient |
| `m₀` | 50 m | Initial mixed layer depth |
| `T₀` | 5 °C | Surface temperature |

## Diagnostics

The `Diagnostics` submodule provides tools for loading, analyzing, and saving simulation data.

### Loading Data

```julia
using LESStudySetup.Diagnostics

# Single-file time series
snapshots = load_snapshots("snapshots.jld2"; architecture = CPU())

# Distributed (MPI rank-partitioned) checkpoint
snapshot = load_distributed_checkpoint("checkpoint_prefix", iteration)

# Load subdomain from distributed data
subdomain = load_distributed_checkpoint_subdomain("prefix", iteration;
    xlims = (2000.0, 7000.0),
    ylims = (2000.0, 7000.0),
    zlims = (-80.0, -30.0))
```

### Saving/Loading Subdomains

```julia
# Save extracted subdomain
save_subdomain_snapshot("output.jld2", subdomain; iteration = 100)

# Load subdomain
loaded = load_subdomain_snapshot("output.jld2")
```

### Computed Diagnostics

| Function | Description |
|----------|-------------|
| `ζ(snapshots)` | Relative vorticity |
| `ub(snapshots)`, `vb(snapshots)`, `wb(snapshots)` | Buoyancy fluxes |
| `uw(snapshots)`, `vw(snapshots)` | Momentum fluxes |
| `KE(snapshots)` | Kinetic energy |
| `MLD(snapshots)` | Mixed layer depth |
| `BLD1D(snapshots)` | Boundary layer depth (1D) |
| `PV(snapshots)` | Potential vorticity |

## Running Experiments

### Local (Single CPU/GPU)

```julia
using LESStudySetup

set_value!(; Δh = 100.0, Δz = 5.0, Q = 40.0, τw = 0.1)
simulation = idealized_setup(GPU(); stop_time = 1days)
run!(simulation)
```

### Distributed (MPI + Multi-GPU)

```julia
using MPI
MPI.Init()
using LESStudySetup

arch = Distributed(GPU(), partition = Partition(4, 4))  # 16 GPUs
set_value!(; Δh = 5.0, Δz = 1.125)
simulation = idealized_setup(arch; stop_time = 10days)

# Attach per-rank output writers
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
# Run all tests
julia --project -e 'using Pkg; Pkg.test()'

# Run tests directly
julia --project test/runtests.jl
```

Tests include:
- Subdomain extraction and round-trip save/load
- Distributed checkpoint loading across MPI ranks

## Visualization

Analysis scripts in `visualize_results/` demonstrate post-processing workflows:

- `visualize_front.jl` – Compare simulations to analytical frontal solutions
- `visualize_nonhydro.jl` – Nonhydrostatic experiment diagnostics
- `visualize_spectra.jl` – Spectral analysis
- `visualize_fluxes.jl` – Flux diagnostics

## Authors

- Simone Silvestri
- Shirui Peng
- Abigail Bodner

## License

See [LICENSE](LICENSE) for details.
