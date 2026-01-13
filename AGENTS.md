# AGENTS.md

This file provides repo-specific guidance for agentic coding tools.
Scope: entire repository unless overridden by a deeper AGENTS.md.

## Cursor/Copilot rules
- No `.cursor/rules/`, `.cursorrules`, or `.github/copilot-instructions.md` present in this repo.

## Build / Run / Test Commands
### Environment setup
- Instantiate dependencies:
  - `julia --project -e 'using Pkg; Pkg.instantiate()'`
- CI uses Julia build step (equivalent to `Pkg.build()`):
  - `julia --project -e 'using Pkg; Pkg.build()'`

### Tests
- Run full test suite:
  - `julia --project -e 'using Pkg; Pkg.test()'`
  - `julia --project test/runtests.jl`
- Run a single test file directly:
  - `julia --project test/loadsave_nhyles.jl`
  - `julia --project test/<other_test>.jl`

### Experiments (local)
- Hydrostatic experiment:
  - `julia --project experiments/hydrostatic_experiment.jl`
- Nonhydrostatic experiment:
  - `julia --project experiments/nonhydrostatic_experiment.jl`
- Benchmark run:
  - `julia --project experiments/benchmark_nonhydrostatic.jl`

### HPC batch scripts (Perlmutter)
- Environment setup:
  - `source setup_perlmutter.sh`
- Production run:
  - `sbatch job.sh` (executes `experiments/nonhydrostatic_experiment.jl`)
- Benchmark run:
  - `sbatch prof.sh` (executes `experiments/benchmark_nonhydrostatic.jl`)

### Lint / Format
- No linting or formatting configuration in repo.
- If formatting is required, use JuliaFormatter explicitly and do not add configs unless requested:
  - `julia --project -e 'using JuliaFormatter; format(".")'`

## Code Style Guidelines
### General Julia style
- Functions use `snake_case` (e.g., `idealized_setup`, `set_value!`).
- Types/structs use `PascalCase` (e.g., `ProblemConstants`).
- Use `!` suffix for mutating functions (e.g., `compute_v_from_continuity!`).
- Use concise docstrings with triple quotes (`"""`) for public functions.

### Imports and module layout
- Prefer explicit `using`/`import` statements at top of files.
- Keep module organization consistent with `src/LESStudySetup.jl` re-export pattern.
- For submodules (e.g., `src/Diagnostics/`), keep related helpers colocated.

### Parameters and configuration
- Use the singleton `parameters` (`ProblemConstants`) from `src/parameters.jl`.
- Update settings via `set_value!` or `set!(parameters; ...)` rather than new globals.
- Use `@kwdef` for configuration structs and keyword-heavy APIs.

### Naming and scientific notation
- Unicode variable names are common and acceptable (e.g., `Δh`, `ρ₀`, `σ²`).
- Grid dimensions use `Nx`, `Ny`, `Nz`, and lengths use `Lx`, `Ly`, `Lz`.
- Keep physical parameters named after Oceananigans conventions.

### Error handling and logging
- Prefer Julia logging macros (`@info`, `@warn`, `@error`) over silent failures.
- Let exceptions propagate unless there is a clear recovery path.
- Avoid empty `catch` blocks or swallowing errors.

### Performance and kernels
- Use `KernelAbstractions.@kernel` with `launch!` for compute kernels.
- Use `architecture(grid)` to decide CPU vs GPU execution paths.
- Apply `@inbounds` inside tight loops when safe.
- Call `fill_halo_regions!` after mutating fields that require halos.

### Oceananigans patterns
- Model configuration uses dispatch on `Val` (e.g., `model_type(Val(...))`).
- Boundary conditions use `FluxBoundaryCondition` and `FieldBoundaryConditions`.
- Fields and grids rely on Oceananigans types (`RectilinearGrid`, `Field`).

### Tests
- Tests use `Test.jl` and live under `test/`.
- Use `@testset` and keep tests deterministic.
- For integration tests, follow `test/loadsave_nhyles.jl` as a template.

### Scripts
- Experiments live in `experiments/` and should be runnable as standalone scripts.
- Visualization scripts live in `visualize_results/` and use Makie.

## Notes for agents
- Keep changes minimal and consistent with existing scientific style.
- Avoid introducing new dependencies unless explicitly requested.
- When modifying diagnostics, ensure GPU/CPU compatibility is preserved.
- If you add new tests, wire them through `test/runtests.jl`.

## Repository layout
- `src/` is the package code.
- `src/LESStudySetup.jl` re-exports Oceananigans and includes submodules.
- `src/Diagnostics/` contains analysis utilities and I/O helpers.
- `experiments/` are runnable scripts with side effects.
- `visualize_results/` scripts may assume local data files.
- `test/` contains regression/integration tests.

## Diagnostics guidelines
- Keep diagnostics pure where possible; return fields or `NamedTuple`s.
- When reconstructing data, respect grid topology and halo regions.
- Use `load_distributed_snapshot` helpers for MPI outputs.
- Favor `Field` objects over raw arrays for compatibility.

## Distributed / MPI
- Use `Oceananigans.DistributedComputations` utilities for ranks.
- Synchronize distributed kernels with `MPI.Barrier` when required.
- Avoid writing MPI-specific code into shared kernels.

## Field utilities
- Use `interior(field)` for array access in tests.
- Use `location(field)` when matching grid staggering.
- Use `XFaceField`, `YFaceField`, `CenterField` as appropriate.
- Initialize fields with `set!` and then fill halos.

## Data I/O
- Snapshot data uses JLD2 in tests and diagnostics.
- Use `save_subdomain_snapshot` and `load_subdomain_snapshot` patterns.
- Keep file paths explicit; avoid writing into repo without user request.

## CI expectations
- CI currently runs on `ubuntu-latest` with Julia 1.9/1.10.
- Keep tests deterministic and avoid reliance on GPUs in CI.

## Safety
- Do not delete user data or large result directories.
- Avoid long-running simulations in tests.
- Prefer small grids for unit tests.

## Editing guidelines
- Keep edits minimal and scoped.
- Avoid refactors unrelated to the requested change.
- Follow existing naming and layout patterns.

## Simulation setup tips
- Use `TimeStepWizard` for CFL-based timesteps.
- Use `Simulation(model; Δt, stop_time)` pattern.
- Register callbacks via `simulation.callbacks`.
