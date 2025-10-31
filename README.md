# LESStudySetup

LESStudySetup is a Julia package that orchestrates large-eddy simulations (LES) of submesoscale ocean fronts on top of the [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) modeling framework. The package bundles parameter definitions, model setup utilities, diagnostics, experiment scripts, and visualization recipes that streamline running idealized hydrostatic or nonhydrostatic studies and post-processing their output.

## Repository layout

The repository is organized around a Julia package located under `src/`, supported by experiment drivers, diagnostics, and plotting scripts:

| Path | Purpose |
| --- | --- |
| `Project.toml`, `Manifest.toml` | Julia environment that pins Oceananigans and plotting dependencies. |
| `src/` | Package entry point (`LESStudySetup.jl`), physical parameters (`parameters.jl`), initial and background states, model builders, and diagnostics. |
| `src/Diagnostics/` | A submodule with utilities for loading distributed outputs, computing mixed-layer metrics, spectra, and exporting pointwise diagnostics. |
| `experiments/` | Standalone drivers for hydrostatic, nonhydrostatic, benchmark, and wind-mixing experiments that build models and run simulations. |
| `test/` | Regression and integration checks that exercise diagnostic loaders and snapshot saving. |
| `visualize_results/` | Analysis notebooks written as Julia scripts that load saved data products and produce publication-style figures. |
| `print_results.jl` | Utility helpers for summarizing model output. |
| `job.sh`, `setup_perlmutter.sh`, `prof.sh` | Batch scripts and profiling helpers for HPC environments. |

## Core package components

* `LESStudySetup.jl` ties everything together: it re-exports Oceananigans, exposes the mutable [`parameters`](src/parameters.jl) singleton for global constants, and includes initial-condition generators and model setup helpers. 【F:src/LESStudySetup.jl†L1-L19】【F:src/parameters.jl†L1-L97】
* `model_setup.jl` builds either hydrostatic or nonhydrostatic models, configures advection and closures, and provides utilities to impose background flows that satisfy continuity. 【F:src/model_setup.jl†L1-L94】
* `Diagnostics/Diagnostics.jl` exposes loaders for time-series data, writers for pointwise diagnostics, and includes submodules for mixed-layer depth calculations, spectra, boundary layer diagnostics, and filtering. 【F:src/Diagnostics/Diagnostics.jl†L1-L98】

Together, these files let you define a simulation by tweaking `parameters`, constructing a grid and `NonhydrostaticModel`/`HydrostaticFreeSurfaceModel`, then running one of the experiment scripts.

## Detailed file guide

### `src/Diagnostics/load_distributed_snapshot.jl`

This file contains the heavy-lifting utilities for reconstructing model snapshots that were written in a distributed, rank-per-tile format:

* `load_distributed_checkpoint` and `load_distributed_snapshot` stitch together velocity and temperature fields from rank-specific JLD2 files into global `Oceananigans.Field` objects. They infer the domain decomposition (`Px × Py` tiles), rebuild an appropriate `RectilinearGrid`, and populate halo regions before returning a dictionary of fields. 【F:src/Diagnostics/load_distributed_snapshot.jl†L1-L121】
* `load_distributed_checkpoint_subdomain` extends that logic to extract an arbitrary subdomain (possibly wrapping around periodic boundaries), optionally restricting to vertical levels and computing mixed-layer diagnostics such as `MLD`, `MLD2/3`, and depth-averaged vertical kinetic energy (`Ew`). It normalizes requested coordinates, determines which ranks overlap the desired region, copies the relevant interiors, and fills halos in the output fields. 【F:src/Diagnostics/load_distributed_snapshot.jl†L123-L332】【F:src/Diagnostics/load_distributed_snapshot.jl†L333-L463】
* `load_subdomain_snapshot` reloads previously saved subdomain data (produced by the test script below) by reading the stored grid, metadata, and field data, then reconstructing typed `Field` objects and filling halos. It supports selecting a single vertical level by creating a thin grid when needed. 【F:src/Diagnostics/load_distributed_snapshot.jl†L465-L569】

Understanding these loaders is key when working with large Oceananigans simulations that ran on many MPI ranks: you rarely read the raw arrays directly; instead, you reconstruct fields through these helpers to preserve geometry, halo layout, and architecture compatibility.

### `test/loadsave_nhyles.jl`

This integration test demonstrates how to call `load_distributed_checkpoint_subdomain` on a production-scale nonhydrostatic run and save its output for later analysis:

1. It selects a target subdomain in physical coordinates (here a 25 km wide along-front band and specific vertical levels) and loops over representative iterations. 【F:test/loadsave_nhyles.jl†L1-L31】
2. For each iteration it requests the subdomain, enabling vertical kinetic energy (`getEw = true`) and multiple mixed-layer depth thresholds (`getMLD = 3`). 【F:test/loadsave_nhyles.jl†L33-L47】
3. It then writes the resulting fields, along with the grid and metadata, to a compact JLD2 file whose structure matches what `load_subdomain_snapshot` expects. 【F:test/loadsave_nhyles.jl†L49-L78】

Use this script as a template for carving out manageable chunks of a large simulation so you can analyze them locally or share them with collaborators without shipping the full 3-D dataset.

### `visualize_results/visualize_front.jl`

This plotting script reproduces analytical frontal solutions and compares them to simulation data:

* It imports Oceananigans diagnostics, CUDA transfer helpers, and several Makie plotting utilities, then defines constants (Rossby, Burger, Froude numbers) and domain grids used for the analytic solution. 【F:visualize_results/visualize_front.jl†L1-L73】
* `solve_X_vec` and `fields_at_time` evaluate semi-analytical solutions for buoyancy, velocity, and streamfunction fields, returning them on the chosen `(x, z)` grid for a given time. 【F:visualize_results/visualize_front.jl†L75-L149】
* Utility functions such as `coarsen_binned_vectorized` and `resample_contour_periodic` help reduce high-resolution curves and handle periodic contours before plotting. 【F:visualize_results/visualize_front.jl†L151-L226】
* The remainder of the file (not shown above) loads simulation output, coarsens or interpolates it, and produces Makie figures saved under `results/`. Working through the script is a good way to learn how LESStudySetup uses analytical baselines to interpret numerical runs.

### `visualize_results/visualize_nonhydro.jl`

This complementary script focuses on nonhydrostatic experiment output and more elaborate diagnostics:

* It gathers a broad suite of tools—mixed-layer depth calculators, coarse-grained flux diagnostics, and MathTeX rendering—before defining helpers to enumerate saved iterations across raw rank output or precomputed subdomains. 【F:visualize_results/visualize_nonhydro.jl†L1-L34】
* The script then loads hydrostatic reference snapshots, computes mixed-layer depths, relative vorticity, and strain tensors, and prepares shifted/periodic views of the domain for centered plotting. 【F:visualize_results/visualize_nonhydro.jl†L36-L84】
* Large commented sections outline figure layouts for temperature, velocity, and diagnostic fields, showcasing how to combine Oceananigans field data with Makie’s declarative plotting API. These blocks serve as worked examples you can adapt when designing new diagnostics for the project.

## Working with the package

1. **Instantiate the environment**
   ```julia
   julia --project
   using Pkg
   Pkg.instantiate()
   ```
2. **Explore experiments** – Run scripts in `experiments/` to launch canonical setups. They all import `LESStudySetup`, mutate `parameters` as needed, call the model constructors, and advance a simulation.
3. **Load existing data** – Use `Diagnostics.load_snapshots` for single-file time series or the distributed loaders described above when working with rank-partitioned outputs.
4. **Visualize results** – The Makie scripts in `visualize_results/` double as tutorials for post-processing. Start by adapting the smaller `visualize_front.jl` before diving into the more complex nonhydrostatic workflow.

## Next steps for newcomers

* **Understand parameterization** – Read through `parameters.jl` to learn what physical constants are available and how to adjust them via `set_value!`. 【F:src/parameters.jl†L1-L120】
* **Customize diagnostics** – Explore the rest of `src/Diagnostics/` to see how mixed-layer depth or spectral diagnostics are computed, then extend them to suit your analysis.
* **Reproduce figures** – Run the visualization scripts on saved subdomain snapshots; this teaches you how the repository expects data to be organized and highlights common plotting patterns.
* **Add tests** – Mirror `test/loadsave_nhyles.jl` when adding new loaders or diagnostics to ensure they can read real-world outputs and save portable subsets for collaboration.

With this structure in mind, you can navigate the codebase, hook into Oceananigans workflows, and extend LESStudySetup with new experiments or diagnostics tailored to your LES research.