# ===============================================================================
#                    SUBDOMAIN EXTRACTION AND COARSE-GRAINING WORKFLOW
# ===============================================================================
#
# This script demonstrates a complete workflow for:
#   1. Extracting checkpoint metadata (simulation time, iteration)
#   2. Dividing a 100km x 100km domain into 10km x 10km tiles with halo bands
#   3. Saving each tile to a separate file
#   4. Loading w and T fields from a specified tile
#   5. Computing coarse-grained fields and residuals
#   6. Performing quadrant analysis of w' and b'
#
# The halo band ensures boundary artifacts from filtering don't affect the
# core region. After filtering, only the core region is retained.
#
# ===============================================================================

using Oceananigans
using Oceananigans: location
using JLD2
using LESStudySetup
using LESStudySetup.Diagnostics
using LESStudySetup.Diagnostics: load_distributed_checkpoint_subdomain
using LESStudySetup.Diagnostics: load_checkpoint_clock
using LESStudySetup.Diagnostics: compute_subdomain_tiles
using LESStudySetup.Diagnostics: save_subdomain_with_halo
using LESStudySetup.Diagnostics: coarse_graining!
using StatsBase: fit, Histogram

# ===============================================================================
# SECTION 1: USER CONFIGURATION
# ===============================================================================
#
# Modify these parameters to customize the workflow for your data.
#

# --- Data Paths ---
const CHECKPOINT_DIR = "/orcd/data/abodner/002/shared_datasets/nhyles_output/"
const CHECKPOINT_PREFIX = CHECKPOINT_DIR * "iteration16x/nonhydrostatic_checkpoint_"
const OUTPUT_DIR = CHECKPOINT_DIR * "subdomains/"

# --- Checkpoint Selection ---
const ITERATION = 164410

# --- Domain Tiling Parameters ---
const DOMAIN_LX = 100e3           # Full domain x-extent (m)
const DOMAIN_LY = 100e3           # Full domain y-extent (m)  
const TILE_SIZE = 10e3            # Core tile size (m): 10km x 10km
const Z_LIMITS = (-81.0, 0.0)     # Vertical extent (m): 72 cells at dz=1.125m

# --- Coarse-Graining Parameters ---
const KERNEL = :gaussian          # Filter kernel: :gaussian, :tophat, or :lanczos
const CUTOFF = 300.0              # Filter cutoff scale (m) used when saving tiles
const BORDER = :reflect           # Boundary handling: :reflect, :circular

# Halo width = 2x cutoff for Gaussian (captures >95% of kernel weight)
const HALO_WIDTH = 2 * CUTOFF     # 600m halo on each side

# Filter cutoff for coarse-graining (can differ from CUTOFF used for tile saving)
# If FILTER_CUTOFF > CUTOFF, tile will be reloaded from checkpoint with larger halo
const FILTER_CUTOFF = 300.0       # Filter cutoff scale (m) for coarse-graining

# Whether to save reloaded tiles (when FILTER_CUTOFF > CUTOFF requires reload)
const SAVE_RELOADED_TILE = false  # Set true to save tile with larger halo

# --- Physical Constants ---
const alpha = 2e-4                # Thermal expansion coefficient (1/K)
const g = 9.81                    # Gravitational acceleration (m/s^2)

# --- Quadrant Analysis Parameters ---
const N_BINS = 50                 # Number of histogram bins

# --- Processing Options ---
const SAVE_ALL_TILES = true       # Set false to skip tile extraction step
const TARGET_TILE = 1             # Which tile to process for coarse-graining

# ===============================================================================
# SECTION 2: EXTRACT CHECKPOINT METADATA
# ===============================================================================

println("\n" * "="^70)
println("STEP 1: Loading Checkpoint Metadata")
println("="^70)

clock_info = load_checkpoint_clock(CHECKPOINT_PREFIX, ITERATION)

println("+---------------------------------------------------------------------+")
println("| Checkpoint Information                                              |")
println("+---------------------------------------------------------------------+")
println("|  Iteration:        $(lpad(clock_info.iteration, 10))                              |")
println("|  Simulation time:  $(lpad(round(clock_info.time_days, digits=3), 10)) days                        |")
println("|                    $(lpad(round(clock_info.time, digits=1), 10)) seconds                     |")
println("+---------------------------------------------------------------------+")

# ===============================================================================
# SECTION 3: DIVIDE DOMAIN INTO TILES
# ===============================================================================

println("\n" * "="^70)
println("STEP 2: Computing Domain Tiles")
println("="^70)

tiles = compute_subdomain_tiles(;
    Lx = DOMAIN_LX,
    Ly = DOMAIN_LY,
    tile_size = TILE_SIZE,
    halo_width = HALO_WIDTH
)

n_tiles = length(tiles)
println("\nDomain configuration:")
println("  * Full domain: $(DOMAIN_LX/1e3) km x $(DOMAIN_LY/1e3) km")
println("  * Tile grid: $(round(Int, sqrt(n_tiles))) x $(round(Int, sqrt(n_tiles))) = $n_tiles tiles")
println("  * Core tile: $(TILE_SIZE/1e3) km x $(TILE_SIZE/1e3) km")
println("  * With halo: $((TILE_SIZE + 2*HALO_WIDTH)/1e3) km x $((TILE_SIZE + 2*HALO_WIDTH)/1e3) km")
println("  * Vertical: $(Z_LIMITS) m")

# ===============================================================================
# SECTION 4: SAVE TILES (OPTIONAL)
# ===============================================================================

if SAVE_ALL_TILES
    println("\n" * "="^70)
    println("STEP 3: Extracting and Saving Tiles")
    println("="^70)
    
    mkpath(OUTPUT_DIR)
    println("Output directory: $OUTPUT_DIR\n")
    
    for (i, tile) in enumerate(tiles)
        output_file = OUTPUT_DIR * "subdomain$(tile.tile_id)_iter$(ITERATION).jld2"
        
        # Skip if output file already exists
        if isfile(output_file)
            println("  Tile $(lpad(tile.tile_id, 3))/$(n_tiles): Already exists, skipping")
            continue
        end
        
        # Progress indicator
        print("  Tile $(lpad(tile.tile_id, 3))/$(n_tiles): ")
        print("x=[$(Int(tile.core_xlims[1]/1e3)),$(Int(tile.core_xlims[2]/1e3))]km, ")
        print("y=[$(Int(tile.core_ylims[1]/1e3)),$(Int(tile.core_ylims[2]/1e3))]km ... ")
        
        # Load subdomain with halo (using periodic boundary wrapping)
        snapshot = load_distributed_checkpoint_subdomain(CHECKPOINT_PREFIX, ITERATION;
            xlims = tile.full_xlims,
            ylims = tile.full_ylims,
            zlims = Z_LIMITS,
            getEw = false,
            getMLD = 0
        )
        
        # Save with halo metadata
        save_subdomain_with_halo(output_file, snapshot;
            core_xlims = tile.core_xlims,
            core_ylims = tile.core_ylims,
            halo_width = HALO_WIDTH,
            zlims = Z_LIMITS,
            iteration = ITERATION,
            clock_time = clock_info.time,
            clock_time_days = clock_info.time_days
        )
        
        println("Done")
    end
    
    println("\nAll tiles saved successfully!")
else
    println("\n[Skipping tile extraction - SAVE_ALL_TILES = false]")
end

# ===============================================================================
# SECTION 5: LOAD TILE AND PERFORM COARSE-GRAINING
# ===============================================================================

println("\n" * "="^70)
println("STEP 4: Coarse-Graining Tile $TARGET_TILE")
println("="^70)

input_file = OUTPUT_DIR * "subdomain$(TARGET_TILE)_iter$(ITERATION).jld2"
println("Loading: $input_file")

snapshot = load_subdomain_snapshot(input_file; variables = ("w", "T"))

# Warn if using a larger-than-default filter cutoff
if FILTER_CUTOFF > 300.0
    println("⚠️  WARNING: Filter cutoff ($FILTER_CUTOFF m) is larger than default (300 m).")
end

# Determine if reload is needed due to larger filter cutoff
required_halo = 2 * FILTER_CUTOFF
active_halo_width = HALO_WIDTH  # Default: use saved halo

if required_halo > HALO_WIDTH
    @warn "FILTER_CUTOFF=$FILTER_CUTOFF requires halo=$(required_halo)m, " *
          "but tile was saved with halo=$(HALO_WIDTH)m."
    
    # Get tile info for reload
    tile = tiles[TARGET_TILE]
    
    # Check if a previously saved tile with larger halo exists
    reloaded_file = OUTPUT_DIR * "subdomain$(TARGET_TILE)_iter$(ITERATION)_halo$(Int(required_halo)).jld2"
    
    if isfile(reloaded_file)
        # Load from existing file with larger halo
        println("Loading existing tile with larger halo: $reloaded_file")
        snapshot = load_subdomain_snapshot(reloaded_file; variables = ("w", "T"))
    else
        # Reload from checkpoint
        println("Reloading tile $(TARGET_TILE) from checkpoint with halo=$(required_halo)m...")
        
        # Compute expanded limits with larger halo
        expanded_xlims = (tile.core_xlims[1] - required_halo, tile.core_xlims[2] + required_halo)
        expanded_ylims = (tile.core_ylims[1] - required_halo, tile.core_ylims[2] + required_halo)
        
        snapshot = load_distributed_checkpoint_subdomain(CHECKPOINT_PREFIX, ITERATION;
            xlims = expanded_xlims,
            ylims = expanded_ylims,
            zlims = Z_LIMITS,
            getEw = false,
            getMLD = 0
        )
        
        # Optionally save the reloaded tile with larger halo
        if SAVE_RELOADED_TILE
            println("Saving reloaded tile to: $reloaded_file")
            save_subdomain_with_halo(reloaded_file, snapshot;
                core_xlims = tile.core_xlims,
                core_ylims = tile.core_ylims,
                halo_width = required_halo,
                zlims = Z_LIMITS,
                iteration = ITERATION,
                clock_time = clock_info.time,
                clock_time_days = clock_info.time_days
            )
        end
    end
    
    # Update active halo for cropping
    active_halo_width = required_halo
end

# Display loaded metadata
grid = snapshot[:grid]
println("\nLoaded subdomain:")
println("  * Grid size: $(grid.Nx) x $(grid.Ny) x $(grid.Nz) cells")
println("  * Resolution: dx=$(grid.Δxᶜᵃᵃ)m, dz=$(grid.Lz / grid.Nz)m")

if haskey(snapshot, :clock_time_days)
    println("  * Simulation time: $(round(snapshot[:clock_time_days], digits=3)) days")
end
if haskey(snapshot, :core_xlims)
    println("  * Core region (valid after filtering):")
    println("      x: $(snapshot[:core_xlims]) m")
    println("      y: $(snapshot[:core_ylims]) m")
end
if haskey(snapshot, :halo_width)
    println("  * Halo width: $(snapshot[:halo_width]) m")
end

# -----------------------------------------------------------------------------
# Compute buoyancy from temperature
# -----------------------------------------------------------------------------

w = snapshot[:w]
T = snapshot[:T]

# Buoyancy: b = alpha * g * T
b = compute!(Field(alpha * g * T))

println("\nInput fields:")
println("  * w range: $(extrema(interior(w)))")
println("  * T range: $(extrema(interior(T)))")
println("  * b range: $(extrema(interior(b)))")

# -----------------------------------------------------------------------------
# Apply coarse-graining filter
# -----------------------------------------------------------------------------

println("\nApplying coarse-graining filter:")
println("  * Kernel: $KERNEL")
println("  * Cutoff: $CUTOFF m")
println("  * Border: $BORDER")

# Allocate output fields for filtered quantities
w_bar = ZFaceField(grid, Float32)
b_bar = CenterField(grid, Float32)

# Apply filter
t_start = time()
coarse_graining!(w, w_bar; kernel=KERNEL, cutoff=CUTOFF, border=BORDER)
coarse_graining!(b, b_bar; kernel=KERNEL, cutoff=CUTOFF, border=BORDER)
t_filter = time() - t_start

println("  * Filtering completed in $(round(t_filter, digits=2)) seconds")

# -----------------------------------------------------------------------------
# Compute residuals (fine-scale fluctuations)
# -----------------------------------------------------------------------------

# w' = w - w_bar (vertical velocity fluctuation)
# b' = b - b_bar (buoyancy fluctuation)
wp_full = interior(w) .- interior(w_bar)
bp_full = interior(b) .- interior(b_bar)

# ===============================================================================
# SECTION 6: CROP TO CORE REGION
# ===============================================================================

println("\n" * "="^70)
println("STEP 5: Extracting Core Region (Discarding Halo)")
println("="^70)

# Compute indices for the valid core region
# Use active_halo_width (may differ from HALO_WIDTH if tile was reloaded)
dh = grid.Δxᶜᵃᵃ
halo_cells = ceil(Int, active_halo_width / dh)

# Get dimensions for ZFaceField (w) and CenterField (b) separately
# ZFaceField has Nz+1 z-faces, CenterField has Nz z-cells
Nx_w, Ny_w, Nz_w = size(interior(w))  # w is ZFaceField: Nz+1 z-faces
Nx_b, Ny_b, Nz_b = size(interior(b))  # b is CenterField: Nz z-cells

# Core region indices - same for x/y, different for z
core_x_range = (halo_cells + 1):(Nx_w - halo_cells)
core_y_range = (halo_cells + 1):(Ny_w - halo_cells)
core_z_faces = 1:Nz_w      # For ZFaceField (w, w_bar): all faces
core_z_centers = 1:Nz_b    # For CenterField (b, b_bar): all cells

println("Full subdomain:")
println("  w (ZFaceField):   $Nx_w x $Ny_w x $Nz_w")
println("  b (CenterField):  $Nx_b x $Ny_b x $Nz_b")
println("Halo cells: $halo_cells on each horizontal side")
println("Core region: $(length(core_x_range)) x $(length(core_y_range)) cells")
println("  w z-range: 1:$Nz_w ($(Nz_w) faces)")
println("  b z-range: 1:$Nz_b ($(Nz_b) cells)")

# Extract core region arrays - use correct z-range for each field type
w_bar_core = interior(w_bar)[core_x_range, core_y_range, core_z_faces]
b_bar_core = interior(b_bar)[core_x_range, core_y_range, core_z_centers]
wp_core = wp_full[core_x_range, core_y_range, core_z_faces]
bp_core = bp_full[core_x_range, core_y_range, core_z_centers]

println("\nCore region statistics:")
println("  * w_bar range: $(extrema(w_bar_core))")
println("  * b_bar range: $(extrema(b_bar_core))")
println("  * w' range: $(extrema(wp_core))")
println("  * b' range: $(extrema(bp_core))")

# ===============================================================================
# SECTION 7: QUADRANT ANALYSIS
# ===============================================================================

println("\n" * "="^70)
println("STEP 6: Quadrant Analysis (w' vs b')")
println("="^70)

# Interpolate w' to cell centers (average adjacent z-faces)
# w' is on z-faces (Nz+1), b' is at cell centers (Nz)
# Average w'[k] and w'[k+1] to get w' at cell center k
wp_centered = (wp_core[:, :, 1:end-1] .+ wp_core[:, :, 2:end]) ./ 2  # Now (Nx, Ny, Nz)
bp_centered = bp_core  # Already at cell centers, same Nz as wp_centered

# Flatten for histogram
w_vec = vec(wp_centered)
b_vec = vec(bp_centered)

# Compute 2D histogram
h = fit(Histogram, (w_vec, b_vec), nbins=N_BINS)
counts = h.weights
w_edges = h.edges[1]
b_edges = h.edges[2]

# Quadrant statistics
# Q1: w'>0, b'>0 (warm updrafts - buoyancy-driven convection)
# Q2: w'<0, b'>0 (warm downdrafts - counter-gradient)
# Q3: w'<0, b'<0 (cold downdrafts - convective plumes)
# Q4: w'>0, b'<0 (cold updrafts - counter-gradient)

n_total = length(w_vec)
n_Q1 = count((w_vec .> 0) .& (b_vec .> 0))
n_Q2 = count((w_vec .< 0) .& (b_vec .> 0))
n_Q3 = count((w_vec .< 0) .& (b_vec .< 0))
n_Q4 = count((w_vec .> 0) .& (b_vec .< 0))

# Compute w'b' (vertical buoyancy flux)
wb_flux = w_vec .* b_vec
mean_wb = sum(wb_flux) / n_total

println("\nQuadrant distribution:")
println("  +---------------+---------------+")
println("  |    Q2 (w'<0)  |    Q1 (w'>0)  |  b' > 0")
println("  |    $(lpad(round(100*n_Q2/n_total, digits=1), 5))%     |    $(lpad(round(100*n_Q1/n_total, digits=1), 5))%     |")
println("  +---------------+---------------+")
println("  |    Q3 (w'<0)  |    Q4 (w'>0)  |  b' < 0")
println("  |    $(lpad(round(100*n_Q3/n_total, digits=1), 5))%     |    $(lpad(round(100*n_Q4/n_total, digits=1), 5))%     |")
println("  +---------------+---------------+")

# Gradient-consistent quadrants: Q1 + Q3 (expected for convection)
# Counter-gradient quadrants: Q2 + Q4
gradient_frac = (n_Q1 + n_Q3) / n_total
counter_frac = (n_Q2 + n_Q4) / n_total

println("\nFlux analysis:")
println("  * Mean w'b': $(mean_wb) m^2/s^3")
println("  * Gradient-consistent (Q1+Q3): $(round(100*gradient_frac, digits=1))%")
println("  * Counter-gradient (Q2+Q4): $(round(100*counter_frac, digits=1))%")

println("\n2D Histogram:")
println("  * w' bins: $(length(w_edges)-1), range $(extrema(w_edges))")
println("  * b' bins: $(length(b_edges)-1), range $(extrema(b_edges))")
println("  * Max counts per bin: $(maximum(counts))")

# ===============================================================================
# SECTION 8: SUMMARY AND OUTPUT
# ===============================================================================

println("\n" * "="^70)
println("WORKFLOW COMPLETE")
println("="^70)

println("\nOutputs available in memory:")
println("  * w_bar_core: Coarse-grained vertical velocity (core region)")
println("  * b_bar_core: Coarse-grained buoyancy (core region)")
println("  * wp_core: Fine-scale w fluctuation (core region)")
println("  * bp_core: Fine-scale b fluctuation (core region)")
println("  * h: 2D histogram of (w', b') for quadrant analysis")
println("  * counts, w_edges, b_edges: Histogram data")

println("\nFor visualization, consider:")
println("  using CairoMakie")
println("  heatmap(w_edges, b_edges, counts')")
println("  # or")
println("  heatmap(w_bar_core[:,:,end])")
