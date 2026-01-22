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
using Statistics: std, quantile
using CairoMakie
using CairoMakie.Makie.Colors: RGB, RGBA, red, green, blue

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
const N_BINS = 60                       # Number of histogram bins
const COMPUTE_DEPTH_PROFILES = true     # Compute depth-resolved quadrant analysis
const SAVE_FIGURES = true               # Save quadrant analysis figures

# Fixed depth regions for aggregation (meters, negative down)
const DEPTH_SURFACE = (-10.0, 0.0)      # Near-surface layer
const DEPTH_MIXED = (-60.0, -10.0)      # Mixed layer interior  
const DEPTH_PYCNOCLINE = (-81.0, -60.0) # Pycnocline/entrainment zone

# Threshold for masking weak fluctuations based on |w'b'| magnitude
# Points with |w'b'| below the THRESHOLD_PERCENTILE of the distribution are masked
# This filters noise near origin AND along both axes (where flux is negligible)
const THRESHOLD_PERCENTILE = 0.50       # Mask bottom 50% of |w'b'| distribution

# Y-slices for x-z quadrant visualization (fraction of Ny_core)
const Y_SLICE_FRACS = (0.25, 0.5, 0.75) # 3 slices at 25%, 50%, 75% of domain

# Z-levels for x-y visualization (layer centers, meters)
const Z_LEVEL_FULL = -40.0              # Representative depth for full column
const Z_LEVEL_SURFACE = -5.0            # Center of surface layer
const Z_LEVEL_MIXED = -35.0             # Center of mixed layer interior
const Z_LEVEL_DEEP = -70.5              # Center of pycnocline/deep layer

# Quadrant analysis visualization settings
const QUADRANT_NAMES = ["Q1: w'>0, b'>0", "Q2: w'<0, b'>0", 
                        "Q3: w'<0, b'<0", "Q4: w'>0, b'<0"]
const QUADRANT_COLORS = [RGB(0.894, 0.102, 0.110),   # Q1: red - warm updrafts
                         RGB(0.216, 0.494, 0.722),   # Q2: blue - warm downdrafts
                         RGB(0.302, 0.686, 0.290),   # Q3: green - cold downdrafts
                         RGB(0.596, 0.306, 0.639)]   # Q4: purple - cold updrafts

# --- Processing Options ---
const SAVE_ALL_TILES = false      # Set false to skip tile extraction step
const TARGET_TILE = 4             # Which tile to process for coarse-graining

# ===============================================================================
# HELPER FUNCTIONS FOR QUADRANT ANALYSIS
# ===============================================================================

"""
    assign_quadrant(w, b) -> Int

Assign quadrant index based on signs of w' and b':
  - 1 = Q1: w'>0, b'>0 (warm updrafts - buoyancy-driven convection)
  - 2 = Q2: w'<0, b'>0 (warm downdrafts - counter-gradient)
  - 3 = Q3: w'<0, b'<0 (cold downdrafts - convective plumes)
  - 4 = Q4: w'>0, b'<0 (cold updrafts - counter-gradient)
"""
function assign_quadrant(w, b)
    if w > 0 && b > 0
        return 1  # Q1
    elseif w < 0 && b > 0
        return 2  # Q2
    elseif w < 0 && b < 0
        return 3  # Q3
    else
        return 4  # Q4
    end
end

"""
    assign_quadrants(w_arr, b_arr, mask_3d) -> Array{Int}

Vectorized quadrant assignment with masking. Returns 0 for masked points.
"""
function assign_quadrants(w_arr, b_arr, mask_3d)
    Q = similar(w_arr, Int)
    for i in eachindex(w_arr, b_arr, mask_3d)
        Q[i] = mask_3d[i] ? assign_quadrant(w_arr[i], b_arr[i]) : 0
    end
    return Q
end

"""
    quadrant_stats_in_region(wp, bp, mask, depth_name) -> NamedTuple

Compute quadrant statistics for a depth region, including:
- Quadrant fractions (Q1-Q4)
- Gradient-consistent vs counter-gradient fractions
- Mean w'b' flux decomposed by quadrant type
"""
function quadrant_stats_in_region(wp, bp, mask, depth_name)
    n_sig = count(mask)
    if n_sig == 0
        return (name=depth_name, n_total=0, 
                Q1_frac=0.0, Q2_frac=0.0, Q3_frac=0.0, Q4_frac=0.0,
                mean_wb=0.0, gradient_frac=0.0, counter_frac=0.0,
                wb_gradient=0.0, wb_counter=0.0)
    end
    
    # Quadrant masks (using significant points only)
    Q1_mask = (wp .> 0) .& (bp .> 0) .& mask  # warm updrafts
    Q2_mask = (wp .< 0) .& (bp .> 0) .& mask  # warm downdrafts  
    Q3_mask = (wp .< 0) .& (bp .< 0) .& mask  # cold downdrafts
    Q4_mask = (wp .> 0) .& (bp .< 0) .& mask  # cold updrafts
    
    n_Q1, n_Q2, n_Q3, n_Q4 = count(Q1_mask), count(Q2_mask), count(Q3_mask), count(Q4_mask)
    
    # Flux contributions by quadrant
    wb_Q1 = sum(wp[Q1_mask] .* bp[Q1_mask])  # positive (upward buoyancy)
    wb_Q2 = sum(wp[Q2_mask] .* bp[Q2_mask])  # negative
    wb_Q3 = sum(wp[Q3_mask] .* bp[Q3_mask])  # positive (downward cold)
    wb_Q4 = sum(wp[Q4_mask] .* bp[Q4_mask])  # negative
    
    # Gradient-consistent: Q1 + Q3 (w'b' > 0, drives convection)
    # Counter-gradient: Q2 + Q4 (w'b' < 0, restratifying)
    gradient_frac = (n_Q1 + n_Q3) / n_sig
    counter_frac = (n_Q2 + n_Q4) / n_sig
    wb_gradient = (wb_Q1 + wb_Q3) / n_sig
    wb_counter = (wb_Q2 + wb_Q4) / n_sig
    
    mean_wb = (wb_Q1 + wb_Q2 + wb_Q3 + wb_Q4) / n_sig
    
    return (name=depth_name, n_total=n_sig,
            Q1_frac=n_Q1/n_sig, Q2_frac=n_Q2/n_sig,
            Q3_frac=n_Q3/n_sig, Q4_frac=n_Q4/n_sig,
            mean_wb=mean_wb,
            gradient_frac=gradient_frac, counter_frac=counter_frac,
            wb_gradient=wb_gradient, wb_counter=wb_counter)
end

"""
    print_quadrant_summary(stats)

Print formatted quadrant statistics for a depth region.
"""
function print_quadrant_summary(stats)
    println("\n  $(stats.name) (N = $(stats.n_total) significant points):")
    println("    Quadrant fractions:")
    println("      Q1 (w'>0,b'>0): $(round(100*stats.Q1_frac, digits=1))%")
    println("      Q2 (w'<0,b'>0): $(round(100*stats.Q2_frac, digits=1))%")
    println("      Q3 (w'<0,b'<0): $(round(100*stats.Q3_frac, digits=1))%")
    println("      Q4 (w'>0,b'<0): $(round(100*stats.Q4_frac, digits=1))%")
    println("    Flux decomposition:")
    println("      Gradient-consistent (Q1+Q3): $(round(100*stats.gradient_frac, digits=1))%")
    println("      Counter-gradient (Q2+Q4):    $(round(100*stats.counter_frac, digits=1))%")
    println("      Mean w'b' (total):    $(round(stats.mean_wb, sigdigits=3)) m²/s³")
    println("      Mean w'b' (gradient): $(round(stats.wb_gradient, sigdigits=3)) m²/s³")
    println("      Mean w'b' (counter):  $(round(stats.wb_counter, sigdigits=3)) m²/s³")
end

"""
    compute_mean_depth_per_bin(wp, bp, z_arr, w_edges, b_edges) -> Matrix

Compute the mean depth for each histogram bin in the (w', b') plane.
"""
function compute_mean_depth_per_bin(wp, bp, z_arr, w_edges, b_edges)
    nw, nb = length(w_edges)-1, length(b_edges)-1
    depth_sum = zeros(nw, nb)
    depth_count = zeros(Int, nw, nb)
    
    for k in axes(wp, 3), j in axes(wp, 2), i in axes(wp, 1)
        wi = searchsortedfirst(w_edges, wp[i,j,k]) - 1
        bi = searchsortedfirst(b_edges, bp[i,j,k]) - 1
        if 1 <= wi <= nw && 1 <= bi <= nb
            depth_sum[wi, bi] += z_arr[k]
            depth_count[wi, bi] += 1
        end
    end
    mean_depth = depth_sum ./ max.(depth_count, 1)
    mean_depth[depth_count .== 0] .= NaN
    return mean_depth
end

"""
    depth_indices(z_arr, depth_limits) -> Vector{Int}

Find indices in z_arr that fall within depth_limits = (z_min, z_max).
"""
function depth_indices(z_arr, depth_limits)
    return findall(depth_limits[1] .<= z_arr .<= depth_limits[2])
end

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

# Get core dimensions
Nx_core, Ny_core, Nz_core = size(wp_centered)

# Flatten for histogram
w_vec = vec(wp_centered)
b_vec = vec(bp_centered)

# -----------------------------------------------------------------------------
# Mask weak fluctuations (noise filtering based on |w'b'| magnitude)
# -----------------------------------------------------------------------------
# Compute |w'b'| for each point
wb_magnitude_vec = abs.(w_vec .* b_vec)
wb_magnitude_3d = abs.(wp_centered .* bp_centered)

# Use percentile-based threshold: dynamically adapts to data distribution
# This filters: (1) near origin, (2) along w'-axis, (3) along b'-axis
wb_threshold = quantile(wb_magnitude_vec, THRESHOLD_PERCENTILE)

# Create mask: true for significant points (|w'b'| above threshold)
significant_mask = wb_magnitude_vec .>= wb_threshold
sig_mask_3d = wb_magnitude_3d .>= wb_threshold

# Statistics for output
w_std = std(w_vec)
b_std = std(b_vec)

println("\nMasking weak fluctuations (|w'b'| percentile method):")
println("  * w' std: $(round(w_std, sigdigits=3)) m/s")
println("  * b' std: $(round(b_std, sigdigits=3)) m/s²")
println("  * |w'b'| threshold ($(Int(THRESHOLD_PERCENTILE*100))th percentile): $(round(wb_threshold, sigdigits=3)) m²/s³")
println("  * Significant points: $(count(significant_mask))/$(length(significant_mask)) ($(round(100*count(significant_mask)/length(significant_mask), digits=1))%)")

# -----------------------------------------------------------------------------
# Z-coordinates for cell centers
# -----------------------------------------------------------------------------
dz_core = abs(Z_LIMITS[2] - Z_LIMITS[1]) / Nz_core
z_centers = collect(range(Z_LIMITS[1] + dz_core/2, Z_LIMITS[2] - dz_core/2, length=Nz_core))

# Find depth indices for each region
k_surface = depth_indices(z_centers, DEPTH_SURFACE)
k_mixed = depth_indices(z_centers, DEPTH_MIXED)
k_pycnocline = depth_indices(z_centers, DEPTH_PYCNOCLINE)

println("\nDepth regions:")
println("  * Surface $(DEPTH_SURFACE): z-indices $(first(k_surface)):$(last(k_surface)) ($(length(k_surface)) levels)")
println("  * Mixed layer $(DEPTH_MIXED): z-indices $(first(k_mixed)):$(last(k_mixed)) ($(length(k_mixed)) levels)")
println("  * Pycnocline $(DEPTH_PYCNOCLINE): z-indices $(first(k_pycnocline)):$(last(k_pycnocline)) ($(length(k_pycnocline)) levels)")

# -----------------------------------------------------------------------------
# Compute depth-resolved quadrant statistics
# -----------------------------------------------------------------------------
stats_global = quadrant_stats_in_region(wp_centered, bp_centered, sig_mask_3d, "Global")
stats_surface = quadrant_stats_in_region(
    wp_centered[:,:,k_surface], bp_centered[:,:,k_surface], 
    sig_mask_3d[:,:,k_surface], "Surface $(DEPTH_SURFACE) m")
stats_mixed = quadrant_stats_in_region(
    wp_centered[:,:,k_mixed], bp_centered[:,:,k_mixed],
    sig_mask_3d[:,:,k_mixed], "Mixed Layer $(DEPTH_MIXED) m")
stats_pycnocline = quadrant_stats_in_region(
    wp_centered[:,:,k_pycnocline], bp_centered[:,:,k_pycnocline],
    sig_mask_3d[:,:,k_pycnocline], "Pycnocline $(DEPTH_PYCNOCLINE) m")

println("\n" * "="^70)
println("QUADRANT ANALYSIS BY DEPTH REGION")
println("="^70)
for stats in [stats_global, stats_surface, stats_mixed, stats_pycnocline]
    print_quadrant_summary(stats)
end

# Compute 2D histogram (significant points only)
h = fit(Histogram, (w_vec[significant_mask], b_vec[significant_mask]), nbins=N_BINS)
counts = h.weights
w_edges = collect(h.edges[1])
b_edges = collect(h.edges[2])

println("\n2D Histogram (significant points only):")
println("  * w' bins: $(length(w_edges)-1), range $(extrema(w_edges))")
println("  * b' bins: $(length(b_edges)-1), range $(extrema(b_edges))")
println("  * Max counts per bin: $(maximum(counts))")

# =============================================================================
# SECTION 7A: FIGURE 1 - 2×2 Histogram Layout
# =============================================================================

fig1_path = OUTPUT_DIR * "quadrant_histograms_tile$(TARGET_TILE)_iter$(ITERATION).pdf"
if SAVE_FIGURES && !isfile(fig1_path)
    println("\n" * "="^70)
    println("STEP 7: Generating Quadrant Analysis Figures")
    println("="^70)
    
    set_theme!(theme_latexfonts(), fontsize=12, figure_padding = 10)
    
    # Create masks for each depth region (significant points only)
    sig_surface = sig_mask_3d[:,:,k_surface]
    sig_mixed = sig_mask_3d[:,:,k_mixed]
    
    # Extract significant points for histograms
    wp_surface_sig = vec(wp_centered[:,:,k_surface])[vec(sig_surface)]
    bp_surface_sig = vec(bp_centered[:,:,k_surface])[vec(sig_surface)]
    wp_mixed_sig = vec(wp_centered[:,:,k_mixed])[vec(sig_mixed)]
    bp_mixed_sig = vec(bp_centered[:,:,k_mixed])[vec(sig_mixed)]
    
    # Compute histograms (all use significant points only)
    h_global = fit(Histogram, (w_vec[significant_mask], b_vec[significant_mask]), nbins=N_BINS)
    h_surface = fit(Histogram, (wp_surface_sig, bp_surface_sig), nbins=N_BINS)
    h_mixed = fit(Histogram, (wp_mixed_sig, bp_mixed_sig), nbins=N_BINS)
    
    # Get axis limits from global histogram edges
    w_lim = (first(h_global.edges[1]), last(h_global.edges[1]))
    b_lim = (first(h_global.edges[2]), last(h_global.edges[2]))
    
    # Helper to shade the masked region |w'b'| < wb_threshold using band!
    # The hyperbola b = ±wb_thresh/w defines the boundary
    function add_wb_threshold_mask!(ax, wb_thresh, w_lim, b_lim; n_pts=200, alpha=0.7)
        # For w > 0: shade between b = -wb_thresh/w and b = +wb_thresh/w
        w_pos = range(wb_thresh / abs(b_lim[2]), w_lim[2], length=n_pts)
        w_pos = collect(filter(w -> w > 1e-12, w_pos))
        if length(w_pos) > 1
            b_upper = clamp.(wb_thresh ./ w_pos, b_lim[1], b_lim[2])
            b_lower = clamp.(-wb_thresh ./ w_pos, b_lim[1], b_lim[2])
            band!(ax, w_pos, b_lower, b_upper; color = (:white, alpha))
        end
        
        # For w < 0: shade between b = -wb_thresh/w and b = +wb_thresh/w
        w_neg = range(w_lim[1], -wb_thresh / abs(b_lim[2]), length=n_pts)
        w_neg = collect(filter(w -> w < -1e-12, w_neg))
        if length(w_neg) > 1
            b_upper = clamp.(-wb_thresh ./ w_neg, b_lim[1], b_lim[2])
            b_lower = clamp.(wb_thresh ./ w_neg, b_lim[1], b_lim[2])
            band!(ax, w_neg, b_lower, b_upper; color = (:white, alpha))
        end
    end
    
    # Mean depth per bin
    mean_depth = compute_mean_depth_per_bin(wp_centered, bp_centered, z_centers, 
                                             collect(h_global.edges[1]), collect(h_global.edges[2]))
    
    fig1 = Figure(size = (720, 560))
    
    # [1,1] Global histogram with log counts
    ax11 = Axis(fig1[1,1]; xlabel=L"w^\prime~\text{(m s^{-1})}", ylabel=L"b^\prime~\text{(m s^{-2})}",
                title=L"\text{(a) Global, log counts}", limits=(w_lim, b_lim))
    hm11 = heatmap!(ax11, h_global.edges[1], h_global.edges[2], log10.(1 .+ h_global.weights); 
                    rasterize=true, colormap=Reverse(:grays))
    add_wb_threshold_mask!(ax11, wb_threshold, w_lim, b_lim)
    Colorbar(fig1[1,2], hm11, label=L"\log_{10}(1+N)")
    
    # [1,2] Global histogram colored by mean depth
    ax12 = Axis(fig1[1,3]; xlabel=L"w^\prime~\text{(m s^{-1})}", 
                title=L"\text{(b) Global, mean depth}", limits=(w_lim, b_lim))
    hm12 = heatmap!(ax12, h_global.edges[1], h_global.edges[2], mean_depth; 
                    rasterize=true, colormap=Reverse(:deep), colorrange=(Z_LIMITS[1], 0))
    add_wb_threshold_mask!(ax12, wb_threshold, w_lim, b_lim)
    Colorbar(fig1[1,4], hm12, label=L"\bar{z}~\text{(m)}")
    hideydecorations!(ax12, ticks = false)
    
    # [2,1] Surface region histogram  
    ax21 = Axis(fig1[2,1]; xlabel=L"w^\prime~\text{(m s^{-1})}", ylabel=L"b^\prime~\text{(m s^{-2})}",
                title=L"\text{(c) Surface }z\in[-10,0]~\text{m}", limits=(w_lim, b_lim))
    hm21 = heatmap!(ax21, h_surface.edges[1], h_surface.edges[2], log10.(1 .+ h_surface.weights); 
                    rasterize=true, colormap=Reverse(:grays))
    add_wb_threshold_mask!(ax21, wb_threshold, w_lim, b_lim)
    Colorbar(fig1[2,2], hm21, label=L"\log_{10}(1+N)")
    
    # [2,2] Mixed layer histogram
    ax22 = Axis(fig1[2,3]; xlabel=L"w^\prime~\text{(m s^{-1})}", 
                title=L"\text{(d) Mixed layer }z\in[-60,-10]~\text{m}", limits=(w_lim, b_lim))
    hm22 = heatmap!(ax22, h_mixed.edges[1], h_mixed.edges[2], log10.(1 .+ h_mixed.weights); 
                    rasterize=true, colormap=Reverse(:grays))
    add_wb_threshold_mask!(ax22, wb_threshold, w_lim, b_lim)
    Colorbar(fig1[2,4], hm22, label=L"\log_{10}(1+N)")
    hideydecorations!(ax22, ticks = false)
    
    colgap!(fig1.layout, 2, 15)
    rowgap!(fig1.layout, 1, 10)
    resize_to_layout!(fig1)
    
    save(fig1_path, fig1; pt_per_unit=1)
    println("  Saved: $fig1_path")
elseif SAVE_FIGURES
    println("  Skipping Figure 1: $fig1_path already exists")
end

# =============================================================================
# SECTION 7B: FIGURE 2 - 3×1 x-z Slices with Quadrant Spatial Distribution
# =============================================================================

fig2_path = OUTPUT_DIR * "quadrant_xz_slices_tile$(TARGET_TILE)_iter$(ITERATION).pdf"
if SAVE_FIGURES && !isfile(fig2_path)
    # Assign quadrant categories to 3D field
    Q_field = assign_quadrants(wp_centered, bp_centered, sig_mask_3d)
    
    # Compute |w'b'| for alpha/intensity scaling
    wb_magnitude = abs.(wp_centered .* bp_centered)
    wb_ref = quantile(vec(wb_magnitude[sig_mask_3d]), 0.99)  # 99th percentile for scaling
    alpha_field = clamp.(wb_magnitude ./ wb_ref, 0, 1)
    
    # Get tile's absolute position in the full domain (from tiles array)
    tile_info = tiles[TARGET_TILE]
    x_start, x_end = tile_info.core_xlims[1] / 1e3, tile_info.core_xlims[2] / 1e3  # km
    y_start, y_end = tile_info.core_ylims[1] / 1e3, tile_info.core_ylims[2] / 1e3  # km
    z_start, z_end = Z_LIMITS[1], Z_LIMITS[2]
    
    # y-slice indices and their absolute y-positions
    j_slices = [max(1, round(Int, f * Ny_core)) for f in Y_SLICE_FRACS]
    y_positions = [round(y_start + f * (y_end - y_start), digits=1) for f in Y_SLICE_FRACS]
    
    # Create categorical colormap
    QCMAP = cgrad(QUADRANT_COLORS, 4, categorical=true)
    
    fig2 = Figure(size = (540, 480))
    
    panel_labels = ["(a)", "(b)", "(c)"]
    local ax_first  # Declare local to avoid scope ambiguity
    for (row, (j_slice, y_pos)) in enumerate(zip(j_slices, y_positions))
        ax = Axis(fig2[row, 1]; 
                  xlabel = row == 3 ? L"x~\text{(km)}" : "",
                  ylabel = L"z~\text{(m)}",
                  title = L"\text{%$(panel_labels[row])}~y = %$(y_pos)~\text{km}",
                  limits = ((x_start, x_end), (z_start, z_end)))
        
        if row == 1
            ax_first = ax
        end
        
        # Extract slice
        Q_slice = Q_field[:, j_slice, :]
        alpha_slice = alpha_field[:, j_slice, :]
        
        # Create RGBA image for quadrant visualization with alpha
        rgba_data = fill(RGBA(1.0, 1.0, 1.0, 0.0), Nx_core, Nz_core)
        for k in 1:Nz_core, i in 1:Nx_core
            q = Q_slice[i, k]
            if q > 0
                c = QUADRANT_COLORS[q]
                rgba_data[i, k] = RGBA(red(c), green(c), blue(c), alpha_slice[i, k])
            end
        end
        
        image!(ax, (x_start, x_end), (z_start, z_end), rgba_data)
        
        if row < 3
            hidexdecorations!(ax, ticks = false)
        end
    end
    
    # Legend in first subplot using PolyElements
    legend_elements = [PolyElement(color=c) for c in QUADRANT_COLORS]
    axislegend(ax_first, legend_elements, QUADRANT_NAMES, position = :rt, 
               labelsize=10, patchsize = (15, 10), framevisible = false, 
               padding = (0f0, 0f0, 0f0, 0f0), patchlabelgap = 3, rowgap = 1)
    
    resize_to_layout!(fig2)
    
    save(fig2_path, fig2; pt_per_unit=1)
    println("  Saved: $fig2_path")
elseif SAVE_FIGURES
    println("  Skipping Figure 2: $fig2_path already exists")
end

# =============================================================================
# SECTION 7C: FIGURE 3 - 2×2 x-y Slices at Different Z-levels
# =============================================================================

fig3_path = OUTPUT_DIR * "quadrant_xy_slices_tile$(TARGET_TILE)_iter$(ITERATION).pdf"
if SAVE_FIGURES && !isfile(fig3_path)
    # Find z-indices for each level
    find_z_index(z_level) = argmin(abs.(z_centers .- z_level))
    
    z_levels = [Z_LEVEL_FULL, Z_LEVEL_SURFACE, Z_LEVEL_MIXED, Z_LEVEL_DEEP]
    z_titles = [L"\text{(a) Full depth repr. }z=%$(Int(Z_LEVEL_FULL))~\text{m}",
                L"\text{(b) Surface }z=%$(Int(Z_LEVEL_SURFACE))~\text{m}",
                L"\text{(c) Mixed layer }z=%$(Int(Z_LEVEL_MIXED))~\text{m}",
                L"\text{(d) Below ML }z=%$((Z_LEVEL_DEEP))~\text{m}"]
    
    fig3 = Figure(size = (560, 560))
    
    local ax_first  # Declare local to avoid scope ambiguity
    for (idx, (z_lev, ztitle)) in enumerate(zip(z_levels, z_titles))
        row = (idx - 1) ÷ 2 + 1
        col = (idx - 1) % 2 + 1
        
        k = find_z_index(z_lev)
        
        ax = Axis(fig3[row, col];
                  xlabel = row == 2 ? L"x~\text{(km)}" : "",
                  ylabel = col == 1 ? L"y~\text{(km)}" : "",
                  title = ztitle,
                  aspect = 1,
                  limits = ((x_start, x_end), (y_start, y_end)))
        
        if idx == 1
            ax_first = ax
        end
        
        Q_slice = Q_field[:, :, k]
        alpha_slice = alpha_field[:, :, k]
        
        # Create RGBA image
        rgba_data = fill(RGBA(1.0, 1.0, 1.0, 0.0), Nx_core, Ny_core)
        for j in 1:Ny_core, i in 1:Nx_core
            q = Q_slice[i, j]
            if q > 0
                c = QUADRANT_COLORS[q]
                rgba_data[i, j] = RGBA(red(c), green(c), blue(c), alpha_slice[i, j])
            end
        end
        
        image!(ax, (x_start, x_end), (y_start, y_end), rgba_data; rasterize=true)
        
        if row == 1
            hidexdecorations!(ax, ticks = false)
        end
        if col == 2
            hideydecorations!(ax, ticks = false)
        end
    end
    
    # Legend in first subplot using PolyElements
    legend_elements = [PolyElement(color=c) for c in QUADRANT_COLORS]
    axislegend(ax_first, legend_elements, QUADRANT_NAMES, position = :rb, 
               labelsize=10, patchsize = (15, 10), framevisible = false, 
               padding = (0f0, 0f0, 0f0, 0f0), patchlabelgap = 3, rowgap = 1)
    
    colgap!(fig3.layout, 1, 10)
    rowgap!(fig3.layout, 1, 10)
    resize_to_layout!(fig3)
    
    save(fig3_path, fig3; pt_per_unit=1)
    println("  Saved: $fig3_path")
elseif SAVE_FIGURES
    println("  Skipping Figure 3: $fig3_path already exists")
end

# ===============================================================================
# SECTION 8: SUMMARY AND OUTPUT
# ===============================================================================

println("\n" * "="^70)
println("WORKFLOW COMPLETE")
println("="^70)

println("\nOutputs available in memory:")
println("  * w_bar_core: Coarse-grained vertical velocity (core region)")
println("  * b_bar_core: Coarse-grained buoyancy (core region)")
println("  * wp_centered: Fine-scale w' at cell centers (core region)")
println("  * bp_centered: Fine-scale b' at cell centers (core region)")
println("  * Q_field: Quadrant assignment (1-4, 0=masked)")
println("  * sig_mask_3d: Significant points mask")
println("  * z_centers: Z-coordinates for cell centers")

if SAVE_FIGURES
    println("\nSaved figures:")
    println("  * $(OUTPUT_DIR)quadrant_histograms_tile$(TARGET_TILE)_iter$(ITERATION).pdf")
    println("  * $(OUTPUT_DIR)quadrant_xz_slices_tile$(TARGET_TILE)_iter$(ITERATION).pdf")
    println("  * $(OUTPUT_DIR)quadrant_xy_slices_tile$(TARGET_TILE)_iter$(ITERATION).pdf")
end
