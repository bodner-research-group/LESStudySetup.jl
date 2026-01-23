# ===============================================================================
#                    PUBLICATION QUADRANT ANALYSIS FIGURES
# ===============================================================================
#
# This script generates two publication-quality figures for quadrant analysis
# of w' and b' fluctuations from LES data:
#
#   Figure 1 (3×2): X-Y analysis at three Z-levels (-5m, -35m, -70.5m)
#     - Left column: 2D histograms of w' vs b'
#     - Right column: Quadrant spatial patterns in x-y planes
#
#   Figure 2 (3×2): X-Z analysis at three Y-slices (30km, 50km, 70km)
#     - Left column: 2D histograms of w' vs b'
#     - Right column: Quadrant spatial patterns in x-z planes
#
# Pipeline:
#   1. Load data with 20km halo for 10km filter
#   2. Coarse-grain w and b fields (10km Gaussian)
#   3. Compute residuals w' = w - w̄, b' = b - b̄
#   4. Extract core region (remove halo)
#   5. Mask bottom 50% of |w'b'| distribution
#   6. Generate publication figures
#
# ===============================================================================

using Oceananigans
using Oceananigans: location, architecture, fill_halo_regions!
using JLD2
using LESStudySetup
using LESStudySetup.Diagnostics
using LESStudySetup.Diagnostics: load_distributed_checkpoint_subdomain
using LESStudySetup.Diagnostics: load_checkpoint_clock
using LESStudySetup.Diagnostics: coarse_graining!
using LESStudySetup.Diagnostics: save_subdomain_with_halo
using LESStudySetup.Diagnostics: load_subdomain_snapshot
using StatsBase: fit, Histogram
using Statistics: std, quantile
using CairoMakie
using CairoMakie.Makie.Colors: RGB, RGBA, red, green, blue

# ===============================================================================
# SECTION 1: CONFIGURATION
# ===============================================================================

# --- Data Paths ---
const CHECKPOINT_DIR = "/orcd/data/abodner/002/shared_datasets/nhyles_output/"
const CHECKPOINT_PREFIX = CHECKPOINT_DIR * "iteration16x/nonhydrostatic_checkpoint_"
const OUTPUT_DIR = CHECKPOINT_DIR * "publication_figures/"

# --- Checkpoint Selection ---
const ITERATION = 164410

# --- Grid Parameters (from simulation) ---
const Δh = 4.8828125              # Horizontal grid spacing (m)
const Δz = 1.125                  # Vertical grid spacing (m)
const Lz = 252.0                  # Domain depth (m)
const Nz = 72                     # Number of vertical levels

# --- Figure 1: X-Y Analysis at Z-levels ---
const XY_CORE_XLIMS = (-40e3, 10e3)    # Core region: 50km in x
const XY_CORE_YLIMS = (25e3, 75e3)     # Core region: 50km in y
const XY_HALO = 20e3                    # Halo width for 10km filter
const XY_FULL_XLIMS = (XY_CORE_XLIMS[1] - XY_HALO, XY_CORE_XLIMS[2] + XY_HALO)  # (-60km, 30km)
const XY_FULL_YLIMS = (XY_CORE_YLIMS[1] - XY_HALO, XY_CORE_YLIMS[2] + XY_HALO)  # (5km, 95km)
const Z_TARGETS = [-5.0, -35.0, -70.5]  # Target depths (m)

# --- Figure 2: X-Z Analysis at Y-slices ---
const XZ_CORE_XLIMS = (-40e3, 10e3)    # Same x-range as Figure 1
const XZ_ZLIMS = (-81.0, 0.0)          # Full depth
const Y_TARGETS = [30e3, 50e3, 70e3]   # Target y-positions (m)
const Y_BAND_HALF_WIDTH = 50.0         # Load ±50m around each y-slice

# --- Coarse-Graining Parameters ---
const FILTER_CUTOFF = 10e3             # 10km filter cutoff
const FILTER_KERNEL = :gaussian        # Gaussian kernel
const FILTER_BORDER = :reflect         # Boundary handling for x-y
const FILTER_BORDER_XZ = :ycircular    # For narrow y-bands, treat y as circular

# --- Physical Constants ---
const α_thermal = 2e-4                 # Thermal expansion coefficient (1/K)
const g_gravity = 9.81                 # Gravitational acceleration (m/s²)

# --- Quadrant Analysis Parameters ---
const N_BINS = 60                      # Number of histogram bins
const THRESHOLD_PERCENTILE = 0.50      # Mask bottom 50% of |w'b'|

# --- Visualization Settings ---
const QUADRANT_COLORS = [
    RGB(0.894, 0.102, 0.110),   # Q1: red - warm updrafts
    RGB(0.216, 0.494, 0.722),   # Q2: blue - warm downdrafts
    RGB(0.302, 0.686, 0.290),   # Q3: green - cold downdrafts
    RGB(0.596, 0.306, 0.639),   # Q4: purple - cold updrafts
]
const QUADRANT_NAMES = [L"Q1: w'>0, b'>0", L"Q2: w'<0, b'>0", 
                        L"Q3: w'<0, b'<0", L"Q4: w'>0, b'<0"]

# --- Data Saving Options ---
const SAVE_LOADED_DATA = true          # Save loaded data before coarse-graining
const SUBDOMAIN_DIR = OUTPUT_DIR * "subdomains/"  # Directory for saved subdomains

# ===============================================================================
# SECTION 2: HELPER FUNCTIONS
# ===============================================================================

"""
    compute_z_indices(z_targets, Δz, Lz) -> Vector{Int}

Find grid indices closest to target z-depths.
Grid formula: z = -Lz + (k-1) * Δz, so k = (z + Lz) / Δz + 1
"""
function compute_z_indices(z_targets, Δz, Lz)
    return [round(Int, (z + Lz) / Δz) + 1 for z in z_targets]
end

"""
    compute_actual_z(k, Δz, Lz) -> Float64

Compute actual z-depth for grid index k (cell center).
"""
function compute_actual_z(k, Δz, Lz)
    return -Lz + (k - 1) * Δz + Δz/2
end

"""
    compute_y_band_limits(y_target, half_width) -> Tuple

Return (y_min, y_max) for loading a narrow y-band.
"""
function compute_y_band_limits(y_target, half_width)
    return (y_target - half_width, y_target + half_width)
end

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
    assign_quadrants(w_arr, b_arr, mask) -> Array{Int}

Vectorized quadrant assignment with masking. Returns 0 for masked points.
"""
function assign_quadrants(w_arr, b_arr, mask)
    Q = similar(w_arr, Int)
    for i in eachindex(w_arr, b_arr, mask)
        Q[i] = mask[i] ? assign_quadrant(w_arr[i], b_arr[i]) : 0
    end
    return Q
end

"""
    create_quadrant_rgba(Q_field, alpha_field, colors) -> Array{RGBA}

Create RGBA image array from quadrant field with alpha based on |w'b'| magnitude.
"""
function create_quadrant_rgba(Q_field, alpha_field, colors)
    Nx, Ny = size(Q_field)
    rgba_data = fill(RGBA(1.0, 1.0, 1.0, 0.0), Nx, Ny)
    for j in 1:Ny, i in 1:Nx
        q = Q_field[i, j]
        if q > 0
            c = colors[q]
            rgba_data[i, j] = RGBA(red(c), green(c), blue(c), alpha_field[i, j])
        end
    end
    return rgba_data
end

"""
    add_wb_threshold_mask!(ax, wb_thresh, w_lim, b_lim; n_pts=200, alpha=0.7)

Shade the masked region |w'b'| < wb_threshold using hyperbolic bands.
"""
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

"""
    extract_core_xy(data, halo_cells) -> Array

Remove halo from array in both horizontal dimensions (for x-y analysis).
"""
function extract_core_xy(data, halo_cells)
    return data[(halo_cells+1):(end-halo_cells), (halo_cells+1):(end-halo_cells), :]
end

"""
    extract_core_x(data, halo_cells) -> Array

Remove halo from array in x-dimension only (for x-z analysis).
"""
function extract_core_x(data, halo_cells)
    return data[(halo_cells+1):(end-halo_cells), :, :]
end

"""
    create_compact_snapshot_from_levels(snapshot, z_indices, Δz, Lz) -> Dict

Create a compact snapshot with a grid matching the discrete z-levels.

When `load_distributed_checkpoint_subdomain()` is called with `levels`,
it creates a windowed field on the full-depth grid. This function creates
a new grid with Nz = length(z_indices) and copies the field data, ensuring
the grid dimensions match the actual data shape for consistent save/load.

Arguments:
- `snapshot`: Dict from `load_distributed_checkpoint_subdomain()` with `:grid`, `:u`, `:v`, `:w`, `:T`
- `z_indices`: Vector of z-level indices that were loaded
- `Δz`: Vertical grid spacing (m)
- `Lz`: Total domain depth (m)

Returns:
- New Dict with compact grid and fields with matching dimensions
"""
function create_compact_snapshot_from_levels(snapshot, z_indices, Δz, Lz)
    old_grid = snapshot[:grid]
    Nz_compact = length(z_indices)
    
    # Compute z-bounds for the discrete levels
    z_centers = [-Lz + (k - 1) * Δz + Δz/2 for k in z_indices]
    z_min = minimum(z_centers) - Δz/2
    z_max = maximum(z_centers) + Δz/2
    
    # Extract x and y extents from original grid
    x_min = old_grid.xᶜᵃᵃ[1] - old_grid.Δxᶜᵃᵃ/2
    x_max = old_grid.xᶜᵃᵃ[old_grid.Nx] + old_grid.Δxᶜᵃᵃ/2
    y_min = old_grid.yᵃᶜᵃ[1] - old_grid.Δyᵃᶜᵃ/2
    y_max = old_grid.yᵃᶜᵃ[old_grid.Ny] + old_grid.Δyᵃᶜᵃ/2
    
    compact_grid = RectilinearGrid(architecture(old_grid);
        size = (old_grid.Nx, old_grid.Ny, Nz_compact),
        x = (x_min, x_max),
        y = (y_min, y_max),
        z = (z_min, z_max),
        topology = (Bounded, Bounded, Bounded))
    
    compact_snapshot = Dict{Symbol, Any}(:grid => compact_grid)
    
    for (name, FieldType) in [(:u, XFaceField), (:v, YFaceField), (:w, ZFaceField), (:T, CenterField)]
        if haskey(snapshot, name)
            new_field = FieldType(compact_grid)
            interior(new_field) .= interior(snapshot[name])
            fill_halo_regions!(new_field)
            compact_snapshot[name] = new_field
        end
    end
    
    println("  Created compact snapshot: $(old_grid.Nx) × $(old_grid.Ny) × $(Nz_compact)")
    return compact_snapshot
end

# ===============================================================================
# SECTION 3: FIGURE 1 - X-Y ANALYSIS AT Z-LEVELS
# ===============================================================================

println("\n" * "="^70)
println("FIGURE 1: X-Y Quadrant Analysis at Z-Levels")
println("="^70)

# Ensure output directory exists
mkpath(OUTPUT_DIR)

# 3.1 COMPUTE Z-INDICES
z_indices = compute_z_indices(Z_TARGETS, Δz, Lz)
actual_z = [compute_actual_z(k, Δz, Lz) for k in z_indices]

println("\nTarget vs Actual Z-levels:")
for (i, (target, actual, k)) in enumerate(zip(Z_TARGETS, actual_z, z_indices))
    println("  Level $i: target=$(target)m → index=$k → actual=$(round(actual, digits=2))m")
end

# 3.2 LOAD DATA
println("\nLoading data for Figure 1...")
println("  X range (with halo): $(XY_FULL_XLIMS)")
println("  Y range (with halo): $(XY_FULL_YLIMS)")
println("  Z indices: $z_indices")

# Check if saved subdomain exists first
mkpath(SUBDOMAIN_DIR)
xy_save_file = SUBDOMAIN_DIR * "subdomain_xy_levels_iter$(ITERATION).jld2"

if isfile(xy_save_file)
    println("\nLoading from saved subdomain: $xy_save_file")
    snapshot_xy = load_subdomain_snapshot(xy_save_file; variables=("u", "v", "w", "T"))
else
    println("\nLoading from distributed checkpoint...")
    snapshot_xy = load_distributed_checkpoint_subdomain(CHECKPOINT_PREFIX, ITERATION;
        xlims = XY_FULL_XLIMS,
        ylims = XY_FULL_YLIMS,
        levels = z_indices,
        getEw = false,
        getMLD = 0
    )
    
    # Save for future runs
    if SAVE_LOADED_DATA
        println("\nCompacting snapshot for storage...")
        snapshot_xy_compact = create_compact_snapshot_from_levels(snapshot_xy, z_indices, Δz, Lz)
        
        println("Saving compacted subdomain to: $xy_save_file")
        save_subdomain_with_halo(xy_save_file, snapshot_xy_compact;
            core_xlims = XY_CORE_XLIMS,
            core_ylims = XY_CORE_YLIMS,
            halo_width = XY_HALO,
            levels = z_indices,
            iteration = ITERATION
        )
    end
end

grid_xy = snapshot_xy[:grid]
println("\nLoaded subdomain grid: $(grid_xy.Nx) × $(grid_xy.Ny) × $(grid_xy.Nz)")

# 3.3 COMPUTE BUOYANCY
println("\nComputing buoyancy from temperature...")
w_xy = snapshot_xy[:w]
T_xy = snapshot_xy[:T]
b_xy = compute!(Field(α_thermal * g_gravity * T_xy))

println("  w range: $(extrema(interior(w_xy)))")
println("  T range: $(extrema(interior(T_xy)))")
println("  b range: $(extrema(interior(b_xy)))")

# 3.4 COARSE-GRAIN
println("\nApplying coarse-graining filter...")
println("  Kernel: $FILTER_KERNEL")
println("  Cutoff: $(FILTER_CUTOFF/1e3) km")
println("  Border: $FILTER_BORDER")

w_bar_xy = ZFaceField(grid_xy, Float32)
b_bar_xy = CenterField(grid_xy, Float32)

t_start = time()
coarse_graining!(w_xy, w_bar_xy; kernel=FILTER_KERNEL, cutoff=FILTER_CUTOFF, border=FILTER_BORDER)
coarse_graining!(b_xy, b_bar_xy; kernel=FILTER_KERNEL, cutoff=FILTER_CUTOFF, border=FILTER_BORDER)
t_filter = time() - t_start
println("  Filtering completed in $(round(t_filter, digits=2)) seconds")

# 3.5 COMPUTE RESIDUALS
println("\nComputing residuals w' and b'...")
wp_full_xy = interior(w_xy) .- interior(w_bar_xy)
bp_full_xy = interior(b_xy) .- interior(b_bar_xy)

# 3.6 EXTRACT CORE REGION
halo_cells_xy = ceil(Int, XY_HALO / Δh)
println("\nExtracting core region (removing $halo_cells_xy halo cells per side)...")

# For ZFaceField, we have Nz+1 faces; for CenterField, Nz cells
# Since we loaded specific levels, the z-dimension is already small
wp_core_xy = extract_core_xy(wp_full_xy, halo_cells_xy)
bp_core_xy = extract_core_xy(bp_full_xy, halo_cells_xy)

# Get w' from all but the lowest level
# For discrete levels loaded via `levels`, w is set only in upper levels
# If w has one extra z-level, average them
wp_centered_xy = wp_core_xy[:, :, 2:end]
bp_centered_xy = bp_core_xy

Nx_core, Ny_core, Nz_core = size(wp_centered_xy)
println("  Core region size: $Nx_core × $Ny_core × $Nz_core")
println("  w' range: $(extrema(wp_centered_xy))")
println("  b' range: $(extrema(bp_centered_xy))")

# 3.7 BUILD FIGURE 1
println("\nGenerating Figure 1...")

set_theme!(theme_latexfonts(), fontsize=10, figure_padding=8)
fig1 = Figure(size = (700, 850))

# Panel labels
panel_labels_left = ["(a)", "(c)", "(e)"]
panel_labels_right = ["(b)", "(d)", "(f)"]

# Coordinate vectors for spatial plots (in km)
x_km = range(XY_CORE_XLIMS[1]/1e3, XY_CORE_XLIMS[2]/1e3, length=Nx_core)
y_km = range(XY_CORE_YLIMS[1]/1e3, XY_CORE_YLIMS[2]/1e3, length=Ny_core)

for (row, (z_idx, z_target, z_actual)) in enumerate(zip(z_indices, Z_TARGETS, actual_z))
    # Extract this level
    wp_level = wp_centered_xy[:, :, row]
    bp_level = bp_centered_xy[:, :, row]
    
    # Compute mask based on |w'b'| magnitude
    wb_mag = abs.(wp_level .* bp_level)
    wb_thresh = quantile(vec(wb_mag), THRESHOLD_PERCENTILE)
    mask = wb_mag .>= wb_thresh
    
    # Quadrant assignment
    Q_field = assign_quadrants(wp_level, bp_level, mask)
    
    # --- Left column: Histogram ---
    wp_sig = vec(wp_level)[vec(mask)]
    bp_sig = vec(bp_level)[vec(mask)]
    h = fit(Histogram, (wp_sig, bp_sig), nbins=N_BINS)
    
    w_lim = (first(h.edges[1]), last(h.edges[1]))
    b_lim = (first(h.edges[2]), last(h.edges[2]))
    
    ax_hist = Axis(fig1[row, 1]; 
        xlabel = row == 3 ? L"w'~\text{(m s}^{-1}\text{)}" : "",
        ylabel = L"b'~\text{(m s}^{-2}\text{)}",
        title = L"%$(panel_labels_left[row])~z = %$((z_target))~\text{m}",
        limits = (w_lim, b_lim))
    
    hm = heatmap!(ax_hist, h.edges[1], h.edges[2], log10.(1 .+ h.weights); 
                  colormap=Reverse(:grays), rasterize=true)
    add_wb_threshold_mask!(ax_hist, wb_thresh, w_lim, b_lim)
    Colorbar(fig1[row, 2], hm; label=L"\log_{10}(1+N)")
    
    # --- Right column: X-Y Pattern ---
    ax_xy = Axis(fig1[row, 3];
        xlabel = row == 3 ? L"x~\text{(km)}" : "",
        ylabel = L"y~\text{(km)}",
        title = L"%$(panel_labels_right[row])~\text{Quadrant pattern}",
        aspect = 1,
        limits = ((x_km[1], x_km[end]), (y_km[1], y_km[end])))
    
    # Alpha based on |w'b'| magnitude
    wb_ref = quantile(vec(wb_mag[mask]), 0.99)
    alpha_field = clamp.(wb_mag ./ wb_ref, 0, 1)
    rgba_data = create_quadrant_rgba(Q_field, alpha_field, QUADRANT_COLORS)
    
    image!(ax_xy, (x_km[1], x_km[end]), (y_km[1], y_km[end]), rgba_data; rasterize=true)
    
    # Hide x-decorations for non-bottom rows
    if row < 3
        hidexdecorations!(ax_hist, ticks=false)
        hidexdecorations!(ax_xy, ticks=false)
    end
    
    # Statistics
    n_sig = count(mask)
    n_total = length(mask)
    println("  Level $row (z=$((z_target))m): $(n_sig)/$(n_total) significant points ($(round(100*n_sig/n_total, digits=1))%)")
end

# Legend at top
legend_elements = [PolyElement(color=c) for c in QUADRANT_COLORS]
Legend(fig1[0, 3], legend_elements, QUADRANT_NAMES; 
       orientation=:horizontal, framevisible=false, 
       labelsize=9, patchsize=(12, 8), padding=(0, 0, 0, 0))

# Adjust layout
colgap!(fig1.layout, 2, 5)
colgap!(fig1.layout, 1, 10)
rowgap!(fig1.layout, 1, 8)
rowgap!(fig1.layout, 2, 8)
resize_to_layout!(fig1)

# Save Figure 1
fig1_pdf = OUTPUT_DIR * "figure1_xy_quadrants_iter$(ITERATION).pdf"
fig1_png = OUTPUT_DIR * "figure1_xy_quadrants_iter$(ITERATION).png"
save(fig1_pdf, fig1; pt_per_unit=1)
save(fig1_png, fig1; px_per_unit=4)
println("\nFigure 1 saved:")
println("  $fig1_pdf")
println("  $fig1_png")

# Clear memory before Figure 2
snapshot_xy = nothing
w_xy = nothing
T_xy = nothing
b_xy = nothing
w_bar_xy = nothing
b_bar_xy = nothing
wp_full_xy = nothing
bp_full_xy = nothing
GC.gc()

# ===============================================================================
# SECTION 4: FIGURE 2 - X-Z ANALYSIS AT Y-SLICES
# ===============================================================================

println("\n" * "="^70)
println("FIGURE 2: X-Z Quadrant Analysis at Y-Slices")
println("="^70)

# Storage for processed data from each y-slice
xz_data = []

for (i, y_target) in enumerate(Y_TARGETS)
    println("\n--- Processing Y-slice $i: y = $(y_target/1e3) km ---")
    
    # 4.1 LOAD narrow y-band
    y_lims = compute_y_band_limits(y_target, Y_BAND_HALF_WIDTH)
    x_lims_with_halo = (XZ_CORE_XLIMS[1] - XY_HALO, XZ_CORE_XLIMS[2] + XY_HALO)
    
    println("  Loading y-band: $y_lims")
    println("  X range (with halo): $x_lims_with_halo")
    println("  Z range: $XZ_ZLIMS")
    
    # Check if saved subdomain exists first
    xz_save_file = SUBDOMAIN_DIR * "subdomain_xz_y$(Int(y_target/1e3))km_iter$(ITERATION).jld2"
    
    if isfile(xz_save_file)
        println("  Loading from saved subdomain: $xz_save_file")
        snapshot_xz = load_subdomain_snapshot(xz_save_file; variables=("u", "v", "w", "T"))
    else
        println("  Loading from distributed checkpoint...")
        snapshot_xz = load_distributed_checkpoint_subdomain(CHECKPOINT_PREFIX, ITERATION;
            xlims = x_lims_with_halo,
            ylims = y_lims,
            zlims = XZ_ZLIMS,
            getEw = false,
            getMLD = 0
        )
        
        # Save for future runs
        if SAVE_LOADED_DATA
            println("  Saving loaded subdomain to: $xz_save_file")
            save_subdomain_with_halo(xz_save_file, snapshot_xz;
                core_xlims = XZ_CORE_XLIMS,
                core_ylims = y_lims,
                halo_width = XY_HALO,
                zlims = XZ_ZLIMS,
                iteration = ITERATION
            )
        end
    end
    
    grid_xz = snapshot_xz[:grid]
    println("  Loaded grid: $(grid_xz.Nx) × $(grid_xz.Ny) × $(grid_xz.Nz)")
    
    # 4.2 COMPUTE BUOYANCY
    w_xz = snapshot_xz[:w]
    T_xz = snapshot_xz[:T]
    b_xz = compute!(Field(α_thermal * g_gravity * T_xz))
    
    # 4.3 COARSE-GRAIN
    w_bar_xz = ZFaceField(grid_xz, Float32)
    b_bar_xz = CenterField(grid_xz, Float32)
    
    # Use ycircular border for narrow y-band
    coarse_graining!(w_xz, w_bar_xz; kernel=FILTER_KERNEL, cutoff=FILTER_CUTOFF, border=FILTER_BORDER_XZ)
    coarse_graining!(b_xz, b_bar_xz; kernel=FILTER_KERNEL, cutoff=FILTER_CUTOFF, border=FILTER_BORDER_XZ)
    
    # 4.4 COMPUTE RESIDUALS
    wp_full_xz = interior(w_xz) .- interior(w_bar_xz)
    bp_full_xz = interior(b_xz) .- interior(b_bar_xz)
    
    # 4.5 EXTRACT CORE (x-dimension only) and middle y-slice
    halo_cells_xz = ceil(Int, XY_HALO / Δh)
    j_mid = size(wp_full_xz, 2) ÷ 2 + 1  # Middle of the narrow y-band
    
    wp_slice = extract_core_x(wp_full_xz, halo_cells_xz)[:, j_mid, :]
    bp_slice = extract_core_x(bp_full_xz, halo_cells_xz)[:, j_mid, :]
    
    # Interpolate w to cell centers in z (w has Nz+1 faces, b has Nz cells)
    if size(wp_slice, 2) > size(bp_slice, 2)
        wp_centered = (wp_slice[:, 1:end-1] .+ wp_slice[:, 2:end]) ./ 2
    else
        wp_centered = wp_slice
    end
    bp_centered = bp_slice
    
    println("  Core slice size: $(size(wp_centered))")
    println("  w' range: $(extrema(wp_centered))")
    println("  b' range: $(extrema(bp_centered))")
    
    push!(xz_data, (wp=wp_centered, bp=bp_centered, y=y_target))
    
    # Clear this iteration's large arrays
    snapshot_xz = nothing
    w_xz = nothing
    T_xz = nothing
    b_xz = nothing
    w_bar_xz = nothing
    b_bar_xz = nothing
    GC.gc()
end

# 4.6 BUILD FIGURE 2
println("\nGenerating Figure 2...")

set_theme!(theme_latexfonts(), fontsize=10, figure_padding=8)
fig2 = Figure(size = (720, 700))

# Coordinate vectors
Nx_xz = size(xz_data[1].wp, 1)
Nz_xz = size(xz_data[1].wp, 2)
x_km_xz = range(XZ_CORE_XLIMS[1]/1e3, XZ_CORE_XLIMS[2]/1e3, length=Nx_xz)
z_m = range(XZ_ZLIMS[1] + Δz/2, XZ_ZLIMS[2] - Δz/2, length=Nz_xz)  # Cell centers

for (row, data) in enumerate(xz_data)
    wp_slice = data.wp
    bp_slice = data.bp
    y_km = data.y / 1e3
    
    # Compute mask
    wb_mag = abs.(wp_slice .* bp_slice)
    wb_thresh = quantile(vec(wb_mag), THRESHOLD_PERCENTILE)
    mask = wb_mag .>= wb_thresh
    
    # Quadrant assignment
    Q_field = assign_quadrants(wp_slice, bp_slice, mask)
    
    # --- Left column: Histogram ---
    wp_sig = vec(wp_slice)[vec(mask)]
    bp_sig = vec(bp_slice)[vec(mask)]
    h = fit(Histogram, (wp_sig, bp_sig), nbins=N_BINS)
    
    w_lim = (first(h.edges[1]), last(h.edges[1]))
    b_lim = (first(h.edges[2]), last(h.edges[2]))
    
    ax_hist = Axis(fig2[row, 1]; 
        xlabel = row == 3 ? L"w'~\text{(m s}^{-1}\text{)}" : "",
        ylabel = L"b'~\text{(m s}^{-2}\text{)}",
        title = L"%$(panel_labels_left[row])~y = %$(Int(y_km))~\text{km}",
        limits = (w_lim, b_lim))
    
    hm = heatmap!(ax_hist, h.edges[1], h.edges[2], log10.(1 .+ h.weights); 
                  colormap=Reverse(:grays), rasterize=true)
    add_wb_threshold_mask!(ax_hist, wb_thresh, w_lim, b_lim)
    Colorbar(fig2[row, 2], hm; label=L"\log_{10}(1+N)")
    
    # --- Right column: X-Z Pattern ---
    # Non-square panel with some vertical exaggeration
    # Actual ratio: 50km / 81m ≈ 617:1
    # With ~100x vertical exaggeration, display aspect ~6:1
    ax_xz = Axis(fig2[row, 3];
        xlabel = row == 3 ? L"x~\text{(km)}" : "",
        ylabel = L"z~\text{(m)}",
        title = L"%$(panel_labels_right[row])~\text{Quadrant pattern}",
        limits = ((x_km_xz[1], x_km_xz[end]), (z_m[1], z_m[end])))
    
    # Alpha based on |w'b'| magnitude
    wb_ref = quantile(vec(wb_mag[mask]), 0.99)
    alpha_field = clamp.(wb_mag ./ wb_ref, 0, 1)
    rgba_data = create_quadrant_rgba(Q_field, alpha_field, QUADRANT_COLORS)
    
    image!(ax_xz, (x_km_xz[1], x_km_xz[end]), (z_m[1], z_m[end]), rgba_data; rasterize=true)
    
    # Hide x-decorations for non-bottom rows
    if row < 3
        hidexdecorations!(ax_hist, ticks=false)
        hidexdecorations!(ax_xz, ticks=false)
    end
    
    # Statistics
    n_sig = count(mask)
    n_total = length(mask)
    println("  Slice $row (y=$(Int(y_km))km): $(n_sig)/$(n_total) significant points ($(round(100*n_sig/n_total, digits=1))%)")
end

# Legend at top
Legend(fig2[0, 3], legend_elements, QUADRANT_NAMES; 
       orientation=:horizontal, framevisible=false, 
       labelsize=9, patchsize=(12, 8), padding=(0, 0, 0, 0))

# Adjust layout
colgap!(fig2.layout, 2, 5)
colgap!(fig2.layout, 1, 10)
rowgap!(fig2.layout, 1, 8)
rowgap!(fig2.layout, 2, 8)
resize_to_layout!(fig2)

# Save Figure 2
fig2_pdf = OUTPUT_DIR * "figure2_xz_quadrants_iter$(ITERATION).pdf"
fig2_png = OUTPUT_DIR * "figure2_xz_quadrants_iter$(ITERATION).png"
save(fig2_pdf, fig2; pt_per_unit=1)
save(fig2_png, fig2; px_per_unit=4)
println("\nFigure 2 saved:")
println("  $fig2_pdf")
println("  $fig2_png")

# ===============================================================================
# SECTION 5: SUMMARY
# ===============================================================================

println("\n" * "="^70)
println("WORKFLOW COMPLETE")
println("="^70)

println("\nConfiguration used:")
println("  Checkpoint: iteration $ITERATION")
println("  Filter: $(FILTER_KERNEL) kernel, $(FILTER_CUTOFF/1e3) km cutoff")
println("  Masking: bottom $(Int(THRESHOLD_PERCENTILE*100))% of |w'b'| masked")

println("\nFigure 1 (X-Y at Z-levels):")
println("  Core domain: x ∈ $(XY_CORE_XLIMS./1e3) km, y ∈ $(XY_CORE_YLIMS./1e3) km")
println("  Z-levels: $Z_TARGETS m")

println("\nFigure 2 (X-Z at Y-slices):")
println("  Core domain: x ∈ $(XZ_CORE_XLIMS./1e3) km, z ∈ $XZ_ZLIMS m")
println("  Y-slices: $(Y_TARGETS./1e3) km")

println("\nOutput files:")
println("  $fig1_pdf")
println("  $fig1_png")
println("  $fig2_pdf")
println("  $fig2_png")

if SAVE_LOADED_DATA
    println("\nSaved subdomain files:")
    println("  $(SUBDOMAIN_DIR)subdomain_xy_levels_iter$(ITERATION).jld2")
    for y_target in Y_TARGETS
        println("  $(SUBDOMAIN_DIR)subdomain_xz_y$(Int(y_target/1e3))km_iter$(ITERATION).jld2")
    end
end
