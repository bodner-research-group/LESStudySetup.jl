# ===============================================================================
#                    SUBDOMAIN XY LEVELS CORRECTION SCRIPT
# ===============================================================================
#
# One-time correction script to fix the subdomain_xy_levels file that was
# saved with zlims instead of levels metadata.
#
# The original file has:
#   - grid with Nz=72 (full depth)
#   - field data with shape (Nx, Ny, 3) (windowed to 3 levels)
#   - metadata/zlims = (-5.0, -70.5) ← WRONG
#
# This script creates a corrected file with:
#   - grid with Nz=3 (compact, matching data)
#   - field data with shape (Nx, Ny, 3)
#   - metadata/levels = [220, 193, 163] ← CORRECT
#
# Usage:
#   julia --project test/correct_subdomain_xy_levels.jl
#
# ===============================================================================

using JLD2
using Oceananigans
using Oceananigans: location

# --- Configuration (same as quadrant_analysis_publication.jl) ---
const CHECKPOINT_DIR = "/orcd/data/abodner/002/shared_datasets/nhyles_output/"
const SUBDOMAIN_DIR = CHECKPOINT_DIR * "publication_figures/subdomains/"
const ITERATION = 164410

const Δz = 1.125
const Lz = 252.0
const Z_TARGETS = [-5.0, -35.0, -70.5]

# Compute z-indices (same formula as quadrant_analysis_publication.jl)
z_indices = [round(Int, (z + Lz) / Δz) + 1 for z in Z_TARGETS]

# File paths
xy_original_file = SUBDOMAIN_DIR * "subdomain_xy_levels_iter$(ITERATION).jld2"
xy_corrected_file = SUBDOMAIN_DIR * "subdomain_xy_levels_iter$(ITERATION)_corrected.jld2"

println("=" ^ 70)
println("SUBDOMAIN XY LEVELS CORRECTION SCRIPT")
println("=" ^ 70)
println("\nConfiguration:")
println("  Original file: $xy_original_file")
println("  Corrected file: $xy_corrected_file")
println("  Z-targets: $Z_TARGETS")
println("  Z-indices: $z_indices")

# Verify original file exists
if !isfile(xy_original_file)
    error("Original file not found: $xy_original_file")
end

println("\n" * "-" ^ 70)
println("Reading original file...")
println("-" ^ 70)

# Open original and create corrected
jldopen(xy_original_file, "r") do old_file
    # Read original grid (has Nz=72)
    old_grid = old_file["grid"]
    println("\nOriginal grid: $(old_grid.Nx) × $(old_grid.Ny) × $(old_grid.Nz)")
    
    # Read field data (correct shape: Nx × Ny × 3)
    u_data = old_file["fields/u/data"]
    v_data = old_file["fields/v/data"]
    w_data = old_file["fields/w/data"]
    T_data = old_file["fields/T/data"]
    
    println("\nField data shapes:")
    println("  u: $(size(u_data))")
    println("  v: $(size(v_data))")
    println("  w: $(size(w_data))")
    println("  T: $(size(T_data))")
    
    # Determine Nz from data (should be 3)
    Nz_compact = size(T_data, 3)
    println("\nCompact Nz from data: $Nz_compact")
    
    if Nz_compact != length(z_indices)
        @warn "Data Nz ($Nz_compact) does not match expected levels ($(length(z_indices)))"
    end
    
    # Compute z-bounds for the discrete levels
    z_centers = [-Lz + (k - 1) * Δz + Δz/2 for k in z_indices]
    z_min = minimum(z_centers) - Δz/2
    z_max = maximum(z_centers) + Δz/2
    println("\nComputed z-range: ($z_min, $z_max)")
    println("  Z-centers: $z_centers")
    
    # Extract x and y extents from original grid
    x_min = old_grid.xᶜᵃᵃ[1] - old_grid.Δxᶜᵃᵃ/2
    x_max = old_grid.xᶜᵃᵃ[old_grid.Nx] + old_grid.Δxᶜᵃᵃ/2
    y_min = old_grid.yᵃᶜᵃ[1] - old_grid.Δyᵃᶜᵃ/2
    y_max = old_grid.yᵃᶜᵃ[old_grid.Ny] + old_grid.Δyᵃᶜᵃ/2
    
    println("\nX-extent: ($x_min, $x_max)")
    println("Y-extent: ($y_min, $y_max)")
    
    # Create compact grid with Nz = 3
    compact_grid = RectilinearGrid(CPU();
        size = (old_grid.Nx, old_grid.Ny, Nz_compact),
        x = (x_min, x_max),
        y = (y_min, y_max),
        z = (z_min, z_max),
        topology = (Bounded, Bounded, Bounded))
    
    println("\n" * "-" ^ 70)
    println("Creating compact grid and fields...")
    println("-" ^ 70)
    println("\nCompact grid: $(compact_grid.Nx) × $(compact_grid.Ny) × $(compact_grid.Nz)")
    
    # Create fields on compact grid
    u = XFaceField(compact_grid)
    v = YFaceField(compact_grid)
    w = ZFaceField(compact_grid)
    T = CenterField(compact_grid)
    
    interior(u) .= u_data
    interior(v) .= v_data
    interior(w)[:,:,2:end] .= w_data
    interior(T) .= T_data
    
    fill_halo_regions!(u)
    fill_halo_regions!(v)
    fill_halo_regions!(w)
    fill_halo_regions!(T)
    
    println("\nFields created and halos filled.")
    
    # Write corrected file
    println("\n" * "-" ^ 70)
    println("Writing corrected file...")
    println("-" ^ 70)
    
    jldopen(xy_corrected_file, "w") do new_file
        # Save compact grid
        new_file["grid"] = compact_grid
        
        # Save fields
        new_file["fields/u/data"] = Array(interior(u))
        new_file["fields/u/location"] = location(u)
        new_file["fields/v/data"] = Array(interior(v))
        new_file["fields/v/location"] = location(v)
        new_file["fields/w/data"] = Array(interior(w))
        new_file["fields/w/location"] = location(w)
        new_file["fields/T/data"] = Array(interior(T))
        new_file["fields/T/location"] = location(T)
        
        # Copy existing metadata (except zlims)
        if haskey(old_file, "metadata/core_xlims")
            new_file["metadata/core_xlims"] = old_file["metadata/core_xlims"]
            println("  Copied core_xlims: $(old_file["metadata/core_xlims"])")
        end
        if haskey(old_file, "metadata/core_ylims")
            new_file["metadata/core_ylims"] = old_file["metadata/core_ylims"]
            println("  Copied core_ylims: $(old_file["metadata/core_ylims"])")
        end
        if haskey(old_file, "metadata/halo_width")
            new_file["metadata/halo_width"] = old_file["metadata/halo_width"]
            println("  Copied halo_width: $(old_file["metadata/halo_width"])")
        end
        if haskey(old_file, "metadata/xlims")
            new_file["metadata/xlims"] = old_file["metadata/xlims"]
            println("  Copied xlims: $(old_file["metadata/xlims"])")
        end
        if haskey(old_file, "metadata/ylims")
            new_file["metadata/ylims"] = old_file["metadata/ylims"]
            println("  Copied ylims: $(old_file["metadata/ylims"])")
        end
        if haskey(old_file, "metadata/iteration")
            new_file["metadata/iteration"] = old_file["metadata/iteration"]
            println("  Copied iteration: $(old_file["metadata/iteration"])")
        end
        if haskey(old_file, "metadata/clock_time")
            new_file["metadata/clock_time"] = old_file["metadata/clock_time"]
            println("  Copied clock_time: $(old_file["metadata/clock_time"])")
        end
        if haskey(old_file, "metadata/clock_time_days")
            new_file["metadata/clock_time_days"] = old_file["metadata/clock_time_days"]
            println("  Copied clock_time_days: $(old_file["metadata/clock_time_days"])")
        end
        
        # NOTE: Intentionally NOT copying zlims since we're replacing with levels
        if haskey(old_file, "metadata/zlims")
            println("  SKIPPED zlims (replacing with levels): $(old_file["metadata/zlims"])")
        end
        
        # Save CORRECT levels metadata (NOT zlims)
        new_file["metadata/levels"] = z_indices
        println("  Added levels: $z_indices")
    end
end

println("\n" * "=" ^ 70)
println("CORRECTION COMPLETE")
println("=" ^ 70)
println("\nCorrected file saved to:")
println("  $xy_corrected_file")
println("\nTo use the corrected file, either:")
println("  1. Update quadrant_analysis_publication.jl to load from the corrected file")
println("  2. Or replace the original (backup recommended first):")
println("     mv $xy_corrected_file $xy_original_file")
println()
