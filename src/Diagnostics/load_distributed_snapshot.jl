
function load_distributed_checkpoint(filename, iteration;
                                     architecture = CPU(),
                                     partition,
                                     metadata = nothing,
                                     level = nothing)

    # Recent Oceananigans checkpoints store only the prognostic state, not the grid, so the
    # global grid is reconstructed from `parameters` exactly as in `idealized_setup`. When a
    # metadata file is provided its parameters are applied first, before the grid is built.
    if !isnothing(metadata)
        params = jldopen(metadata)["parameters"]
        set_value!(params)
    end

    snapshot = Dict()

    Px = partition.x
    Py = partition.y

    Nx = ceil(Int, parameters.Lx / parameters.Δh)
    Ny = ceil(Int, parameters.Ly / parameters.Δh)
    Nz = ceil(Int, parameters.Lz / parameters.Δz)

    grid = RectilinearGrid(architecture;
                           size = (Nx, Ny, Nz),
                           x = (0, parameters.Lx),
                           y = (0, parameters.Ly),
                           z = (-parameters.Lz, 0),
                           halo = (6, 6, 6))

    nx = Nx ÷ Px
    ny = Ny ÷ Py

    indices = isnothing(level) ? (Colon(), Colon(), Colon()) : (Colon(), Colon(), UnitRange(level, level))

    u = XFaceField(grid; indices)
    v = YFaceField(grid; indices)
    w = ZFaceField(grid; indices)
    T = CenterField(grid; indices)

    for rank in 0 : (Px * Py - 1)
        @info "loading rank $rank of $(Px * Py - 1)"

        Rx = div(rank, Py) + 1
        Ry = mod(rank, Py) + 1

        file = jldopen(filename * "_rank$(rank)_iteration$(iteration).jld2")

        # Halos are inferred from the center field, whose interior is exactly (nx, ny, Nz).
        Tdata = file["simulation/model/tracers/T/data"]
        Hx = (size(Tdata, 1) - nx) ÷ 2
        Hy = (size(Tdata, 2) - ny) ÷ 2
        Hz = (size(Tdata, 3) - Nz) ÷ 2

        udata = file["simulation/model/velocities/u/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
        vdata = file["simulation/model/velocities/v/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
        wdata = file["simulation/model/velocities/w/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
        Tdata = Tdata[Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]

        close(file)

        irange = 1 + (Rx - 1) * nx : Rx * nx
        jrange = 1 + (Ry - 1) * ny : Ry * ny

        interior(u, irange, jrange, :) .= udata[:, :, indices[3]]
        interior(v, irange, jrange, :) .= vdata[:, :, indices[3]]
        interior(w, irange, jrange, :) .= wdata[:, :, indices[3]]
        interior(T, irange, jrange, :) .= Tdata[:, :, indices[3]]
    end

    snapshot[:u] = u
    snapshot[:v] = v
    snapshot[:w] = w
    snapshot[:T] = T

    return snapshot
end

"""
    load_checkpoint_clock(filename_prefix, iteration)

Extract clock information (time and iteration) from a checkpoint file.

This function reads the simulation clock from a distributed checkpoint,
providing the simulation time and iteration number. Useful for displaying
metadata about when a snapshot was saved.

Arguments:
- `filename_prefix`: Path prefix to checkpoint files (without rank suffix).
                     Example: "/path/to/checkpoint" for files like "checkpoint_rank0_iteration1000.jld2"
- `iteration`: Checkpoint iteration number

Returns:
- NamedTuple with fields:
  - `time`: Simulation time in seconds
  - `iteration`: Model iteration number
  - `time_days`: Simulation time converted to days

Example:
```julia
clock = load_checkpoint_clock("/path/to/checkpoint", 32207)
@info "Simulation time: \$(clock.time_days) days (iteration \$(clock.iteration))"
```
"""
function load_checkpoint_clock(filename_prefix, iteration)
    filepath = filename_prefix * "_rank0_iteration$(iteration).jld2"
    
    jldopen(filepath, "r") do file
        clock = file["simulation/model/clock"]
        time_seconds = clock.time
        iter = clock.iteration
        time_days = time_seconds / 86400.0  # Convert to days
        return (time = time_seconds, iteration = iter, time_days = time_days)
    end
end

"""
    compute_subdomain_tiles(; Lx=100e3, Ly=100e3, tile_size=10e3, halo_width=600.0)

Compute tile boundaries for dividing a domain into subdomains with halo overlap.

The domain is divided into a grid of tiles. Each tile has:
- A "core" region: the actual tile (e.g., 10km × 10km)
- A "full" region: core + halo band on all sides (e.g., 11.2km × 11.2km)

The halo band ensures that coarse-graining operations near tile boundaries
don't suffer from edge artifacts. After filtering, only the core region
contains valid data.

Arguments:
- `Lx`: Full domain x-extent in meters (default: 100km)
- `Ly`: Full domain y-extent in meters (default: 100km)
- `tile_size`: Core tile size in meters (default: 10km)
- `halo_width`: Halo band width in meters on each side (default: 600m)

Returns:
- Vector of NamedTuples, each containing:
  - `tile_id`: Integer tile identifier (1-indexed, row-major order)
  - `tile_ix`, `tile_iy`: Tile indices in x and y directions
  - `core_xlims`, `core_ylims`: Core region coordinate bounds
  - `full_xlims`, `full_ylims`: Full region bounds (core + halo)

Example:
```julia
tiles = compute_subdomain_tiles(; Lx=100e3, Ly=100e3, tile_size=10e3, halo_width=600.0)
# Returns 100 tiles (10×10 grid)
# tiles[1] = (tile_id=1, tile_ix=1, tile_iy=1, 
#             core_xlims=(0.0, 10000.0), core_ylims=(0.0, 10000.0),
#             full_xlims=(-600.0, 10600.0), full_ylims=(-600.0, 10600.0))
```

Note: Negative coordinates and coordinates exceeding domain size are handled
by periodic wrapping in `load_distributed_checkpoint_subdomain`.
"""
function compute_subdomain_tiles(;
    Lx::Real = 100e3,
    Ly::Real = 100e3,
    tile_size::Real = 10e3,
    halo_width::Real = 600.0)
    
    # Compute number of tiles in each direction
    n_tiles_x = round(Int, Lx / tile_size)
    n_tiles_y = round(Int, Ly / tile_size)
    
    # Validate that domain divides evenly
    if abs(n_tiles_x * tile_size - Lx) > 1e-6
        @warn "Domain Lx=$Lx does not divide evenly by tile_size=$tile_size"
    end
    if abs(n_tiles_y * tile_size - Ly) > 1e-6
        @warn "Domain Ly=$Ly does not divide evenly by tile_size=$tile_size"
    end
    
    tiles = Vector{NamedTuple}()
    tile_id = 0
    
    for ix in 1:n_tiles_x
        for iy in 1:n_tiles_y
            tile_id += 1
            
            # Core region bounds
            x0 = (ix - 1) * tile_size
            x1 = ix * tile_size
            y0 = (iy - 1) * tile_size
            y1 = iy * tile_size
            
            core_xlims = (x0, x1)
            core_ylims = (y0, y1)
            
            # Full region bounds (halo extends into neighboring tiles)
            # Note: negative coordinates and coordinates > domain size
            # are handled by periodic wrapping in load_distributed_checkpoint_subdomain
            full_xlims = (x0 - halo_width, x1 + halo_width)
            full_ylims = (y0 - halo_width, y1 + halo_width)
            
            push!(tiles, (;
                tile_id,
                tile_ix = ix,
                tile_iy = iy,
                core_xlims,
                core_ylims,
                full_xlims,
                full_ylims
            ))
        end
    end
    
    @info "Created $(length(tiles)) tiles: $(n_tiles_x)×$(n_tiles_y) grid"
    @info "  Core tile size: $(tile_size/1e3) km × $(tile_size/1e3) km"
    @info "  Halo width: $(halo_width) m"
    @info "  Full tile size: $((tile_size + 2*halo_width)/1e3) km × $((tile_size + 2*halo_width)/1e3) km"
    
    return tiles
end

function load_distributed_snapshot(filename, iteration; 
                                   architecture = CPU(),
                                   metadata = nothing,
                                   level = nothing)

    snapshot = Dict()

    file = jldopen(filename * "_rank0.jld2")

    Px = file["grid/architecture/partition/x"]
    Py = file["grid/architecture/partition/y"]

    nx = file["grid/Nx"]
    ny = file["grid/Ny"]
    Nz = file["grid/Nz"]

    Nx = nx * Px
    Ny = ny * Py

    Lx = file["grid/Lx"] * Px
    Ly = file["grid/Ly"] * Py
    Lz = file["grid/Lz"]

    grid = RectilinearGrid(architecture; size = (Nx, Ny, Nz), extent = (Lx, Ly, Lz))

    indices = isnothing(level) ? (Colon(), Colon(), Colon()) : (Colon(), Colon(), UnitRange(level, level))
    
    u = XFaceField(grid; indices)
    v = YFaceField(grid; indices)
    w = ZFaceField(grid; indices)
    T = CenterField(grid; indices)

    close(file)

    for rank in 0 : (Px * Py - 1)
        @info "loading rank $rank of $(Px * Py - 1)"

        file = jldopen(filename * "_rank$(rank).jld2")

        Rx = file["grid/architecture/local_index/1"]
        Ry = file["grid/architecture/local_index/2"]

        udata = file["timeseries/u/" * iteration]
        vdata = file["timeseries/v/" * iteration]
        wdata = file["timeseries/w/" * iteration]
        Tdata = file["timeseries/T/" * iteration]
        
        irange = 1 + (Rx - 1) * nx : Rx * nx
        jrange = 1 + (Ry - 1) * ny : Ry * ny

	interior(u, irange, jrange, :) .= udata[:, :, indices[3]] 
        interior(v, irange, jrange, :) .= vdata[:, :, indices[3]]
        interior(w, irange, jrange, :) .= wdata[:, :, indices[3]]
        interior(T, irange, jrange, :) .= Tdata[:, :, indices[3]]
    end

    snapshot[:u] = u
    snapshot[:v] = v
    snapshot[:w] = w
    snapshot[:T] = T

    if !isnothing(metadata)
        params = jldopen(metadata)["parameters"]
        set_value!(params)
    end
    
    return snapshot
end

function load_distributed_checkpoint_subdomain(filename, iteration;
                                               architecture = CPU(),
                                               partition,
                                               metadata = nothing,
                                               xlims = nothing,
                                               ylims = nothing,
                                               zlims = nothing,
                                               levels = nothing,
                                               getEw = false,
                                               getMLD = 0, Δρ = 0.03)

    # Helper function to handle periodic coordinate normalization
    function normalize_periodic_coords(coord_min, coord_max, domain_size)
        if coord_min > -1 && coord_max < domain_size+1
            # No wraparound: single segment
            return [(coord_min, coord_max)]
        end
        # Normalize coordinates to [0, domain_size) range
        coord_min_norm = mod(coord_min, domain_size)
        coord_max_norm = mod(coord_max, domain_size)
        
        # Handle wraparound case
        if coord_min_norm+1 > coord_max_norm
            # Domain wraps around: split into two segments
            return [(coord_min_norm, domain_size), (0.0, coord_max_norm)]
        else
            # Normal case: single segment
            return [(coord_min_norm, coord_max_norm)]
        end
    end

    # Recent Oceananigans checkpoints store only the prognostic state, so the full-domain
    # layout is reconstructed from `parameters` (as in `idealized_setup`) plus the supplied
    # `partition`. Per-rank halos are inferred from each rank's center field inside the loop.
    if !isnothing(metadata)
        params = jldopen(metadata)["parameters"]
        set_value!(params)
    end

    Px = partition.x
    Py = partition.y

    # Full domain parameters
    Nx_full = ceil(Int, parameters.Lx / parameters.Δh)
    Ny_full = ceil(Int, parameters.Ly / parameters.Δh)
    Nz = ceil(Int, parameters.Lz / parameters.Δz) # total points in z
    Lx_full = parameters.Lx
    Ly_full = parameters.Ly
    Lz_full = parameters.Lz

    nx = Nx_full ÷ Px # points per rank in x
    ny = Ny_full ÷ Py # points per rank in y

    # Grid spacing
    Δx = Lx_full / Nx_full
    Δy = Ly_full / Ny_full
    Δz = Lz_full / Nz

    # Handle default limits (full domain)
    if isnothing(xlims)
        xlims = (0.0, Lx_full)
    end
    if isnothing(ylims)
        ylims = (0.0, Ly_full)
    end
    if isnothing(zlims)
        zlims = (-Lz_full, 0.0)
    end

    # NEW: Handle periodic coordinates - get segments that may wrap around
    x_segments = normalize_periodic_coords(xlims[1], xlims[2], Lx_full)
    y_segments = normalize_periodic_coords(ylims[1], ylims[2], Ly_full)

    # Calculate total subdomain size (based on requested range, not normalized)
    Nx_sub = round(Int, (xlims[2] - xlims[1]) / Δx)
    Ny_sub = round(Int, (ylims[2] - ylims[1]) / Δy)
    
    # Handle z-limits (unchanged)
    z_start_idx = max(1, round(Int, (zlims[1] + Lz_full) / Δz) + 1)
    z_end_idx = min(Nz, round(Int, (zlims[2] + Lz_full) / Δz))
    Nz_sub = z_end_idx - z_start_idx + 1

    # Create subdomain grid with REQUESTED coordinates (not normalized)
    x_topoloty = (xlims[2] - xlims[1]) > Lx_full - Δx/2 ? Oceananigans.Grids.Periodic : Bounded
    y_topoloty = (ylims[2] - ylims[1]) > Ly_full - Δy/2 ? Oceananigans.Grids.Periodic : Bounded
    grid_topology = (x_topoloty, y_topoloty, Bounded)
    grid = RectilinearGrid(architecture;
                           size = (Nx_sub, Ny_sub, Nz_sub),
                           x = xlims,
                           y = ylims,
                           z = zlims,
                           topology = grid_topology)

    @info "Created subdomain grid:"
    @info " Size: ($Nx_sub, $Ny_sub, $Nz_sub)"
    @info " X extent: $(xlims)"
    @info " Y extent: $(ylims)"
    @info " Z extent: $(zlims)"

    # Determine vertical indices for field creation and data extraction
    if isnothing(levels)
        field_indices = (Colon(), Colon(), Colon())
        data_z_range = z_start_idx:z_end_idx
        data_z_range_w = z_start_idx:(z_end_idx+1)
    else
        field_indices = (Colon(), Colon(), UnitRange(1, length(levels)))
        data_z_range = levels
        data_z_range_w = levels
    end

    # Create fields on subdomain grid
    u =  XFaceField(grid; indices=field_indices)
    v =  YFaceField(grid; indices=field_indices)
    w =  ZFaceField(grid; indices=field_indices)
    T = CenterField(grid; indices=field_indices)
    if getMLD >= 1
        MLD = Field{Center, Center, Nothing}(grid)
        if getEw
            Ew = Field{Center, Center, Nothing}(grid)
        end
    end
    if getMLD >= 2
        MLD2 = Field{Center, Center, Nothing}(grid)
        if getEw
            Ew2 = Field{Center, Center, Nothing}(grid)
        end
    end
    if getMLD >= 3
        MLD3 = Field{Center, Center, Nothing}(grid)
        if getEw
            Ew3 = Field{Center, Center, Nothing}(grid)
        end
    end

    # NEW: Load data from all segments (handling periodic wraparound)
    @info "Loading data from $(length(x_segments)) x-segments and $(length(y_segments)) y-segments"
    
    for (seg_x_idx, (x_min_seg, x_max_seg)) in enumerate(x_segments)
        for (seg_y_idx, (y_min_seg, y_max_seg)) in enumerate(y_segments)
            
            @info "Processing segment: x ∈ ($x_min_seg, $x_max_seg), y ∈ ($y_min_seg, $y_max_seg)"
            
            # Convert segment coordinates to grid indices
            x_start_idx = max(1, round(Int, x_min_seg / Δx) + 1)
            x_end_idx = min(Nx_full, round(Int, x_max_seg / Δx))
            y_start_idx = max(1, round(Int, y_min_seg / Δy) + 1)
            y_end_idx = min(Ny_full, round(Int, y_max_seg / Δy))
            
            # Calculate which ranks overlap with requested subdomain
            rank_x_start = div(x_start_idx - 1, nx) + 1
            rank_x_end = div(x_end_idx - 1, nx) + 1
            rank_y_start = div(y_start_idx - 1, ny) + 1
            rank_y_end = div(y_end_idx - 1, ny) + 1
            
            @info "Segment indices: x=$x_start_idx:$x_end_idx, y=$y_start_idx:$y_end_idx"
            @info "Loading from ranks: Rx=$rank_x_start:$rank_x_end, Ry=$rank_y_start:$rank_y_end"
            
            # Calculate where this segment maps to in the OUTPUT grid
            # Key insight: map normalized segment coordinates to requested coordinate space
            if seg_x_idx == 1
                # First x-segment starts at beginning of output
                out_x_start = 1
                out_x_end = round(Int, (x_max_seg - x_min_seg) / Δx)
            else
                # Second x-segment (wraparound case)
                prev_seg_width = round(Int, (x_segments[1][2] - x_segments[1][1]) / Δx)
                out_x_start = prev_seg_width + 1
                out_x_end = prev_seg_width + round(Int, (x_max_seg - x_min_seg) / Δx)
            end
            
            if seg_y_idx == 1
                # First y-segment starts at beginning of output
                out_y_start = 1
                out_y_end = round(Int, (y_max_seg - y_min_seg) / Δy)
            else
                # Second y-segment (wraparound case)
                prev_seg_width = round(Int, (y_segments[1][2] - y_segments[1][1]) / Δy)
                out_y_start = prev_seg_width + 1
                out_y_end = prev_seg_width + round(Int, (y_max_seg - y_min_seg) / Δy)
            end
            
            @info "Output mapping: x=$out_x_start:$out_x_end, y=$out_y_start:$out_y_end"
            
            # Load data from necessary ranks for this segment
            for Rx in rank_x_start:rank_x_end
                for Ry in rank_y_start:rank_y_end
                    # Rank number from the (Rx, Ry) tile position (0-indexed, y varies fastest)
                    rank = (Rx - 1) * Py + (Ry - 1)
                    file_Rx = Rx
                    file_Ry = Ry
                    @info "Loading from rank $rank (Rx=$Rx, Ry=$Ry)"

                    file = jldopen(filename * "_rank$(rank)_iteration$(iteration).jld2")

                    # Load data from this rank (excluding halos); halos inferred from the center field
                    Tdata = file["simulation/model/tracers/T/data"]
                    Hx = (size(Tdata, 1) - nx) ÷ 2
                    Hy = (size(Tdata, 2) - ny) ÷ 2
                    Hz = (size(Tdata, 3) - Nz) ÷ 2

                    udata = file["simulation/model/velocities/u/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
                    vdata = file["simulation/model/velocities/v/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
                    wdata = file["simulation/model/velocities/w/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
                    Tdata = Tdata[Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]

                    close(file)

                    if getMLD >= 1
                        grid_r = RectilinearGrid(architecture;
                                                 size = (nx, ny, Nz),
                                                 x = (0, Δx*nx),
                                                 y = (0, Δy*ny),
                                                 z = (-Lz_full, 0.0),
                                                 topology = (Bounded, Bounded, Bounded))
                        T_r = CenterField(grid_r)
                        interior(T_r) .= Tdata
                        α, ρ₀ = parameters.α, parameters.ρ₀
                        if !getEw
                            MLDdata = interior(MixedLayerDepth(grid_r, (; T=T_r); ΔT = abs(Δρ / ρ₀ / α)))
                        else
                            w_r = ZFaceField(grid_r)
                            interior(w_r) .= wdata
                            (MLDdata, Ewdata) = MixedLayerDepth(grid_r, (; T=T_r); w=w_r, ΔT = abs(Δρ / ρ₀ / α))
                            MLDdata, Ewdata = interior(MLDdata), interior(Ewdata)
                        end
                    end
                    if getMLD >= 2
                        if !getEw
                            MLD2data = interior(MixedLayerDepth(grid_r, (; T=T_r); ΔT = abs(2Δρ / ρ₀ / α)))
                        else
                            (MLD2data, Ew2data) = MixedLayerDepth(grid_r, (; T=T_r); w=w_r, ΔT = abs(2Δρ / ρ₀ / α))
                            MLD2data, Ew2data = interior(MLD2data), interior(Ew2data)
                        end
                    end
                    if getMLD >= 3
                        if !getEw
                            MLD3data = interior(MixedLayerDepth(grid_r, (; T=T_r); ΔT = abs(3Δρ / ρ₀ / α)))
                        else
                            (MLD3data, Ew3data) = MixedLayerDepth(grid_r, (; T=T_r); w=w_r, ΔT = abs(3Δρ / ρ₀ / α))
                            MLD3data, Ew3data = interior(MLD3data), interior(Ew3data)
                        end
                    end
                    
                    # Calculate global index ranges for this rank's data
                    rank_x_global_start = 1 + (file_Rx - 1) * nx
                    rank_x_global_end = file_Rx * nx
                    rank_y_global_start = 1 + (file_Ry - 1) * ny
                    rank_y_global_end = file_Ry * ny
                    
                    # Calculate overlap with our target segment
                    x_overlap_start = max(rank_x_global_start, x_start_idx)
                    x_overlap_end = min(rank_x_global_end, x_end_idx)
                    y_overlap_start = max(rank_y_global_start, y_start_idx)
                    y_overlap_end = min(rank_y_global_end, y_end_idx)
                    
                    # Skip if no overlap
                    if x_overlap_start > x_overlap_end || y_overlap_start > y_overlap_end
                        @info "No overlap with global x=$rank_x_global_start:$rank_x_global_end, y=$rank_y_global_start:$rank_y_global_end"
                        continue
                    end
                    
                    # Calculate indices within this rank's data array
                    rank_x_start_local = x_overlap_start - rank_x_global_start + 1
                    rank_x_end_local = x_overlap_end - rank_x_global_start + 1
                    rank_y_start_local = y_overlap_start - rank_y_global_start + 1
                    rank_y_end_local = y_overlap_end - rank_y_global_start + 1
                    
                    # Calculate indices within THIS SEGMENT of the output
                    seg_x_start = x_overlap_start - x_start_idx + 1
                    seg_x_end = x_overlap_end - x_start_idx + 1
                    seg_y_start = y_overlap_start - y_start_idx + 1
                    seg_y_end = y_overlap_end - y_start_idx + 1
                    
                    # Map segment indices to final output indices
                    final_x_start = out_x_start + seg_x_start - 1
                    final_x_end = out_x_start + seg_x_end - 1
                    final_y_start = out_y_start + seg_y_start - 1
                    final_y_end = out_y_start + seg_y_end - 1
                    
                    @info "Copying data: rank local ($rank_x_start_local:$rank_x_end_local, $rank_y_start_local:$rank_y_end_local) → output ($final_x_start:$final_x_end, $final_y_start:$final_y_end)"
                    
                    # Copy data with proper indexing
                    # if isnothing(levels)
                    # Full vertical range
                    interior(u, final_x_start:final_x_end, final_y_start:final_y_end, :) .=
                        udata[rank_x_start_local:rank_x_end_local, rank_y_start_local:rank_y_end_local, data_z_range]
                    interior(v, final_x_start:final_x_end, final_y_start:final_y_end, :) .=
                        vdata[rank_x_start_local:rank_x_end_local, rank_y_start_local:rank_y_end_local, data_z_range]
                    interior(w, final_x_start:final_x_end, final_y_start:final_y_end, :) .=
                        wdata[rank_x_start_local:rank_x_end_local, rank_y_start_local:rank_y_end_local, data_z_range_w]
                    interior(T, final_x_start:final_x_end, final_y_start:final_y_end, :) .=
                        Tdata[rank_x_start_local:rank_x_end_local, rank_y_start_local:rank_y_end_local, data_z_range]
                    if getMLD >= 1
                        interior(MLD, final_x_start:final_x_end, final_y_start:final_y_end, 1) .=
                            MLDdata[rank_x_start_local:rank_x_end_local, rank_y_start_local:rank_y_end_local, 1]
                        if getEw
                            interior(Ew, final_x_start:final_x_end, final_y_start:final_y_end, 1) .=
                                Ewdata[rank_x_start_local:rank_x_end_local, rank_y_start_local:rank_y_end_local, 1]
                        end
                    end
                    if getMLD >= 2
                        interior(MLD2, final_x_start:final_x_end, final_y_start:final_y_end, 1) .=
                            MLD2data[rank_x_start_local:rank_x_end_local, rank_y_start_local:rank_y_end_local, 1]
                        if getEw
                            interior(Ew2, final_x_start:final_x_end, final_y_start:final_y_end, 1) .=
                                Ew2data[rank_x_start_local:rank_x_end_local, rank_y_start_local:rank_y_end_local, 1]
                        end
                    end
                    if getMLD >= 3
                        interior(MLD3, final_x_start:final_x_end, final_y_start:final_y_end, 1) .=
                            MLD3data[rank_x_start_local:rank_x_end_local, rank_y_start_local:rank_y_end_local, 1]
                        if getEw
                            interior(Ew3, final_x_start:final_x_end, final_y_start:final_y_end, 1) .=
                                Ew3data[rank_x_start_local:rank_x_end_local, rank_y_start_local:rank_y_end_local, 1]
                        end
                    end
                    # else
                    #     # Vertical levels
                    #     interior(u, final_x_start:final_x_end, final_y_start:final_y_end, levels) .=
                    #         udata[rank_x_start_local:rank_x_end_local, rank_y_start_local:rank_y_end_local, levels]
                    #     interior(v, final_x_start:final_x_end, final_y_start:final_y_end, levels) .=
                    #         vdata[rank_x_start_local:rank_x_end_local, rank_y_start_local:rank_y_end_local, levels]
                    #     interior(w, final_x_start:final_x_end, final_y_start:final_y_end, levels) .=
                    #         wdata[rank_x_start_local:rank_x_end_local, rank_y_start_local:rank_y_end_local, levels]
                    #     interior(T, final_x_start:final_x_end, final_y_start:final_y_end, levels) .=
                    #         Tdata[rank_x_start_local:rank_x_end_local, rank_y_start_local:rank_y_end_local, levels]
                    # end
                end
            end
        end
    end

    fill_halo_regions!(u)
    fill_halo_regions!(v)
    fill_halo_regions!(w)
    fill_halo_regions!(T)
    snapshot = Dict()
    snapshot[:u] = u
    snapshot[:v] = v
    snapshot[:w] = w
    snapshot[:T] = T
    snapshot[:grid] = grid
    if getMLD >= 1
        fill_halo_regions!(MLD)
        snapshot[:MLD] = MLD
        if getEw
            fill_halo_regions!(Ew)
            snapshot[:Ew] = Ew
        end
    end
    if getMLD >= 2
        fill_halo_regions!(MLD2)
        snapshot[:MLD2] = MLD2
        if getEw
            fill_halo_regions!(Ew2)
            snapshot[:Ew2] = Ew2
        end
    end
    if getMLD >= 3
        fill_halo_regions!(MLD3)
        snapshot[:MLD3] = MLD3
        if getEw
            fill_halo_regions!(Ew3)
            snapshot[:Ew3] = Ew3
        end
    end

    return snapshot
end

"""
    load_subdomain_snapshot(filename)

Load a snapshot of subdomain fields from a JLD2 file that was previously
saved by the `load_subdomain` workflow.

This function reconstructs the `Oceananigans.Field` objects from the saved
grid, data arrays, and location metadata, making them ready for analysis
and visualization.

Arguments:
- `filename`: The path to the JLD2 file containing the subdomain snapshot.

Returns:
- A `Dict` containing the reconstructed fields (e.g., `:u`, `:v`, `:T`) and the `grid`.
"""
function load_subdomain_snapshot(filename; T=Float32, variables = ("u", "v", "w", "T"), level = nothing)
    # Create an empty dictionary to store the results
    snapshot = Dict{Symbol, Any}()

    # Open the JLD2 file in read-only mode
    jldopen(filename, "r") do file
        # 1. Load the grid object. This is essential for reconstructing fields.
        grid = file["grid"]

        @info "Loaded grid: $grid"

        if "levels" in keys(file["metadata"])
            levels = file["metadata/levels"]
            
            @info "Found levels: $levels"
        end
        
        # Load halo-aware metadata if present (from save_subdomain_with_halo)
        if haskey(file, "metadata/core_xlims")
            snapshot[:core_xlims] = Tuple(file["metadata/core_xlims"])
            @info "Found core_xlims: $(snapshot[:core_xlims])"
        end
        if haskey(file, "metadata/core_ylims")
            snapshot[:core_ylims] = Tuple(file["metadata/core_ylims"])
            @info "Found core_ylims: $(snapshot[:core_ylims])"
        end
        if haskey(file, "metadata/halo_width")
            snapshot[:halo_width] = file["metadata/halo_width"]
            @info "Found halo_width: $(snapshot[:halo_width])"
        end
        if haskey(file, "metadata/clock_time")
            snapshot[:clock_time] = file["metadata/clock_time"]
            @info "Found clock_time: $(snapshot[:clock_time]) seconds"
        end
        if haskey(file, "metadata/clock_time_days")
            snapshot[:clock_time_days] = file["metadata/clock_time_days"]
            @info "Found clock_time_days: $(snapshot[:clock_time_days]) days"
        end
        if haskey(file, "metadata/iteration")
            snapshot[:iteration] = file["metadata/iteration"]
            @info "Found iteration: $(snapshot[:iteration])"
        end
        
        # 2. Iterate through the saved fields group to find all field names.
        field_names = keys(file["fields"])
        @info "Found fields: $field_names"
 
        for var in variables
            if var in field_names
                @info "Loading field $var."
                field_symbol = Symbol(var)
                field_group = file["fields/$field_symbol"]

                # 3. For each field, read its raw data array and its location.
                # 4. Reconstruct the Field object on the grid at the correct location.
                if isnothing(level)
                    snapshot[:grid] = grid
                    data = T.(field_group["data"])
                    loc = field_group["location"]
                    field = Field{loc[1], loc[2], loc[3]}(grid, T)
                else
                    i = "levels" in keys(file["metadata"]) ? findfirst(levels .== level) : level

                    if isnothing(i)
                        @warn "Level $level not found in $filename."
                        continue
                    else
                        gridl = RectilinearGrid(grid.architecture,T;
                                                size = (grid.Nx, grid.Ny, 1),
                                                x = (-grid.Lx/2,grid.Lx/2),
                                                y = (0,grid.Ly),
                                                z = (-grid.Lz/grid.Nz,0),
                                                topology = (grid.Lx>parameters.Lx-1 ? Oceananigans.Grids.Periodic : Bounded,grid.Ly>parameters.Ly-1 ? Oceananigans.Grids.Periodic : Bounded,Bounded))
                        @info "Loading level $level with one-layer grid $gridl."
                        snapshot[:grid] = gridl
                        snapshot[:level] = level
                        data = var in ("u", "v", "w", "T") ? T.(field_group["data"][:,:,i]) : T.(field_group["data"][:,:,1])
                        loc = field_group["location"]
                        loc3 = var=="w" ? Nothing : loc[3]
                        ind3 = var in ("u", "v", "T") ? Colon() : UnitRange(1, 1)
                        field = Field{loc[1], loc[2], loc3}(gridl, T; indices=(Colon(), Colon(), ind3))
                    end
                end

                # 5. Fill the interior of the newly created field with the loaded data.
                interior(field) .= data
                fill_halo_regions!(field)

                # 6. Store the fully reconstructed, usable field in the snapshot dictionary.
                snapshot[field_symbol] = field
            else
                @warn "Field $var not found in $filename."
            end
        end
    end

    @info "Snapshot successfully loaded."
    return snapshot
end

"""
    extract_subdomain(snapshot; architecture=CPU(), xlims=nothing, ylims=nothing, zlims=nothing)

Extract a subdomain from a full-domain snapshot based on coordinate limits.

This function takes a snapshot dictionary containing fields (`:u`, `:v`, `:w`, `:T`)
and a `:grid`, and extracts a subdomain defined by the specified coordinate limits.

Arguments:
- `snapshot`: A Dict containing fields and grid from a full-domain simulation.
- `architecture`: The architecture for the new subdomain grid (default: CPU()).
- `xlims`: Tuple (xmin, xmax) specifying x-coordinate limits. Default: full x-extent.
- `ylims`: Tuple (ymin, ymax) specifying y-coordinate limits. Default: full y-extent.
- `zlims`: Tuple (zmin, zmax) specifying z-coordinate limits. Default: full z-extent.

Returns:
- A new Dict containing subdomain fields (`:u`, `:v`, `:w`, `:T`) and `:grid`.
"""
function extract_subdomain(snapshot;
                           architecture = CPU(),
                           xlims = nothing,
                           ylims = nothing,
                           zlims = nothing)

    # Get source grid
    source_grid = snapshot[:grid]

    # Get full domain extents
    Lx_full = source_grid.Lx
    Ly_full = source_grid.Ly
    Lz_full = source_grid.Lz
    Nx_full = source_grid.Nx
    Ny_full = source_grid.Ny
    Nz_full = source_grid.Nz

    # Compute grid spacing
    Δx = Lx_full / Nx_full
    Δy = Ly_full / Ny_full
    Δz = Lz_full / Nz_full

    # Get source domain origin (for grids that don't start at 0)
    x_origin = source_grid.xᶜᵃᵃ[1] - Δx/2
    y_origin = source_grid.yᵃᶜᵃ[1] - Δy/2
    z_origin = source_grid.z.cᵃᵃᶜ[1] - Δz/2

    # Handle default limits (full domain)
    if isnothing(xlims)
        xlims = (x_origin, x_origin + Lx_full)
    end
    if isnothing(ylims)
        ylims = (y_origin, y_origin + Ly_full)
    end
    if isnothing(zlims)
        zlims = (z_origin, z_origin + Lz_full)
    end

    # Convert coordinate limits to index ranges (1-indexed)
    # For cell-centered fields, find indices where cell centers fall within limits
    x_start_idx = max(1, floor(Int, (xlims[1] - x_origin) / Δx) + 1)
    x_end_idx = min(Nx_full, floor(Int, (xlims[2] - x_origin) / Δx))
    y_start_idx = max(1, floor(Int, (ylims[1] - y_origin) / Δy) + 1)
    y_end_idx = min(Ny_full, floor(Int, (ylims[2] - y_origin) / Δy))
    z_start_idx = max(1, floor(Int, (zlims[1] - z_origin) / Δz) + 1)
    z_end_idx = min(Nz_full, floor(Int, (zlims[2] - z_origin) / Δz))

    # Calculate subdomain size
    Nx_sub = x_end_idx - x_start_idx + 1
    Ny_sub = y_end_idx - y_start_idx + 1
    Nz_sub = z_end_idx - z_start_idx + 1

    # Calculate actual subdomain coordinates
    x_sub_min = x_origin + (x_start_idx - 1) * Δx
    x_sub_max = x_origin + x_end_idx * Δx
    y_sub_min = y_origin + (y_start_idx - 1) * Δy
    y_sub_max = y_origin + y_end_idx * Δy
    z_sub_min = z_origin + (z_start_idx - 1) * Δz
    z_sub_max = z_origin + z_end_idx * Δz

    # Determine topology for subdomain
    # Use Bounded topology for subdomain (it's a slice of the original domain)
    sub_topology = (Bounded, Bounded, Bounded)

    # Create subdomain grid
    sub_grid = RectilinearGrid(architecture;
                               size = (Nx_sub, Ny_sub, Nz_sub),
                               x = (x_sub_min, x_sub_max),
                               y = (y_sub_min, y_sub_max),
                               z = (z_sub_min, z_sub_max),
                               topology = sub_topology)

    @info "Extracting subdomain:"
    @info "  Source grid: ($Nx_full, $Ny_full, $Nz_full)"
    @info "  Subdomain grid: ($Nx_sub, $Ny_sub, $Nz_sub)"
    @info "  Index ranges: x=$x_start_idx:$x_end_idx, y=$y_start_idx:$y_end_idx, z=$z_start_idx:$z_end_idx"

    # Create output snapshot
    sub_snapshot = Dict{Symbol, Any}()
    sub_snapshot[:grid] = sub_grid

    # Extract each field
    field_configs = [
        (:u, XFaceField),
        (:v, YFaceField),
        (:w, ZFaceField),
        (:T, CenterField)
    ]

    for (field_name, FieldType) in field_configs
        if haskey(snapshot, field_name)
            source_field = snapshot[field_name]

            # Create new field on subdomain grid
            sub_field = FieldType(sub_grid)

            # Determine ranges for this field type
            # Face-centered fields with Bounded topology have one extra point in their dimension
            # XFaceField: +1 in x, YFaceField: +1 in y, ZFaceField: +1 in z
            if field_name == :u
                # XFaceField: need one extra x-point for Bounded topology
                x_range_src = x_start_idx:(x_end_idx + 1)
                y_range_src = y_start_idx:y_end_idx
                z_range_src = z_start_idx:z_end_idx
            elseif field_name == :v
                # YFaceField: need one extra y-point for Bounded topology
                x_range_src = x_start_idx:x_end_idx
                y_range_src = y_start_idx:(y_end_idx + 1)
                z_range_src = z_start_idx:z_end_idx
            elseif field_name == :w
                # ZFaceField: need one extra z-point for Bounded topology
                x_range_src = x_start_idx:x_end_idx
                y_range_src = y_start_idx:y_end_idx
                z_range_src = z_start_idx:(z_end_idx + 1)
            else
                # CenterField: no extra points
                x_range_src = x_start_idx:x_end_idx
                y_range_src = y_start_idx:y_end_idx
                z_range_src = z_start_idx:z_end_idx
            end

            # Copy interior data - interior(sub_field) already has correct size
            interior(sub_field) .=
                interior(source_field, x_range_src, y_range_src, z_range_src)

            fill_halo_regions!(sub_field)
            sub_snapshot[field_name] = sub_field
        end
    end

    @info "Subdomain extraction complete."
    return sub_snapshot
end

"""
    save_subdomain_snapshot(filename, snapshot; iteration=nothing, xlims=nothing, ylims=nothing, zlims=nothing, levels=nothing)

Save a subdomain snapshot to a JLD2 file in a format compatible with `load_subdomain_snapshot`.

This function saves a snapshot dictionary containing fields and grid to a JLD2 file,
along with optional metadata about the iteration and coordinate limits.

Arguments:
- `filename`: Path to the output JLD2 file.
- `snapshot`: A Dict containing fields (`:u`, `:v`, `:w`, `:T`, etc.) and `:grid`.
- `iteration`: (Optional) Iteration number to store in metadata.
- `xlims`: (Optional) Tuple of x-coordinate limits to store in metadata.
- `ylims`: (Optional) Tuple of y-coordinate limits to store in metadata.
- `zlims`: (Optional) Tuple of z-coordinate limits to store in metadata.
- `levels`: (Optional) Array of vertical levels to store in metadata.

The file structure matches what `load_subdomain_snapshot` expects:
- `grid`: The grid object
- `fields/<name>/data`: Interior data as Array
- `fields/<name>/location`: Field location tuple
- `metadata/...`: Optional metadata (iteration, coordinate limits)
"""
function save_subdomain_snapshot(filename, snapshot;
                                 iteration = nothing,
                                 xlims = nothing,
                                 ylims = nothing,
                                 zlims = nothing,
                                 levels = nothing)

    # Ensure the snapshot has a grid
    if !haskey(snapshot, :grid)
        error("Snapshot must contain a :grid key.")
    end

    sub_grid = snapshot[:grid]

    @info "Saving subdomain to $filename..."

    jldopen(filename, "w") do file
        # Save the grid
        file["grid"] = sub_grid

        # Save each field's interior data and location metadata
        for field_name in keys(snapshot)
            if field_name != :grid
                field = snapshot[field_name]

                # Convert to standard Array on the CPU before saving
                field_data = Array(interior(field))

                file["fields/$field_name/data"] = field_data
                file["fields/$field_name/location"] = location(field)
            end
        end

        # Save metadata
        if !isnothing(iteration)
            file["metadata/iteration"] = iteration
        end
        if !isnothing(xlims)
            file["metadata/xlims"] = xlims
        end
        if !isnothing(ylims)
            file["metadata/ylims"] = ylims
        end
        if !isnothing(zlims)
            file["metadata/zlims"] = zlims
        end
        if !isnothing(levels)
            file["metadata/levels"] = levels
        end
    end

    @info "Successfully saved subdomain data."
    return nothing
end

"""
    save_subdomain_with_halo(filename, snapshot; kwargs...)

Save a subdomain snapshot with halo band metadata for coarse-graining workflows.

This function extends `save_subdomain_snapshot` by recording both the core region
(where filtered output is valid) and the full region (including halo buffer).
The halo band provides padding for spatial filtering operations, ensuring that
edge artifacts don't contaminate the core region.

Arguments:
- `filename`: Output JLD2 file path
- `snapshot`: Dict containing fields (`:u`, `:v`, `:w`, `:T`, etc.) and `:grid`

Keyword Arguments:
- `core_xlims`: Tuple (xmin, xmax) for core region x-bounds (required)
- `core_ylims`: Tuple (ymin, ymax) for core region y-bounds (required)
- `halo_width`: Width of halo band in meters (required)
- `zlims`: Tuple (zmin, zmax) for vertical bounds (optional)
- `iteration`: Checkpoint iteration number (optional)
- `clock_time`: Simulation time in seconds (optional)
- `clock_time_days`: Simulation time in days (optional)

File Structure:
```
├── grid                    # Oceananigans grid object
├── fields/
│   ├── u/data, u/location
│   ├── v/data, v/location  
│   ├── w/data, w/location
│   └── T/data, T/location
└── metadata/
    ├── core_xlims          # Valid region after filtering
    ├── core_ylims
    ├── halo_width          # Halo band width in meters
    ├── xlims               # Full region (core + halo)
    ├── ylims
    ├── zlims
    ├── iteration
    ├── clock_time          # Time in seconds
    └── clock_time_days     # Time in days
```

Example:
```julia
# Save a tile with 600m halo for 300m Gaussian filter
save_subdomain_with_halo("tile_1.jld2", snapshot;
    core_xlims = (0.0, 10000.0),
    core_ylims = (0.0, 10000.0),
    halo_width = 600.0,
    zlims = (-81.0, 0.0),
    iteration = 32207,
    clock_time_days = 3.5)
```
"""
function save_subdomain_with_halo(filename, snapshot;
    core_xlims::Tuple{Real,Real},
    core_ylims::Tuple{Real,Real},
    halo_width::Real,
    zlims::Union{Nothing, Tuple{Real,Real}} = nothing,
    levels::Union{Nothing, AbstractVector{<:Integer}} = nothing,
    iteration::Union{Nothing, Integer} = nothing,
    clock_time::Union{Nothing, Real} = nothing,
    clock_time_days::Union{Nothing, Real} = nothing)
    
    # Mutual exclusivity check
    if !isnothing(zlims) && !isnothing(levels)
        error("Cannot specify both `zlims` and `levels`. Use one or the other.")
    end
    
    # Warning if neither specified
    if isnothing(zlims) && isnothing(levels)
        @warn "Neither `zlims` nor `levels` specified. Z-dimension metadata will not be saved."
    end
    
    # Validate snapshot has required keys
    if !haskey(snapshot, :grid)
        error("Snapshot must contain a :grid key.")
    end
    
    sub_grid = snapshot[:grid]
    
    # Compute full region from core + halo
    full_xlims = (core_xlims[1] - halo_width, core_xlims[2] + halo_width)
    full_ylims = (core_ylims[1] - halo_width, core_ylims[2] + halo_width)
    
    @info "Saving subdomain with halo to $filename..."
    @info "  Core region: x=$(core_xlims), y=$(core_ylims)"
    @info "  Full region: x=$(full_xlims), y=$(full_ylims)"
    @info "  Halo width: $(halo_width) m"
    
    jldopen(filename, "w") do file
        # Save grid
        file["grid"] = sub_grid
        
        # Save each field's interior data and location
        for field_name in keys(snapshot)
            if field_name != :grid && field_name != :clock
                field = snapshot[field_name]
                field_data = Array(interior(field))
                
                file["fields/$field_name/data"] = field_data
                file["fields/$field_name/location"] = location(field)
            end
        end
        
        # Save halo-aware metadata
        file["metadata/core_xlims"] = core_xlims
        file["metadata/core_ylims"] = core_ylims
        file["metadata/halo_width"] = halo_width
        file["metadata/xlims"] = full_xlims
        file["metadata/ylims"] = full_ylims
        
        # Optional metadata
        if !isnothing(zlims)
            file["metadata/zlims"] = zlims
        end
        if !isnothing(levels)
            file["metadata/levels"] = collect(levels)
        end
        if !isnothing(iteration)
            file["metadata/iteration"] = iteration
        end
        if !isnothing(clock_time)
            file["metadata/clock_time"] = clock_time
        end
        if !isnothing(clock_time_days)
            file["metadata/clock_time_days"] = clock_time_days
        end
    end
    
    @info "Successfully saved subdomain with halo."
    return nothing
end
