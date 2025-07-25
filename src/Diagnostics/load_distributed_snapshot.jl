
function load_distributed_checkpoint(filename, iteration; 
                                     architecture = CPU(),
                                     metadata = nothing,
                                     level = nothing)

    snapshot = Dict()

    file = jldopen(filename * "0_iteration$(iteration).jld2")

    Px = file["NonhydrostaticModel/grid"].architecture.partition.x
    Py = file["NonhydrostaticModel/grid"].architecture.partition.y

    nx = file["NonhydrostaticModel/grid"].Nx
    ny = file["NonhydrostaticModel/grid"].Ny
    Nz = file["NonhydrostaticModel/grid"].Nz

    Hx = file["NonhydrostaticModel/grid"].Hx
    Hy = file["NonhydrostaticModel/grid"].Hy
    Hz = file["NonhydrostaticModel/grid"].Hz

    Nx = nx * Px
    Ny = ny * Py

    Lx = file["NonhydrostaticModel/grid"].Lx * Px
    Ly = file["NonhydrostaticModel/grid"].Ly * Py
    Lz = file["NonhydrostaticModel/grid"].Lz

    grid = RectilinearGrid(architecture; size = (Nx, Ny, Nz), extent = (Lx, Ly, Lz))

    indices = isnothing(level) ? (Colon(), Colon(), Colon()) : (Colon(), Colon(), UnitRange(level, level))

    u = XFaceField(grid; indices)
    v = YFaceField(grid; indices)
    w = ZFaceField(grid; indices)
    T = CenterField(grid; indices)

    close(file)

    for rank in 0 : (Px * Py-1)
        @info "loading rank $rank of $(Px * Py - 1)"

        file = jldopen(filename * "$(rank)_iteration$(iteration).jld2")

        Rx = file["NonhydrostaticModel/grid"].architecture.local_index[1]
        Ry = file["NonhydrostaticModel/grid"].architecture.local_index[2]

        udata = file["NonhydrostaticModel/u/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
        vdata = file["NonhydrostaticModel/v/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
        wdata = file["NonhydrostaticModel/w/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
        Tdata = file["NonhydrostaticModel/T/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
        
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

function load_distributed_snapshot(filename, iteration; 
                                   architecture = CPU(),
                                   metadata = nothing,
                                   level = nothing)

    snapshot = Dict()

    file = jldopen(filename * "0.jld2")

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

        file = jldopen(filename * "$(rank).jld2")

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
                                            metadata = nothing,
                                            xlims = nothing,
                                            ylims = nothing,
                                            zlims = nothing,
                                            levels = nothing,
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

    # Read metadata from rank 0 to understand full grid layout
    file = jldopen(filename * "0_iteration$(iteration).jld2")
    Px = file["NonhydrostaticModel/grid"].architecture.partition.x
    Py = file["NonhydrostaticModel/grid"].architecture.partition.y
    nx = file["NonhydrostaticModel/grid"].Nx # points per rank in x
    ny = file["NonhydrostaticModel/grid"].Ny # points per rank in y
    Nz = file["NonhydrostaticModel/grid"].Nz # total points in z
    Hx = file["NonhydrostaticModel/grid"].Hx
    Hy = file["NonhydrostaticModel/grid"].Hy
    Hz = file["NonhydrostaticModel/grid"].Hz
    
    # Full domain parameters
    Nx_full = nx * Px
    Ny_full = ny * Py
    Lx_full = file["NonhydrostaticModel/grid"].Lx * Px
    Ly_full = file["NonhydrostaticModel/grid"].Ly * Py
    Lz_full = file["NonhydrostaticModel/grid"].Lz
    
    # Grid spacing
    Δx = Lx_full / Nx_full
    Δy = Ly_full / Ny_full
    Δz = Lz_full / Nz
    
    close(file)

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
    x_topoloty = (xlims[2] - xlims[1]) > Lx_full - Δx/2 ? Periodic : Bounded
    y_topoloty = (ylims[2] - ylims[1]) > Ly_full - Δy/2 ? Periodic : Bounded
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
    u = XFaceField(grid; indices=field_indices)
    v = YFaceField(grid; indices=field_indices)
    w = ZFaceField(grid; indices=field_indices)
    T = CenterField(grid; indices=field_indices)
    if getMLD >= 1
        MLD = Field{Center, Center, Nothing}(grid; indices=(Colon(), Colon(), UnitRange(1, 1)))
    end
    if getMLD >= 2
        MLD2 = Field{Center, Center, Nothing}(grid; indices=(Colon(), Colon(), UnitRange(1, 1)))
    end
    if getMLD >= 3
        MLD3 = Field{Center, Center, Nothing}(grid; indices=(Colon(), Colon(), UnitRange(1, 1)))
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
            
            # Adjust to align with rank boundaries for efficient loading
            rank_x_start = div(x_start_idx - 1, nx) + 1
            rank_x_end = div(x_end_idx - 1, nx) + 1
            rank_y_start = div(y_start_idx - 1, ny) + 1
            rank_y_end = div(y_end_idx - 1, ny) + 1
            
            # Adjust grid indices to align with rank boundaries
            x_start_idx = (rank_x_start - 1) * nx + 1
            x_end_idx = rank_x_end * nx
            y_start_idx = (rank_y_start - 1) * ny + 1
            y_end_idx = rank_y_end * ny
            
            @info "Segment aligned indices: x=$x_start_idx:$x_end_idx, y=$y_start_idx:$y_end_idx"
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
                    # Calculate rank number (0-indexed)
                    rank = (Rx - 1) * Px + (Ry - 1)
                    @info "Loading from rank $rank (Rx=$Rx, Ry=$Ry)"
                    
                    file = jldopen(filename * "$(rank)_iteration$(iteration).jld2")
                    
                    # Get the actual local indices for this rank from the file
                    file_Rx = file["NonhydrostaticModel/grid"].architecture.local_index[1]
                    file_Ry = file["NonhydrostaticModel/grid"].architecture.local_index[2]
                    
                    # Load data from this rank (excluding halos)
                    udata = file["NonhydrostaticModel/u/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
                    vdata = file["NonhydrostaticModel/v/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
                    wdata = file["NonhydrostaticModel/w/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
                    Tdata = file["NonhydrostaticModel/T/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
                    
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
                        MLDdata = interior(compute!(MixedLayerDepth(grid_r, (; T=T_r); ΔT = abs(Δρ / ρ₀ / α))))
                    end
                    if getMLD >= 2
                        MLD2data = interior(compute!(MixedLayerDepth(grid_r, (; T=T_r); ΔT = abs(2Δρ / ρ₀ / α))))
                    end
                    if getMLD >= 3
                        MLD3data = interior(compute!(MixedLayerDepth(grid_r, (; T=T_r); ΔT = abs(3Δρ / ρ₀ / α))))
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
                    end
                    if getMLD >= 2
                        interior(MLD2, final_x_start:final_x_end, final_y_start:final_y_end, 1) .=
                            MLD2data[rank_x_start_local:rank_x_end_local, rank_y_start_local:rank_y_end_local, 1]
                    end
                    if getMLD >= 3
                        interior(MLD3, final_x_start:final_x_end, final_y_start:final_y_end, 1) .=
                            MLD3data[rank_x_start_local:rank_x_end_local, rank_y_start_local:rank_y_end_local, 1]
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
    end
    if getMLD >= 2
        fill_halo_regions!(MLD2)
        snapshot[:MLD2] = MLD2
    end
    if getMLD >= 3
        fill_halo_regions!(MLD3)
        snapshot[:MLD3] = MLD3
    end

    if !isnothing(metadata)
        params = jldopen(metadata)["parameters"]
        set_value!(params)
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

                    i = findfirst(levels .== level)
                    if isnothing(i)
                        @warn "Level $level not found in $filename."
                        continue
                    else
                        gridl = RectilinearGrid(grid.architecture;
                                                size = (grid.Nx, grid.Ny, 1),
                                                x = (-grid.Lx/2,grid.Lx/2),
                                                y = (0,grid.Ly),
                                                z = (-grid.Lz/grid.Nz,0),
                                                topology = (grid.Lx>parameters.Lx-1 ? Periodic : Bounded,grid.Ly>parameters.Ly-1 ? Periodic : Bounded,Bounded))
                        @info "Loading level $level with one-layer grid $gridl."
                        snapshot[:grid] = gridl
                        snapshot[:level] = level
                        data = T.(field_group["data"][:,:,i])
                        loc = field_group["location"]
                        loc3 = var=="w" ? Nothing : loc[3]
                        ind3 = var=="w" ? UnitRange(1, 1) : Colon()
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