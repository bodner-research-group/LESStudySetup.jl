
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

using Oceananigans, JLD2
function load_distributed_checkpoint_subdomain(filename, iteration;
                                             architecture = CPU(),
                                             metadata = nothing,
                                             xlims = nothing,
                                             ylims = nothing, 
                                             zlims = nothing,
                                             level = nothing)
    
    # Read metadata from rank 0 to understand full grid layout
    file = jldopen(filename * "0_iteration$(iteration).jld2")
    
    Px = file["NonhydrostaticModel/grid"].architecture.partition.x
    Py = file["NonhydrostaticModel/grid"].architecture.partition.y
    
    nx = file["NonhydrostaticModel/grid"].Nx  # points per rank in x
    ny = file["NonhydrostaticModel/grid"].Ny  # points per rank in y
    Nz = file["NonhydrostaticModel/grid"].Nz  # total points in z
    
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
    
    # Convert physical limits to grid indices
    x_start_idx = max(1, round(Int, xlims[1] / Δx) + 1)
    x_end_idx = min(Nx_full, round(Int, xlims[2] / Δx))
    y_start_idx = max(1, round(Int, ylims[1] / Δy) + 1)
    y_end_idx = min(Ny_full, round(Int, ylims[2] / Δy))
    z_start_idx = max(1, round(Int, (zlims[1] + Lz_full) / Δz) + 1)
    z_end_idx = min(Nz, round(Int, (zlims[2] + Lz_full) / Δz))
    
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
    
    # Calculate subdomain grid size
    Nx_sub = x_end_idx - x_start_idx + 1
    Ny_sub = y_end_idx - y_start_idx + 1
    Nz_sub = isnothing(level) ? (z_end_idx - z_start_idx + 1) : 1
    
    # Calculate actual physical coordinates for subdomain boundaries
    x_min_actual = (x_start_idx - 1) * Δx
    x_max_actual = x_end_idx * Δx  
    y_min_actual = (y_start_idx - 1) * Δy
    y_max_actual = y_end_idx * Δy
    z_min_actual = -Lz_full + (z_start_idx - 1) * Δz
    z_max_actual = isnothing(level) ? (-Lz_full + z_end_idx * Δz) : (-Lz_full + z_start_idx * Δz)
    
    # Create subdomain grid - use Bounded topology for subdomains
    grid_topology = (Bounded, Bounded, Bounded)
    
    grid = RectilinearGrid(architecture; 
                          size = (Nx_sub, Ny_sub, Nz_sub),
                          x = (x_min_actual, x_max_actual),
                          y = (y_min_actual, y_max_actual), 
                          z = (z_min_actual, z_max_actual),
                          topology = grid_topology)
    
    @info "Created subdomain grid:"
    @info "  Size: ($Nx_sub, $Ny_sub, $Nz_sub)"
    @info "  X extent: ($x_min_actual, $x_max_actual)"
    @info "  Y extent: ($y_min_actual, $y_max_actual)" 
    @info "  Z extent: ($z_min_actual, $z_max_actual)"
    
    # Determine vertical indices for field creation and data extraction
    if isnothing(level)
        field_indices = (Colon(), Colon(), Colon())
        data_z_range = z_start_idx:z_end_idx
        data_z_range_w = z_start_idx:(z_end_idx+1)
    else
        field_indices = (Colon(), Colon(), UnitRange(1, 1))
        data_z_range = level:level
    end
    
    # Create fields on subdomain grid
    u = XFaceField(grid; indices=field_indices)
    v = YFaceField(grid; indices=field_indices)
    w = ZFaceField(grid; indices=field_indices)
    T = CenterField(grid; indices=field_indices)
    
    # Load data from necessary ranks only
    @info "Loading subdomain from ranks covering x-ranks $rank_x_start:$rank_x_end, y-ranks $rank_y_start:$rank_y_end"
    
    for Rx in rank_x_start:rank_x_end
        for Ry in rank_y_start:rank_y_end
            # Calculate rank number (0-indexed) - matches original function logic
            rank = (Ry - 1) * Px + (Rx - 1)
            
            @info "Loading from rank $rank (Rx=$Rx, Ry=$Ry)"
            
            file = jldopen(filename * "$(rank)_iteration$(iteration).jld2")
            
            # Get the actual local indices for this rank from the file (like original function)
            file_Rx = file["NonhydrostaticModel/grid"].architecture.local_index[1]
            file_Ry = file["NonhydrostaticModel/grid"].architecture.local_index[2]
            
            # Load data from this rank (excluding halos like in original)
            udata = file["NonhydrostaticModel/u/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
            vdata = file["NonhydrostaticModel/v/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
            wdata = file["NonhydrostaticModel/w/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
            Tdata = file["NonhydrostaticModel/T/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
            
            close(file)
            
            # Calculate global index ranges for this rank's data (like original function)
            rank_x_global_start = 1 + (file_Rx - 1) * nx
            rank_x_global_end = file_Rx * nx
            rank_y_global_start = 1 + (file_Ry - 1) * ny  
            rank_y_global_end = file_Ry * ny
            
            # Calculate overlap with our target subdomain
            x_overlap_start = max(rank_x_global_start, x_start_idx)
            x_overlap_end = min(rank_x_global_end, x_end_idx)
            y_overlap_start = max(rank_y_global_start, y_start_idx)
            y_overlap_end = min(rank_y_global_end, y_end_idx)
            
            # Skip if no overlap
            if x_overlap_start > x_overlap_end || y_overlap_start > y_overlap_end
                continue
            end
            
            # Calculate indices within this rank's data array
            rank_x_start_local = x_overlap_start - rank_x_global_start + 1
            rank_x_end_local = x_overlap_end - rank_x_global_start + 1  
            rank_y_start_local = y_overlap_start - rank_y_global_start + 1
            rank_y_end_local = y_overlap_end - rank_y_global_start + 1
            
            # Calculate indices within subdomain 
            sub_x_start = x_overlap_start - x_start_idx + 1
            sub_x_end = x_overlap_end - x_start_idx + 1
            sub_y_start = y_overlap_start - y_start_idx + 1  
            sub_y_end = y_overlap_end - y_start_idx + 1
            
            @info "Rank $rank: global range ($rank_x_global_start:$rank_x_global_end, $rank_y_global_start:$rank_y_global_end)"
            @info "  -> local slice ($rank_x_start_local:$rank_x_end_local, $rank_y_start_local:$rank_y_end_local)" 
            @info "  -> subdomain slice ($sub_x_start:$sub_x_end, $sub_y_start:$sub_y_end)"
            
            # Copy data with proper indexing
            if isnothing(level)
                # Full vertical range
                interior(u, sub_x_start:sub_x_end, sub_y_start:sub_y_end, :) .= 
                    udata[rank_x_start_local:rank_x_end_local, rank_y_start_local:rank_y_end_local, data_z_range .- Hz]
                interior(v, sub_x_start:sub_x_end, sub_y_start:sub_y_end, :) .= 
                    vdata[rank_x_start_local:rank_x_end_local, rank_y_start_local:rank_y_end_local, data_z_range .- Hz]
                interior(w, sub_x_start:sub_x_end, sub_y_start:sub_y_end, :) .= 
                    wdata[rank_x_start_local:rank_x_end_local, rank_y_start_local:rank_y_end_local, data_z_range_w .- Hz]
                interior(T, sub_x_start:sub_x_end, sub_y_start:sub_y_end, :) .= 
                    Tdata[rank_x_start_local:rank_x_end_local, rank_y_start_local:rank_y_end_local, data_z_range .- Hz]
            else
                # Single vertical level (like original function)
                z_level_local = level - Hz
                interior(u, sub_x_start:sub_x_end, sub_y_start:sub_y_end, 1) .= 
                    udata[rank_x_start_local:rank_x_end_local, rank_y_start_local:rank_y_end_local, z_level_local]
                interior(v, sub_x_start:sub_x_end, sub_y_start:sub_y_end, 1) .= 
                    vdata[rank_x_start_local:rank_x_end_local, rank_y_start_local:rank_y_end_local, z_level_local]
                interior(w, sub_x_start:sub_x_end, sub_y_start:sub_y_end, 1) .= 
                    wdata[rank_x_start_local:rank_x_end_local, rank_y_start_local:rank_y_end_local, z_level_local]
                interior(T, sub_x_start:sub_x_end, sub_y_start:sub_y_end, 1) .= 
                    Tdata[rank_x_start_local:rank_x_end_local, rank_y_start_local:rank_y_end_local, z_level_local]
            end
        end
    end
    
    snapshot = Dict()
    snapshot[:u] = u
    snapshot[:v] = v  
    snapshot[:w] = w
    snapshot[:T] = T
    snapshot[:grid] = grid
    
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
function load_subdomain_snapshot(filename)
    # Create an empty dictionary to store the results
    snapshot = Dict{Symbol, Any}()

    # Open the JLD2 file in read-only mode
    jldopen(filename, "r") do file
        # 1. Load the grid object. This is essential for reconstructing fields.
        grid = file["grid"]
        snapshot[:grid] = grid

        @info "Loaded grid: $grid"

        # 2. Iterate through the saved fields group to find all field names.
        field_names = keys(file["fields"])
        @info "Found fields: $field_names"

        for name in field_names
            field_symbol = Symbol(name)
            field_group = file["fields/$name"]

            # 3. For each field, read its raw data array and its location.
            data = field_group["data"]
            loc = field_group["location"]

            # 4. Reconstruct the Field object on the grid at the correct location.
            field = Field(loc, grid)

            # 5. Fill the interior of the newly created field with the loaded data.
            interior(field) .= data

            # 6. Store the fully reconstructed, usable field in the snapshot dictionary.
            snapshot[field_symbol] = field
        end
    end

    @info "Snapshot successfully loaded."
    return snapshot
end