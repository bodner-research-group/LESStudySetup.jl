
function load_distributed_checkpoint(filename, iteration; 
                                     architecture = CPU(),
                                     metadata = nothing,
                                     level = nothing)

    snapshot = Dict()

    file = jldopen(filename * "0_iteration$(iteration).jld2")

    Px = file["grid"].architecture.partition.x
    Py = file["grid"].architecture.partition.y

    nx = file["grid"].Nx
    ny = file["grid"].Ny
    Nz = file["grid"].Nz

    Hx = file["grid"].Hx
    Hy = file["grid"].Hy
    Hz = file["grid"].Hz

    Nx = nx * Px
    Ny = ny * Py

    Lx = file["grid"].Lx * Px
    Ly = file["grid"].Ly * Py
    Lz = file["grid"].Lz

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

        Rx = file["grid"].architecture.local_index[1]
        Ry = file["grid"].architecture.local_index[2]

        udata = file["u/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
        vdata = file["v/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
        wdata = file["w/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
        Tdata = file["T/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
        
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
                                   metadata = nothing)

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
        
    u = XFaceField(grid)
    v = YFaceField(grid)
    w = ZFaceField(grid)
    T = CenterField(grid)

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

        interior(u, irange, jrange, :) .= udata
        interior(v, irange, jrange, :) .= vdata
        interior(w, irange, jrange, :) .= wdata
        interior(T, irange, jrange, :) .= Tdata
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