
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

    for rank in 0 : (Px * Py-1)
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
    
    return snapshots
end