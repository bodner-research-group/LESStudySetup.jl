using FFTW
using Oceananigans.Utils
using Oceananigans.Units
using Oceananigans.BoundaryConditions
using Oceananigans.Operators
using Oceananigans.Fields: instantiated_location
using KernelAbstractions: @kernel, @index
using ImageFiltering

@kernel function _horizontal_box_filter!(new_field, grid, field)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        new_field[i, j, k] = field[i, j, k]
        nn = (field[i, j, k], field[i + 1, j, k], field[i - 1, j, k], field[i, j + 1, k], field[i, j - 1, k], field[i + 1, j + 1, k], field[i - 1, j + 1, k], field[i + 1, j - 1, k], field[i - 1, j - 1, k])    
        new_field[i, j, k] = sum(nn) / 9
    end
end

@kernel function _horizontal_gauss_filter!(new_field, grid, field)
    i, j, k = @index(Global, NTuple)
    loc  = instantiated_location(field)
    ℑxy₁ = get_first_interpolation(loc)
    ℑxy₂ = get_second_interpolation(loc)

    @inbounds new_field[i, j, k] = ℑxy₂(i, j, k, grid, ℑxy₁, field)
end

@kernel function _top_hat_filter!(new_field, field, G, Ng)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        nn = field[i-Ng÷2:i+Ng÷2, j-Ng÷2:j+Ng÷2, k].*G    
        new_field[i, j, k] = sum(nn)
    end
end

@inline get_first_interpolation(::Tuple{<:Face,   <:Face,   <:Any}) = ℑxyᶜᶜᵃ
@inline get_first_interpolation(::Tuple{<:Center, <:Face,   <:Any}) = ℑxyᶠᶜᵃ
@inline get_first_interpolation(::Tuple{<:Face,   <:Center, <:Any}) = ℑxyᶜᶠᵃ
@inline get_first_interpolation(::Tuple{<:Center, <:Center, <:Any}) = ℑxyᶠᶠᵃ

@inline get_second_interpolation(::Tuple{<:Face,   <:Face,   <:Any}) = ℑxyᶠᶠᵃ
@inline get_second_interpolation(::Tuple{<:Center, <:Face,   <:Any}) = ℑxyᶜᶠᵃ
@inline get_second_interpolation(::Tuple{<:Face,   <:Center, <:Any}) = ℑxyᶠᶜᵃ
@inline get_second_interpolation(::Tuple{<:Center, <:Center, <:Any}) = ℑxyᶜᶜᵃ

function spatial_filtering(u::Field; 
                           smoothing_range = 20kilometer, 
                           kernel! = _horizontal_box_filter!)

    Δx = u.grid.Δxᶜᵃᵃ
    Δy = u.grid.Δyᵃᶜᵃ

    Δs = sqrt(Δy * Δx)

    u̅₁ = deepcopy(u)
    u̅₂ = deepcopy(u)

    iterations = ceil(Int, smoothing_range / Δs) ÷ 2

    arch = architecture(u)
    grid = u.grid

    for iter in 1:iterations
        if isodd(iter)
            launch!(arch, grid, :xyz, kernel!, u̅₁, grid, u̅₂)
            fill_halo_regions!(u̅₁)
        else
            launch!(arch, grid, :xyz, kernel!, u̅₂, grid, u̅₁)
            fill_halo_regions!(u̅₂)
        end
    end

    return ifelse(isodd(iterations), u̅₁, u̅₂)
end

function spatial_convolution!(u::Field, u̅₁::Field; 
    smoothing_range = 1kilometer, 
    kernel! = _top_hat_filter!)

    Δx = u.grid.Δxᶜᵃᵃ
    Δy = u.grid.Δyᵃᶜᵃ

    Δs = sqrt(Δy * Δx)
    Ng = ceil(Int, smoothing_range / Δs)
    if Ng % 2 == 0
    Ng += 1
    end

    xg, yg = -Ng÷2:Ng÷2, -Ng÷2:Ng÷2
    xg, yg = Δx*repeat(xg, 1, Ng), Δy*repeat(yg', Ng, 1)
    rg = sqrt.(xg.^2 + yg.^2)
    G = zeros(Float32, Ng, Ng)
    G[rg .< smoothing_range/2] .= 1
    G ./= sum(G)

    arch = architecture(u)
    grid = u.grid

    launch!(arch, grid, :xyz, kernel!, u̅₁, u, G, Ng)
    fill_halo_regions!(u̅₁)

    return nothing
end

@kernel function _corse_grain!(d̂, nx, ny, nx2, ny2, Nx, Ny)
    i, j, k = @index(Global, NTuple)
    i′ = ifelse(i > Nx ÷ 2, i - Nx ÷ 2, Nx ÷ 2 - i)
    j′ = ifelse(j > Ny ÷ 2, j - Ny ÷ 2, Ny ÷ 2 - j)

    # remove the outside frame till nx and ny
    outside_frame   = (i′ <= nx) | (j′ <= ny)
    smooth_region_x = (i′ > nx)  & (i′ <= nx2)
    smooth_region_y = (j′ > ny)  & (j′ <= ny2)

    only_smooth_x = smooth_region_x & !smooth_region_y
    only_smooth_y = smooth_region_y & !smooth_region_x
    smooth_both   = smooth_region_x & smooth_region_y

    scaling_x = (i′ - nx) / (nx2 - nx)
    scaling_y = (j′ - ny) / (ny2 - ny)
    scaling_b = scaling_x * scaling_y
    
    scaling = ifelse(outside_frame, 0, 
              ifelse(only_smooth_x, scaling_x,
              ifelse(only_smooth_y, scaling_y, 
              ifelse(smooth_both, scaling_b, 1))))

    @inbounds d̂[i, j, k] = scaling * d̂[i, j, k]
end

function spectral_filtering(u::Field; xcutoff = 20kilometer, ycutoff = 20kilometer)

    Δx = u.grid.Δxᶜᵃᵃ
    Δy = u.grid.Δyᵃᶜᵃ
    
    nx  = ceil(Int, xcutoff / Δx)
    ny  = ceil(Int, ycutoff / Δy)

    nx2 = ceil(Int, xcutoff / Δx * π / 2) 
    ny2 = ceil(Int, xcutoff / Δx * π / 2) 

    u̅ = deepcopy(u)
    d = interior(u)

    d̂ = fft(d)

    Nx, Ny, _ = size(u)

    arch = architecture(u)
    grid = u.grid

    launch!(arch, grid, :xyz, _corse_grain!, d̂, nx, ny, nx2, ny2, Nx, Ny)

    d = ifft(d̂)
    
    set!(u̅, real.(d))
    fill_halo_regions!(u̅)

    return u̅
end

function symmetric_filtering(u::Field; cutoff = 20kilometer)

    Δx = u.grid.Δxᶜᵃᵃ
    Δy = u.grid.Δyᵃᶜᵃ

    u̅l = deepcopy(u)
    u̅h = deepcopy(u)
    d = interior(u)
    dl = deepcopy(d)
    dh = deepcopy(d)

    Nx, Ny, Nz = size(d)
    Nfx, Nfy = Int(Int64(Nx)/2), Int(Int64(Ny)/2)

    # frequencies and wavenumbers
    kx = (fftfreq(Nx)[1:Nfx])/Δx
    ky = (fftfreq(Ny)[1:Nfy])/Δy
    kx, ky = repeat(kx, 1, Nfy), repeat(ky', Nfx, 1)
    k = 2π * sqrt.(kx.^2 + ky.^2)
    kc = 2π/cutoff

    # Fourier transform
    for iz = 1:Nz
        d̂ = (rfft(d[:,:,iz]))

        d̂l = deepcopy(d̂)
        d̂h = deepcopy(d̂)
        d̂l[k.>kc] .= 0
        d̂h[k.<=kc] .= 0
        
        # Inverse Fourier transform
        dl[:, :, iz] = irfft(d̂l, Nx) 
        dh[:, :, iz] = irfft(d̂h, Nx) 
    end
    
    set!(u̅l, dl)
    set!(u̅h, dh)
    fill_halo_regions!(u̅l)
    fill_halo_regions!(u̅h)

    return u̅l, u̅h
end

using FFTW
using ImageFiltering
using Base.Threads

# --- Helper Functions ---

# Computes the periodic distance between two coordinates.
function periodic_distance(x, y, Lx, Ly)
    dx = abs(x)
    dx = dx > Lx/2 ? Lx - dx : dx
    dy = abs(y)
    dy = dy > Ly/2 ? Ly - dy : dy
    return sqrt(dx^2 + dy^2)
end

"""
    build_tophat_kernel(grid_info, cutoff; Lx, Ly, method=:spectral)

Builds a normalized tophat kernel based on a given grid and cutoff length.
- If `method == :spectral`, the kernel is built on the full grid using vectorized operations  
  (suitable for FFT-based convolution).
- If `method == :physical`, the kernel is built only over its nonzero support (which is more  
  efficient for direct physical-space convolution).

The returned kernel is normalized to sum to one.
"""
function build_tophat_kernel(grid_info, cutoff; Lx=1e5, Ly=1e5, method=:physical)
    if method == :spectral
        # Extract grid dimensions and construct node vectors.
        Nx, Ny = grid_info.Nx, grid_info.Ny
        # Here we assume that the grid spans [0, Lx] and [0, Ly].
        xu = range(0, Lx, length=Nx)
        yu = range(0, Ly, length=Ny)
        # Allocate the full-grid kernel.
        kernel = zeros(Float32, Nx, Ny)
        # Use broadcasting to assign 1.0 where the periodic distance is less than cutoff/2.
        @. kernel = Float32(periodic_distance(xu, yu', Lx, Ly) < cutoff/2)
        return kernel / sum(kernel)
    elseif method == :physical
        # For a compact kernel, compute grid spacings.
        Nx, Ny = grid_info.Nx, grid_info.Ny
        dx = Lx / Nx
        dy = Ly / Ny
        # Determine half-widths in grid cells.
        half_width_x = floor(Int, (cutoff/2) / dx)
        half_width_y = floor(Int, (cutoff/2) / dy)
        # Allocate a kernel array covering only the nonzero support.
        kernel = zeros(Float32, 2 * half_width_x + 1, 2 * half_width_y + 1)
        for j in 1:size(kernel, 2), i in 1:size(kernel, 1)
            # Compute the physical distance from the kernel center.
            x = (i - half_width_x - 1) * dx
            y = (j - half_width_y - 1) * dy
            kernel[i, j] = (sqrt(x^2 + y^2) <= cutoff/2) ? 1.0f0 : 0.0f0
        end
        return kernel / sum(kernel)
    else
        throw(ArgumentError("Invalid method: $method. Use :spectral or :physical"))
    end
end

"""
    _sinc(x::Real)

Computes sinc(x) = sin(πx) / (πx), handling the singularity at x=0.
Note: This definition uses πx, common in signal processing. 
The user formula uses sinc(k_c ζ) = sin(k_c ζ) / (k_c ζ). 
We will implement the user's version directly.
"""
function _custom_sinc(x::T) where T <: AbstractFloat
    if iszero(x)
        return T(1.0)
    else
        # Use sinpi/pi for better precision near multiples of π if needed,
        # but here the argument is k_c*ζ, not directly related to π initially.
        return sin(x) / x
    end
end

"""
    lanczos1D(ζ::Real, kc::Real, a::Integer)

Computes the 1D Lanczos filter kernel G_1D(ζ).
G_1D(ζ) = sinc(k_c ζ) * sinc(k_c ζ / a)
where sinc(x) = sin(x)/x.
Requires `kc > 0` and `a` to be a non-zero positive integer.
"""
function lanczos1D(ζ::T, kc::Real, a::Integer) where T <: AbstractFloat
    if kc <= 0
        throw(ArgumentError("Cutoff wavenumber kc must be positive."))
    end
    if a <= 0
        throw(ArgumentError("Parameter 'a' must be a positive integer."))
    end
    
    # Calculate arguments for the sinc functions
    arg1 = T(kc * ζ)
    arg2 = T(arg1 / a) # Or T(kc * ζ / a)
    
    # Compute the two sinc terms
    sinc1 = _custom_sinc(arg1)
    sinc2 = _custom_sinc(arg2)
    
    return sinc1 * sinc2
end

# Computes the periodic distance along one dimension.
function periodic_distance_1d(x_coord::Real, L_dim::Real)
    # Distance from 0 on a periodic domain [0, L_dim]
    # Ensure x_coord is within [0, L_dim] if necessary, though formula works generally
    dx = abs(x_coord) 
    return dx > L_dim / 2 ? L_dim - dx : dx
end


"""
    build_lanczos_kernel(grid_info, kc, a; Lx, Ly, lobes::Integer=a, method=:physical)

Builds a normalized 2D Lanczos filter kernel.

The 1D kernel is G_1D(ζ) = sinc(k_c ζ) * sinc(k_c ζ / a), with sinc(x) = sin(x)/x.
The 2D kernel is G_2D(ζ, η) = G_1D(ζ) * G_1D(η).

Arguments:
- `grid_info`: An object containing grid dimensions (e.g., `Nx`, `Ny`). Needs to have fields `Nx` and `Ny`.
- `kc`: The cutoff wavenumber (must be positive).
- `a`: The Lanczos parameter (must be a positive integer).
- `Lx`, `Ly`: Physical dimensions of the domain.
- `lobes` (optional, default=`a`): Determines the spatial extent of the kernel in the `:physical` method. The kernel extends roughly `lobes * π / kc` distance from the center along each axis. Corresponds to including 'lobes' number of zero-crossings of the sinc(k_c ζ / a) term.
- `method` (optional, default=`:physical`):
    - `:spectral`: Builds the kernel on the full grid, centered at index (1, 1) suitable for FFT-based convolution. Uses periodic distances.
    - `:physical`: Builds a compact kernel covering only the significant non-zero region, suitable for direct convolution. Uses standard Euclidean distances from the kernel center.

Returns:
- `kernel`: A 2D array (Float32) containing the filter kernel, normalized to sum to one.
"""
function build_lanczos_kernel(grid_info, kc::Real, a::Integer; 
                              Lx::Real=1e5, Ly::Real=1e5, 
                              lobes::Integer=a, 
                              method::Symbol=:physical)
    
    # --- Input validation ---
    if kc <= 0
        throw(ArgumentError("Cutoff wavenumber kc must be positive."))
    end
    if a <= 0
        throw(ArgumentError("Parameter 'a' must be a positive integer."))
    end
     if lobes <= 0
        throw(ArgumentError("Parameter 'lobes' must be a positive integer."))
    end

    Nx, Ny = grid_info.Nx, grid_info.Ny
    T = Float32 # Use Float32 as in the example

    if method == :spectral
        # --- Spectral Method (Full Grid, Periodic) ---
        dx = Lx / Nx
        dy = Ly / Ny
        
        # Create coordinate vectors representing distances from origin (0,0) or index (1,1)
        # considering periodicity.
        x_coords = T[periodic_distance_1d(i * dx, Lx) for i in 0:(Nx-1)]
        y_coords = T[periodic_distance_1d(j * dy, Ly) for j in 0:(Ny-1)]

        # Allocate the full-grid kernel
        kernel = zeros(T, Nx, Ny)

        # Compute kernel values using broadcasting
        # G(ζ, η) = G1D(ζ) * G1D(η)
        # Note: x_coords correspond to ζ, y_coords correspond to η
        # We need to evaluate lanczos1D for each x_coord and each y_coord
        lanczos_x = lanczos1D.(x_coords, kc, a)
        lanczos_y = lanczos1D.(y_coords, kc, a)

        # kernel[i, j] = lanczos_x[i] * lanczos_y[j]
        kernel .= lanczos_x * lanczos_y' # Outer product

        # Normalize the kernel
        sum_kernel = sum(kernel)
        if isapprox(sum_kernel, 0.0; atol=sqrt(eps(T)))
             @warn "Sum of spectral Lanczos kernel is close to zero. Normalization might fail or be inaccurate."
             # Avoid division by zero, return unnormalized or handle as error?
             # Let's return the unnormalized kernel with a warning.
             return kernel 
        else
             return kernel ./ sum_kernel
        end
        
    elseif method == :physical
        # --- Physical Method (Compact Kernel, Euclidean Distance) ---
        dx = Lx / Nx
        dy = Ly / Ny

        # Determine kernel extent based on 'lobes' parameter
        # The zero-crossings of sinc(k_c*ζ/a) occur at k_c*ζ/a = n*π, so ζ = n*a*π/kc
        # We extend the kernel to the 'lobes'-th zero-crossing distance.
        half_width_dist_x = T(lobes * π / kc)
        half_width_dist_y = T(lobes * π / kc) # Assuming isotropic cutoff for extent

        # Determine half-widths in grid cells (ensure it includes the center point)
        # Use ceil to ensure the extent reaches the target distance.
        half_width_nx = ceil(Int, half_width_dist_x / dx)
        half_width_ny = ceil(Int, half_width_dist_y / dy)

        # Kernel dimensions
        kernel_nx = 2 * half_width_nx + 1
        kernel_ny = 2 * half_width_ny + 1

        # Allocate the compact kernel array
        kernel = zeros(T, kernel_nx, kernel_ny)

        # Center indices of the kernel array
        center_ix = half_width_nx + 1
        center_iy = half_width_ny + 1

        # Fill the kernel array
        for j in 1:kernel_ny
            # Physical distance in y from the center
            η = T((j - center_iy) * dy)
            lanczos_y = lanczos1D(η, kc, a) # Pre-calculate y-component for the row
            
            # If η is outside the cutoff radius, G1D(η) might be zero or negligible.
            # However, the loop structure covers the calculated width.
            # If lanczos_y is effectively zero, the whole row might be zero. Optimization possible.

            for i in 1:kernel_nx
                # Physical distance in x from the center
                ζ = T((i - center_ix) * dx)
                
                # Calculate the 2D kernel value: G(ζ, η) = G1D(ζ) * G1D(η)
                kernel[i, j] = lanczos1D(ζ, kc, a) * lanczos_y
            end
        end

        # Normalize the kernel
        sum_kernel = sum(kernel)
         if isapprox(sum_kernel, 0.0; atol=sqrt(eps(T)))
             @warn "Sum of physical Lanczos kernel is close to zero. Normalization might fail or be inaccurate."
             # Return unnormalized kernel
             return kernel
         else
             return kernel ./ sum_kernel
         end

    else
        throw(ArgumentError("Invalid method: $method. Use :spectral or :physical"))
    end
end

# --- Main Coarse-Graining Function ---

"""
    coarse_graining!(u::Field, u̅l::Field; kernel=:tophat, cutoff=20kilometer,
                     method=:fft, kernel_build_method=:spectral)

Computes a coarse-grained (filtered) field `u̅l` from the field `u` using a specified 
kernel and cutoff. Two convolution methods are available:
  - `method = :fft`: FFT-based convolution (assumes periodic boundaries).
  - `method = :physical`: Direct physical-space convolution via ImageFiltering.jl 
    (with periodic/circular padding).

A separate helper (`build_tophat_kernel`) builds the tophat kernel. The parameter 
`kernel_build_method` selects whether to build the full-grid (:spectral) or compact (:physical)
kernel. If the cutoff is no greater than the grid spacing, the original field is returned.
"""
function coarse_graining!(u::Field, u̅l::Field; kernel=:tophat, cutoff=20kilometer,
                           method=:physical)
    # Extract interior data and grid nodes.
    d = interior(u)
    xu, yu, _ = nodes(u)
    Nx, Ny, Nz = size(d)
    Lx = parameters.Lx
    Ly = parameters.Ly
    Δh = parameters.Δh
    dx = Lx / Nx
    dy = Ly / Ny
    
    # Early exit: if cutoff is no greater than the grid spacing, return the original field.
    if cutoff < 2Δh
        println("Return the original field due to small cutoff!")
        set!(u̅l, d)
        fill_halo_regions!(u̅l)
        return nothing
    end

    # Allocate the output array.
    dl = similar(d)

    # Currently, only :tophat, :gaussian, or lanczos kernel is implemented.
    if kernel == :tophat
        grid_info = (; Nx=Nx, Ny=Ny)
        Gl = build_tophat_kernel(grid_info, cutoff; Lx=Lx, Ly=Ly, method=method)
    elseif kernel == :gaussian
        Gl = Kernel.gaussian((floor(Int, (cutoff/2) / dx), floor(Int, (cutoff/2) / dy)))
    elseif kernel == :lanczos
        grid_info = (; Nx=Nx, Ny=Ny)
        Gl = build_lanczos_kernel(grid_info, 1/cutoff, 2; Lx=Lx, Ly=Ly, method=method)
    else
        error("Kernel $(kernel) not implemented.")
    end

    # Decide on the convolution method based on the kernel's nonzero fraction relative to the full grid.
    nonzero_fraction = count(!iszero, Gl) / prod(size(d))
    if method == :physical && nonzero_fraction > 0.2
        @info "Kernel is large (nonzero fraction = $(nonzero_fraction)); switching to FFT method."
        method = :spectral
        if kernel == :tophat
            Gl = build_tophat_kernel(grid_info, cutoff; Lx=Lx, Ly=Ly, method=method)
        elseif kernel == :gaussian
            Gl = Kernel.gaussian((floor(Int, (cutoff/2) / dx), floor(Int, (cutoff/2) / dy)))
        elseif kernel == :lanczos
            Gl = build_lanczos_kernel(grid_info, 1/cutoff, 2; Lx=Lx, Ly=Ly, method=method)
        end
    end

    if method != :physical && method != :spectral
        error("Method $(method) not recognized. Use :spectral or :physical.")
    end

    if method == :physical
        t0 = time()
        for iz in 1:Nz
            dl[:, :, iz] .= imfilter(d[:, :, iz], centered(Gl), Pad(:circular))
            if time()-t0 > 10
                println(":physical is slow (>10s); switching to FFT method.")
                method = :spectral
                if kernel == :tophat
                    Gl = build_tophat_kernel(grid_info, cutoff; Lx=Lx, Ly=Ly, method=method)
                elseif kernel == :gaussian
                    Gl = Kernel.gaussian((floor(Int, (cutoff/2) / dx), floor(Int, (cutoff/2) / dy)))
                elseif kernel == :lanczos
                    Gl = build_lanczos_kernel(grid_info, 1/cutoff, 2; Lx=Lx, Ly=Ly, method=method)
                end
            end
            #println("$(time()-t0)s")
        end
    end

    if method == :spectral
        # Precompute the FFT of the kernel.
        Ĝl = rfft(Gl)
        for iz in 1:Nz
            d_slice = d[:, :, iz]
            d_hat = rfft(d_slice)
            filtered_hat = Ĝl .* d_hat
            dl[:, :, iz] .= irfft(filtered_hat, Nx)
        end
    end

    set!(u̅l, dl)
    fill_halo_regions!(u̅l)
    return nothing
end
