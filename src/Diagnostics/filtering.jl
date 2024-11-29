using FFTW
using Oceananigans.Utils
using Oceananigans.Units
using Oceananigans.BoundaryConditions
using Oceananigans.Operators
using Oceananigans.Fields: instantiated_location
using KernelAbstractions: @kernel, @index

@kernel function _horizontal_box_filter!(new_field, grid, field)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        new_field[i, j, k] = field[i, j, k]
        nn = (field[i, j, k], field[i + 1, j, k], field[i - 1, j, k], field[i, j + 1, k], field[i, j - 1, k])    
        new_field[i, j, k] = sum(nn) / 5
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

function coarse_graining!(u::Field, u̅l::Field; kernel = :tophat, cutoff = 20kilometer)

    d = interior(u)

    xu, yu, _ = nodes(u)
    Nx, Ny, Nz = size(d)
    Lx = parameters.Lx 
    Ly = parameters.Ly 
    
    Gl = zeros(Float32, Nx, Ny)
    dl = zeros(Float32, Nx, Ny, Nz)
    if kernel == :tophat
        xG, yG = repeat(xu, 1, Ny), repeat(yu', Nx, 1)
        rG1 = sqrt.(xG.^2 + yG.^2)
        rG2 = sqrt.((xG.-Lx).^2 + yG.^2)   
        rG3 = sqrt.(xG.^2 + (yG.-Ly).^2)
        rG4 = sqrt.((xG.-Lx).^2 + (yG.-Ly).^2)
        Gl[rG1 .< cutoff/2] .= 1
        Gl[rG2 .< cutoff/2] .= 1
        Gl[rG3 .< cutoff/2] .= 1
        Gl[rG4 .< cutoff/2] .= 1
    end
    Ĝl = rfft(Gl/sum(Gl))

    # Fourier transform for convolution
    for iz = 1:Nz
        
        # Inverse Fourier transform
        dl[:, :, iz] = irfft(Ĝl .* (rfft(d[:,:,iz])), Nx) 
    end
    
    set!(u̅l, dl)
    fill_halo_regions!(u̅l)

    return nothing
end