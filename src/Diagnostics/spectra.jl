using FFTW
using Oceananigans.Grids: φnode
using Statistics: mean
using CUDA, Tullio

struct Spectrum{S, F}
    spec :: S
    freq :: F
end

import Base

Base.:(+)(s::Spectrum, t::Spectrum) = Spectrum(s.spec .+ t.spec, s.freq)
Base.:(*)(s::Spectrum, t::Spectrum) = Spectrum(s.spec .* t.spec, s.freq)
Base.:(/)(s::Spectrum, t::Int)      = Spectrum(s.spec ./ t, s.freq)

Base.real(s::Spectrum) = Spectrum(real.(s.spec), s.freq)
Base.abs(s::Spectrum)  = Spectrum( abs.(s.spec), s.freq)

function power_cospectrum_1d(var1, var2, x)

    Nx  = length(x)
    Nfx = Int64(Nx)
    
    spectra = zeros(ComplexF64, Int(Nfx/2))
    
    dx = x[2] - x[1]

    freqs = fftfreq(Nfx, 1.0 / dx) # 0, +ve freq,-ve freqs (lowest to highest)
    freqs = freqs[1:Int(Nfx/2)] .* 2.0 .* π
    
    fourier1   = fft(var1) / Nfx
    fourier2   = fft(var2) / Nfx
    spectra[1] += fourier1[1] .* conj(fourier2[1]) .+ fourier2[1] .* conj(fourier1[1])

    for m in 2:Int(Nfx/2)
        spectra[m] += fourier1[m] .* conj(fourier2[m]) .+ fourier2[m] .* conj(fourier1[m])
    end
    return Spectrum(spectra, freqs)
end

"""
Calculates the 1D isotropic cross-power spectrum from two 2D fields.

Includes optional Hann windowing with correct normalization and improved binning efficiency.
Normalization ensures sum(spectra .* Δfs) ≈ sum(v1_prime .* v2_prime) / (Nx * Ny) / window_correction_factor.

Args:
    var1 (Matrix): First 2D data field.
    var2 (Matrix): Second 2D data field (can be the same as var1 for auto-spectrum).
    window (Union{Nothing, Symbol}): Specify `:hann` for Hann window, or nothing.
    L_filter (Union{Real, Nothing}): Filter scale for Gaussian filter (in grid units).
    use_gpu (Bool): Whether to use GPU acceleration.

Returns:
    Tuple: A tuple containing:
        - freqs: Wavenumber bin centers (k).
        - spectra: Binned spectral density values ([var1][var2]/k).
        - Δfs: Wavenumber width of each bin (Δk).
        - v1_prime: Mean-subtracted (and potentially windowed) first input field.
        - v2_prime: Mean-subtracted (and potentially windowed) second input field.
        - Δx: Grid spacing in the first dimension.
        - Δy: Grid spacing in the second dimension.
        - window_correction_factor: The factor applied to correct for window energy loss (1.0 if no window).
"""
function isotropic_powerspectrum(var1::AbstractMatrix{T1}, var2::AbstractMatrix{T2}; T::DataType=Float32,
                                 L_filter::Union{Real, Nothing}=nothing, use_gpu::Bool=true,
                                 window::Union{Nothing, Symbol} = nothing) where {T1<:Real, T2<:Real}

    if size(var1) != size(var2)
        error("Input variables var1 and var2 must have the same dimensions.")
    end
    CT = T == Float32 ? ComplexF32 : ComplexF64

    # Check for GPU and select the device/array type
    can_use_gpu = use_gpu && CUDA.functional()
    device = can_use_gpu ? CuArray : Array
    println(can_use_gpu ? "✅ GPU detected, using CUDA.jl." : "🖥️  No functional GPU, using CPU.")

    # --- Grid parameters ---
    Nx, Ny = size(var1)
    # Correctly calculate grid spacing assuming uniform grid
    Δx = parameters.Δh
    Δy = parameters.Δh

    # --- Preprocessing ---
    v1_mean = mean(var1)
    v2_mean = mean(var2)
    v1_prime = var1 .- v1_mean
    v2_prime = var2 .- v2_mean

    # --- Windowing (Optional) ---
    window_correction_factor = 1.0
    local v1_fft_input, v2_fft_input
    if isnothing(window)
        v1_fft_input = v1_prime
        v2_fft_input = v2_prime
    elseif window == :hann
        # Calculate Hann window
        wx = sin.(π .* (0:Nx-1) ./ (Nx-1)).^2
        wy = sin.(π .* (0:Ny-1) ./ (Ny-1)).^2
        w = reshape(wx, Nx, 1) .* reshape(wy, 1, Ny)
        # Apply window
        v1_fft_input = w .* v1_prime
        v2_fft_input = w .* v2_prime
        # Calculate window correction factor (mean square of the window)
        window_correction_factor = mean(w .^ 2)
    elseif window == :xhann
        # Calculate xHann window
        wx = sin.(π .* (0:Nx-1) ./ (Nx-1)).^2
        w = reshape(wx, Nx, 1) 
        # Apply window
        v1_fft_input = w .* v1_prime
        v2_fft_input = w .* v2_prime
        # Calculate window correction factor (mean square of the window)
        window_correction_factor = mean(w .^ 2)
    else
        error("Unsupported window type: $window")
    end

    # --- Fourier transform ---
    # Note: v1_prime and v2_prime returned are the *original* mean-subtracted fields,
    #       *before* potential windowing, as this is often more useful for reference.
    #       If windowed fields are needed, return v1_fft_input, v2_fft_input instead.
    v̂1 = rfft(device(v1_fft_input))
    v̂2 = rfft(device(v2_fft_input))

    Nkx, Nky = size(v̂1) # Nkx = div(Nx,2)+1, Nky = Ny

    # --- Compute raw cross power spectrum (unweighted, unnormalized) ---
    S_raw = Array(v̂1 .* conj(v̂2)) # Raw coefficients product
    @info "Raw cross-power spectrum computed."

    # Compute filtered spectrum if L_filter is specified
    S_filtered = nothing
    if !isnothing(L_filter)
        # --- Create Gaussian filter G_hat in spectral space ---
        # The standard deviation `sigma` is derived from the filter scale `L`.
        σ = L_filter / parameters.Δh# / sqrt(12.0)

        # Create centered coordinate grids for the physical kernel
        x = fftshift(-Nx/2 : Nx/2 - 1)
        y = fftshift(-Ny/2 : Ny/2 - 1)

        # Generate the physical-space Gaussian kernel, normalized
        @tullio G[i, j] := T(exp(-(x[i]^2 + y[j]^2) / (2.0 * σ^2)))
        G ./= sum(G)

        # Get the transfer function G_hat and move to device
        # Note: We create the kernel on the CPU as it's a small, one-time cost.
        Ĝ = device(rfft(G))

        # --- Apply filter and compute filtered spectrum ---
        v̂1_filt = Ĝ .* v̂1
        v̂2_filt = Ĝ .* v̂2
        S_filtered = Array(v̂1_filt .* conj(v̂2_filt))
    end

    # --- Wavenumbers ---
    # Correct calculation using spacing Δx, Δy
    kx_vec = 2π .* rfftfreq(Nx, 1/Δx) # Frequencies for rfft output
    ky_vec = 2π .* fftfreq(Ny, 1/Δy)  # Frequencies for standard fft dimension
    Δkx = length(kx_vec) > 1 ? kx_vec[2] - kx_vec[1] : 0.0
    Δky = length(ky_vec) > 1 ? abs(ky_vec[2] - ky_vec[1]) : 0.0
    Δk = Δkx ≈ 0.0 || Δky ≈ 0.0 ? max(Δkx, Δky) : sqrt(Δkx * Δky)

    # Compute wavenumber magnitude grid
    k_mag_sq = zeros(T, Nkx, Nky)
    # Outer product equivalent using broadcasting for potentially better efficiency
    # Ensure ky_vec is treated as a row vector for broadcasting
    k_mag_sq .= (kx_vec.^2) .+ reshape(ky_vec.^2, 1, Nky)
    # Avoid sqrt of negative numbers due to potential floating point issues near zero
    k = sqrt.(max.(0.0, k_mag_sq))

    # --- Binning parameters ---
    # Check if Δx or Δy are zero or very small
    kmin, kmax = 0.0, maximum(k)
    klen = kmax - kmin

    # Avoid division by zero if Δk is effectively zero
    Nk = Δk > 1e-12 ? max(1, ceil(Int, klen / Δk)) : 1
    # Define bin edges
    kbins = range(kmin, stop=kmax + Δk, length=Nk+1) # Extend slightly to ensure max k included

    # --- Optimized Weight Calculation ---
    is_nx_even = iseven(Nx)
    nyquist_i_index = div(Nx, 2) + 1
    # Weights depend only on the kx index (i). Create a column vector.
    weights_col = map(i -> (i == 1 || (is_nx_even && i == nyquist_i_index)) ? 1.0 : 2.0, 1:Nkx)

    # Broadcast the column to the full 2D grid size
    weights = device(weights_col .* ones(T, 1, Nky))

    # --- Vectorized Bin Assignment and Scatter-Add ---
    # Move data to target device
    k_device = device(k)
    S_raw_device = device(S_raw)
    S_filt_device = isnothing(S_filtered) ? nothing : device(S_filtered)

    # A. Find bin index for EVERY wavenumber point in a single, vectorized step
    # `searchsortedfirst` is efficient for this. `kbins` is small, so keep it on CPU.
    bin_indices = searchsortedfirst.(Ref(kbins), k_device) .- 1
    # Clamp indices to be within 1:Nk, as `searchsortedfirst` can return Nk+1
    CUDA.allowscalar() do
        @tullio bin_indices[i] = clamp(bin_indices[i], 1, Nk)
    end

    # B. Sum spectra into bins using a high-performance scatter-add
    spectra_summed = zeros(CT, Nk) |> device
    weighted_S_raw = S_raw_device .* weights

    if can_use_gpu
        # GPU Path: Use a custom kernel with atomic operations for race-free summation
        CUDA.@sync begin
            @cuda threads=256 blocks=cld(length(weighted_S_raw), 256) scatter_add_kernel!(spectra_summed, vec(weighted_S_raw), vec(bin_indices))
        end
    else
        # CPU Path: Use Tullio.jl for a highly optimized, multithreaded scatter-add
        @tullio spectra_summed[bin_indices[i]] += weighted_S_raw[i]
    end

    # Repeat for filtered spectrum if it exists
    spectra_summed_filtered = nothing
    if !isnothing(S_filtered)
        spectra_summed_filtered = zeros(CT, Nk) |> device
        weighted_S_filt = S_filt_device .* weights
        if can_use_gpu
            CUDA.@sync begin
                @cuda threads=256 blocks=cld(length(weighted_S_filt), 256) scatter_add_kernel!(spectra_summed_filtered, vec(weighted_S_filt), vec(bin_indices))
            end
        else
            @tullio spectra_summed_filtered[bin_indices[i]] += weighted_S_filt[i]
        end
    end

    # --- Final Normalization ---
    freqs = (kbins[1:Nk] .+ kbins[2:Nk+1]) ./ 2
    Δfs = kbins[2:Nk+1] .- kbins[1:Nk]
    # Factor includes:
    # 1 / (Nx * Ny)^2 : From FFTW normalization convention (forward and backward)
    # 1 / window_correction_factor : To correct for energy loss from windowing
    norm_factor = 1 / (Nx * Ny)^2 / window_correction_factor
    # Apply normalization, avoiding division by zero if Δfs is zero
    # Handle division by zero for empty bins
    inv_Δfs = map(df -> df > 1e-12 ? 1.0 / df : 0.0, Δfs)
    spectra = Array(spectra_summed) .* (norm_factor .* inv_Δfs)
    if !isnothing(S_filtered)
        spectra_filtered = Array(spectra_summed_filtered) .* (norm_factor .* inv_Δfs)
        return (spec=spectra, specf=spectra_filtered, freq=freqs)
    else
        return Spectrum(spectra, freqs)
    end

    # Return relevant quantities
    # return Spectrum(spectra, freqs)#freqs, spectra, Δfs, v1_prime, v2_prime, Δx, Δy, window_correction_factor
end

"""
Custom CUDA kernel for a high-performance, race-free scatter-add operation.
"""
function scatter_add_kernel!(dest::CuDeviceVector, src::CuDeviceVector, indices::CuDeviceVector{Int})
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if i <= length(src)
        # Atomically add the source value to the destination bin
        CUDA.@atomic dest[indices[i]] += src[i]
    end
    return
end