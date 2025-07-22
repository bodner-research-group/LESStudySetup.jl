using FFTW
using Oceananigans.Grids: φnode
using Statistics: mean

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
    x (Vector): Coordinates for the first dimension.
    y (Vector): Coordinates for the second dimension.
    window (Union{Nothing, Symbol}): Specify `:hann` for Hann window, or nothing.

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
function isotropic_powerspectrum(var1::AbstractMatrix{T1}, var2::AbstractMatrix{T2},
                                 x::AbstractVector{T3}, y::AbstractVector{T4};
                                 window::Union{Nothing, Symbol} = nothing) where {T1<:Real, T2<:Real, T3<:Real, T4<:Real}

    if size(var1) != size(var2)
        error("Input variables var1 and var2 must have the same dimensions.")
    end
    if size(var1) != (length(x), length(y))
        error("Input variable dimensions must match length of x and y coordinates.")
    end

    Nx, Ny = size(var1)
    # Correctly calculate grid spacing assuming uniform grid
    Δx = x[2] - x[1]
    Δy = y[2] - y[1]
    # Optional: Add checks for uniform spacing if needed
    # if !all(diff(x) .≈ Δx) || !all(diff(y) .≈ Δy)
    #     @warn "Grid spacing may not be uniform. Using first difference."
    # end

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
    v̂1 = rfft(v1_fft_input)
    v̂2 = rfft(v2_fft_input)

    Nkx, Nky = size(v̂1) # Nkx = div(Nx,2)+1, Nky = Ny

    # --- Compute raw cross power spectrum (unweighted, unnormalized) ---
    S_raw = v̂1 .* conj(v̂2) # Raw coefficients product

    # --- Wavenumbers ---
    # Correct calculation using spacing Δx, Δy
    kx_vec = 2π .* rfftfreq(Nx, 1/Δx) # Frequencies for rfft output
    ky_vec = 2π .* fftfreq(Ny, 1/Δy)  # Frequencies for standard fft dimension

    # Compute wavenumber magnitude grid
    k_mag_sq = zeros(Float64, Nkx, Nky)
    # Outer product equivalent using broadcasting for potentially better efficiency
    # Ensure ky_vec is treated as a row vector for broadcasting
    k_mag_sq .= (kx_vec.^2) .+ reshape(ky_vec.^2, 1, Nky)
    # Avoid sqrt of negative numbers due to potential floating point issues near zero
    k = sqrt.(max.(0.0, k_mag_sq))

    # --- Binning parameters ---
    # Check if Δx or Δy are zero or very small
    Δkx = if length(kx_vec) > 1 kx_vec[2] - kx_vec[1] else 0.0 end # 2π / (Nx*Δx)
    Δky = if length(ky_vec) > 1 abs(ky_vec[2] - ky_vec[1]) else 0.0 end # 2π / (Ny*Δy)
    Δk = if Δkx ≈ 0.0 || Δky ≈ 0.0
             max(Δkx, Δky) # Handle 1D cases
         else
             sqrt(Δkx * Δky) # Geometric mean for 2D
         end
    kmin, kmax = 0.0, maximum(k)
    klen = kmax - kmin

    # Avoid division by zero if Δk is effectively zero
    if Δk < 1e-12 && klen > 1e-9 # Check if range is non-trivial but spacing is zero
         @warn "Wavenumber spacing Δk is near zero ($Δk) but range is $klen. Check grid setup. Defaulting Nk=1."
         Nk = 1
    elseif Δk < 1e-12
         Nk = 1 # Only one bin if spacing is zero (e.g., single point in k-space)
    else
         Nk = max(1, ceil(Int, klen / Δk)) # Ensure at least one bin
    end
    # Define bin edges
    kbins = range(kmin, stop=kmax + Δk, length=Nk+1) # Extend slightly to ensure max k included

    # --- Prepare for efficient binning ---
    k_flat = vec(k)
    S_flat = vec(S_raw)
    i_indices = repeat(1:Nkx, Nky) # Get the 'i' index (rfft dimension) for each flattened element

    # Pre-calculate weights based on the 'i' index
    is_nx_even = iseven(Nx)
    nyquist_i_index = div(Nx, 2) + 1
    # Use 1.0 and 2.0 for type stability with complex S_flat multiplication
    weights = map(i -> (i == 1 || (is_nx_even && i == nyquist_i_index)) ? 1.0 : 2.0, i_indices)

    # --- Bin and sum using vectorized indexing (more efficient) ---
    spectra_summed = zeros(ComplexF64, Nk) # Store the summed values per bin
    freqs = zeros(Float64, Nk)
    Δfs = zeros(Float64, Nk) # Store bin widths

    for bin_i in 1:Nk
        bin_start = kbins[bin_i]
        bin_end = kbins[bin_i+1]

        # Create boolean mask for elements in the current bin
        # Handle k=0 inclusion carefully for the first bin
        if bin_i == 1
            # Include lower bound for the first bin (k=0)
            idx = (bin_start .≤ k_flat) .& (k_flat .≤ bin_end)
            # Special case: If kmax is exactly 0, ensure the single k=0 point is captured
            if kmax == 0.0 && bin_start == 0.0
                idx = (k_flat .== 0.0)
            end
        else
            # Exclude lower bound for subsequent bins
            idx = (bin_start .< k_flat) .& (k_flat .≤ bin_end)
        end

        if any(idx)
            # Sum the weighted raw spectral components within this bin
            spectra_summed[bin_i] = sum(S_flat[idx] .* weights[idx])
        # else: bin_sum remains zero
        end

        # Store bin center frequency and bin width
        freqs[bin_i] = (bin_start + bin_end) / 2
        Δfs[bin_i] = bin_end - bin_start
    end

    # --- Final Normalization ---
    # Factor includes:
    # 1 / (Nx * Ny)^2 : From FFTW normalization convention (forward and backward)
    # 1 / window_correction_factor : To correct for energy loss from windowing
    norm_factor = 1 / (Nx * Ny)^2 / window_correction_factor
    # Apply normalization, avoiding division by zero if Δfs is zero
    spectra = map( (sum_val, df) -> df > 1e-12 ? sum_val * norm_factor / df : zero(ComplexF64),
                   spectra_summed, Δfs)

    # Return relevant quantities
    return Spectrum(spectra, freqs)#freqs, spectra, Δfs, v1_prime, v2_prime, Δx, Δy, window_correction_factor
end