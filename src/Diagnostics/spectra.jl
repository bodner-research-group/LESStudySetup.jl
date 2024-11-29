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

function isotropic_powerspectrum(var1, var2, x, y)

    Nx, Ny = size(var1)
    
    Δx, Δy = x[2] - x[1], y[2] - y[1]

    # Fourier transform
    v̂1 = (rfft(var1 .- mean(var1))) * Δx * Δy
    v̂2 = (rfft(var2 .- mean(var2))) * Δx * Δy
    Nfx, Nfy = size(v̂1, 1), size(v̂1, 1)
    v̂1 = v̂1[:,1:Nfx]
    v̂2 = v̂2[:,1:Nfy]

    # Compute the power spectrum
    S = v̂1 .* conj(v̂2)

    # wavenumbers
    kx = (fftfreq(Nx)[1:Nfx])/Δx
    ky = (fftfreq(Ny)[1:Nfy])/Δy
    k = sqrt.(kx.^2 .+ ky'.^2)

    # group the spectra by wavenumber bins and compute the mean
    Δkx, Δky = kx[2] - kx[1], ky[2] - ky[1]
    Δk = sqrt(Δkx * Δky)
    kmin, kmax = minimum(k), maximum(k)
    klen = kmax - kmin
    Nk = ceil(Int, klen/Δk)
    kbins = range(0, Nk*Δk, length = Nk+1)
    spectra = []
    freqs = []
    for i = 1:Nk
        idx = kbins[i] .< k .≤ kbins[i+1]
        if !isnan(mean(S[idx]))
            push!(spectra, Δk * sum(S[idx]))
            push!(freqs, (kbins[i]+kbins[i+1])/2)
        end
    end
    return Spectrum(spectra, 2π*freqs)
end
