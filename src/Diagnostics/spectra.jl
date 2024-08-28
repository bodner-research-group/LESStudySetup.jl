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

function isotropic_powerspectrum(var1, var2, x, y; nfactor=100)

    Nx, Ny = size(var1)
    Nfx, Nfy = Int(Int64(Nx)/2), Int(Int64(Ny)/2)

    # Hann window
    wx = sin.(π*(1:Nx)/Nx).^2
    wy = sin.(π*(1:Ny)/Ny).^2
    w = wx.*wy'
    
    Δx, Δy = x[2] - x[1], y[2] - y[1]

    # Fourier transform
    v̂1 = (rfft(w.*var1))[2:Nfx+1,2:Nfy+1]
    v̂2 = (rfft(w.*var2))[2:Nfx+1,2:Nfy+1]

    # Compute the power spectrum
    S = vec(v̂1 .* conj(v̂2))

    # frequencies and wavenumbers
    kx = (fftfreq(Nx)[2:Nfx+1])/Δx
    ky = (fftfreq(Ny)[2:Nfy+1])/Δy
    kx, ky = repeat(kx, 1, Nfy), repeat(ky', Nfx, 1)
    k = vec(sqrt.(kx.^2 + ky.^2))

    # group the spectra by wavenumber bins and compute the mean
    kmin, kmax = min(k...), max(k...)
    kdis, klen = kmax - kmin, length(k)
    Nk = Int(Int64(klen/nfactor))
    kbins = range(kmin-1e-3kdis, kmax+1e-3kdis, length = Nk+1)
    spectra = []
    freqs = []
    for i = 1:Nk
        idx = kbins[i] .< k .≤ kbins[i+1]
        if !isnan(mean(S[idx]))
            push!(spectra, mean(S[idx]))
            push!(freqs, mean(k[idx]))
        end
    end
    return Spectrum(spectra.*freqs, 2π*freqs)
end
