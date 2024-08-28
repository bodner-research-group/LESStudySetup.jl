@inline function transformX(x′, p)

    x = x′ - 25e3
    x = ifelse(x < 0, x + 100e3, x)

    return ifelse(x <= p.Lx / 2, 
                 2x / p.Lx * π * (1 + p.Lf) - π/2 * p.Lf,
                 2(p.Lx - x) / p.Lx * π * (1 + p.Lf) - π/2 * p.Lf)
end

@inline transformR(r, p) = 2(p.R - r) / p.R * π * p.Le - π/2 * (p.Le - 1)

@inline function uᵢ(x, y, z)
    Lf = parameters.Lf
    Le = parameters.Le
    R  = parameters.Lx / 4

    return eddy_tangential_velocity(x, y, z, R, Lf, Le, sin)
end

@inline minus_cos(θ) = - cos(θ)

@inline function vᵢ(x, y, z)
    Lf = parameters.Lf
    Le = parameters.Le
    R  = parameters.Lx / 4

    return eddy_tangential_velocity(x, y, z, R, Lf, Le, minus_cos)
end

@inline function eddy_tangential_velocity(x, y, z, R, Lf, Le, trig)
    # divide into 4 regions

    # if x < 50e3 && y < 50e3 # Region 1: warm eddy!
        x′ = x - R
        y′ = y - R
        
        r  = sqrt(x′^2 + y′^2)
        ξ  = transformR(r, (; R, Le))
        uθ = warm_eddy_velocity(ξ, z, r, R, Lf)
        θ  = atan(y′, x′)
        u1 = trig(θ) * uθ
    
    # elseif x < 50e3 && y >= 50e3 # Region 2: cold eddy!
        x′ = x - R
        y′ = y - 3R

        r  = sqrt(x′^2 + y′^2)
        ξ  = transformR(r, (; R, Le))
        uθ = cold_eddy_velocity(ξ, z, r, R, Lf)
        θ  = atan(y′, x′)
        u2 = trig(θ) * uθ

    # elseif x >= 50e3 && y < 50e3 # Region 3: cold eddy!
        x′ = x - 3R
        y′ = y - R

        r  = sqrt(x′^2 + y′^2)
        ξ  = transformR(r, (; R, Le))
        uθ = cold_eddy_velocity(ξ, z, r, R, Lf)
        θ  = atan(y′, x′)
        u3 = trig(θ) * uθ

    # else # Region 4: warm eddy!
        x′ = x - 3R
        y′ = y - 3R

        r = sqrt(x′^2 + y′^2)
        ξ = transformR(r, (; R, Le))
        uθ = warm_eddy_velocity(ξ, z, r, R, Lf)
        θ  = atan(y′, x′)
        u4 = trig(θ) * uθ
    
    return u1 + u2 + u3 + u4
end

@inline function Tᵢ(x, y, z)

    Le = parameters.Le
    R = 25e3

    # divide into 4 regions
    if x < 50e3 && y < 50e3 # Region 1: warm eddy!
        x′ = x - R
        y′ = y - R
        
        r = sqrt(x′^2 + y′^2)
        ξ = transformR(r, (; R, Le))
        return warm_eddy(ξ, x, z)
    
    elseif x < 50e3 && y >= 50e3 # Region 2: cold eddy!
        x′ = x - R
        y′ = y - 3R

        r = sqrt(x′^2 + y′^2)
        ξ = transformR(r, (; R, Le))
        return cold_eddy(ξ, x, z)

    elseif x >= 50e3 && y < 50e3 # Region 3: cold eddy!
        x′ = x - 3R
        y′ = y - R

        r = sqrt(x′^2 + y′^2)
        ξ = transformR(r, (; R, Le))
        return cold_eddy(ξ, x, z)

    else # Region 4: warm eddy!
        x′ = x - 3R
        y′ = y - 3R

        r = sqrt(x′^2 + y′^2)
        ξ = transformR(r, (; R, Le))
        return warm_eddy(ξ, x, z)
    end
end

@inline function T̅(χ) 
    ΔT = parameters.ΔTᶠ
    T₀ = parameters.T₀
    
    T₃ = 1 - (π - χ - sin(π - χ) * cos(π - χ)) / π
    Tₘ = Int(χ > 3.1415926535897) + Int(0 < χ < 3.1415926535897) * T₃

    return ΔT * Tₘ + T₀
end

# Mixed layer depth profile
@inline function h̅⁺(ξ)
    Δh = parameters.Δm
    h₀ = parameters.m₀

    hᵪ = 1 - (π - ξ - sin(π - ξ) * cos(π - ξ)) / π
    hₘ = Int(ξ > 3.1415926535897) + Int(0 < ξ < 3.1415926535897) * hᵪ

    return Δh * hₘ + h₀
end

@inline function h̅⁻(ξ)
    Δh = parameters.Δm
    h₀ = parameters.m₀

    hᵪ = 1 - (π - ξ - sin(π - ξ) * cos(π - ξ)) / π
    hₘ = Int(ξ > 3.1415926535897) + Int(0 < ξ < 3.1415926535897) * hᵪ

    return h₀ - Δh * hₘ
end

""" eddy with isopycnals pushed up """
@inline function cold_eddy(ξ, x, z)

    Lz = parameters.Lz
    T₀ = parameters.T₀
    ΔT = parameters.ΔTᵉ
    Lf = parameters.Lf
    Lx = parameters.Lx
    a  = parameters.a

    χ  = transformX(x, (; Lf, Lx))

    Tˢ = T̅(χ)
    h  = h̅⁻(ξ)

    if z > - h
        return Tˢ
    else
        return (Tˢ - T₀ + a*ΔT) / (Lz - h)^2 * (Lz + z)^2 + T₀ - a*ΔT
    end
end

""" eddy with isopycnals pushed down """
@inline function warm_eddy(ξ, x, z)

    Lz = parameters.Lz
    T₀ = parameters.T₀
    ΔT = parameters.ΔTᵉ
    Lf = parameters.Lf
    Lx = parameters.Lx
    a  = parameters.a

    χ  = transformX(x, (; Lf, Lx))
    Tˢ = T̅(χ)
    h  = h̅⁺(ξ)

    if z > - h
        return Tˢ
    else
        return (Tˢ - T₀ + a*ΔT) / (Lz - h)^2 * (Lz + z)^2 + T₀ - a*ΔT
    end
end

# Free surface as a function of the eddy radius 
@inline  η(r, p) = exp(-(r / p.R)^2 / p.σ²) * p.Φ
@inline ∂η(r, p) = - 2r / p.R^2 / p.σ² * exp(-(r / p.R)^2 / p.σ²) * p.Φ

@inline function ηᵢ(x, y, z) 
    Lf = parameters.Lf
    σ² = parameters.σ²
    Φ  = parameters.Φ
    R  = 25e3

    if x < 50e3 && y < 50e3 # Region 1: warm eddy!
        x′  = x - R
        y′  = y - R
        sng = 1

    elseif x < 50e3 && y >= 50e3 # Region 2: cold eddy!
        x′  = x - R
        y′  = y - 3R
        sng = - 1

    elseif x >= 50e3 && y < 50e3 # Region 3: cold eddy!
        x′  = x - 3R
        y′  = y - R
        sng = - 1

    else # Region 4: warm eddy!
        x′  = x - 3R
        y′  = y - 3R
        sng = 1

    end

    r  = sqrt(x′^2 + y′^2)
    return sng * η(r, (; R, Lf, σ², Φ))
end

@inline function warm_eddy_velocity(ξ, z, r, R, Lf)

    Lz = parameters.Lz
    ΔT = parameters.ΔTᵉ
    f  = parameters.f
    α  = parameters.α
    g  = parameters.g
    σ² = parameters.σ²
    Φ  = parameters.Φ

    ∂b∂ξ = - g * α * ΔT * (sin(ξ)^2 - cos(ξ)^2 + 1) / π
    ∂b∂ξ = Int(0 < ξ < 3.1415926535897) * ∂b∂ξ
    ∂ξ∂r = - 2π / R * Lf

    uθᴮ = - g / f * ∂η(r, (; R, Lf, σ², Φ))

    h = h̅⁺(ξ)
    if z > - h
        return ∂ξ∂r * ∂b∂ξ / f * (Lz - h) / 3 + uθᴮ
    else
        return ∂ξ∂r * ∂b∂ξ / f * (Lz + z)^3 / (Lz - h)^2 / 3 + uθᴮ
    end
end

@inline function cold_eddy_velocity(ξ, z, r, R, Lf)

    Lz = parameters.Lz
    ΔT = parameters.ΔTᵉ
    f  = parameters.f
    α  = parameters.α
    g  = parameters.g
    σ² = parameters.σ²
    Φ  = parameters.Φ

    ∂b∂ξ = g * α * ΔT * (sin(ξ)^2 - cos(ξ)^2 + 1) / π
    ∂b∂ξ = Int(0 < ξ < 3.1415926535897) * ∂b∂ξ
    ∂ξ∂r = - 2π / R * Lf

    uθᴮ = g / f * ∂η(r, (; R, Lf, σ², Φ))

    h = h̅⁻(ξ)
    if z > - h
        return ∂ξ∂r * ∂b∂ξ / f * (Lz - h) / 3 + uθᴮ
    else
        return ∂ξ∂r * ∂b∂ξ / f * (Lz + z)^3 / (Lz - h)^2 / 3 + uθᴮ
    end
end

""" eddy with isopycnals pushed up """
@inline function Tᶻ(x, y, z)

    Lz = parameters.Lz
    T₀ = parameters.T₀
    ΔT = parameters.ΔTᵉ
    h₀ = parameters.m₀
    Δh = parameters.Δm
    a  = parameters.a

    Tˢ = T̅(1)
    h  = h₀ - Δh / 2
    
    if z > - h
        return Tˢ
    else
        return (Tˢ - T₀ + a * ΔT) / (Lz - h)^2 * (Lz + z)^2 + T₀ - a * ΔT
    end
end