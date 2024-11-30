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

@inline minus_sin(θ) = - sin(θ)
@inline minus_cos(θ) = - cos(θ)

@inline function vᵢ(x, y, z)
    Lf = parameters.Lf
    Le = parameters.Le
    R  = parameters.Lx / 4

    return eddy_tangential_velocity(x, y, z, R, Lf, Le, minus_cos)
end

@inline function eddy_tangential_velocity(x, y, z, R, Lf, Le, trig)
    if abs(x - 50e3) > 25e3 && abs(y - 50e3) > 25e3
        x += x < 25e3 ? 50e3 : -50e3
        y += y < 25e3 ? 50e3 : -50e3
    elseif x < 25e3
        x, y = 100e3 - y, 50e3 + x
        trig = trig == sin ? minus_cos : minus_sin
    elseif y < 25e3
        x, y = 50e3 - y, x
        trig = trig == sin ? minus_cos : minus_sin
    elseif x > 75e3
        x, y = 100e3 - y, x - 50e3
        trig = trig == sin ? minus_cos : minus_sin
    elseif y > 75e3
        x, y = 150e3 - y, x
        trig = trig == sin ? minus_cos : minus_sin
    end

    # divide into 4 regions

    # Region 1: warm eddy! (x < 50e3 && y < 50e3)
    x′ = x - R
    y′ = y - R
    
    r  = sqrt(x′^2 + y′^2)
    ξ  = transformR(r, (; R, Le))
    uθ = warm_eddy_velocity(ξ, z, r, R, Lf)
    θ  = atan(y′, x′)
    u1 = trig(θ) * uθ

    for (i, j) in zip([1, 5, 5], [5, 1, 5])
        x′ = x - i * R
        y′ = y - j * R
        
        r  = sqrt(x′^2 + y′^2)
        ξ  = transformR(r, (; R, Le))
        uθ = warm_eddy_velocity(ξ, z, r, R, Lf)
        θ  = atan(y′, x′)
        u1 += trig(θ) * uθ
    end

    # Region 2: cold eddy! (x < 50e3 && y >= 50e3)
    x′ = x - R
    y′ = y - 3R

    r  = sqrt(x′^2 + y′^2)
    ξ  = transformR(r, (; R, Le))
    uθ = cold_eddy_velocity(ξ, z, r, R, Lf)
    θ  = atan(y′, x′)
    u2 = trig(θ) * uθ

    for (i, j) in zip([1, 5, 5], [-1, 3, -1])
        x′ = x - i * R
        y′ = y - j * R
        
        r  = sqrt(x′^2 + y′^2)
        ξ  = transformR(r, (; R, Le))
        uθ = cold_eddy_velocity(ξ, z, r, R, Lf)
        θ  = atan(y′, x′)
        u2 += trig(θ) * uθ
    end

    # Region 3: cold eddy! (x >= 50e3 && y < 50e3)
    x′ = x - 3R
    y′ = y - R

    r  = sqrt(x′^2 + y′^2)
    ξ  = transformR(r, (; R, Le))
    uθ = cold_eddy_velocity(ξ, z, r, R, Lf)
    θ  = atan(y′, x′)
    u3 = trig(θ) * uθ

    for (i, j) in zip([3, -1, -1], [5, 1, 5])
        x′ = x - i * R
        y′ = y - j * R
        
        r  = sqrt(x′^2 + y′^2)
        ξ  = transformR(r, (; R, Le))
        uθ = cold_eddy_velocity(ξ, z, r, R, Lf)
        θ  = atan(y′, x′)
        u3 += trig(θ) * uθ
    end

    # Region 4: warm eddy! 
    x′ = x - 3R
    y′ = y - 3R

    r = sqrt(x′^2 + y′^2)
    ξ = transformR(r, (; R, Le))
    uθ = warm_eddy_velocity(ξ, z, r, R, Lf)
    θ  = atan(y′, x′)
    u4 = trig(θ) * uθ
    
    for (i, j) in zip([-1, 3, -1], [3, -1, -1])
        x′ = x - i * R
        y′ = y - j * R
        
        r  = sqrt(x′^2 + y′^2)
        ξ  = transformR(r, (; R, Le))
        uθ = warm_eddy_velocity(ξ, z, r, R, Lf)
        θ  = atan(y′, x′)
        u4 += trig(θ) * uθ
    end

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
        return warm_eddy(ξ, x, z) + Tᶠ(x, z)
    
    elseif x < 50e3 && y >= 50e3 # Region 2: cold eddy!
        x′ = x - R
        y′ = y - 3R

        r = sqrt(x′^2 + y′^2)
        ξ = transformR(r, (; R, Le))
        return cold_eddy(ξ, x, z) + Tᶠ(x, z)

    elseif x >= 50e3 && y < 50e3 # Region 3: cold eddy!
        x′ = x - 3R
        y′ = y - R

        r = sqrt(x′^2 + y′^2)
        ξ = transformR(r, (; R, Le))
        return cold_eddy(ξ, x, z) + Tᶠ(x, z)

    else # Region 4: warm eddy!
        x′ = x - 3R
        y′ = y - 3R

        r = sqrt(x′^2 + y′^2)
        ξ = transformR(r, (; R, Le))
        return warm_eddy(ξ, x, z) + Tᶠ(x, z)
    end
end

@inline T̅(χ) = parameters.T₀

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

    if abs(x - 50e3) > 25e3 && abs(y - 50e3) > 25e3
        x += x < 25e3 ? 50e3 : -50e3
        y += y < 25e3 ? 50e3 : -50e3
    elseif x < 25e3
        x, y = 100e3 - y, 50e3 + x
    elseif y < 25e3
        x, y = 50e3 - y, x
    elseif x > 75e3
        x, y = 100e3 - y, x - 50e3
    elseif y > 75e3
        x, y = 150e3 - y, x
    end

    # divide into 4 regions

    #if x < 50e3 && y < 50e3 # Region 1: warm eddy!
    x′  = x - R
    y′  = y - R
    r  = sqrt(x′^2 + y′^2)
    η1 = η(r, (; R, Lf, σ², Φ))

    for (i, j) in zip([1, 5, 5], [5, 1, 5])
        x′ = x - i * R
        y′ = y - j * R
        
        r  = sqrt(x′^2 + y′^2)
        η1 += η(r, (; R, Lf, σ², Φ))
    end

    #elseif x < 50e3 && y >= 50e3 # Region 2: cold eddy!
    x′  = x - R
    y′  = y - 3R
    r  = sqrt(x′^2 + y′^2)
    η2 = - η(r, (; R, Lf, σ², Φ))

    for (i, j) in zip([1, 5, 5], [-1, 3, -1])
        x′ = x - i * R
        y′ = y - j * R
        
        r  = sqrt(x′^2 + y′^2)
        η2 += - η(r, (; R, Lf, σ², Φ))
    end

    #elseif x >= 50e3 && y < 50e3 # Region 3: cold eddy!
    x′  = x - 3R
    y′  = y - R
    r  = sqrt(x′^2 + y′^2)
    η3 = - η(r, (; R, Lf, σ², Φ))

    for (i, j) in zip([3, -1, -1], [5, 1, 5])
        x′ = x - i * R
        y′ = y - j * R
        
        r  = sqrt(x′^2 + y′^2)
        η3 += - η(r, (; R, Lf, σ², Φ))
    end

    #else # Region 4: warm eddy!
    x′  = x - 3R
    y′  = y - 3R
    r  = sqrt(x′^2 + y′^2)
    η4 = η(r, (; R, Lf, σ², Φ))

    for (i, j) in zip([-1, 3, -1], [3, -1, -1])
        x′ = x - i * R
        y′ = y - j * R
        
        r  = sqrt(x′^2 + y′^2)
        η4 += η(r, (; R, Lf, σ², Φ))
    end
    #end

    
    return η1 + η2 + η3 + η4
end

@inline function warm_eddy_velocity(ξ, z, r, R, Lf)

    Lz = parameters.Lz
    ΔT = parameters.ΔTᵉ
    f  = parameters.f
    α  = parameters.α
    g  = parameters.g
    σ² = parameters.σ²
    Φ  = parameters.Φ
    a  = parameters.a
    Δh = parameters.Δm

    ∂b∂ξ = - 2g * α * (sin(ξ)^2 - cos(ξ)^2 + 1) / π * a * ΔT * Δh / 3
    ∂b∂ξ = Int(0 < ξ < 3.1415926535897) * ∂b∂ξ
    ∂ξ∂r = - 2π / R * Lf

    uθᴮ = - g / f * ∂η(r, (; R, Lf, σ², Φ))

    h = h̅⁺(ξ)
    if z > - h
        return ∂ξ∂r * ∂b∂ξ / f + uθᴮ
    else
        return ∂ξ∂r * ∂b∂ξ / f * (Lz + z)^3 / (Lz - h)^3 + uθᴮ
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
    a  = parameters.a
    Δh = parameters.Δm

    ∂b∂ξ = 2g * α * (sin(ξ)^2 - cos(ξ)^2 + 1) / π * a * ΔT * Δh / 3
    ∂b∂ξ = Int(0 < ξ < 3.1415926535897) * ∂b∂ξ
    ∂ξ∂r = - 2π / R * Lf

    uθᴮ = g / f * ∂η(r, (; R, Lf, σ², Φ))

    h = h̅⁻(ξ)
    if z > - h
        return ∂ξ∂r * ∂b∂ξ / f + uθᴮ
    else
        return ∂ξ∂r * ∂b∂ξ / f * (Lz + z)^3 / (Lz - h)^3 + uθᴮ
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

""" temperature for pure fronts """
@inline function Tᶠ(x, z)

    Lx  = parameters.Lx
    Lz  = parameters.Lz
    Δz  = parameters.Δz
    Lf  = parameters.Lf * Lx / 45
    N²s = parameters.N²s
    N²T = parameters.N²T
    M²₀ = parameters.M²₀
    h₀  = parameters.m₀
    Δh  = parameters.Δmᶠ
    α   = parameters.α
    g   = parameters.g
    
    ΔTₒ = M²₀ * Lf / (α * g)
    ΓT = 0.5 / (α * g) * ((N²s+0.1*N²T)*z+Δh*((N²s-N²T)*log(cosh((z+h₀)/Δh)/cosh(h₀/Δh))+0.9*N²T*log(cosh((z+1.5h₀)/Δh)/cosh(1.5h₀/Δh))))
    if z <= - (Lz - Δz)
        return ΓT
    elseif x <= Lx/2
        return ΔTₒ * 0.25 * (1-tanh(x/(0.5*Lf))+tanh((x-Lx/2)/(0.5*Lf)))*(tanh((z+h₀)/Δh)+1) + ΓT
    else
        return ΔTₒ * 0.25 * (tanh((x-Lx/2)/(0.5*Lf))-tanh((x-Lx)/(0.5*Lf))-1)*(tanh((z+h₀)/Δh)+1) + ΓT
    end
end

""" velocity for pure fronts """
@inline function vᶠ(x, y, z)

    Lx  = parameters.Lx
    Lf  = parameters.Lf * Lx / 45
    M²₀ = parameters.M²₀
    h₀  = parameters.m₀
    Δh  = parameters.Δmᶠ
    f   = parameters.f

    if z <= -h₀
        return 0.0
    elseif x <= Lx/2
        return M²₀ * 0.5 / f * (sech((x-Lx/2)/(0.5*Lf))^2-sech(x/(0.5*Lf))^2)*(Δh*log(cosh((z+h₀)/Δh))+z+h₀) 
    else
        return M²₀ * 0.5 / f * (sech((x-Lx/2)/(0.5*Lf))^2-sech((x-Lx)/(0.5*Lf))^2)*(Δh*log(cosh((z+h₀)/Δh))+z+h₀)
    end
end

""" velocity for eddies and fronts """
@inline function vᵢᶠ(x, y, z)
    return vᵢ(x, y, z) + vᶠ(x, y, z)
end