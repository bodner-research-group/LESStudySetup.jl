import Oceananigans.Fields: set!

using Base
using Base: getproperty
using Adapt 

"""
    struct ProblemConstants

A mutable struct representing the constants used in the LES study setup.

# Fields
- `Δρ::Float64`: Density difference at the fronts.
- `ρ₀::Float64`: Reference density of the fluid.
- `Δh::Float64`: Horizontal grid spacing.
- `Δz::Float64`: Vertical grid spacing.
- `Lx::Float64`: Zonal domain size.
- `Ly::Float64`: Meridional domain size.
- `Lz::Float64`: Vertical domain depth.
- `Lf::Float64`: Width of the front.
- `θ::Float64`:  Vertical gradient of potential temperature.
- `f::Float64`:  Coriolis parameter.
- `τw::Float64`: Surface stress.
- `α::Float64`:  Thermal expansion coefficient.
- `Jᵀ::Float64`: Heat flux at the top.
"""
@kwdef mutable struct ProblemConstants
    N²s :: Float64 = 5e-7
    N²T :: Float64 = 1e-4
    M²₀ :: Float64 = 5e-7
    ΔTᵉ :: Float64 = 0.5
    a   :: Float64 = 1.2
    ρ₀  :: Float64 = 1020
    T₀  :: Float64 = 5
    cp  :: Float64 = 3995
    Δh  :: Float64 = 1kilometers
    m₀  :: Float64 = 50
    Δmᶠ :: Float64 = 10
    Δm  :: Float64 = 30
    Δz  :: Float64 = 4meters
    Lx  :: Float64 = 100kilometers
    Ly  :: Float64 = 100kilometers
    Lz  :: Float64 = 256meters
     f  :: Float64 = 1e-4
    τw  :: Float64 = 0.1
     θ  :: Float64 = 30
     Q  :: Float64 = 10
     α  :: Float64 = 2e-4
     Φ  :: Float64 = 0.075
    Lf  :: Float64 = 0.9
    Le  :: Float64 = 0.9
    σ²  :: Float64 = 0.15
     g  :: Float64 = Oceananigans.BuoyancyFormulations.g_Earth
end

Base.show(io::IO, c::ProblemConstants) =
    print(io, "├── surface stratification:       N²s = ", c.N²s, "\n",
              "├── pycnocline stratification:    N²T = ", c.N²T, "\n",
              "├── frontal density gradient:     M²₀ = ", c.M²₀, "\n",
              "├── eddy temperature difference:  ΔTᵉ = ", c.ΔTᵉ, "\n",
              "├── eddy temperature amplitude:     a = ", c.a, "\n",
              "├── reference density:             ρ₀ = ", c.ρ₀, "\n",
              "├── surface temperature:           T₀ = ", c.T₀, "\n",
              "├── heat capacity:                 cp = ", c.cp, "\n",
              "├── initial mixed layer            m₀ = ", c.m₀, "\n",
              "├── frontal mld difference        Δmᶠ = ", c.Δmᶠ, "\n",
              "├── mld difference                 Δm = ", c.Δm, "\n",
              "├── horizontal spacing:            Δh = ", c.Δh, "\n",
              "├── vertical spacing:              Δz = ", c.Δz, "\n",
              "├── x-domain size:                 Lx = ", c.Lx, "\n",
              "├── y-domain size:                 Ly = ", c.Ly, "\n",
              "├── z-domain size:                 Lz = ", c.Lz, "\n",
              "├── Coriolis parameter:             f = ", c.f,  "\n",
              "├── wind stress:                   τw = ", c.τw, "\n",
              "├── wind angle                      θ = ", c.θ,  "\n",
              "├── heat flux:                      Q = ", c.Q,  "\n", 
              "├── thermal expansion:              α = ", c.α,  "\n", 
              "├── barotropic vortex:              Φ = ", c.Φ,  "\n",
              "├── gravity:                        g = ", c.g,  "\n",
              "├── Initial vortex spread:         σ² = ", c.σ², "\n",
              "├── Eddy frontal width:            Le = ", c.Le, "\n",
              "└── Frontal width:                 Lf = ", c.Lf, "\n")

# The constants of the idealized setup
const parameters = ProblemConstants()

# Setting problem constants
set_value!(; kwargs...)      = set!(parameters; kwargs...)
set_value!(var::Symbol, val) = parameters[var] = val

function set_value!(params::ProblemConstants) 
    for name in propertynames(params)
        set_value!(name, getproperty(params, name))
    end
end

Base.getindex(c::ProblemConstants, var::Symbol)     = @eval $c.$var
Base.setindex!(c::ProblemConstants, v, var::Symbol) = @eval $c.$var = $v

function set!(c::ProblemConstants; kwargs...)
    for (fldname, value) in kwargs
        if fldname ∈ propertynames(c)
            c[fldname] = value
        end
    end

    return nothing
end

tuplify(p::ProblemConstants) = NamedTuple{propertynames(parameters)}(Tuple(getproperty(parameters, prop) for prop in propertynames(parameters)))
