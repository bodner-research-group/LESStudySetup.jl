using Oceananigans.Advection: div_𝐯u, div_𝐯v

using Oceananigans.Advection: 
            AbstractAdvectionScheme,
            _advective_tracer_flux_x,
            _advective_tracer_flux_y,
            _advective_tracer_flux_z,
            required_halo_size

using Oceananigans.Fields: ZeroField
using Oceananigans.Utils: SumOfArrays
using Oceananigans.Operators

using Adapt 

import Base
import Oceananigans.Advection: div_Uc, U_dot_∇u, U_dot_∇v

"""
    struct ForcedAdvection{N, FT, A, U, V, W} <: AbstractAdvectionScheme{N, FT}

A structure representing advection from the prognostic velocities plus additional 
background velocities.
"""
struct ForcedAdvection{N, FT, A, U, V, W} <: AbstractAdvectionScheme{N, FT}
    scheme :: A
    u_background :: U
    v_background :: V
    w_background :: W

    ForcedAdvection{N, FT}(s::A, u::U, v::V, w::W) where {N, FT, A, U, V, W} = new{N, FT, A, U, V, W}(s, u, v, w)
end

Adapt.adapt_structure(to, s::ForcedAdvection{N, FT}) where {N, FT} = 
    ForcedAdvection{N, FT}(Adapt.adapt(to, s.scheme), 
                           Adapt.adapt(to, s.u_background), 
                           Adapt.adapt(to, s.v_background), 
                           Adapt.adapt(to, s.w_background))

function ForcedAdvection(; scheme,
                         u_background = ZeroField(eltype(scheme)),
                         v_background = ZeroField(eltype(scheme)),
                         w_background = ZeroField(eltype(scheme)))

    N  = required_halo_size(scheme)
    FT = eltype(scheme)

    return ForcedAdvection{N, FT}(scheme, u_background, v_background, w_background)
end

Base.show(io::IO, f::ForcedAdvection) = 
    print(io, "ForcedAdvection with:", '\n',
              "├── scheme: ", summary(f.scheme), '\n',
              "└── background velocities: ", '\n',
              "    ├── u: ", summary(f.u_background), '\n',
              "    ├── v: ", summary(f.v_background), '\n',
              "    └── w: ", summary(f.w_background))

@inline function U_dot_∇u(i, j, k, grid::RectilinearGrid, advection::ForcedAdvection, U) 

    scheme = advection.scheme

    u = SumOfArrays{2}(U.u, advection.u_background)
    v = SumOfArrays{2}(U.v, advection.v_background)
    w = SumOfArrays{2}(U.w, advection.w_background)

    total_velocities = (; u, v, w)

    return div_𝐯u(i, j, k, grid, scheme, total_velocities, U.u)
    #return div_𝐯v(i, j, k, grid, scheme, total_velocities, total_velocitie.u)-div_𝐯v(i, j, k, grid, scheme, advection, advection.u_background)
end

@inline function U_dot_∇v(i, j, k, grid::RectilinearGrid, advection::ForcedAdvection, U) 

    scheme = advection.scheme

    u = SumOfArrays{2}(U.u, advection.u_background)
    v = SumOfArrays{2}(U.v, advection.v_background)
    w = SumOfArrays{2}(U.w, advection.w_background)

    total_velocities = (; u, v, w)

    return div_𝐯v(i, j, k, grid, scheme, total_velocities, U.v)
    #return div_𝐯v(i, j, k, grid, scheme, total_velocities, total_velocitie.v)-div_𝐯v(i, j, k, grid, scheme, advection, advection.v_background)
end

@inline function div_Uc(i, j, k, grid, advection::ForcedAdvection, U, c)

    scheme = advection.scheme

    u = SumOfArrays{2}(U.u, advection.u_background)
    v = SumOfArrays{2}(U.v, advection.v_background)
    w = SumOfArrays{2}(U.w, advection.w_background)

    return 1/Vᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, _advective_tracer_flux_x, scheme, u, c) +
                                    δyᵃᶜᵃ(i, j, k, grid, _advective_tracer_flux_y, scheme, v, c) +
                                    δzᵃᵃᶜ(i, j, k, grid, _advective_tracer_flux_z, scheme, w, c))
end
