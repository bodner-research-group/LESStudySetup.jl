using Oceananigans.Advection: 
            AbstractAdvectionScheme,
            _advective_tracer_flux_x,
            _advective_tracer_flux_y,
            _advective_tracer_flux_z,
            ZeroU,
            required_halo_size_x,
            required_halo_size_y,
            required_halo_size_z
            

using Oceananigans.Fields: ZeroField
using Oceananigans.Utils: SumOfArrays
using Oceananigans.Operators

using Adapt 

import Base
import Oceananigans.Advection: div_Uc, U_dot_∇u, U_dot_∇v
import Oceananigans.Advection: div_𝐯u, div_𝐯v, div_𝐯w
import Oceananigans.Advection: materialize_advection
import Oceananigans.TimeSteppers: time_discretization

"""
    struct ForcedAdvection{N, FT, TD, A, U, V, W, B} <: AbstractAdvectionScheme{N, FT, TD}

A structure representing advection from the prognostic velocities plus additional
background velocities. The advection term is calculated as:
```math
    (U +  u′)⋅ ∇u′
```
When `advect_background` is `true`, the advection of the background flow by the prognostic
velocity is also included, in flux form ``∇⋅(u′U) = u′⋅∇U``:
```math
    (U +  u′)⋅ ∇u′ + u′ ⋅ ∇U
```
The background self-advection ``U ⋅ ∇U`` is always neglected.
"""
struct ForcedAdvection{N, FT, TD, A, U, V, W, B} <: AbstractAdvectionScheme{N, FT, TD}
    scheme :: A
    u_background :: U
    v_background :: V
    w_background :: W
    advect_background :: B

    ForcedAdvection{N, FT, TD}(s::A, u::U, v::V, w::W, b::B) where {N, FT, TD, A, U, V, W, B} =
        new{N, FT, TD, A, U, V, W, B}(s, u, v, w, b)
end

time_discretization(advection::ForcedAdvection) = time_discretization(advection.scheme)

materialize_advection(advection::ForcedAdvection, grid) =
    ForcedAdvection(; scheme = materialize_advection(advection.scheme, grid),
                      u_background = advection.u_background,
                      v_background = advection.v_background,
                      w_background = advection.w_background,
                      advect_background = advection.advect_background)

Adapt.adapt_structure(to, s::ForcedAdvection{N, FT, TD}) where {N, FT, TD} =
    ForcedAdvection{N, FT, TD}(Adapt.adapt(to, s.scheme),
                               Adapt.adapt(to, s.u_background),
                               Adapt.adapt(to, s.v_background),
                               Adapt.adapt(to, s.w_background),
                               Adapt.adapt(to, s.advect_background))

function ForcedAdvection(; scheme,
                         u_background = ZeroField(eltype(scheme)),
                         v_background = ZeroField(eltype(scheme)),
                         w_background = ZeroField(eltype(scheme)),
                         advect_background = false)

    N  = max(required_halo_size_x(scheme),
             required_halo_size_y(scheme),
             required_halo_size_z(scheme))

    FT = eltype(scheme)
    TD = typeof(time_discretization(scheme))

    b = advect_background isa Val ? advect_background : Val(advect_background)

    return ForcedAdvection{N, FT, TD}(scheme, u_background, v_background, w_background, b)
end

Base.show(io::IO, f::ForcedAdvection) =
    print(io, "ForcedAdvection with:", '\n',
              "├── scheme: ", summary(f.scheme), '\n',
              "├── advect background: ", f.advect_background isa Val{true}, '\n',
              "└── background velocities: ", '\n',
              "    ├── u: ", summary(f.u_background), '\n',
              "    ├── v: ", summary(f.v_background), '\n',
              "    └── w: ", summary(f.w_background))

# Advection of the background velocity by the prognostic velocity, u′⋅∇U = ∇⋅(u′U),
# enabled only when the scheme carries `Val(true)`.
@inline advect_background_u(i, j, k, grid, scheme, ::Val{false}, U, u_background) = zero(grid)
@inline advect_background_u(i, j, k, grid, scheme, ::Val{true},  U, u_background) =
    div_𝐯u(i, j, k, grid, scheme, U, u_background)

@inline advect_background_v(i, j, k, grid, scheme, ::Val{false}, U, v_background) = zero(grid)
@inline advect_background_v(i, j, k, grid, scheme, ::Val{true},  U, v_background) =
    div_𝐯v(i, j, k, grid, scheme, U, v_background)

@inline function U_dot_∇u(i, j, k, grid::RectilinearGrid, advection::ForcedAdvection, U)

    scheme = advection.scheme

    u = SumOfArrays{2}(U.u, advection.u_background)
    v = SumOfArrays{2}(U.v, advection.v_background)
    w = SumOfArrays{2}(U.w, advection.w_background)

    total_velocities = (; u, v, w)

    return div_𝐯u(i, j, k, grid, scheme, total_velocities, U.u) +
           advect_background_u(i, j, k, grid, scheme, advection.advect_background, U, advection.u_background)
end

@inline function U_dot_∇v(i, j, k, grid::RectilinearGrid, advection::ForcedAdvection, U)

    scheme = advection.scheme

    u = SumOfArrays{2}(U.u, advection.u_background)
    v = SumOfArrays{2}(U.v, advection.v_background)
    w = SumOfArrays{2}(U.w, advection.w_background)

    total_velocities = (; u, v, w)

    return div_𝐯v(i, j, k, grid, scheme, total_velocities, U.v) +
           advect_background_v(i, j, k, grid, scheme, advection.advect_background, U, advection.v_background)
end

@inline function div_𝐯u(i, j, k, grid, advection::ForcedAdvection, U, u)

    scheme = advection.scheme

    tu = SumOfArrays{2}(U.u, advection.u_background)
    tv = SumOfArrays{2}(U.v, advection.v_background)
    tw = SumOfArrays{2}(U.w, advection.w_background)

    total_velocities = (; u = tu, v = tv, w = tw)

    return div_𝐯u(i, j, k, grid, scheme, total_velocities, u) +
           advect_background_u(i, j, k, grid, scheme, advection.advect_background, U, advection.u_background)
end

@inline function div_𝐯v(i, j, k, grid, advection::ForcedAdvection, U, v)

    scheme = advection.scheme

    tu = SumOfArrays{2}(U.u, advection.u_background)
    tv = SumOfArrays{2}(U.v, advection.v_background)
    tw = SumOfArrays{2}(U.w, advection.w_background)

    total_velocities = (; u = tu, v = tv, w = tw)

    return div_𝐯v(i, j, k, grid, scheme, total_velocities, v) +
           advect_background_v(i, j, k, grid, scheme, advection.advect_background, U, advection.v_background)
end

@inline function div_𝐯w(i, j, k, grid, advection::ForcedAdvection, U, w) 

    scheme = advection.scheme

    tu = SumOfArrays{2}(U.u, advection.u_background)
    tv = SumOfArrays{2}(U.v, advection.v_background)
    tw = SumOfArrays{2}(U.w, advection.w_background)

    total_velocities = (; u = tu, v = tv, w = tw)

    return div_𝐯w(i, j, k, grid, scheme, total_velocities, w)
end

@inline div_𝐯u(i, j, k, grid, ::ForcedAdvection, ::ZeroU, u) = zero(grid)
@inline div_𝐯v(i, j, k, grid, ::ForcedAdvection, ::ZeroU, v) = zero(grid)
@inline div_𝐯w(i, j, k, grid, ::ForcedAdvection, ::ZeroU, w) = zero(grid)

@inline div_𝐯u(i, j, k, grid, ::ForcedAdvection, U, ::ZeroField) = zero(grid)
@inline div_𝐯v(i, j, k, grid, ::ForcedAdvection, U, ::ZeroField) = zero(grid)
@inline div_𝐯w(i, j, k, grid, ::ForcedAdvection, U, ::ZeroField) = zero(grid)

@inline function div_Uc(i, j, k, grid, advection::ForcedAdvection, U, c)

    scheme = advection.scheme

    u = SumOfArrays{2}(U.u, advection.u_background)
    v = SumOfArrays{2}(U.v, advection.v_background)
    w = SumOfArrays{2}(U.w, advection.w_background)

    return 1/Vᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, _advective_tracer_flux_x, scheme, u, c) +
                                    δyᵃᶜᵃ(i, j, k, grid, _advective_tracer_flux_y, scheme, v, c) +
                                    δzᵃᵃᶜ(i, j, k, grid, _advective_tracer_flux_z, scheme, w, c))
end

@inline div_Uc(i, j, k, grid, ::ForcedAdvection, ::ZeroU, c) = zero(grid)
@inline div_Uc(i, j, k, grid, ::ForcedAdvection, U, ::ZeroField) = zero(grid)