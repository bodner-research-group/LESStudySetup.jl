#=
sponge_layer.jl
    Create a forcing function that relaxes the fields to zero in the bottom σ of the domain
    at a rate c * internal wave freq that decays quadratically away from the boundary
=#

using Oceananigans

@inline function sponge_layer_damping_profile(; σ=0.2)
    @inline one_dim(x) = (-1 + σ) < x < (1 - σ) ? 0 : ( (abs(x) - (1 - σ)) / σ )^2
    (x, y, z) -> one_dim(z / parameters.Lz) + one_dim(x / parameters.Lx)
    #(x, y, z) -> (z / parameters.Lz) > (σ-1) ? 0 : (((z / parameters.Lz) + 1 - σ) / σ)^2
end

@inline function sponge_layer_damping(; σ=0.2, c=0.2, target=0)
    rate = c * parameters.N₀ / 2π
    mask=sponge_layer_damping_profile(; σ)
    return Relaxation(; rate, mask, target)
end

@inline function get_sponge_layer_forcing(; σ=0.2, c=0.2)
    #(b, v) = get_filament_state(; verbose=false)
    damping = sponge_layer_damping(; σ, c)
    T_damping = sponge_layer_damping(; σ, c, target=(x, y, z, t)->Tᵢ(x, y, z))
    return (; u=damping, v=damping, w=damping, T=T_damping)
end