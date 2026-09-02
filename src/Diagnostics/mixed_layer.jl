const c = Center()
const f = Face()

#####
##### MixedLayerDepthField
#####

# b can be temperature (T) or density (ρ)
@kernel function _compute_mixed_layer_properties!(h, Ew, grid, b, Δb, surface, stratification, w)
    i, j = @index(Global, NTuple)

    Nz = grid.Nz

    k_start = Nz - 1
    z_ij = znode(i, j, k_start+1, grid, c, c, f)
    if !isnothing(w) 
       w_ij = @inbounds w[i, j, k_start+1]
       Ew_integrated = w_ij^2/2
    end
    if !surface
        while z_ij > -10.0
            k_start = k_start-1
            z_ij = znode(i, j, k_start+1, grid, c, c, f)
            if !isnothing(w) 
               w_ij = @inbounds w[i, j, k_start+1]
               Ew_integrated += w_ij^2/2 * Δzᶜᶜᶠ(i, j, k_start+1, grid)
            end
        end
    end

    if stratification
        α = parameters.α
        g = parameters.g
        Nh² = 0.0
    end

    b_surface = @inbounds b[i, j, k_start+1]

    @unroll for k in k_start : -1 : 1 # scroll from point just below surface

        b⁺ = @inbounds b[i, j, k+1]
        bᵏ = @inbounds b[i, j, k]

        # If temperature decreases downwards, both are > 0
        # If density increases downwards, both are < 0
        Δb⁺ = b_surface - b⁺
        Δbᵏ = b_surface - bᵏ

        zᵏ = znode(i, j, k, grid, c, c, c)
        Δz⁺ = Δzᶜᶜᶠ(i, j, k+1, grid)

        # Assuming temperature decreases downwards and density increases upwards
        # Linearly interpolate to find mixed layer height
        inside_mixed_layer = (Δb⁺ < Δb) & (Δbᵏ < Δb)
        just_below_mixed_layer = (Δb⁺ < Δb) & (Δbᵏ >= Δb)
        new_z_ij = zᵏ + (Δb - Δbᵏ) / (Δb⁺ - Δbᵏ) * Δz⁺
        
        # Replace z_ij if we found a new mixed layer depth
        replace_z = (just_below_mixed_layer | inside_mixed_layer) & !inactive_node(i, j, k, grid, c, c, c)
        z_ij = ifelse(replace_z, new_z_ij, z_ij)
        if just_below_mixed_layer 
            if stratification
                if zᵏ < z_ij * 1.1
                    Nh² = α * g * z_ij * 0.1 / Δz⁺ * (Δb⁺ - Δbᵏ)
                else
                    b_ij = bᵏ + (z_ij - zᵏ) / Δz⁺ * (Δbᵏ - Δb⁺)
                    for l in k : -1 : 2
                        Δzˡ = Δzᶜᶜᶠ(i, j, l, grid)
                        zˡ = znode(i, j, l-1, grid, c, c, c)
                        if zˡ < z_ij * 1.1
                            b⁻ = @inbounds b[i, j, l-1]
                            bˡ = @inbounds b[i, j, l]
                            b_new = (z_ij * 1.1 - zˡ) / Δzˡ * (b⁻ - bˡ) + bˡ
                            Nh² = α * g * (b_new - b_ij) / (0.1 * z_ij)
                            break
                        end
                    end
                end
            end
            if !isnothing(w) 
                w_ij = @inbounds w[i, j, k]
                Ew_integrated += w_ij^2/2 * (Δb - Δbᵏ) / (Δb⁺ - Δbᵏ) * Δz⁺
            end
            break
        elseif !isnothing(w)
            w_ij = @inbounds w[i, j, k]
            Ew_integrated += w_ij^2/2 * Δz⁺
        end
    end

    # Note "-" since `h` is supposed to be "depth" rather than "height"
    if stratification
        @inbounds h[i, j, 1] = Nh²
    else
        @inbounds h[i, j, 1] = - z_ij
    end
    if !isnothing(w) 
        @inbounds Ew[i, j, 1] = z_ij < 0 ? Ew_integrated / (-z_ij) : 0.0
    end
end

"""
    MixedLayerDepth(model_or_grid; tracers, w=nothing, ΔT=0.2)

Computes the mixed layer depth based on a temperature difference criterion.

Arguments:
- `grid`: An Oceananigans `grid` object.
- `tracers`: A `NamedTuple` containing the temperature field, e.g., `(; T)`.
- `w`: (Optional) The vertical velocity field `w`. If provided, the mean vertical
       kinetic energy (w^2/2) within the mixed layer is also computed.
- `ΔT`: The temperature difference criterion for the mixed layer depth.
- `surface`: Whether to consider the surface layer (default: `false`).
- `stratification`: Whether to consider stratification (default: `false`).

Returns:
- If `w` is not provided, returns a `Field{Center, Center, Nothing}` containing the mixed layer depth.
- If `w` is provided, returns a `NamedTuple` `(h, Ew)` containing the mixed layer depth `h` and the
  mean vertical kinetic energy `Ew`.
"""
function MixedLayerDepth(grid, tracers; w=nothing, ΔT=0.2, surface = false, stratification = false)

    arch = architecture(grid)
    
    # Create output fields
    h = Field{Center, Center, Nothing}(grid)
    
    Ew = isnothing(w) ? nothing : Field{Center, Center, Nothing}(grid)

    launch!(arch, h.grid, :xy, 
            _compute_mixed_layer_properties!, 
            h, Ew, h.grid, tracers.T, ΔT, surface, stratification, w)

    fill_halo_regions!(h)
    if !isnothing(Ew)
        fill_halo_regions!(Ew)
        return (h, Ew)
    end

    return h
end


# """
#     MixedLayerDepthOperand

# An operand for calculating the mixed layer depth, defined by a temperature- or
# density-based criterion. This object stores the fields and parameters required
# for the calculation.
# """
# struct MixedLayerDepthOperand{B, FT, G, S, N, W}
#     temperature_or_buoyancy :: B
#     criterion :: FT
#     grid :: G
#     surface :: S
#     stratification :: N
#     vertical_velocity :: W
# end

# # Show a concise summary
# Base.summary(op::MixedLayerDepthOperand) = "MixedLayerDepthOperand"

# """
#     MixedLayerDepth(grid; T, w=nothing, ΔT=0.2, kw...)

# Constructs a `Field{Center, Center, Nothing}` to store the mixed layer depth.
# The calculation is based on a temperature difference criterion `ΔT`.

# The `compute!` function for this field can also compute the depth-averaged
# vertical kinetic energy if the vertical velocity field `w` is provided.

# Arguments
# =========
# - `grid`: The grid on which to compute the mixed layer depth.
# - `tracers`: The tracer fields (`NamedTuple`) on which to base the calculation.
# - `w`: (Optional) The vertical velocity field (`Field`). Provide this to enable kinetic energy calculation.
# - `ΔT`: The temperature difference criterion for the mixed layer depth, i.e., h where T(z=0) - T(z=-h) = ΔT.
# - `surface`: Whether to consider the surface layer (default: `false`).
# - `stratification`: Whether to consider stratification (default: `false`).
# - `kw...`: Additional keywords passed to the `Field` constructor.
# """
# function MixedLayerDepth(grid, tracers; w=nothing, ΔT = 0.2, surface = false, stratification = false, kw...)
#     operand = MixedLayerDepthOperand(tracers.T, abs(ΔT), grid, surface, stratification, w)
#     return Field{Center, Center, Nothing}(grid; operand, kw...)
# end

# # Define a new type alias for clarity
# const MixedLayerDepthField = Field{Center, Center, Nothing, <:MixedLayerDepthOperand}

# """
#     compute!(h::MixedLayerDepthField, Ew=nothing)

# Compute the mixed layer depth `h`. If an additional field `Ew` is provided,
# also compute the depth-averaged vertical kinetic energy integrated within the mixed layer.

# The results are stored in `h` and `Ew`.
# """
# function compute!(h::MixedLayerDepthField, time=nothing)
#     arch = architecture(h)

#     # Unpack operand fields
#     b    = h.operand.temperature_or_buoyancy
#     Δb   = h.operand.criterion
#     surface = h.operand.surface
#     stratification = h.operand.stratification
#     w    = h.operand.vertical_velocity

#     Ew = isnothing(w) ? nothing : Field{Center, Center, Nothing}(h.grid)

#     launch!(arch, h.grid, :xy, _compute_mixed_layer_properties!, h, Ew, h.grid, b, Δb, surface, stratification, w)

#     fill_halo_regions!(h)
#     if !isnothing(Ew)
#         fill_halo_regions!(Ew)
#         return (h=h, Ew=Ew)
#     end

#     return h
# end

# function MixedLayerN²(grid, tracers; w=nothing,ΔT = 0.2, surface = false, stratification = true, kw...)
#     operand = MixedLayerDepthOperand(tracers.T, abs(ΔT), grid, surface, stratification, w)
#     return Field{Center, Center, Nothing}(grid; operand, kw...)
# end

# const MixedLayerN²Field = Field{Center, Center, Nothing, <:MixedLayerDepthOperand}

# function computeNh!(Nh::MixedLayerN²Field, time=nothing)
#     arch = architecture(Nh)
#     b    = Nh.operand.temperature_or_buoyancy
#     Δb   = Nh.operand.criterion
#     surface = Nh.operand.surface
#     stratification = Nh.operand.stratification
#     w    = Nh.operand.vertical_velocity
#     Ew = isnothing(w) ? nothing : Field{Center, Center, Nothing}(Nh.grid)
#     launch!(arch, Nh.grid, :xy, compute_mixed_layer_properties!, Nh, Ew, Nh.grid, b, Δb, surface, stratification, w)
#     fill_halo_regions!(Nh)
#     if !isnothing(Ew)
#         fill_halo_regions!(Ew)
#         return (Nh=Nh, Ew=Ew)
#     end
#     return Nh
# end