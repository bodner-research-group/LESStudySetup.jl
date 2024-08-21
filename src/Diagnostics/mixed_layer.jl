const c = Center()
const f = Face()

@inline z_bottom(i, j, grid) = znode(i, j, 1, grid, c, c, f)
@inline bottom(i, j, grid)   = znode(i, j, 1, grid, c, c, f)

#####
##### MixedLayerDepthField
#####

# b can be temperature (T) or density (ρ)
@kernel function compute_mld!(h, grid, b, Δb)
    i, j = @index(Global, NTuple)

    Nz = grid.Nz
    
    k_start   = Nz - 2
    z_ij = znode(i, j, k_start, grid, c, c, f)
    while z_ij > -9.99
        k_start = k_start-1
        z_ij = znode(i, j, k_start, grid, c, c, f)
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
            break
        end
    end

    # Note "-" since `h` is supposed to be "depth" rather than "height"
    @inbounds h[i, j, 1] = - z_ij
end

struct MixedLayerDepthOperand{B, FT, G}
    temperature_operation :: B
    mixed_layer_temperature_differential :: FT
    grid :: G
end

Base.summary(op::MixedLayerDepthOperand) = "MixedLayerDepthOperand"

function MixedLayerDepth(grid, tracers; ΔT = 0.2, kw...)
    operand = MixedLayerDepthOperand(tracers.T, abs(ΔT), grid)
    return Field{Center, Center, Nothing}(grid; operand, kw...)
end

const MixedLayerDepthField = Field{Center, Center, Nothing, <:MixedLayerDepthOperand}

function compute!(h::MixedLayerDepthField, time=nothing)
    arch = architecture(h)
    b    = h.operand.temperature_operation
    Δb   = h.operand.mixed_layer_temperature_differential
    launch!(arch, h.grid, :xy, compute_mld!, h, h.grid, b, Δb)
    fill_halo_regions!(h)
    return h
end