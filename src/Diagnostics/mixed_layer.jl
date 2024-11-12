const c = Center()
const f = Face()

@inline z_bottom(i, j, grid) = znode(i, j, 1, grid, c, c, f)
@inline bottom(i, j, grid)   = znode(i, j, 1, grid, c, c, f)

#####
##### MixedLayerDepthField
#####

# b can be temperature (T) or density (ρ)
@kernel function compute_mld!(h, grid, b, Δb, surface, stratification)
    i, j = @index(Global, NTuple)

    Nz = grid.Nz
    
    if surface
        k_start = Nz - 1
        z_ij = znode(i, j, k_start, grid, c, c, f)
    else    
        k_start   = Nz - 2
        z_ij = znode(i, j, k_start, grid, c, c, f)
        while z_ij > -9.99
            k_start = k_start-1
            z_ij = znode(i, j, k_start, grid, c, c, f)
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
            break
        end
    end

    # Note "-" since `h` is supposed to be "depth" rather than "height"
    if stratification
        @inbounds h[i, j, 1] = Nh²
    else
        @inbounds h[i, j, 1] = - z_ij
    end
end

struct MixedLayerDepthOperand{B, FT, G, S, N}
    temperature_operation :: B
    mixed_layer_temperature_differential :: FT
    grid :: G
    surface :: S
    stratification :: N
end

Base.summary(op::MixedLayerDepthOperand) = "MixedLayerDepthOperand"

function MixedLayerDepth(grid, tracers; ΔT = 0.2, surface = false, stratification = false, kw...)
    operand = MixedLayerDepthOperand(tracers.T, abs(ΔT), grid, surface, stratification)
    return Field{Center, Center, Nothing}(grid; operand, kw...)
end

const MixedLayerDepthField = Field{Center, Center, Nothing, <:MixedLayerDepthOperand}

function compute!(h::MixedLayerDepthField, time=nothing)
    arch = architecture(h)
    b    = h.operand.temperature_operation
    Δb   = h.operand.mixed_layer_temperature_differential
    surface = h.operand.surface
    stratification = h.operand.stratification
    launch!(arch, h.grid, :xy, compute_mld!, h, h.grid, b, Δb, surface, stratification)
    fill_halo_regions!(h)
    return h
end

function MixedLayerN²(grid, tracers; ΔT = 0.2, surface = false, stratification = true, kw...)
    operand = MixedLayerDepthOperand(tracers.T, abs(ΔT), grid, surface, stratification)
    return Field{Center, Center, Nothing}(grid; operand, kw...)
end

const MixedLayerN²Field = Field{Center, Center, Nothing, <:MixedLayerDepthOperand}

function computeNh!(Nh::MixedLayerN²Field, time=nothing)
    arch = architecture(Nh)
    b    = Nh.operand.temperature_operation
    Δb   = Nh.operand.mixed_layer_temperature_differential
    surface = Nh.operand.surface
    stratification = Nh.operand.stratification
    launch!(arch, h.grid, :xy, compute_mld!, Nh, Nh.grid, b, Δb, surface, stratification)
    fill_halo_regions!(h)
    return h
end