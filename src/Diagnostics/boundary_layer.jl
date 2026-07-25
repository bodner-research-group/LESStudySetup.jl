const c = Center()
const f = Face()

#####
##### BoundaryLayerDepthField
#####

# κ is diffusivity, and κc is diffusivity at the boundary layer
@kernel function compute_bld!(h, grid, κ, κc)
    i, j = @index(Global, NTuple)

    Nz = grid.Nz
    
    k_start   = Nz - 1
    z_ij = znode(i, j, k_start, grid, c, c, f)

    @unroll for k in k_start : -1 : 1 # scroll from point just below surface

        κᵏ = @inbounds κ[i, j, k]
        κ⁺ = @inbounds κ[i, j, k+1]

        zᵏ = znode(i, j, k, grid, c, c, c)
        Δz⁺ = Δzᶜᶜᶠ(i, j, k+1, grid)

        # Assuming κ increases upwards
        # Linearly interpolate to find boundary layer height
        inside_boundary_layer = (κ⁺ > κc) & (κᵏ > κc)
        just_below_boundary_layer = (κ⁺ > κc) & (κᵏ <= κc)
        new_z_ij = zᵏ + (κc - κᵏ) / (κ⁺ - κᵏ) * Δz⁺
        
        # Replace z_ij if we found a new boundary layer depth
        replace_z = (just_below_boundary_layer | inside_boundary_layer) & !inactive_node(i, j, k, grid, c, c, c)
        z_ij = ifelse(replace_z, new_z_ij, z_ij)
        if just_below_boundary_layer 
            break
        end
    end

    # Note "-" since `h` is supposed to be "depth" rather than "height"
    @inbounds h[i, j, 1] = - z_ij
end

struct BoundaryLayerDepthOperand{B, FT, G}
    diffusivity_operation :: B
    boundary_layer_criterion :: FT
    grid :: G
end

Base.summary(op::BoundaryLayerDepthOperand) = "BoundaryLayerDepthOperand"

function BoundaryLayerDepth(grid, fields; κc = 1e-5, kw...)
    operand = BoundaryLayerDepthOperand(fields.κ, κc, grid)
    return Field{Center, Center, Nothing}(grid; operand, kw...)
end

const BoundaryLayerDepthField = Field{Center, Center, Nothing, <:BoundaryLayerDepthOperand}

function compute!(h::BoundaryLayerDepthField, time=nothing)
    arch = architecture(h)
    κ    = h.operand.diffusivity_operation
    κc  = h.operand.boundary_layer_criterion
    launch!(arch, h.grid, :xy, compute_bld!, h, h.grid, κ, κc)
    fill_halo_regions!(h)
    return h
end
