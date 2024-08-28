const c = Center()
const f = Face()

#####
##### BoundaryLayerDepthField
#####

# b is buoyancy, and (u,v,w) is velocity vector
@kernel function compute_bld!(h, grid, b, u, v, w, Uₜ²,Ric)
    i, j = @index(Global, NTuple)

    Nz = grid.Nz
    
    k_start   = Nz - 1
    z_ij = znode(i, j, k_start, grid, c, c, f)
   
    k_ref = k_start + 1
    z_ref = znode(i, j, k_ref, grid, c, c, f)
    b_ref = @inbounds b[i, j, k_ref]
    u_ref = @inbounds u[i, j, k_ref]
    v_ref = @inbounds v[i, j, k_ref]
    w_ref = @inbounds w[i, j, k_ref]
    b_ref⁻ = @inbounds b[i, j, k_ref-1]
    u_ref⁻ = @inbounds u[i, j, k_ref-1]
    v_ref⁻ = @inbounds v[i, j, k_ref-1]
    w_ref⁻ = @inbounds w[i, j, k_ref-1]

    @unroll for k in k_start : -1 : 1 # scroll from point just below surface

        bᵏ = @inbounds b[i, j, k]
        uᵏ = @inbounds u[i, j, k]
        vᵏ = @inbounds v[i, j, k]
        wᵏ = @inbounds w[i, j, k]
        b⁺ = @inbounds b[i, j, k+1]
        u⁺ = @inbounds u[i, j, k+1]
        v⁺ = @inbounds v[i, j, k+1]
        w⁺ = @inbounds w[i, j, k+1]

        # Compute bulk Richardson number
        Δb = b_ref - bᵏ
        Δu² = ((u_ref - uᵏ)^2 + (v_ref - vᵏ)^2 + (w_ref - wᵏ)^2)/3
        zᵏ = znode(i, j, k, grid, c, c, f)
        Riᵏ = Δu² + Uₜ²[k] > 0 ? -zᵏ * Δb / (Δu² + Uₜ²[k]) : 0.0 # avoid division by zero
        Δb⁺ = b_ref - b⁺
        Δu²⁺ = ((u_ref - u⁺)^2 + (v_ref - v⁺)^2 + (w_ref - w⁺)^2)/3
        z⁺ = znode(i, j, k+1, grid, c, c, f)
        Ri⁺ = Δu²⁺ + Uₜ²[k+1] > 0 ? -z⁺ * Δb⁺ / (Δu²⁺ + Uₜ²[k+1]) : 0.0 # avoid division by zero

        Δz⁺ = Δzᶜᶜᶠ(i, j, k+1, grid)

        # Assuming Ri decreases upwards
        # Linearly interpolate to find boundary layer height
        inside_boundary_layer = (Ri⁺ < Ric) & (Riᵏ < Ric)
        just_below_boundary_layer = (Ri⁺ < Ric)  & (Riᵏ >= Ric)
        new_z_ij = zᵏ + (Ric - Riᵏ) / (Ri⁺ - Riᵏ) * Δz⁺
        
        # Replace z_ij if we found a new boundary layer depth
        replace_z = (just_below_boundary_layer | inside_boundary_layer) & !inactive_node(i, j, k, grid, c, c, c)
        z_ij = ifelse(replace_z, new_z_ij, z_ij)

        # Update reference bouyancy and velocity
        if 0.1 * zᵏ < z_ref
            new_z_ref = 0.1 * zᵏ
            w1,w2 = z_ref/new_z_ref, (new_z_ref - z_ref)/new_z_ref
            z_ref⁻ = znode(i, j, k_ref-1, grid, c, c, f)
            if z_ref⁻ > new_z_ref
                k_ref -= 1
                b_ref⁻ = @inbounds b[i, j, k_ref-1]
                u_ref⁻ = @inbounds u[i, j, k_ref-1]
                v_ref⁻ = @inbounds v[i, j, k_ref-1]
                w_ref⁻ = @inbounds w[i, j, k_ref-1]
            end
            b_ref = w1 * b_ref + w2 * b_ref⁻
            u_ref = w1 * u_ref + w2 * u_ref⁻
            v_ref = w1 * v_ref + w2 * v_ref⁻
            w_ref = w1 * w_ref + w2 * w_ref⁻
            z_ref = new_z_ref
        end

        if just_below_boundary_layer 
            break
        end
    end

    # Note "-" since `h` is supposed to be "depth" rather than "height"
    @inbounds h[i, j, 1] = - z_ij
end

struct BoundaryLayerDepthOperand{B, U, V, W, Uₜ², FT, G}
    buoyancy_operation :: B
    xvelocity_operation :: U
    yvelocity_operation :: V
    zvelocity_operation :: W
    turbulence_operation :: Uₜ²
    boundary_layer_ri_criterion :: FT
    grid :: G
end

Base.summary(op::BoundaryLayerDepthOperand) = "BoundaryLayerDepthOperand"

function BoundaryLayerDepth(grid, fields; Ric = 0.3, kw...)
    operand = BoundaryLayerDepthOperand(fields.B, fields.uc, fields.vc, fields.wc, fields.Uₜ², abs(Ric), grid)
    return Field{Center, Center, Nothing}(grid; operand, kw...)
end

const BoundaryLayerDepthField = Field{Center, Center, Nothing, <:BoundaryLayerDepthOperand}

function compute!(h::BoundaryLayerDepthField, time=nothing)
    arch = architecture(h)
    b    = h.operand.buoyancy_operation
    u    = h.operand.xvelocity_operation
    v    = h.operand.yvelocity_operation
    w    = h.operand.zvelocity_operation
    Uₜ²  = h.operand.turbulence_operation
    Ric  = h.operand.boundary_layer_ri_criterion
    launch!(arch, h.grid, :xy, compute_bld!, h, h.grid, b, u, v, w, Uₜ², Ric)
    fill_halo_regions!(h)
    return h
end