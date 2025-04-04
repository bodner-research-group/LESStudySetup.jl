using LESStudySetup
using CairoMakie
using SixelTerm
using Printf, Dates
using Oceananigans: compute!
using Oceananigans.Grids: xnodes, ynodes, znodes
using Oceananigans.Units: kilometer
using Oceananigans.BoundaryConditions
using Statistics: mean
using LESStudySetup.Diagnostics
using LESStudySetup.Diagnostics: load_snapshots
using LESStudySetup.Diagnostics: N², M²
set_theme!(theme_latexfonts(), fontsize=12,figure_padding = 5)

function cumtrapz(X, Y) 
    # Check matching vector length
    @assert length(X) == length(Y)
    # Initialize Output
    out = similar(X)
    out[1] = 0
    # Iterate over arrays
    for i in 2:length(X)
      out[i] = out[i-1] + 0.5*(X[i] - X[i-1])*(Y[i] + Y[i-1])
    end
    # Return output
    out
end

cooling, wind, dTf,a = 50, 0.1, -1,1.0
# Examples! (fill in the correct filename and metadata filename)
# cooling, wind, dTf = 25, 0.02, -1
cooling = @sprintf("%03d", cooling)
wind = replace("$(wind)","." => "" )
a = replace("$(a)","." => "" )
if dTf < 0
    fileparams = "hydrostatic_twin_simulation"
else
    if length(wind) < 2
        wind = "0" * wind
    end
    dTf = @sprintf("%1d", dTf)
    fileparams = "free_surface_short_test_$(cooling)_wind_$(wind)_dTf_$(dTf)_a_$(a)"
end
filehead = "./"
filename = filehead * "hydrostatic_snapshots_" * fileparams * ".jld2"
metadata = filehead * "experiment_" * fileparams * "_metadata.jld2"
filesave = "./results/"
initfile = filehead * "hydrostatic_snapshots_hydrostatic_background.jld2"

# load all the data!!
println("Loading data from $filename...")
snapshots = load_snapshots(filename; metadata, variables = (:u, :v, :w, :T, :κc, :κu, :pHY′))
initsnaps = load_snapshots(initfile)
u0 = initsnaps[:u][1]
v0 = initsnaps[:v][1]

# Let's pick the last snapshot!
times = snapshots[:T].times
snapshot_number = length(times)÷2 + 1
nday = @sprintf("%3.1f", (times[snapshot_number])/60^2/24)
println("Start on day $(nday)...")

function kebudgets(snapshot_number; nobudget=false, u0 = 0, v0 = 0)
    t0 = time()
    α = parameters.α
    g = parameters.g
    T = snapshots[:T][snapshot_number]
    _, _, zT = nodes(T)
    u = snapshots[:u][snapshot_number]
    v = snapshots[:v][snapshot_number]
    w = snapshots[:w][snapshot_number]
    grid = u.grid
    b = compute!(Field(α * g * T))
    κu = snapshots[:κu][snapshot_number]
    p = snapshots[:pHY′][snapshot_number]
    px, py, pz = compute!(Field(∂x(p))), compute!(Field(∂y(p))), compute!(Field(∂z(p)))
    fx, fy, fz = compute!(Field(∂z(κu * ∂z(u)))), compute!(Field(∂z(κu * ∂z(v)))), compute!(Field(∂z(κu * ∂z(w))))
    println("computed variables with time $(time()-t0) seconds...")

    U, V, W, B = mean(u, dims = (1,2)), mean(v, dims = (1,2)), mean(w, dims = (1,2)), mean(b, dims = (1,2))
    uᵖ, vᵖ, wᵖ, bᵖ = compute!(Field(u - U)), compute!(Field(v - V)), compute!(Field(w - W)), compute!(Field(b - B))
    Px, Py, Pz = mean(px, dims = (1,2)), mean(py, dims = (1,2)), mean(pz, dims = (1,2))
    pxᵖ, pyᵖ, pzᵖ = compute!(Field(px - Px)), compute!(Field(py - Py)), compute!(Field(pz - Pz))
    Fx, Fy, Fz = mean(fx, dims = (1,2)), mean(fy, dims = (1,2)), mean(fz, dims = (1,2))
    fxᵖ, fyᵖ, fzᵖ = compute!(Field(fx - Fx)), compute!(Field(fy - Fy)), compute!(Field(fz - Fz))
    println("computed mean and anoamalous terms with time $(time()-t0) seconds...")

    mke = compute!(Field(0.5 * (U^2 + V^2)))
    ke = compute!(Field(0.5 * (uᵖ^2 + vᵖ^2)))
    kef = compute!(Field(ke + 0.5 * wᵖ^2))

    if nobudget
        return mke, mean(ke, dims = (1,2)), mean(kef, dims = (1,2))
    else
        sfilename = "./hydrostatic_free_surface_" * fileparams * ".jld2"
        η = FieldTimeSeries(sfilename, string(:η); architecture=CPU(), backend = OnDisk())[snapshot_number]
        η3d = ZFaceField(grid)
        set!(η3d, repeat(interior(η,:,:,1), outer=[1,1,length(zT)+1]))
        ηxᵖ, ηyᵖ = compute!(Field(∂x(η3d))), compute!(Field(∂y(η3d)))

        mke1 = compute!(Field(-U * mean(compute!(Field(∂z(w * uᵖ))), dims = (1,2)) - mean(compute!(Field(∂z(w * vᵖ))), dims = (1,2)) * V))
        mke2 = compute!(Field(U * (Fx) + (Fy) * V))
        mke3 = compute!(Field(U * (- Px) + (- Py) * V))
        mke4 = compute!(Field(B * W))
        mke5 = compute!(Field(@at (Nothing, Nothing, Center) W * ∂z(mke)))
        println("computed mean KE terms with time $(time()-t0) seconds...")
        mkeb = [mke1, mke2, mke3, mke4, mke5]
        
        ke2 = compute!(Field(mean(compute!(Field(-uᵖ * wᵖ)), dims = (1,2)) * ∂z(U) - mean(compute!(Field(wᵖ * vᵖ)), dims = (1,2)) * ∂z(V)))
        println("compute eddy KE term 2 with time $(time()-t0) seconds...")
        ke3 = mean(compute!(Field(bᵖ * wᵖ)), dims = (1,2))
        println("compute eddy KE term 3 with time $(time()-t0) seconds...")
        ke4 = mean(compute!(Field(uᵖ * (fxᵖ - pxᵖ - g * ηxᵖ) + (fyᵖ - pyᵖ - g * ηyᵖ) * vᵖ)), dims = (1,2))
        ke4f = compute!(Field(ke4 + mean(compute!(Field(wᵖ * (fzᵖ - pzᵖ))), dims = (1,2))))
        println("compute eddy KE term 4 with time $(time()-t0) seconds...")

        ue, ve = compute!(Field(ke * (u + u0))), compute!(Field(ke * (v + v0)))
        D = KernelFunctionOperation{Center, Center, Center}(Oceananigans.Operators.div_xyᶜᶜᶜ, grid, ue,ve)
        ke1 = mean(compute!(Field(D + ∂z(w * ke))), dims = (1,2))
        ue, ve = compute!(Field(kef * (u + u0))), compute!(Field(kef * (v + v0)))
        D = KernelFunctionOperation{Center, Center, Center}(Oceananigans.Operators.div_xyᶜᶜᶜ, grid, ue,ve)
        ke1f = mean(compute!(Field(D + ∂z(w * kef))), dims = (1,2))
        println("compute eddy KE term 1 with time $(time()-t0) seconds...")
        keb = [ke1, ke2, ke3, ke4, ke1f, ke4f]

        return mke, mkeb, mean(ke, dims = (1,2)), mean(kef, dims = (1,2)), keb, zT
    end

end

# Compute the budgets for the last snapshot
mke, mkeb, ke, kef, keb, zT = kebudgets(snapshot_number; nobudget=false, u0 = u0, v0 = v0)

mke2, ke2, kef2 = kebudgets(snapshot_number+1; nobudget=true)
Δt = times[snapshot_number+1] - times[snapshot_number]

# Plot the buoyancy variance terms
fig = Figure(size = (640, 640))
gabc= fig[1, 1] = GridLayout()
axis_kwargs = (ylabel = L"z~\text{(m)}",limits = (nothing, (-200, 0)))
ax_a = Axis(gabc[1, 1]; xlabel = L"\text{(10^8~m^2 s^{-3})}", titlealign = :left, title=L"\text{(a)}~ \text{Mean KE budget}", axis_kwargs...)
lines!(ax_a, 1e8*(vec(mke2)-vec(mke))/Δt, zT; label = L"|\langle \mathbf{u} \rangle|^2_t/2~\text{(1)}")
lines!(ax_a, 1e8*vec(mkeb[5]), zT; label = L"\langle w \rangle  |\langle \mathbf{u} \rangle|^2_z/2~\text{(2)}")
lines!(ax_a, 1e8*vec(mkeb[1]), zT; label = L"-\langle \mathbf{u} \rangle \cdot \langle \nabla \cdot (\mathbf{u} \mathbf{u}^\prime)\rangle~\text{(3)}")
lines!(ax_a, 1e8*vec(mkeb[4]), zT; label = L"\langle W \rangle \langle b \rangle~\text{(4)}")
lines!(ax_a, 1e8*vec(mkeb[2])[1:end-1], zT[1:end-1]; label = L"\langle \mathbf{F} \rangle \cdot \langle \mathbf{u} \rangle~\text{(5)}")
lines!(ax_a, 1e8*vec(mkeb[3]), zT; label = L"-\langle \nabla p \rangle \cdot \langle \mathbf{u} \rangle~\text{(6)}")
lines!(ax_a, 1e8*((vec(mke2)-vec(mke))/Δt+vec(mkeb[5])-vec(mkeb[1])-vec(mkeb[2])-vec(mkeb[3])-vec(mkeb[4]))[1:end-1], zT[1:end-1]; label = L"\text(1)+\text(2)-\text(3)-\text(4)-\text(5)-\text(6)", color = :black)
axislegend(ax_a, labelsize=9,position = :rb)
ax_b = Axis(gabc[1, 2]; xlabel = L"\text{(10^8~m^2 s^{-3})}", titlealign = :left, title=L"\text{(b)}~ \text{2D TKE budget}", axis_kwargs...)
hideydecorations!(ax_b, ticks = false)
lines!(ax_b, 1e8*(vec(ke2)-vec(ke))/Δt, zT; label = L"\langle | \mathbf{u}^\prime|^2 \rangle_t/2~\text{(1)}")
lines!(ax_b, 1e8*vec(keb[1]), zT; label = L"\langle \mathbf{u} \cdot \nabla | \mathbf{u}^\prime|^2 \rangle/2~\text{(2)}")
lines!(ax_b, 1e8*vec(keb[2]), zT; label = L"-\langle \mathbf{u}_z \rangle \cdot \langle w^\prime \mathbf{u}^{\prime} \rangle~\text{(3)}")
lines!(ax_b, 1e8*vec(keb[3]), zT; label = L"\langle w^\prime b^\prime \rangle~\text{(4)}")
lines!(ax_b, 1e8*vec(keb[4]), zT; label = L"\langle [\mathbf{F}^{\prime}- (\nabla p)^{\prime}- g(\nabla_h \eta)^{\prime}] \cdot \mathbf{u}^\prime \rangle~\text{(5)}")
lines!(ax_b, 1e8*((vec(ke2)-vec(ke))/Δt+vec(keb[1])-vec(keb[2])-vec(keb[3])-vec(keb[4])), zT; label = L"\text(1)+\text(2)-\text(3)-\text(4)-\text(5)", color = :black)
axislegend(ax_b, labelsize=9,position = :rb)
ax_c = Axis(gabc[1, 3]; xlabel = L"\text{(10^8~m^2 s^{-3})}", titlealign = :left, title=L"\text{(c)}~ \text{3D TKE budget}", axis_kwargs...)
hideydecorations!(ax_c, ticks = false)
lines!(ax_c, 1e8*(vec(kef2)-vec(kef))/Δt, zT; label = L"\langle | \mathbf{u}^\prime|^2 \rangle_t/2~\text{(1)}")
lines!(ax_c, 1e8*vec(keb[5]), zT; label = L"\langle \mathbf{u} \cdot \nabla | \mathbf{u}^\prime|^2 \rangle/2~\text{(2)}")
lines!(ax_c, 1e8*vec(keb[2]), zT; label = L"-\langle \mathbf{u}_z \rangle \cdot \langle w^\prime \mathbf{u}^{\prime} \rangle~\text{(3)}")
lines!(ax_c, 1e8*vec(keb[3]), zT; label = L"\langle w^\prime b^\prime \rangle~\text{(4)}")
lines!(ax_c, 1e8*vec(keb[6]), zT; label = L"\langle [\mathbf{F}^{\prime}- (\nabla p)^{\prime}- g(\nabla_h \eta)^{\prime}] \cdot \mathbf{u}^\prime \rangle~\text{(5)}")
lines!(ax_c, 1e8*((vec(kef2)-vec(kef))/Δt+vec(keb[5])-vec(keb[2])-vec(keb[3])-vec(keb[6])), zT; label = L"\text(1)+\text(2)-\text(3)-\text(4)-\text(5)", color = :black)
axislegend(ax_c, labelsize=9, position = :rb)
colgap!(gabc, 1, 20)
colgap!(gabc, 2, 20)
save(filesave * "kebudgets_" * fileparams * "_d$(nday).pdf", fig)

####################################
N̅² = mean(compute!(Field(N²(snapshots, snapshot_number))), dims=(1,2))
Γx, Γy = M²(snapshots, snapshot_number)
Γx, Γy = mean(compute!(Field(Γx)), dims=(1,2)), mean(compute!(Field(Γy)), dims=(1,2))

t0 = time()
for i = 0:9
    for j = 1:8
        N̅² = compute!(Field(N̅² + mean(compute!(Field(N²(snapshots, snapshot_number+16i+2j))), dims=(1,2))))
        Γxi, Γyi = M²(snapshots, snapshot_number+16i+2j)
        Γx = compute!(Field(Γx + mean(compute!(Field(Γxi)), dims=(1,2))))
        Γy = compute!(Field(Γy + mean(compute!(Field(Γyi)), dims=(1,2))))
        println("$j/8 of the  $(1+i)-th day added $(time()-t0) seconds used ")
    end
end
N̅² = compute!(Field(N̅²/81))
Γx, Γy = compute!(Field(Γx/81)), compute!(Field(Γy/81))
Γxz = compute!(Field(@at (Nothing, Nothing, Center) ∂z(Γx)))
Γyz = compute!(Field(@at (Nothing, Nothing, Center) ∂z(Γy)))

function bvarainces(snapshot_number; budget = false, u0 = 0, v0 = 0)
    α = parameters.α
    g = parameters.g
    T = snapshots[:T][snapshot_number]
    κc = snapshots[:κc][snapshot_number]
    xT, yT, zT = nodes(T)
    B = compute!(Field(α * g * T))
    xT, yT = xT .- mean(xT), yT' .- mean(yT)
    w = snapshots[:w][snapshot_number]
    _, _, zw = nodes(w)

    bz = deepcopy(mean(B, dims=(1,2)))
    bx = deepcopy(mean(B, dims=(2)))
    by = deepcopy(mean(B, dims=(1)))

    d = interior(bz)
    d[1,1,:] = cumtrapz(zw, vec(interior(N̅²)))[2:end]
    set!(bz, d)
    fill_halo_regions!(bz)
    d = interior(bx)
    set!(bx, xT .* interior(Γx))
    fill_halo_regions!(bx)
    set!(by, reshape(yT,(1,:,1)) .* interior(Γy))
    fill_halo_regions!(by)
    
    b̃ = compute!(Field(B - bz - bx - by))

    b̃ᵖ = compute!(Field(b̃ - mean(b̃, dims = (1,2))))
    b̄ = mean(B, dims = (1,2))
    bᵖ = compute!(Field(B - b̄))

    b̃ᵖsq = compute!(Field(b̃ᵖ^2))
    b̄sq = compute!(Field(b̄^2))
    bᵖsq = compute!(Field(bᵖ^2))

    println("computed variances with time $(time()-t0) seconds...")
    if !budget
        return zT, mean(b̃ᵖsq, dims = (1,2)), mean(b̄sq, dims = (1,2)), mean(bᵖsq, dims = (1,2))
    else
        
        w = compute!(Field(@at (Center, Center, Center) snapshots[:w][snapshot_number]))
        w̄ = mean(w, dims = (1,2))
        wᵖ = compute!(Field(w - w̄))
        wbp = mean(b̃ᵖ * wᵖ, dims = (1,2))

        wx, wy = xT .* interior(w), yT .* interior(w)
        wx, wy = wx .- mean(wx, dims = (1,2)), wy .- mean(wy, dims = (1,2))
        wxbp = mean(wx .* interior(b̃ᵖ), dims = (1,2))
        wybp = mean(wy .* interior(b̃ᵖ), dims = (1,2))

        u = snapshots[:u][snapshot_number] 
        v = snapshots[:v][snapshot_number] 
        uᵖ = compute!(Field(u - mean(u, dims = (1,2))))
        vᵖ = compute!(Field(v - mean(v, dims = (1,2))))

        Q = compute!(Field(∂z(κc * ∂z(B))))
        Q̄ = mean(Q, dims = (1,2))
        Qᵖ = compute!(Field(Q - Q̄))
        println("start compute budget terms with time $(time()-t0) seconds...")
        var4 = mean(compute!(Field(b̃ᵖ * Qᵖ)), dims = (1,2))
        println("computed buoyacy variance term 4 with time $(time()-t0) seconds...")
        var3 = compute!(Field(wbp * (N̅² + mean(∂z(b̃ᵖ), dims = (1,2)))))
        println("computed buoyacy variance term 3 with time $(time()-t0) seconds...")
        var2 = compute!(Field(mean(uᵖ * b̃, dims = (1,2)) * Γx + mean(vᵖ * b̃ᵖ, dims = (1,2)) * Γy))
        var2 = interior(var2) + wxbp .* interior(Γxz) + wybp .* interior(Γyz)
        println("computed buoyacy variance term 2 with time $(time()-t0) seconds...")

        mpe1 = mean(compute!(Field(∂z(w * bᵖ))), dims = (1,2))
        mpe1 = compute!(Field(- b̄ * mpe1))
        mpe2 = compute!(Field(b̄ * Q̄))
        println("computed mean PE terms with time $(time()-t0) seconds...")
        mpeb = [mpe1, mpe2]

        pe2 = mean(compute!(Field(- wᵖ * bᵖ)), dims = (1,2))
        println("compute eddy PE term 2 with time $(time()-t0) seconds...")
        pe3 = compute!(Field(pe2 * mean(∂z(b̃ᵖ), dims = (1,2))))
        println("compute eddy PE term 3 with time $(time()-t0) seconds...")
        pe4 = mean(compute!(Field(bᵖ * Qᵖ)), dims = (1,2))
        println("compute eddy PE term 4 with time $(time()-t0) seconds...")

        uc, vc = compute!(Field(b̃ᵖsq * (u + u0))), compute!(Field(b̃ᵖsq * (v + v0)))
        D = KernelFunctionOperation{Center, Center, Center}(Oceananigans.Operators.div_xyᶜᶜᶜ, grid, uc,vc)
        var1 = mean(compute!(Field((D + ∂z(b̃ᵖsq * w)) / 2)), dims = (1,2))
        println("computed buoyacy variance term 1 with time $(time()-t0) seconds...")
        bvb = [var1, var2, var3, var4]
        uc, vc = compute!(Field(bᵖsq * (u + u0))), compute!(Field(bᵖsq * (v + v0)))
        pe1 = mean(compute!(Field((D + ∂z(bᵖsq * w)) / 2)), dims = (1,2))
        println("compute eddy PE term 1 with time $(time()-t0) seconds...")
        peb = [pe1, pe2, pe3, pe4]

        return mean(b̃ᵖsq, dims = (1,2)), mean(b̄sq, dims = (1,2)), mean(bᵖsq, dims = (1,2)), bvb, mpeb, peb, zT, zw
    end
end

t0 = time()
b̃ᵖsq, b̄sq, bᵖsq, bvb, mpeb, peb, zT, zw = bvarainces(snapshot_number;budget=true,u0=u0,v0=v0)
_, b̃ᵖsq2, b̄sq2, bᵖsq2 = bvarainces(snapshot_number+1)
Δt = times[snapshot_number+1] - times[snapshot_number]

# Plot the buoyancy variance terms
fig = Figure(size = (640, 640))
gabc= fig[1, 1] = GridLayout()
axis_kwargs = (ylabel = L"z~\text{(m)}",limits = ((-4,4), (-200, 0)))
ax_a = Axis(gabc[1, 1]; xlabel = L"\text{(10^{13}~m^2 s^{-5})}", titlealign = :left, title=L"\text{(a)}~ \text{Buoyancy variance budget}", axis_kwargs...)
lines!(ax_a, 1e13*(vec(b̃ᵖsq2)-vec(b̃ᵖsq))/2/Δt, zT; label = L"\langle \tilde{b}^{\prime 2} \rangle_t/2~\text{(1)}")
lines!(ax_a, 1e13*vec(bvb[1]), zT; label = L"\langle \mathbf{u} \cdot \nabla (\tilde{b}^{\prime 2}/2) \rangle~\text{(2)}")
lines!(ax_a, 1e13*vec(bvb[2]), zT; label = L"\langle \mathbf{u}^\prime \tilde{b}^\prime \rangle \cdot \mathbf{\Gamma}+\langle (w\mathbf{x})^\prime \tilde{b}^\prime \rangle \cdot \mathbf{\Gamma}_z~\text{(3)}")
lines!(ax_a, 1e13*vec(bvb[3]), zT; label = L"( N^2 +\langle \tilde{b}_z \rangle) \langle w^\prime \tilde{b}^\prime \rangle~\text{(4)}")
lines!(ax_a, 1e13*vec(bvb[4]), zT; label = L"\langle Q^{\prime} \tilde{b}^{\prime}  \rangle~\text{(5)}")
lines!(ax_a, 1e13*(vec(bvb[4])-(vec(b̃ᵖsq2)-vec(b̃ᵖsq))/2/Δt-vec(bvb[1])-vec(bvb[2])-vec(bvb[3])), zT; label = L"\text(5)-\text(1)-\text(2)-\text(3)-\text(4)", color=:black)
axislegend(ax_a, labelsize=9, position = :rb)
ax_b = Axis(gabc[1, 2]; xlabel = L"\text{(10^{13}~m^2 s^{-5})}", titlealign = :left, title=L"\text{(b)}~ \text{Mean PE budget}",limits = ((-100,100), (-200, 0)))
hideydecorations!(ax_b, ticks = false)
lines!(ax_b, 1e13*vec(interior(compute!(Field((b̄sq2-b̄sq)/2/Δt)))), zT; label = L"\langle b \rangle^2_t/2~\text{(1)}")
lines!(ax_b, 1e13*vec(mpeb[1]), zT; label = L"-\langle b \rangle \langle (w^\prime b^\prime)_z \rangle~\text{(2)}")
lines!(ax_b, 1e13*vec(mpeb[2]), zT; label = L"\langle Q \rangle \langle b \rangle~\text{(3)}")
lines!(ax_b, 1e13*(vec(interior(compute!(Field(((b̄sq2-b̄sq)/2/Δt-mpeb[1]-mpeb[2])))))), zT; label = L"\text{(1)}+\text{(2)}-\text{(3)}", color=:black)
axislegend(ax_b, labelsize=9, position = :rb)
ax_c = Axis(gabc[1, 3]; xlabel = L"\text{(10^{13}~m^2 s^{-5})}", titlealign = :left, title=L"\text{(c)}~ \text{Eddy PE budget}", axis_kwargs...)
hideydecorations!(ax_c, ticks = false)
lines!(ax_c, 1e13*vec(interior(compute!(Field((bᵖsq2-bᵖsq)/2/Δt)))), zT; label = L"\langle b^{\prime 2} \rangle_t/2~\text{(1)}")
lines!(ax_c, 1e13*vec(interior(compute!(Field(peb[1])))), zT; label = L"\langle \mathbf{u} \cdot \nabla b^{\prime 2} \rangle/2~\text{(2)}")
lines!(ax_c, 1e13*vec(interior(compute!(Field(peb[3])))), zT; label = L"-\langle \tilde{b}_z \rangle \langle w^\prime b^{\prime} \rangle~\text{(3)}")
lines!(ax_c, 1e13*vec(interior(compute!(Field(peb[2]*N̅²)))), zT; label = L"-N^2 \langle w^\prime b^\prime \rangle~\text{(4)}")
lines!(ax_c, 1e13*vec(interior(compute!(Field(peb[4])))), zT; label = L"\langle Q^{\prime} b^{\prime} \rangle~\text{(5)}")
lines!(ax_c, 1e13*(vec(interior(compute!(Field((bᵖsq2-bᵖsq)/2/Δt+peb[1]-peb[4]-peb[3]-peb[2]*N̅²))))), zT; label = L"\text{(1)}+\text{(2)}-\text{(3)}-\text{(4)}-\text{(5)}", color=:black)
axislegend(ax_c, labelsize=9, position = :rb)
colgap!(gabc, 1, 20)
colgap!(gabc, 2, 20)
ax_b.limits = ((-200,200), (-200, 0))
save(filesave * "bbudgets_" * fileparams * "_d$(nday).pdf", fig)
