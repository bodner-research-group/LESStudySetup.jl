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
set_theme!(Theme(fontsize = 12))

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
snapshots = load_snapshots(filename; metadata)
initsnaps = load_snapshots(initfile)
u0 = initsnaps[:u][1]
v0 = initsnaps[:v][1]

# Let's pick the last snapshot!
times = snapshots[:T].times
snapshot_number = length(times)÷2 + 4
nday = @sprintf("%2.0f", (times[snapshot_number])/60^2/24)
println("Start on day $(nday)...")

N̅² = mean(compute!(Field(N²(snapshots, snapshot_number))), dims=(1,2))
Γx, Γy = M²(snapshots, snapshot_number)
Γx, Γy = mean(compute!(Field(Γx)), dims=(1,2)), mean(compute!(Field(Γy)), dims=(1,2))

for i = 1:9
    N̅² += mean(compute!(Field(N²(snapshots, snapshot_number+8i))), dims=(1,2))
    Γxi, Γyi = M²(snapshots, snapshot_number+8i)
    Γx += mean(compute!(Field(Γxi)), dims=(1,2))
    Γy += mean(compute!(Field(Γyi)), dims=(1,2))
    println("another $i day(s) added")
end
N̅² = compute!(Field(N̅²/10))
Γx, Γy = compute!(Field(Γx/10)), compute!(Field(Γy/10))

function bvarainces(snapshot_number)
    α = parameters.α
    g = parameters.g
    T = snapshots[:T][snapshot_number]
    u = compute!(Field(snapshots[:u][snapshot_number] + u0))
    v = compute!(Field(snapshots[:v][snapshot_number] + v0))
    w = snapshots[:w][snapshot_number]
    xT, yT, zT = nodes(T)
    _, _, zw = nodes(w)
    B = compute!(Field(α * g * T))
    xT, yT = xT .- mean(xT), yT' .- mean(yT)

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
    uᵖ = compute!(Field(u - mean(u, dims = (1,2))))
    vᵖ = compute!(Field(v - mean(v, dims = (1,2))))

    println("compute buoyacy variance term 1...")
    var1 = mean(compute!(Field(b̃ᵖ * (u * ∂x(b̃ᵖ) + v * ∂y(b̃ᵖ) + w * ∂z(b̃ᵖ)))), dims = (1,2))
    println("compute buoyacy variance term 2...")
    var2 = compute!(Field(mean(uᵖ * b̃ᵖ, dims = (1,2)) * Γx + mean(vᵖ * b̃ᵖ, dims = (1,2)) * Γy))
    w = compute!(Field(@at (Center, Center, Center) snapshots[:w][snapshot_number]))
    wᵖ = compute!(Field(w - mean(w, dims = (1,2))))
    wbp = mean(wᵖ * b̃ᵖ, dims = (1,2))
    println("compute buoyacy variance term 3...")
    var3 = compute!(Field((N̅² + mean(∂z(b̃ᵖ), dims = (1,2))) * wbp))
    return wbp, var1, var2, var3, zT, zw
end

wbp, var1, var2, var3, zT, zw = bvarainces(snapshot_number)
for i = 1:9
    wbi, var1i, var2i, var3i, _, _ = bvarainces(snapshot_number+8i)
    wbp += wbi
    var1 += var1i
    var2 += var2i
    var3 += var3i
    println("another $i day(s) added")
end
wbp = compute!(Field(wbp/10))
var1, var2, var3 = compute!(Field(var1/10)), compute!(Field(var2/10)), compute!(Field(var3/10))

# Plot the buoyancy variance terms
fig = Figure(size = (800, 900))

axis_kwargs = (xlabel = "Wavenumber (m⁻¹)", ylabel = "z (m)",
               xscale = log10, limits = (nothing, (-250, 0)))
ax = Axis(fig[1, 1]; title=L"\overline{\langle w^\prime \tilde{b}^\prime \rangle}", xlabel = "(m²/s³)", limits = (nothing, (-250, 0)))
lines!(ax, vec(wbp), zT)
ax = Axis(fig[1, 2]; title="buoyancy variance terms", xlabel = "(m²/s)", limits = (nothing, (-250, 0)))
lines!(ax, vec(var1), zT; label = L"\overline{\langle \mathbf{u} \cdot \nabla (\tilde{b}^{\prime 2}/2) \rangle}")
lines!(ax, vec(var2), zT; label = L"\overline{\langle \mathbf{u}^\prime \tilde{b}^\prime \rangle \cdot \mathbf{\Gamma}}")
lines!(ax, vec(var1 + var2), zT; label = L"\overline{\langle \mathbf{u} \cdot \nabla (\tilde{b}^{\prime 2}/2) \rangle} +\overline{\langle \mathbf{u}' \tilde{b}' \rangle \cdot \mathbf{\Gamma}}")
lines!(ax, vec(var3), zw; label = L"\overline{( N^2 +\langle \partial_z \tilde{b}^\prime \rangle) \langle w^\prime \tilde{b}^\prime \rangle}")
hideydecorations!(ax, ticks = false)
axislegend(ax, position = :rb)

save(filesave * "bvaraince_" * fileparams * "_d$(nday).pdf", fig)
