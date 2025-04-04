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
using LESStudySetup.Diagnostics: load_snapshots, filtered_budgets
using LESStudySetup.Diagnostics: N², M²
set_theme!(theme_latexfonts(), fontsize=12,figure_padding = 5)

cooling, wind, dTf,a = 50, 0.1, -1,1.0
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

# Let's pick the mid snapshot!
times = snapshots[:T].times
snapshot_number = length(times)÷2 + 1
nday = @sprintf("%3.1f", (times[snapshot_number])/60^2/24)
println("Start on day $(nday)...")

# Compute the budgets for the mid snapshot
E̅ₖ, Eₖᵖ, mke, fke, E̅ₚ, Eₚᵖ, mpe, fpe, Σ̅b, Σbᵖ, mbv, fbv = filtered_budgets(snapshots, snapshot_number; u0 = u0, v0 = v0);

E̅ₖ2, Eₖᵖ2, E̅ₚ2, Eₚᵖ2, Σ̅b2, Σbᵖ2 = filtered_budgets(snapshots, snapshot_number+1; budget=false,u0 = u0, v0 = v0);
Δt = times[snapshot_number+1] - times[snapshot_number]

# Plot the buoyancy variance terms
_, _, zT = nodes(snapshots[:T][snapshot_number])
havg(x) = vec(mean(x, dims = (1,2)))
fig = Figure(size = (640, 840))
gabc= fig[1, 1] = GridLayout()
axis_kwargs = (ylabel = L"z~\text{(m)}",limits = (nothing, (-200, 0)))
ax_a = Axis(gabc[1, 1]; titlealign = :left, title=L"\text{(a)}~\bar{E}_k\text{ budget}", axis_kwargs...)
#hidexdecorations!(ax_a, ticks = false)
#lines!(ax_a, 1e8*(E̅ₖ2-havg(E̅ₖ))/Δt, zT; label = L"\bar{E}_{k,t}~\text{(1)}")
#lines!(ax_a, -1e8*(mke[1]), zT; label = L"-\nabla \cdot (\bar{\mathbf{u}} \bar{E}_k)~\text{(2)}")
#lines!(ax_a, -1e8*(mke[2]), zT; label = L"-\nabla \cdot (\bar{\mathbf{u}} \bar{p})~\text{(3)}")
lines!(ax_a, -1e8*(mke[3]), zT; label = L"-\nabla \cdot (\bar{\mathbf{u}} \cdot \mathbf{\tau})~\text{(4)}")
#lines!(ax_a, -1e8*(mke[4]), zT; label = L"- \overline{\mathbf{U} \cdot \nabla \mathbf{u}} \cdot \bar{\mathbf{u}}~\text{(5)}")
lines!(ax_a, 1e8*(mke[5]), zT; label = L"\bar{b} \bar{w}~\text{(6)}")
lines!(ax_a, 1e8*(mke[6])[1:end-1], zT[1:end-1]; label = L"\bar{\mathbf{F}}\cdot \bar{\mathbf{u}}~\text{(7)}")
lines!(ax_a, -1e8*(mke[7]), zT; label = L"-\Pi_h~\text{(8)}")
lines!(ax_a, -1e8*(mke[8]), zT; label = L"-\Pi_v~\text{(9)}")
RHS = -(mke[1]+mke[2]+mke[3]+mke[4])+mke[5]+mke[6]-mke[7]-mke[8]
#lines!(ax_a, 1e8*((E̅ₖ2-havg(E̅ₖ))/Δt-RHS)[1:end-1], zT[1:end-1]; label = L"\text(1)-\sum_{i=2}^9\text(i)", color = :black)
axislegend(ax_a, labelsize=9,position = :rb)
ax_b = Axis(gabc[1, 2]; titlealign = :left, title=L"\text{(b)}~\bar{E}_p\text{ budget}", axis_kwargs...)
#hidexdecorations!(ax_b, ticks = false)
hideydecorations!(ax_b, ticks = false)
#lines!(ax_b, 1e8*(E̅ₚ2-havg(E̅ₚ))/Δt, zT; label = L"\bar{E}_{p,t}~\text{(1)}")
lines!(ax_b, -1e8*(mpe[1]), zT; label = L"-\nabla \cdot (\bar{\mathbf{u}} \bar{E}_p)~\text{(2)}")
lines!(ax_b, -1e8*(mpe[2]), zT; label = L"z \overline{\mathbf{U} \cdot \nabla b} ~\text{(3)}")
lines!(ax_b, -1e8*(mke[5]), zT; label = L"-\bar{b} \bar{w}~\text{(4)}")
#lines!(ax_b, 1e8*(mpe[3]), zT; label = L"z \nabla \cdot \mathbf{q}~\text{(5)}")
lines!(ax_b, 1e8*(mpe[4])[1:end-1], zT[1:end-1]; label = L"-z \bar{Q}~\text{(6)}")
RHS = -(mpe[1]+mpe[2]+mke[5])+mpe[3]+mpe[4]
#lines!(ax_b, 1e8*((E̅ₚ2-havg(E̅ₚ))/Δt-RHS)[1:end-1], zT[1:end-1]; label = L"\text(1)-\sum_{i=2}^6\text(i)", color = :black)
axislegend(ax_b, labelsize=9,position = :rb)
ax_c = Axis(gabc[1, 3]; titlealign = :left, title=L"\text{(c)}~\bar{b}^2/2\text{ budget}", axis_kwargs...)
#hidexdecorations!(ax_c, ticks = false)
hideydecorations!(ax_c, ticks = false)
#lines!(ax_c, 1e8*(Σ̅b2-havg(Σ̅b))/Δt, zT; label = L"\bar{b}^2_t/2~\text{(1)}")
lines!(ax_c, -1e8*(mbv[1]), zT; label = L"-\nabla \cdot (\bar{\mathbf{u}} \bar{b}^2/2)~\text{(2)}")
lines!(ax_c, -1e8*(mbv[2]), zT; label = L"-\bar{b} \overline{\mathbf{U} \cdot \nabla b} ~\text{(3)}")
lines!(ax_c, 1e8*(mbv[3]), zT; label = L"-\bar{b} \nabla \cdot \mathbf{q}~\text{(4)}")
lines!(ax_c, 1e8*(mbv[4])[1:end-1], zT[1:end-1]; label = L"\bar{b} \bar{Q}\cdot~\text{(5)}")
RHS = -(mbv[1]+mbv[2])+mbv[3]+mbv[4]
#lines!(ax_c, 1e8*((Σbᵖ2-havg(Σbᵖ))/Δt-RHS)[1:end-1], zT[1:end-1]; label = L"\text(1)-\sum_{i=2}^5\text(i)", color = :black)
axislegend(ax_c, labelsize=9, position = :rb)

ax_a = Axis(gabc[2, 1]; xlabel = L"\text{(10^8~m^2 s^{-3})}", titlealign = :left, title=L"\text{(d)}~{E}^\prime_k\text{ budget}", axis_kwargs...)
#lines!(ax_a, 1e8*(Eₖᵖ2-havg(Eₖᵖ))/Δt, zT; label = L"E^\prime_{k,t}~\text{(1)}")
lines!(ax_a, -1e8*(fke[1]), zT; label = L"-\nabla \cdot (\mathbf{u} E^\prime_k)~\text{(2)}")
lines!(ax_a, -1e8*(fke[2]), zT; label = L"-\nabla \cdot ({\mathbf{u}}^\prime {p}^\prime)~\text{(3)}")
#lines!(ax_a, -1e8*(fke[3]), zT; label = L"-\nabla \cdot ({\mathbf{u}}^\prime \cdot \mathbf{\tau})~\text{(4)}")
#lines!(ax_a, -1e8*(fke[4]), zT; label = L"- (\mathbf{U} \cdot \nabla \mathbf{u})^\prime \cdot \mathbf{u}^\prime~\text{(5)}")
lines!(ax_a, 1e8*(fke[6]), zT; label = L"b^\prime w^\prime~\text{(6)}")
lines!(ax_a, 1e8*(fke[5])[1:end-1], zT[1:end-1]; label = L"{\mathbf{F}}^\prime\cdot {\mathbf{u}}^\prime~\text{(7)}")
lines!(ax_a, 1e8*(fke[7]), zT; label = L"\mathrm{Tr}_h~\text{(8)}")
lines!(ax_a, 1e8*(fke[8]), zT; label = L"\mathrm{Tr}_v~\text{(9)}")
#lines!(ax_a, 1e8*(fke[9]), zT; label = L"\mathcal{P}_h~\text{(10)}")
#lines!(ax_a, 1e8*(fke[10]), zT; label = L"\mathcal{P}_v~\text{(11)}")
RHS = -(fke[1]+fke[2]+fke[3]+fke[4])+fke[5]+fke[6]+fke[7]+fke[8]+fke[9]+fke[10]
#lines!(ax_a, 1e8*((Eₖᵖ2-havg(Eₖᵖ))/Δt-RHS)[1:end-1], zT[1:end-1]; label = L"\text(1)-\sum_{i=2}^11\text(i)", color = :black)
axislegend(ax_a, labelsize=9,position = :rb)
ax_b = Axis(gabc[2, 2]; xlabel = L"\text{(10^8~m^2 s^{-3})}", titlealign = :left, title=L"\text{(e)}~{E}^\prime_p\text{ budget}", axis_kwargs...)
hideydecorations!(ax_b, ticks = false)
#lines!(ax_b, 1e8*(Eₚᵖ2-havg(Eₚᵖ))/Δt, zT; label = L"{E}^\prime_{p,t}~\text{(1)}")
lines!(ax_b, -1e8*(fpe[1]), zT; label = L"-\nabla \cdot ({\mathbf{u}} {E}^\prime_p)~\text{(2)}")
lines!(ax_b, -1e8*(fpe[2]), zT; label = L"z (\mathbf{U} \cdot \nabla b)^\prime ~\text{(3)}")
lines!(ax_b, -1e8*(fke[6]), zT; label = L"-b^\prime w^\prime~\text{(4)}")
lines!(ax_b, 1e8*(fpe[3]), zT; label = L"z \mathbf{u}^\prime \cdot \nabla \bar{b} ~\text{(5)}")
#lines!(ax_b, -1e8*(mpe[3]), zT; label = L"-z \nabla \cdot \mathbf{q}~\text{(6)}")
lines!(ax_b, 1e8*(fpe[4])[1:end-1], zT[1:end-1]; label = L"-z {Q}^\prime\cdot~\text{(7)}")
RHS = -(fpe[1]+fpe[2]+fke[6]+mpe[3])+fpe[3]+fpe[4]
#lines!(ax_b, 1e8*((Eₚᵖ2-havg(Eₚᵖ))/Δt-RHS)[1:end-1], zT[1:end-1]; label = L"\text(1)-\sum_{i=2}^7\text(i)", color = :black)
axislegend(ax_b, labelsize=9,position = :rb)
ax_c = Axis(gabc[2, 3]; xlabel = L"\text{(10^8~m^2 s^{-5})}", titlealign = :left, title=L"\text{(f)}~b^{\prime 2}/2\text{ budget}", axis_kwargs...)
hideydecorations!(ax_c, ticks = false)
l#ines!(ax_c, 1e8*(Σbᵖ2-havg(Σbᵖ))/Δt, zT; label = L"{b}^{\prime 2}_t/2~\text{(1)}")
lines!(ax_c, -1e8*(fbv[1]), zT; label = L"-\nabla \cdot ({\mathbf{u}} {b}^{\prime 2}/2)~\text{(2)}")
lines!(ax_c, -1e8*(fbv[2]), zT; label = L"-{b}^\prime ({\mathbf{U} \cdot \nabla b})^\prime ~\text{(3)}")
lines!(ax_c, 1e8*(fbv[3]), zT; label = L"-{b}^\prime \mathbf{u}^\prime \cdot \nabla \bar{b}~\text{(4)}")
lines!(ax_c, 1e8*(fbv[4]), zT; label = L"-{b}^\prime \nabla \cdot \mathbf{q}~\text{(5)}")
lines!(ax_c, 1e8*(fbv[5])[1:end-1], zT[1:end-1]; label = L"b^\prime Q^\prime ~\text{(6)}")
RHS = -(fbv[1]+fbv[2])+fbv[3]+fbv[4]+fbv[5]
#lines!(ax_c, 1e8*((Σbᵖ2-havg(Σbᵖ))/Δt-RHS)[1:end-1], zT[1:end-1]; label = L"\text(1)-\sum_{i=2}^6\text(i)", color = :black)
axislegend(ax_c, labelsize=9, position = :rb)
colgap!(gabc, 1, 20)
colgap!(gabc, 2, 20)
rowgap!(gabc, 1, 5)
save(filesave * "filteredbudgets_" * fileparams * "_d$(nday).pdf", fig)
