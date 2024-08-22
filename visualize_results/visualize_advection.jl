using LESStudySetup
using CairoMakie
using SixelTerm
using Printf, Dates
using Oceananigans: compute!,Units

using LESStudySetup.Diagnostics
using LESStudySetup.Diagnostics: load_snapshots
set_theme!(Theme(fontsize = 12))

function visualize_advection(cooling, wind, dTf, a)
    # Examples! (fill in the correct filename and metadata filename)
    # cooling, wind, dTf = 25, 0.02, -1
    cooling = @sprintf("%03d", cooling)
    wind = replace("$(wind)","." => "" )
    a = replace("$(a)","." => "" )
    if dTf < 0
        fileparams = "four_vortices_cooling_$(cooling)_wind_$(wind)"
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

    # load all the data!!
    snapshots = load_snapshots(filename; metadata)

    simulation = LESStudySetup.idealized_setup(CPU();stop_time = 1Units.days, hydrostatic_approximation = true)
    U,V,W =  simulation.model.velocities
    
    # Let's pick the last snapshot!
    times = snapshots[:T].times
    snapshot_number = length(times)÷2
    nday = @sprintf("%.0f", (times[snapshot_number])/60^2/24)
    println("Plotting snapshot $snapshot_number on day $(nday)...")
    t0 = now()


    u = snapshots[:u][snapshot_number];
    v = snapshots[:v][snapshot_number];
    w = snapshots[:w][snapshot_number];
    println("Loading fields wall time: $((now() - t0).value/1e3) seconds.")

    advection0 = compute!(Field(U*∂x(V) + V*∂y(V) + W*∂z(V)));
    advection1 = compute!(Field(u*∂x(v) + v*∂y(v) + w*∂z(v)));
    advection2 = compute!(Field(U*∂x(v) + V*∂y(v) + W*∂z(v)));
    advection3 = compute!(Field(u*∂x(V) + v*∂y(V) + w*∂z(V)));

    # Coordinate arrays
    xa, ya, za = nodes(advection0)
    k = 113
    
    # Plot advections
    fig = Figure(size = (900, 800));
    vbnd = 1
    ax = Axis(fig[1, 1][1,1]; title = L"\mathbf{U}\cdot \nabla V")
    hm1 = heatmap!(ax, 1e-3xa, 1e-3ya, interior(advection0,:,:,k), rasterize = true,colormap = :balance, colorrange = (0, vbnd))
    Colorbar(fig[1, 1][1, 2], hm1)
    ax = Axis(fig[1, 2][1,1]; title = L"\mathbf{u}'\cdot \nabla v'")
    hm2 = heatmap!(ax, 1e-3xa, 1e-3ya, interior(advection1,:,:,k), rasterize = true,colormap = :balance, colorrange = (0, vbnd))
    Colorbar(fig[1, 2][1, 2], hm2)
    ax = Axis(fig[2, 1][1,1]; title = L"\mathbf{U}\cdot \nabla v'")
    hm3 = heatmap!(ax, 1e-3xa, 1e-3ya,interior(advection2,:,:,k), rasterize = true,colormap = :balance, colorrange = (0, vbnd))
    Colorbar(fig[2, 1][1, 2], hm3)
    ax = Axis(fig[2, 2][1,1]; title = L"\mathbf{u}'\cdot \nabla V")
    hm4 = heatmap!(ax, 1e-3xa, 1e-3ya, interior(advection3,:,:,k), rasterize = true,colormap = :balance, colorrange = (0, vbnd))
    Colorbar(fig[2, 2][1, 2], hm4)
    save("results/advectiony_" * fileparams * "_d$(nday)_h$(-za[k]).pdf", fig)
    
    fig = Figure(size = (900, 800));
    vbnd = 1
    ax = Axis(fig[1, 1][1,1]; title = L"|\mathbf{U}\cdot \nabla V|/| \mathbf{u}'\cdot \nabla v' |")
    hm1 = heatmap!(ax, 1e-3xa, 1e-3ya, abs.(interior(advection0,:,:,k)./interior(advection1,:,:,k)), rasterize = true,colormap = :balance, colorrange = (0, vbnd))
    Colorbar(fig[1, 1][1, 2], hm1)
    ax = Axis(fig[1, 2][1,1]; title = L"|\mathbf{u}'\cdot \nabla V|/| \mathbf{u}'\cdot \nabla v' |")
    hm2 = heatmap!(ax, 1e-3xa, 1e-3ya, abs.(interior(advection3,:,:,k)./interior(advection1,:,:,k)), rasterize = true,colormap = :balance, colorrange = (0, vbnd))
    Colorbar(fig[1, 2][1, 2], hm2)
    ax = Axis(fig[2, 1][1,1]; title = L"|\mathbf{U}\cdot \nabla V|/| \mathbf{U}\cdot \nabla v' |")
    hm3 = heatmap!(ax, 1e-3xa, 1e-3ya,abs.(interior(advection0,:,:,k)./interior(advection2,:,:,k)), rasterize = true,colormap = :balance, colorrange = (0, vbnd))
    Colorbar(fig[2, 1][1, 2], hm3)
    ax = Axis(fig[2, 2][1,1]; title = L"|\mathbf{u}'\cdot \nabla V|/| \mathbf{U}\cdot \nabla v' |")
    hm4 = heatmap!(ax, 1e-3xa, 1e-3ya, abs.(interior(advection3,:,:,k)./interior(advection2,:,:,k)), rasterize = true,colormap = :balance, colorrange = (0, vbnd))
    Colorbar(fig[2, 2][1, 2], hm4)
    save("results/advectiony_" * fileparams * "_d$(nday)_h$(-za[k]).pdf", fig)
    return
end

cooling, wind, dTf, a = 50, 0.1, 1, 1.0
visualize_advection(cooling, wind, dTf,a)