using LESStudySetup
using CairoMakie
using Printf, Dates
using Statistics: mean, std
using LESStudySetup.Diagnostics
using LESStudySetup.Diagnostics: load_distributed_checkpoint,load_subdomain_snapshot
using LESStudySetup.Diagnostics: isotropic_powerspectrum
using LESStudySetup.Diagnostics: coarse_graining!
set_theme!(theme_latexfonts(), fontsize=12, figure_padding = 10)
using Contour, Interpolations, Statistics

"""
    analyze_frontal_properties(T)

Analyzes properties along a frontal outcrop identified from a temperature field.

Arguments
=========
- `T`: The 3D temperature field (`Field` object) from Oceananigans.

Returns
=======
A NamedTuple `(s, h_bl, integrated_w_squared)` containing:
- `clevels`: A vector of temperature levels [°C] for each contour.
- `points_x`: A vector of x-coordinates [m] for each contour.
- `points_y`: A vector of y-coordinates [m] for each contour.
"""
function analyze_frontal_properties(T)
    # ===== 1. Extract grid coordinates and surface buoyancy =====
    x, y, z = nodes(T)
    Nz = length(z)
    T_surface = interior(T, :, :, Nz)
    # ===== 2. Find the frontal outcrop using Contour.jl =====
    clevels = Float64[]
    points_x = Vector{Float64}[]
    points_y = Vector{Float64}[]

    α, g = parameters.α, parameters.g
    fig = Figure(size=(640, 640))
    ax = Axis(fig[1, 1], xlabel=L"x~\text{(km)}", ylabel=L"y~\text{(km)}", title="Surface temperature field")
    hm = heatmap!(1e-3x, 1e-3y, T_surface; rasterize = true,  colormap=:thermal)
    # CairoMakie.lines!(ax, 1e-3points_x, 1e-3points_y, color=:red)
    Colorbar(fig[1, 2], hm)
    for cl in levels(contours(x,y,T_surface))
        lvl = level(cl) # the z-value of this contour level
        for line in Contour.lines(cl)
            xs, ys = coordinates(line) # coordinates of this line segment
            xsmin, xsmax = extrema(xs)
            ysmin, ysmax = extrema(ys)
            if (xsmax - xsmin < 0.5 * (x[end]-x[1])) && (ysmax - ysmin > 0.9 * (y[end]-y[1])) # ignore lines that are too short
                CairoMakie.lines!(ax, 1e-3xs, 1e-3ys; label=@sprintf("%.2f°C", lvl)) 
                push!(clevels, lvl)
                push!(points_x, xs)
                push!(points_y, ys)
            end
        end
    end
    axislegend(ax)
    resize_to_layout!(fig)
    save(filesave * "buoyancy_field_" * fileparam * "_$(iteration).pdf", fig; pt_per_unit = 1)
    
    # ===== 3. Calculate arclength (s) along the front =====
    # dx = diff(points_x)
    # dy = diff(points_y)
    # segment_lengths = sqrt.(dx.^2 .+ dy.^2)
    # s = [0.0; cumsum(segment_lengths)] # Arclength vector [m]

    return (clevels, points_x, points_y)
end


# Example Usage 
filehead = "/orcd/data/abodner/002/nhyles_output/" 
fileparam = "subdomain2"
filesave = "results/"
iteration = 52543
# 1. Define the filename of the saved snapshot
output_filename = filehead * "subdomains/" * fileparam * "_snapshot_iter$(iteration).jld2"

# 2. Load the snapshot using the new function
snapshot = load_subdomain_snapshot(output_filename);

# Call the function to get the analysis results
frontal_data = analyze_frontal_properties(snapshot[:T]);

# Plot the results 
# fig = Figure(size=(640, 320))
# ax1 = Axis(fig[1, 1], xlabel=L"s~\text{(km)}", ylabel=L"h_{bl}~\text{(m)}", title="Boundary Layer Depth")
# ax2 = Axis(fig[2, 1], xlabel=L"s~\text{(km)}", ylabel=L"w^2~\text{(m^2 s^{-2})}", title="Depth-Integrated w²")
# lines!(ax1, 1e-3frontal_data.s, frontal_data.h_bl, color=:blue)
# lines!(ax2, 1e-3frontal_data.s, frontal_data.integrated_w_squared, color=:red)
# resize_to_layout!(fig)
# save(filesave * "frontal_properties_" * fileparam * "_" * iteration *".pdf", fig; pt_per_unit = 1)