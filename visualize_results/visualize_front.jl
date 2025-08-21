using LESStudySetup
using CairoMakie
using Printf, Dates
using Statistics: mean, std
using Oceananigans.Architectures: on_architecture
using LESStudySetup.Diagnostics
using LESStudySetup.Diagnostics: load_distributed_checkpoint,load_subdomain_snapshot
using LESStudySetup.Diagnostics: isotropic_powerspectrum
using LESStudySetup.Diagnostics: coarse_graining!, TKE
set_theme!(theme_latexfonts(), fontsize=12, figure_padding = 10)
using Contour, Interpolations, Statistics
set_value!(; Δh = 4.8828125)
using FFTW, CUDA, JLD2
α, g, f = parameters.α, parameters.g, parameters.f

function compute_curvature(x, y)
    """
    Compute curvature at each point using finite differences.
    Uses the formula: κ = |x'y'' - y'x''| / (x'^2 + y'^2)^(3/2)
    """
    n = length(x)
    
    if n < 3
        return zeros(n)
    end
    
    curvature = zeros(n)
    
    for i in 2:(n-1)
        # First derivatives (central difference)
        dx = (x[i+1] - x[i-1]) / 2
        dy = (y[i+1] - y[i-1]) / 2
        
        # Second derivatives
        d2x = x[i+1] - 2*x[i] + x[i-1]
        d2y = y[i+1] - 2*y[i] + y[i-1]
        
        # Curvature formula
        numerator = abs(dx * d2y - dy * d2x)
        denominator = (dx^2 + dy^2)^1.5
        
        curvature[i] = denominator > 1e-12 ? numerator / denominator : 0.0
    end
    
    # Handle endpoints assuming periodicity
    if n >= 3
        # First point (use periodic boundary)
        dx = (x[2] - x[n-1]) / 2  # Wrap around to last point
        dy = (y[2] - y[n-1]) / 2
        d2x = x[2] - 2*x[1] + x[n-1]  # Periodic second derivative
        d2y = y[2] - 2*y[1] + y[n-1]
        denominator = (dx^2 + dy^2)^1.5
        curvature[1] = denominator > 1e-12 ? abs(dx * d2y - dy * d2x) / denominator : 0.0
        
        # Last point (use periodic boundary)
        dx = (x[2] - x[n-1]) / 2  # Wrap around to first point
        dy = (y[2] - y[n-1]) / 2
        d2x = x[2] - 2*x[n] + x[n-1]  # Periodic second derivative
        d2y = y[2] - 2*y[n] + y[n-1]
        denominator = (dx^2 + dy^2)^1.5
        curvature[n] = denominator > 1e-12 ? abs(dx * d2y - dy * d2x) / denominator : 0.0
    end
    
    return curvature
end

function compute_contour_properties(points_x, points_y, center_idx; T=Float32)
    """
    Compute widths, length, and curvature for contour analysis.
    
    Args:
        points_x: Array of 3 arrays containing x-coordinates for each contour
        points_y: Array of 3 arrays containing y-coordinates for each contour  
        center_idx: Index (1-3) indicating which contour is the center
        
    Returns:
        (widths, length, curvature): Tuple containing:
            - widths: Array of local widths at each center point
            - length: Total arc length of center contour
            - curvature: Array of curvature values at each center point
    """
    
    # Extract contours
    center_x, center_y = points_x[center_idx], points_y[center_idx]
    n_center = length(center_x)
    if length(points_x) > 2
        outer_indices = [i for i in 1:3 if i != center_idx]
        outer1_x, outer1_y = points_x[outer_indices[1]], points_y[outer_indices[1]]
        outer2_x, outer2_y = points_x[outer_indices[2]], points_y[outer_indices[2]]
        @info "contour sizes: $(length(center_x)), $(length(outer1_x)), $(length(outer2_x))"

        # Pre-allocate arrays
        widths = Vector{T}(undef, n_center)

        # Compute widths using broadcasting for better performance
        for i in 1:n_center
            # Vectorized distance calculations
            dists1_sq = (center_x[i] .- outer1_x).^2 .+ (center_y[i] .- outer1_y).^2
            dists2_sq = (center_x[i] .- outer2_x).^2 .+ (center_y[i] .- outer2_y).^2
            
            # Find minimum distances (avoid sqrt until necessary)
            min_dist1 = sqrt(minimum(dists1_sq))
            min_dist2 = sqrt(minimum(dists2_sq))
            
            widths[i] = T(min_dist1 + min_dist2)
        end
    else
        @info "No enough contours to compute widths"
        widths = nothing
    end
    
    # Compute arc length of center contour
    if n_center < 2
        arclength = 0.0
        curvature = 0.0
    else
        # Calculate segment lengths
        dx = diff(center_x)
        dy = diff(center_y)
        segment_lengths = sqrt.(dx.^2 .+ dy.^2)
        arclength = cumsum(segment_lengths)
        
        # Compute curvature using finite differences
        curvature = compute_curvature(center_x, center_y)
    end
    
    return widths, arclength, curvature
end

"""
    analyze_frontal_properties(T)

Analyzes properties along a frontal outcrop identified from a temperature field.

Arguments
=========
- `T`: The 3D temperature field (`Field` object) from Oceananigans.

Returns
=======
A NamedTuple `(s, h_bl, integrated_w_squared)` containing:
- `levels_T`: A vector of temperature levels [°C] for each contour.
- `points_x`: A vector of x-coordinates [m] for each contour.
- `points_y`: A vector of y-coordinates [m] for each contour.
"""
function analyze_frontal_properties(snapshot,iteration;use_gpu=false,plans=nothing)
    # ===== 1. Extract grid coordinates and surface buoyancy =====
    x, y, z = nodes(snapshot[:T])
    Nz = length(z)
    # α, g = parameters.α, parameters.g
    T̅ = CenterField(snapshot[:grid], Float32);
    coarse_graining!(snapshot[:T], T̅; kernel=:gaussian, cutoff=300, border = :ycircular, method = :spectral, use_gpu,plans)
    # ∇b = compute!(Field(α * g * (∂x(T̅)^2+∂y(T̅)^2)^0.5))
    # ∇b_surface = interior(∇b, :, :, Nz)
    # xb, yb, _ = nodes(∇b)
    T_surface = interior(T̅, :, :, Nz)

    # ===== 2. Find the frontal outcrop using Contour.jl =====
    clevels = Float32[]
    points_x = Vector{Float32}[]
    points_y = Vector{Float32}[]

    for cl in levels(contours(x,y,T_surface))
        lvl = level(cl) # the z-value of this contour level
        for line in Contour.lines(cl)
            xs, ys = coordinates(line) # coordinates of this line segment
            xsmin, xsmax = extrema(xs)
            ysmin, ysmax = extrema(ys)
            if (xsmax - xsmin < 0.5 * (x[end]-x[1])) && (ysmax - ysmin > 0.9 * (y[end]-y[1])) # ignore lines that are too short
                push!(clevels, lvl)
            end
        end
    end

    levels_T = Float64[]
    for cl in levels(contours(x,y,T_surface,mean(clevels) .+ [-0.194,0,0.194])) # gives initial 2km width
        lvl = level(cl) # the |∇b|-value of this contour level
        for line in Contour.lines(cl)
            xs, ys = coordinates(line) # coordinates of this line segment
            xsmin, xsmax = extrema(xs)
            ysmin, ysmax = extrema(ys)
            if (xsmax - xsmin < 0.5 * (x[end]-x[1])) && (ysmax - ysmin > 0.9 * (y[end]-y[1])) # ignore lines that are too short
                # CairoMakie.lines!(axT, 1e-3xs, 1e-3ys; label=@sprintf("%.2f°C", lvl), linewidth=0.8) 
                push!(levels_T, lvl)
                push!(points_x, xs)
                push!(points_y, ys)
            end
        end
    end

    idx_front = min(levels_T...) .< T_surface .< max(levels_T...)
    MLD_front = vec(idx_front .* interior(snapshot[:MLD3], :, :, 1))
    Ew_front = vec(idx_front .* interior(snapshot[:Ew3], :, :, 1))
    Ew_front = Ew_front[MLD_front.>1e-6]
    MLD_front = MLD_front[MLD_front.>1e-6]

    # axislegend(axT, labelsize=9, padding = (0f0, 0f0, 0f0, 0f0), patchlabelgap = 3)
    # colgap!(fig.layout, 1, Relative(0)) 
    # colgap!(fig.layout, 3, Relative(0)) 
    # colgap!(fig.layout, 5, Relative(0)) 

    # resize_to_layout!(fig)
    # save(filesave * "T0_MLD3_Ew3_" * fileparam * "_$(iteration).pdf", fig; pt_per_unit = 1)

    return (mean(clevels), levels_T, points_x, points_y, MLD_front, Ew_front)
end

"""
    normal_strain_rate(xs, ys, itpU, itpV)

Calculates the strain rate of a 2D velocity field normal to a contour.

This function is essential for analyzing frontogenesis, where the mesoscale strain 
field acts to sharpen or weaken a front. It is vectorized for efficiency.

# Arguments
- `xs::Vector{<:Real}`: A vector of x-coordinates defining the frontal contour.
- `ys::Vector{<:Real}`: A vector of y-coordinates defining the frontal contour.
- `itpU::AbstractInterpolation`: An interpolation object for the U-component of the 
                                 velocity field. Must support `Interpolations.gradient`.
- `itpV::AbstractInterpolation`: An interpolation object for the V-component of the 
                                 velocity field. Must support `Interpolations.gradient`.

# Returns
- `Vector{<:Real}`: A vector containing the normal strain rate `Sn` at each point 
                    along the input contour. Positive values typically correspond to 
                    frontogenesis.

# Example
# Assuming you have U, V on a grid (x, y) and have found a contour (xs, ys)
# x = 1:100; y = 1:100;
# U_grid = rand(100, 100); V_grid = rand(100, 100);
# itpU = linear_interpolation((x, y), U_grid);
# itpV = linear_interpolation((x, y), V_grid);
# xs, ys = find_contour(...); # Placeholder for a contouring algorithm
# Sn = normal_strain_rate(xs, ys, itpU, itpV);
"""
function normal_strain_rate(xs::Vector{<:Real}, ys::Vector{<:Real}, itpU, itpV)
    
    # 1. Calculate tangent vectors to the front using a central difference scheme.
    # `circshift` elegantly handles the periodic boundary conditions.
    tangent_x = circshift(xs, -1) .- circshift(xs, 1)
    tangent_y = circshift(ys, -1) .- circshift(ys, 1)
    tangent_y[1] += 1e5 
    tangent_y[end] += 1e5 

    # 2. Compute the corresponding normal vectors (un-normalized).
    # A 90-degree rotation of the tangent (tx, ty) gives the normal (-ty, tx).
    normal_x_un = -tangent_y
    normal_y_un =  tangent_x

    # 3. Normalize the normal vectors to get unit vectors n = (nx, ny).
    # We add a small epsilon to prevent division by zero for any degenerate points.
    magnitudes = hypot.(normal_x_un, normal_y_un) .+ 1e-12
    nx = normal_x_un ./ magnitudes
    ny = normal_y_un ./ magnitudes

    # 4. Interpolate the velocity gradients onto the front coordinates.
    # The `Ref()` ensures that the function `gradient` is broadcast over the
    # vectors of coordinates `xs` and `ys`.
    ∇U = Interpolations.gradient.(Ref(itpU), xs, ys)
    ∇V = Interpolations.gradient.(Ref(itpV), xs, ys)

    # Unpack the gradients into their components for clarity.
    # ∇U is a vector of tuples, e.g., [(dU/dx₁, dU/dy₁), (dU/dx₂, dU/dy₂), ...]
    dudx = [g[1] for g in ∇U]
    dudy = [g[2] for g in ∇U]
    dvdx = [g[1] for g in ∇V]
    dvdy = [g[2] for g in ∇V]

    # 5. Compute the normal strain rate using the vectorized formula.
    # Sn = nx² * (∂U/∂x) + ny² * (∂V/∂y) + nx*ny * (∂U/∂y + ∂V/∂x)
    Sn = @. nx^2 * dudx + ny^2 * dvdy + nx * ny * (dudy + dvdx)

    return -Sn
end

function get_iterations_regex(filehead, fileparam, directory="."; subdirparam="subdomains", rank = 1008)
    subdir = joinpath(directory, filehead * subdirparam)
    
    if !isdir(subdir)
        return Int[]
    end
    
    files = readdir(subdir)
    if subdirparam == "subdomains"
        pattern = Regex("^" * escape_string(fileparam) * "_snapshot_iter(\\d+)\\.jld2\$")
    elseif subdirparam[1:4] == "iter"
        pattern = Regex("^" * escape_string(fileparam) * "_$(rank)_iteration(\\d+)\\.jld2\$")
    end
    
    iterations = Int[]
    
    for file in files
        m = match(pattern, file)
        if m !== nothing
            push!(iterations, parse(Int, m.captures[1]))
        end
    end
    
    return sort(iterations)
end

#(a) Mesoscale strain rate normal to the front, S_n.
#(b) Frontal width, d.
#(c) Frontal curvature number, Cu.
#(b) Depth-integrated TKE, ∫TKEdz.

@info "Memory usage: " * (Sys.free_memory() |> Base.format_bytes) * " GB free"
filehead = "/orcd/data/abodner/002/shared_datasets/nhyles_output/" 
fileparam = "xband1sublevels"
filesave = "results/"

use_gpu = true
x=rand(Float32,10240,20480,1);
if use_gpu && CUDA.functional()
    p = CUDA.CUFFT.plan_rfft(CuArray(x), (1,2));
    ip = CUDA.CUFFT.plan_irfft(p * CuArray(x), 10240, (1,2));
else
    p = FFTW.plan_rfft(x, (1,2); flags=FFTW.MEASURE, timelimit=Inf);
    ip = FFTW.plan_irfft(p * x, 10240, (1,2); flags=FFTW.MEASURE, timelimit=Inf);
end

# t = [11; 16; 21; 24:80];
# iterations = get_iterations_regex(filehead, fileparam)
# #iterations = iterations[t .>= 57]
# niter = length(iterations)
# width3 = zeros(niter,3)
# MLD3_front = zeros(niter,3)
# Ew3_front = zeros(niter,3)
# rcu3 = zeros(niter,3)
# Cu4 = zeros(niter,4)
# Ewmax = zeros(niter)
# @info "Iterations: $iterations"
# for (i,iteration) in enumerate(iterations)
#     @info "Iteration: $iteration"
#     # 1. Define the filename of the saved snapshot
#     output_filename = filehead * "subdomains/" * fileparam * "_snapshot_iter$(iteration).jld2"

#     # 2. Load the snapshot using the new function
#     snapshot = load_subdomain_snapshot(output_filename;variables = ("T","MLD3","Ew3"),level=224);

#     # Call the function to get the analysis results
#     clevel, levels_T, points_x, points_y, MLD_front, Ew_front = analyze_frontal_properties(snapshot,iteration;use_gpu,plans=(p,ip));

#     if length(levels_T) < 1
#         break
#     end
#     center_idx = argmin(abs.(levels_T.- clevel))
#     widths, arclength, curvature = compute_contour_properties(points_x, points_y, center_idx)
#     if !isnothing(widths)
#         d̅ = mean(widths)/2e3 # normalized width
#         d_10percent = quantile(widths, 0.1)/2e3 # normalized width
#         d_90percent = quantile(widths, 0.9)/2e3 # normalized width
#     else
#         break
#     end
#     h̅ = mean(MLD_front) # mean mixed layer depth
#     h_10percent = quantile(MLD_front, 0.1) # 10th percentile mixed layer depth
#     h_90percent = quantile(MLD_front, 0.9) # 90th percentile mixed layer depth
#     E̅ = mean(Ew_front) # mean Ew
#     E_10percent = quantile(Ew_front, 0.1) # 10th percentile Ew
#     E_90percent = quantile(Ew_front, 0.9) # 90th percentile Ew
#     Ewmax[i] = maximum(interior(snapshot[:Ew3]))
#     @info "T levels: $levels_T"
#     @info "Total arc length: $(arclength[end]/1e3) km"
#     @info "Normalized width: $d̅ (10th percentile: $d_10percent, 90th percentile: $d_90percent)"
#     @info "Mixed layer depth: $h̅ (10th percentile: $h_10percent, 90th percentile: $h_90percent)"
#     @info "Depth-averaged w²/2: $E̅ (10th percentile: $E_10percent, 90th percentile: $E_90percent)"
#     itp = interpolate((nodes(snapshot[:T])[1],nodes(snapshot[:T])[2]), interior(snapshot[:MLD3], :, :, 1), Gridded(Constant()))
#     MLD_cfront = itp.(points_x[center_idx], points_y[center_idx])
#     Cu = 2α * g * (max(levels_T...) - min(levels_T...)) / f^2 * MLD_cfront ./ widths .* curvature
#     C̅u = mean(Cu)
#     C_10percent = quantile(Cu, 0.1)
#     C_90percent = quantile(Cu, 0.9)
#     @info "C̅u: $C̅u (10th percentile: $C_10percent, 90th percentile: $C_90percent, max: $(maximum(Cu)))"
#     r̅cu = mean(1 ./ curvature)/2e3
#     r_10percent = quantile(1 ./ curvature, 0.1)/2e3
#     r_90percent = quantile(1 ./ curvature, 0.9)/2e3
#     width3[i,:] = [d̅, d_10percent, d_90percent]
#     rcu3[i,:] = [r̅cu, r_10percent, r_90percent]
#     MLD3_front[i,:] = [h̅, h_10percent, h_90percent]
#     Cu4[i,:] = [C̅u, C_10percent, C_90percent, maximum(Cu)]
#     Ew3_front[i,:] = [E̅, E_10percent, E_90percent]
#     @info "Memory usage: " * (Sys.free_memory() |> Base.format_bytes) * " GB free"
# end

# jldopen(filesave * "ts_d_h_Cub_Ewmax_" * fileparam * "_cg3hm.jld2", "w") do file
#     file["widths"] = width3
#     file["MLDs_front"] = MLD3_front
#     file["Cus"] = Cu4
#     file["Ewmax"] = Ewmax 
#     file["t"] = t
# end

# file = jldopen(filesave * "ts_d_h_Cub_Ewmax_" * fileparam * "_cg3hm.jld2", "r")
# width3 = file["widths"]
# MLD3_front = file["MLDs_front"]
# Cu4 = file["Cus"]
# Ewmax = file["Ewmax"]
# close(file)
# idx = (.!(iszero.(Ewmax)))
# t = t[idx]

Q, h₀, ρ₀, cₚ = 40, 60, parameters.ρ₀, parameters.cp 
wₛ = (α * g * Q * h₀ / (ρ₀ * cₚ))^(1/3)
using Makie
wcolors = Makie.wong_colors()
# fig = Figure(size=(640, 640))
# gl = fig[1, 1] = GridLayout()
# axw = Axis(gl[1, 1], titlealign = :left, title=L"\text{(a)}", ylabel=L"d/d_0")
# axh = Axis(gl[2, 1], titlealign = :left, title=L"\text{(b)}", ylabel=L"h/h_0")
# axc = Axis(gl[3, 1], titlealign = :left, title=L"\text{(c)}", ylabel=L"Cu_b")
# #axb = Axis(gl[4, 1], titlealign = :left, title=L"\text{(d)}", ylabel=L"Cu_b")
# axE = Axis(gl[4, 1], titlealign = :left, title=L"\text{(d)}", ylabel=L"\overline{(w^2/2)}^h_{\max}/w_*^2", xlabel=L"t~\text{(h)}")
# hidexdecorations!(axw, ticks=false)
# hidexdecorations!(axc, ticks=false)
# hidexdecorations!(axh, ticks=false)
# lines!(axw, t, width3[idx,1], color=wcolors[1])
# fill_between!(axw, t, width3[idx,2], width3[idx,3]; color = wcolors[1], alpha = 0.25)
# lines!(axh, t, MLD3_front[idx,1]/h₀, color=wcolors[2])
# fill_between!(axh, t, MLD3_front[idx,2]/h₀, MLD3_front[idx,3]/h₀; color = wcolors[2], alpha = 0.25)
# # lines!(axc, t, Cu3[:,1], color=wcolors[3])
# # fill_between!(axc, t, Cu3[:,2], Cu3[:,3]; color = wcolors[3], alpha = 0.25)
# lines!(axc, t, Cu4[idx,4], color=wcolors[3])
# # lines!(axb, t, Cu4[:,1], color=wcolors[4])
# lines!(axE, t, Ewmax[idx]/wₛ^2, color=wcolors[4])
# #fill_between!(axE, t, Ew3_front[:,2]/wₛ^2, Ew3_front[:,3]/wₛ^2; color = wcolors[5], alpha = 0.25)
# kwargs = (; markersize = 10, strokewidth = 1, color = :white)
# scatter!(axc, 44, Cu4[argmin(abs.(t.- 44)),4]; kwargs..., strokecolor = wcolors[3])
# scatter!(axE, 30, Ewmax[argmin(abs.(t.- 30))]/wₛ^2; kwargs..., strokecolor = wcolors[4])
# scatter!(axE, 49, Ewmax[argmin(abs.(t.- 49))]/wₛ^2; kwargs..., strokecolor = wcolors[4])
# rowgap!(gl, 3)
# resize_to_layout!(fig)
# save(filesave * "ts_d_h_Cub_Ewmax_" * fileparam * "_cg3hm.pdf", fig; pt_per_unit = 1)

iteration = 49086
# 1. Define the filename of the saved snapshot
output_filename = filehead * "subdomains/" * fileparam * "_snapshot_iter$(iteration).jld2"

# 2. Load the snapshot using the new function
snapshot = load_subdomain_snapshot(output_filename;variables = ("u","v","T","MLD3","Ew3"),level=224);

grid = snapshot[:u].grid
u̅  = XFaceField(grid);
v̅  = YFaceField(grid);
coarse_graining!(snapshot[:v] , v̅; kernel=:gaussian, cutoff=300, border = :ycircular, method = :spectral, use_gpu,plans=(p,ip));
coarse_graining!(snapshot[:u] , u̅; kernel=:gaussian, cutoff=300, border = :ycircular, method = :spectral, use_gpu,plans=(p,ip));
c2g(field) = on_architecture(GPU(),field)
g2c(field) = on_architecture(CPU(),field)
ζ̅ = g2c(compute!(Field((∂x(c2g(v̅)) - ∂y(c2g(u̅))))));
δ̅ = g2c(compute!(Field((∂x(c2g(u̅)) + ∂y(c2g(v̅))))));

# Call the function to get the analysis results
clevel, levels_T, points_x, points_y, MLD_front, Ew_front = analyze_frontal_properties(snapshot,iteration;use_gpu,plans=(p,ip));

center_idx = argmin(abs.(levels_T.- clevel))
widths, arclength, curvature = compute_contour_properties(points_x, points_y, center_idx)
@info "Memory usage: " * (Sys.free_memory() |> Base.format_bytes) * " GB free"

A=rand(Float32,10240,20480,8);
if use_gpu && CUDA.functional()
    p = CUDA.CUFFT.plan_rfft(CuArray(A), (1,2));
    pk = CUDA.CUFFT.plan_rfft(CuArray(A[:,:,1:1]), (1,2));
    ip = CUDA.CUFFT.plan_irfft(p * CuArray(A), 10240, (1,2));
else
    @error "CUDA not available"
end
fileparam = "xband1surf9"
output_filename = filehead * "subdomains/" * fileparam * "_snapshot_iter$(iteration).jld2"
snapshot2 = load_subdomain_snapshot(output_filename; variables = ("u", "v", "w"));
u̅, v̅, w̅, τuu, τvv, τww = TKE(snapshot2; cutoff=300, border=:ycircular, Lx = snapshot2[:grid].Lx, Ly = snapshot2[:grid].Ly,
                                       method=:spectral,use_gpu,plans=(p,pk,ip));
t0 = time()
TKEₛ = g2c(compute!(Field((τvv + τuu + τww)/2/wₛ^2)));
@info "TKE compute time: $(time() - t0)"
@info "Memory usage: " * (Sys.free_memory() |> Base.format_bytes) * " GB free"
SKEₛ = ((u̅ - mean(u̅;dims=(1,2)))^2 + (v̅ - mean(v̅;dims=(1,2)))^2 + (w̅ - mean(w̅;dims=(1,2)))^2)/2/wₛ^2;
SKEₛ = g2c(compute!(Field(SKEₛ)));
@info "SKE compute time: $(time() - t0)"
@info "Memory usage: " * (Sys.free_memory() |> Base.format_bytes) * " GB free"

x, y, z = nodes(snapshot[:T])
fig = Figure(size = (640, 450))
gab = fig[1, 1] = GridLayout()
aspect = 0.25
axis_kwargs = (titlealign = :left,xlabel = L"x~\text{(km)}", limits = ((-12.5,12.5),(0,100)), aspect=aspect)
axT = Axis(fig[1, 1]; ylabel=L"y~\text{(km)}", title=L"\text{(a)}~T_{\text{surf}}~\text{({^\circ}C)}", axis_kwargs...)
axb = Axis(fig[1, 3]; aspect = 0.25, title=L"\text{(b)}~h/h_0", axis_kwargs...)
axc = Axis(fig[1, 5]; aspect = 0.25, title=L"\text{(c) TKE}/w_*^2", axis_kwargs...)
axd = Axis(fig[1, 7]; aspect = 0.25, title=L"\text{(d) SKE}/w_*^2", axis_kwargs...)
hideydecorations!(axb, ticks=false)
hideydecorations!(axc, ticks=false)
hideydecorations!(axd, ticks=false)
hmT = heatmap!(axT, 1e-3x, 1e-3y, interior(snapshot[:T], :, :, 1); rasterize = true,  colormap=:thermal)
hmb = heatmap!(axb, 1e-3x, 1e-3y, interior(snapshot[:MLD3], :, :, 1)/h₀; rasterize = true, colormap=:deep)
idxTKEm = argmax(interior(TKEₛ))
idxSKEm = argmax(interior(SKEₛ))
xc, yc, zc = nodes(TKEₛ);
xd, yd, zd = nodes(SKEₛ);
@info "TKE maximum at x=$(xc[idxTKEm[1]]/1e3)km, y=$(yc[idxTKEm[2]]/1e3)km, z=$(zc[idxTKEm[3]])m"
@info "SKE maximum at x=$(xd[idxSKEm[1]]/1e3)km, y=$(yd[idxSKEm[2]]/1e3)km, z=$(zd[idxSKEm[3]])m"
hmc = heatmap!(axc, 1e-3xc, 1e-3yc, interior(TKEₛ, :, :, idxTKEm[3]); rasterize = true, colormap = :amp, colorrange = (0, max(interior(TKEₛ, :, :, idxTKEm[3])...)))
hmd = heatmap!(axd, 1e-3xd, 1e-3yd, interior(SKEₛ, :, :, idxSKEm[3]); rasterize = true, colormap = :amp, colorrange = (0, max(interior(SKEₛ, :, :, idxSKEm[3])...)))
scatter!(axc, 1e-3xc[idxTKEm[1]], 1e-3yc[idxTKEm[2]]; marker = :star4, markersize = 10, color = :black)
scatter!(axd, 1e-3xd[idxSKEm[1]], 1e-3yd[idxSKEm[2]]; marker = :star4, markersize = 10, color = :black)
Colorbar(fig[1, 2], hmT)
Colorbar(fig[1, 4], hmb)
Colorbar(fig[1, 6], hmc)
Colorbar(fig[1, 8], hmd)
for i = 1:3
    xs, ys = points_x[i], points_y[i]
    lvl = levels_T[i]
    CairoMakie.lines!(axT, 1e-3xs, 1e-3ys; label=@sprintf("%.2f°C", lvl), linewidth=0.8) 
end
@info "T levels: $(levels_T)"
#axislegend(axT, labelsize=9, padding = (0f0, 0f0, 0f0, 0f0), patchlabelgap = 3)
for i = 1:2:7
    colgap!(fig.layout, i, Relative(0)) 
end
for i = 2:2:6
    colgap!(fig.layout, i, 3) 
end

shift(x) = [x[size(x,1)÷2+1:end, :]; x[1:size(x,1)÷2, :]]
filename0 = "./hydrostatic_snapshots_init.jld2"
snapshots = load_snapshots(filename0);
v0 = snapshots[:v][1];
T0 = snapshots[:T][1];
xT, yT, zT = nodes(T0);
initfile = "./hydrostatic_snapshots_free.jld2"
initsnaps = load_snapshots(initfile)
Ub = initsnaps[:u][1];
Vb = compute!(Field(initsnaps[:v][1] - snapshots[:v][1]));
iarrow = 2:32:640
k = length(zT)
arrows!(axb,xT[iarrow]/1e3.-50, yT[iarrow]/1e3, shift(interior(Ub,iarrow,iarrow,k)), shift(interior(Vb,iarrow,iarrow,k)), arrowsize = 3, lengthscale = 1e2,linecolor = :white, arrowcolor = :white, linewidth = 0.2,alpha=0.5)
arrows!(axc,xT[iarrow]/1e3.-50, yT[iarrow]/1e3, shift(interior(Ub,iarrow,iarrow,k)), shift(interior(Vb,iarrow,iarrow,k)), arrowsize = 3, lengthscale = 1e2,linecolor = :black, arrowcolor = :black, linewidth = 0.2,alpha=0.5)
arrows!(axd,xT[iarrow]/1e3.-50, yT[iarrow]/1e3, shift(interior(Ub,iarrow,iarrow,k)), shift(interior(Vb,iarrow,iarrow,k)), arrowsize = 3, lengthscale = 1e2,linecolor = :black, arrowcolor = :black, linewidth = 0.2,alpha=0.5)
resize_to_layout!(fig)
save(filesave * "T0_MLD3_TKE_SKE_" * fileparam * "_$(iteration).pdf", fig; pt_per_unit = 1)

xU,yU,zU = nodes(Ub);
itpU = interpolate((xU,yU), interior(Ub, :, :, length(zU)), Gridded(Linear(Interpolations.Periodic())))
etpU = extrapolate(itpU, Interpolations.Periodic())
xV,yV,zV = nodes(Vb);
itpV = interpolate((xV,yV), interior(Vb, :, :, length(zV)), Gridded(Linear(Interpolations.Periodic())))
etpV = extrapolate(itpV, Interpolations.Periodic())
Sn = normal_strain_rate(points_x[center_idx], points_y[center_idx], etpU, etpV)
itph = interpolate((x,y), interior(snapshot[:MLD3], :, :, 1), Gridded(Constant()))
MLD_cfront = itph.(points_x[center_idx], points_y[center_idx]);
Cu = 2α * g * (max(levels_T...) - min(levels_T...)) / f^2 * MLD_cfront ./ widths .* curvature;
xζ, yζ, _ = nodes(ζ̅);
itpζ = interpolate((xζ,yζ), interior(ζ̅, :, :, 1), Gridded(Constant(Interpolations.Periodic())))
etpζ = extrapolate(itpζ, Interpolations.Periodic())
ζ_cfront = etpζ.(points_x[center_idx], points_y[center_idx]);
xδ, yδ, _ = nodes(δ̅);
itpδ = interpolate((xδ,yδ), interior(δ̅, :, :, 1), Gridded(Constant(Interpolations.Periodic())))
etpδ = extrapolate(itpδ, Interpolations.Periodic())
δ_cfront = etpδ.(points_x[center_idx], points_y[center_idx]);
itpTKE = interpolate((xc,yc), interior(TKEₛ, :, :, idxTKEm[3]), Gridded(Constant(Interpolations.Periodic())))
etpTKE = extrapolate(itpTKE, Interpolations.Periodic())
TKE_cfront = etpTKE.(points_x[center_idx], points_y[center_idx]);
itpSKE = interpolate((xd,yd), interior(SKEₛ, :, :, idxSKEm[3]), Gridded(Constant(Interpolations.Periodic())))
etpSKE = extrapolate(itpSKE, Interpolations.Periodic())
SKE_cfront = etpSKE.(points_x[center_idx], points_y[center_idx]);
itpEw = interpolate((x,y), interior(snapshot[:Ew3], :, :, 1), Gridded(Constant()))
Ew_cfront = itpEw.(points_x[center_idx], points_y[center_idx]);
fig = Figure(size=(640, 720))
limits = ((0, arclength[end]/1e5),nothing)
axis_kwargs = (titlealign = :left, limits)
axa = Axis(fig[1, 1]; title=L"\text{(a)}", ylabel=L"S_n/f", axis_kwargs...)
axb = Axis(fig[2, 1]; title=L"\text{(b)}", ylabel=L"d/d_0", yticklabelcolor = wcolors[1], axis_kwargs...)
axb2 = Axis(fig[2, 1]; ylabel=L"h/h_0", yticklabelcolor = wcolors[2], yaxisposition = :right, limits)
hidespines!(axb2);hidexdecorations!(axb2);
axc = Axis(fig[3, 1]; title=L"\text{(c)}", ylabel=L"Cu", axis_kwargs...)
axd = Axis(fig[4, 1]; title=L"\text{(d)}", ylabel=L"\zeta/f", yticklabelcolor = wcolors[1], axis_kwargs...)
axd2 = Axis(fig[4, 1]; ylabel=L"\delta/f", yticklabelcolor = wcolors[2], yaxisposition = :right, limits)
hidespines!(axd2);hidexdecorations!(axd2);
axe = Axis(fig[5, 1]; title=L"\text{(e)}", xlabel=L"s/L_y", ylabel=L"TKE/w_*^2", yticklabelcolor = wcolors[1], axis_kwargs...)
axe2 = Axis(fig[5, 1]; ylabel=L"SKE/w_*^2", yticklabelcolor = wcolors[2], yaxisposition = :right,limits)
hidexdecorations!(axa, ticks=false)
hidexdecorations!(axb, ticks=false)
hidexdecorations!(axc, ticks=false)
hidexdecorations!(axd, ticks=false)
arclength = [0;arclength]
lines!(axa, arclength/1e5, Sn/f; linewidth=1, color = wcolors[1])
lines!(axb, arclength/1e5, widths/2e3; linewidth=1, color = wcolors[1])
lines!(axb2, arclength/1e5, MLD_cfront/60; linewidth=1, color = wcolors[2])
lines!(axc, arclength/1e5, Cu; linewidth=1, color = wcolors[1])
lines!(axd, arclength/1e5, ζ_cfront/f; linewidth=1, color = wcolors[1])
lines!(axd2, arclength/1e5, δ_cfront/f; linewidth=1, color = wcolors[2])
lines!(axe, arclength/1e5, TKE_cfront; linewidth=1, color = wcolors[1])
#lines!(axe, arclength/1e5, Ew_cfront/wₛ^2; linewidth=1, color = wcolors[1], linestyle = :dash)
lines!(axe2, arclength/1e5, SKE_cfront; linewidth=1, color = wcolors[2])
rowgap!(fig.layout, 3)
resize_to_layout!(fig)
save(filesave * "along_front_" * fileparam * "_$(iteration).pdf", fig; pt_per_unit = 1)