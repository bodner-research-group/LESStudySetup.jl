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
set_value!(; Δh = 4.8828125)
using FFTW, CUDA

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

function compute_contour_properties(points_x, points_y, center_idx)
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
    outer_indices = [i for i in 1:3 if i != center_idx]
    outer1_x, outer1_y = points_x[outer_indices[1]], points_y[outer_indices[1]]
    outer2_x, outer2_y = points_x[outer_indices[2]], points_y[outer_indices[2]]
    
    n_center = length(center_x)
    
    # Pre-allocate arrays
    widths = Vector{Float32}(undef, n_center)
    
    # Compute widths using broadcasting for better performance
    for i in 1:n_center
        # Vectorized distance calculations
        dists1_sq = (center_x[i] .- outer1_x).^2 .+ (center_y[i] .- outer1_y).^2
        dists2_sq = (center_x[i] .- outer2_x).^2 .+ (center_y[i] .- outer2_y).^2
        
        # Find minimum distances (avoid sqrt until necessary)
        min_dist1 = sqrt(minimum(dists1_sq))
        min_dist2 = sqrt(minimum(dists2_sq))
        
        widths[i] = min_dist1 + min_dist2
    end
    
    # Compute arc length of center contour
    if n_center < 2
        length_total = 0.0
        curvature = Float64[]
    else
        # Calculate segment lengths
        dx = diff(center_x)
        dy = diff(center_y)
        segment_lengths = sqrt.(dx.^2 .+ dy.^2)
        length_total = sum(segment_lengths)
        
        # Compute curvature using finite differences
        curvature = compute_curvature(center_x, center_y)
    end
    
    return widths, length_total, curvature
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
    coarse_graining!(snapshot[:T], T̅; kernel=:gaussian, cutoff=100, border = :ycircular, method = :spectral, use_gpu,plans)
    # ∇b = compute!(Field(α * g * (∂x(T̅)^2+∂y(T̅)^2)^0.5))
    # ∇b_surface = interior(∇b, :, :, Nz)
    # xb, yb, _ = nodes(∇b)
    T_surface = interior(T̅, :, :, Nz)

    # ===== 2. Find the frontal outcrop using Contour.jl =====
    clevels = Float32[]
    points_x = Vector{Float32}[]
    points_y = Vector{Float32}[]

    # fig = Figure(size=(640, 640))
    # axT = Axis(fig[1, 1], xlabel=L"x~\text{(km)}", ylabel=L"y~\text{(km)}", aspect = 0.25, title="Surface temperature (°C)")
    # axb = Axis(fig[1, 3], xlabel=L"x~\text{(km)}", aspect = 0.25, title="Mixed layer depth (m)")
    # axc = Axis(fig[1, 5], xlabel=L"x~\text{(km)}", aspect = 0.25, title="Depth-averaged w²/2 (mm² s⁻²)")
    # hideydecorations!(axb, ticks=false)
    # hideydecorations!(axc, ticks=false)
    # hmT = heatmap!(axT, 1e-3x, 1e-3y, T_surface; rasterize = true,  colormap=:thermal)
    # hmb = heatmap!(axb, 1e-3x, 1e-3y, interior(snapshot[:MLD3], :, :, 1); rasterize = true, colormap=:deep)
    # hmc = heatmap!(axc, 1e-3x, 1e-3y, 1e6*interior(snapshot[:Ew3], :, :, 1); rasterize = true, colorrange=(0,100),colormap=:amp)
    # Colorbar(fig[1, 2], hmT)
    # Colorbar(fig[1, 4], hmb)
    # Colorbar(fig[1, 6], hmc)
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

    return (levels_T, points_x, points_y, MLD_front, Ew_front)
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

iterations = get_iterations_regex(filehead, fileparam)
niter = length(iterations)
width3 = zeros(niter,3)
MLD3_front = zeros(niter,3)
Ew3_front = zeros(niter,3)
rcu3 = zeros(niter,3)
Cu4 = zeros(niter,4)
Ewmax = zeros(niter)
α, g, f = parameters.α, parameters.g, parameters.f
@info "Iterations: $iterations"
for (i,iteration) in enumerate(iterations)
    @info "Iteration: $iteration"
    # 1. Define the filename of the saved snapshot
    output_filename = filehead * "subdomains/" * fileparam * "_snapshot_iter$(iteration).jld2"

    # 2. Load the snapshot using the new function
    snapshot = load_subdomain_snapshot(output_filename;variables = ("T","MLD3","Ew3"),level=224);

    # Call the function to get the analysis results
    levels_T, points_x, points_y, MLD_front, Ew_front = analyze_frontal_properties(snapshot,iteration;use_gpu,plans=(p,ip));

    center_idx = argmin(abs.(levels_T.- mean(levels_T)))
    widths, length_total, curvature = compute_contour_properties(points_x, points_y, center_idx)
    d̅ = mean(widths)/2e3 # normalized width
    d_10percent = quantile(widths, 0.1)/2e3 # normalized width
    d_90percent = quantile(widths, 0.9)/2e3 # normalized width
    h̅ = mean(MLD_front) # mean mixed layer depth
    h_10percent = quantile(MLD_front, 0.1) # 10th percentile mixed layer depth
    h_90percent = quantile(MLD_front, 0.9) # 90th percentile mixed layer depth
    E̅ = mean(Ew_front) # mean Ew
    E_10percent = quantile(Ew_front, 0.1) # 10th percentile Ew
    E_90percent = quantile(Ew_front, 0.9) # 90th percentile Ew
    Ewmax[i] = maximum(interior(snapshot[:Ew3]))
    @info "T levels: $levels_T"
    @info "Total arc length: $(length_total/1e3) km"
    @info "Normalized width: $d̅ (10th percentile: $d_10percent, 90th percentile: $d_90percent)"
    @info "Mixed layer depth: $h̅ (10th percentile: $h_10percent, 90th percentile: $h_90percent)"
    @info "Depth-averaged w²/2: $E̅ (10th percentile: $E_10percent, 90th percentile: $E_90percent)"
    itp = interpolate((nodes(snapshot[:T])[1],nodes(snapshot[:T])[2]), interior(snapshot[:MLD3], :, :, 1), Gridded(Constant()))
    MLD_cfront = itp.(points_x[center_idx], points_y[center_idx])
    Cu = 2α * g * (max(levels_T...) - min(levels_T...)) / f^2 * MLD_cfront ./ widths .* curvature
    C̅u = mean(Cu)
    C_10percent = quantile(Cu, 0.1)
    C_90percent = quantile(Cu, 0.9)
    @info "C̅u: $C̅u (10th percentile: $C_10percent, 90th percentile: $C_90percent, max: $(maximum(Cu)))"
    r̅cu = mean(1 ./ curvature)/2e3
    r_10percent = quantile(1 ./ curvature, 0.1)/2e3
    r_90percent = quantile(1 ./ curvature, 0.9)/2e3
    width3[i,:] = [d̅, d_10percent, d_90percent]
    rcu3[i,:] = [r̅cu, r_10percent, r_90percent]
    MLD3_front[i,:] = [h̅, h_10percent, h_90percent]
    Cu4[i,:] = [C̅u, C_10percent, C_90percent, maximum(Cu)]
    Ew3_front[i,:] = [E̅, E_10percent, E_90percent]
end

Q, h₀, ρ₀, cₚ = 40, 60, parameters.ρ₀, parameters.cp 
wₛ = (α * g * Q * h₀ / (ρ₀ * cₚ))^(1/3)
using Makie
wcolors = Makie.wong_colors()
fig = Figure(size=(640, 640))
g = fig[1, 1] = GridLayout()
axw = Axis(g[1, 1], titlealign = :left, title=L"\text{(a)}", ylabel=L"d/d_0")
axh = Axis(g[2, 1], titlealign = :left, title=L"\text{(b)}", ylabel=L"h/h_0")
axc = Axis(g[3, 1], titlealign = :left, title=L"\text{(c)}", ylabel=L"Cu_b")
#axb = Axis(g[4, 1], titlealign = :left, title=L"\text{(d)}", ylabel=L"Cu_b")
axE = Axis(g[4, 1], titlealign = :left, title=L"\text{(d)}", ylabel=L"{\overline{(w^2/2)}^h_{max}/w_*^2", xlabel=L"t~\text{(h)}")
hidexdecorations!(axw, ticks=false)
hidexdecorations!(axc, ticks=false)
hidexdecorations!(axh, ticks=false)
t = [11; 16; 21; 24:51]
lines!(axw, t, width3[:,1], color=wcolors[1])
fill_between!(axw, t, width3[:,2], width3[:,3]; color = wcolors[1], alpha = 0.25)
lines!(axh, t, MLD3_front[:,1]/h₀, color=wcolors[2])
fill_between!(axh, t, MLD3_front[:,2]/h₀, MLD3_front[:,3]/h₀; color = wcolors[2], alpha = 0.25)
# lines!(axc, t, Cu3[:,1], color=wcolors[3])
# fill_between!(axc, t, Cu3[:,2], Cu3[:,3]; color = wcolors[3], alpha = 0.25)
lines!(axc, t, Cu4[:,4], color=wcolors[3])
# lines!(axb, t, Cu4[:,1], color=wcolors[4])
lines!(axE, t, Ewmax/wₛ^2, color=wcolors[4])
#fill_between!(axE, t, Ew3_front[:,2]/wₛ^2, Ew3_front[:,3]/wₛ^2; color = wcolors[5], alpha = 0.25)
rowgap!(g, 3)
resize_to_layout!(fig)
save(filesave * "ts_d_h_Cu_Ewmax_" * fileparam * ".pdf", fig; pt_per_unit = 1)