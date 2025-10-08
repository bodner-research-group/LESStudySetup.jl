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
using Colors, ColorSchemes
set_value!(; Δh = 4.8828125)
using FFTW, CUDA, JLD2
α, g, f = parameters.α, parameters.g, parameters.f
M²₀ = parameters.M²₀
Q, h₀, ρ₀, cₚ = 40, 60, parameters.ρ₀, parameters.cp 
wₛ = (α * g * Q * h₀ / (ρ₀ * cₚ))^(1/3)
using Makie
wcolors = Makie.wong_colors()
filehead = "/orcd/data/abodner/002/shared_datasets/nhyles_output/" 
filesave = "results/"
c2g(field) = on_architecture(GPU(),field)
g2c(field) = on_architecture(CPU(),field)
println(CUDA.pool_status())

"""
    coarsen_binned_vectorized(f, s, n_target)

Vectorized version of binned averaging for better performance.
"""
function coarsen_binned_vectorized(f::Vector{T}, s::Vector{T}, n_target::Int) where T<:Real
    n_orig = length(f)
    s_total = s[end] - s[1] + (s[2] - s[1])
    bin_size = s_total / n_target
    s_min = s[1]
    
    # Vectorized bin assignment
    s_wrapped = mod.(s .- s_min, s_total) .+ s_min
    bin_indices = min.(floor.(Int, (s_wrapped .- s_min) ./ bin_size) .+ 1, n_target)
    
    # Accumulate using sparse-like approach
    f_coarse = zeros(T, n_target)
    bin_counts = zeros(Int, n_target)
    
    for i in 1:n_orig
        idx = bin_indices[i]
        f_coarse[idx] += f[i]
        bin_counts[idx] += 1
    end
    
    # Vectorized averaging
    mask = bin_counts .> 0
    f_coarse[mask] ./= bin_counts[mask]
    
    s_coarse = s_min .+ (0.5:n_target-0.5) * bin_size
    
    return f_coarse, s_coarse
end

"""
    resample_contour_periodic(points_x, points_y, Ly, n_segments; functions=nothing)

Resample a contour uniformly in length segments, assuming periodic boundary 
conditions in the y-dimension. Optionally resample functions defined on the contour.

# Arguments
- `points_x`: x-coordinates of original contour points
- `points_y`: y-coordinates of original contour points  
- `Ly`: period length in y-dimension
- `n_segments`: number of uniform segments to create
- `functions`: optional tuple/vector of functions to resample (e.g., (f1, f2, f3))

# Returns
- `new_x`: x-coordinates of resampled points
- `new_y`: y-coordinates of resampled points (within [0, Ly))
- `resampled_functions`: resampled function values (only if functions provided)
"""
function resample_contour_periodic(points_x::Vector{T}, points_y::Vector{T}, 
                                 Ly::T, n_segments::Int; 
                                 functions=nothing) where T<:Real
    n_points = length(points_x)
    @assert length(points_y) == n_points "x and y must have same length"
    @assert n_points >= 2 "Need at least 2 points"
    
    # Create extended arrays including the periodic closure
    if points_y[end] < points_y[1]
        points_x, points_y = points_x[end:-1:1], points_y[end:-1:1]  # Reverse order
    end
    x_ext = vcat(points_x, points_x[1])
    y_ext = vcat(points_y, points_y[1] + Ly)  # Periodic closure
    
    # Calculate segment lengths vectorized
    dx = diff(x_ext)
    dy = diff(y_ext)
    segment_lengths = sqrt.(dx.^2 + dy.^2)
    
    # Cumulative arc length
    cumulative_length = vcat(0.0, cumsum(segment_lengths))
    total_length = cumulative_length[end]
    
    # Target arc lengths for uniform resampling
    target_lengths = range(0, total_length, length=n_segments+1)[1:end-1]
    
    # Find which segments contain each target length
    segment_indices = searchsortedlast.(Ref(cumulative_length), target_lengths)
    segment_indices = max.(segment_indices, 1)  # Ensure at least 1
    segment_indices = min.(segment_indices, length(segment_lengths))  # Ensure not beyond
    
    # Calculate interpolation parameters
    s0 = cumulative_length[segment_indices]
    s1 = cumulative_length[segment_indices .+ 1]
    
    # Avoid division by zero for segments with zero length
    segment_lens = s1 - s0
    safe_lens = max.(segment_lens, eps(T))
    t = (target_lengths - s0) ./ safe_lens
    
    # Interpolate coordinates
    x0 = x_ext[segment_indices]
    x1 = x_ext[segment_indices .+ 1]
    y0 = y_ext[segment_indices]
    y1 = y_ext[segment_indices .+ 1]
    
    new_x = x0 + t .* (x1 - x0)
    new_y = y0 + t .* (y1 - y0)
    
    # Apply periodic boundary condition to y-coordinates
    new_y = mod.(new_y, Ly)
    
    # Resample functions if provided
    if functions !== nothing
        # Handle both single function and multiple functions
        if isa(functions, Vector) && length(functions[1]) == 1
            # Single function
            func = functions
            @assert length(func) == n_points "Function must have same length as coordinate arrays"
            
            # Extend function with periodic closure
            func_ext = vcat(func, func[1])
            
            # Interpolate function values
            f0 = func_ext[segment_indices]
            f1 = func_ext[segment_indices .+ 1]
            new_func = f0 + t .* (f1 - f0)
            
            return vec(new_x), vec(new_y), vec(new_func)
        else
            # Multiple functions
            resampled_functions = []
            for func in functions
                @assert length(func) == n_points "Function must have same length as coordinate arrays"
                
                # Extend function with periodic closure
                func_ext = vcat(func, func[1])
                
                # Interpolate function values
                f0 = func_ext[segment_indices]
                f1 = func_ext[segment_indices .+ 1]
                new_func = f0 + t .* (f1 - f0)
                
                push!(resampled_functions, new_func)
            end
            
            # Return as tuple if input was tuple, vector if input was vector
            if isa(functions, Tuple)
                resampled_functions = Tuple(resampled_functions)
            end
            
            return vec(new_x), vec(new_y), resampled_functions
        end
    else
        return vec(new_x), vec(new_y)
    end
end

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
        arclength = T(0.0)
        curvature = T(0.0)
    else
        # Calculate segment lengths
        dx = diff(center_x)
        dy = diff(center_y)
        segment_lengths = sqrt.(dx.^2 .+ dy.^2)
        arclength = T.(cumsum(segment_lengths))
        
        # Compute curvature using finite differences
        curvature = T.(compute_curvature(center_x, center_y))
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
    # MLD_front = vec(idx_front .* interior(snapshot[:MLD3], :, :, 1))
    # Ew_front = vec(idx_front .* interior(snapshot[:Ew3], :, :, 1))
    # Ew_front = Ew_front[MLD_front.>1e-6]
    # MLD_front = MLD_front[MLD_front.>1e-6]

    # axislegend(axT, labelsize=9, padding = (0f0, 0f0, 0f0, 0f0), patchlabelgap = 3)
    # colgap!(fig.layout, 1, Relative(0)) 
    # colgap!(fig.layout, 3, Relative(0)) 
    # colgap!(fig.layout, 5, Relative(0)) 

    # resize_to_layout!(fig)
    # save(filesave * "T0_MLD3_Ew3_" * fileparam * "_$(iteration).pdf", fig; pt_per_unit = 1)

    return mean(clevels), levels_T, points_x, points_y, idx_front #MLD_front, Ew_front
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


"""
    cross_correlation_periodic(f, g, total_length; normalize=true, confidence_level=0.95)

Compute cross-correlation of two periodic functions f and g defined on 
uniformly spaced points along a contour.

Cross-correlation R_fg(s) = Σ f(i) * g(i-k) where g is shifted by distance s.
Positive s means g is shifted forward along the contour.

# Arguments
- `f`: first function values at uniform points
- `g`: second function values at uniform points  
- `total_length`: total length of the contour
- `normalize`: if true, normalize by standard deviations (Pearson correlation)
- `confidence_level`: confidence level for bounds (e.g., 0.95 for 95% confidence)

# Returns
- `xcf`: cross-correlation values for each offset
- `length_offsets`: corresponding length offsets in [-total_length/2, total_length/2]
- `upper_bound`: upper confidence bound (scalar)
- `lower_bound`: lower confidence bound (scalar)
"""
function cross_correlation_periodic(f::Vector{T}, g::Vector{T}, total_length::T; 
                                  normalize::Bool=true, confidence_level::T=Float32(0.95)) where T<:Real
    n = length(f)
    @assert length(g) == n "f and g must have same length"
    
    # Center the signals for correlation calculation
    f_centered = f .- mean(f)
    g_centered = g .- mean(g)
    
    if normalize
        # Compute standard deviations for normalization
        σ_f = std(f_centered)
        σ_g = std(g_centered)
        
        # Handle constant signals
        if σ_f ≈ 0 || σ_g ≈ 0
            return zeros(T, n), collect(0:n-1)
        end
        
        f_centered ./= σ_f
        g_centered ./= σ_g
    end
    
    # Use FFT for efficient circular correlation
    # For cross-correlation R_fg(k) = Σ f(i) * g(i-k), we need:
    # calculate covariance function
    cov = irfft(conj.(rfft(f_centered)).*rfft(g_centered), n)

    # indicator functions
    ind1 = ones(n)
    ind2 = ones(n)

    # calculate normalization factors
    norm1 = irfft(conj.(rfft(f_centered.^2)).*rfft(ind2), n)
    norm2 = irfft(conj.(rfft(ind1)).*rfft(g_centered.^2), n)

    # exclude negative and zero normalization factors (due to roundoff error and non-overlapping segments)
    norm1[norm1.≤0] .= Inf
    norm2[norm2.≤0] .= Inf

  # cross-correlation
    xcf = cov./sqrt.(norm1.*norm2)
    
    # Convert to length-based offsets centered at zero
    segment_length = total_length / n
    index_offsets = collect(0:n-1)
    
    # Map indices to length offsets: [0, 1, ..., n-1] -> [-L/2, ..., L/2]
    # For even n: indices 0 to n/2-1 map to [0, L/2), indices n/2 to n-1 map to [-L/2, 0)
    # For odd n: indices 0 to (n-1)/2 map to [0, L/2], indices (n+1)/2 to n-1 map to [-L/2, 0)
    length_offsets = similar(xcf)
    half_n = n ÷ 2
    
    # First half: positive offsets [0, L/2)
    for i in 1:(half_n + (n % 2))  # Include middle point for odd n
        length_offsets[i] = (i - 1) * segment_length
    end
    
    # Second half: negative offsets [-L/2, 0)
    for i in (half_n + (n % 2) + 1):n
        length_offsets[i] = (i - 1 - n) * segment_length
    end
    
    # Reorder xcf to match length_offsets ordering: [0, +, -, 0)
    # FFT output order: [0, 1, 2, ..., n-1] corresponds to [0, +δs, +2δs, ..., -δs]
    # We want: [-L/2, ..., -δs, 0, +δs, ..., +L/2)
    if n % 2 == 0
        # Even n: reorder to [n/2:n-1, 0:n/2-1]
        reorder_idx = vcat(half_n+1:n, 1:half_n)
    else
        # Odd n: reorder to [(n+1)/2:n-1, 0:(n-1)/2]  
        reorder_idx = vcat(half_n+2:n, 1:half_n+1)
    end
    
    xcf = xcf[reorder_idx]
    length_offsets = sort(length_offsets)  # Ensure proper ordering
    
    # Calculate confidence bounds (Bartlett's formula for large samples)
    # Under null hypothesis of no correlation: σ²(ρ̂) ≈ 1/n
    # Use normal approximation: ρ̂ ± z_{α/2} * σ(ρ̂)
    if confidence_level > 0 && confidence_level < 1
        # Standard error under null hypothesis
        std_error = 1 / sqrt(n)
        
        # Critical value for two-tailed test
        α = 1 - confidence_level
        z_critical = sqrt(-2 * log(α/2))  # Approximation: z ≈ sqrt(-2*ln(α/2)) for α small
        
        # More accurate critical values for common confidence levels
        if abs(confidence_level - 0.95) < 1e-10
            z_critical = 1.96
        elseif abs(confidence_level - 0.99) < 1e-10
            z_critical = 2.576
        elseif abs(confidence_level - 0.90) < 1e-10
            z_critical = 1.645
        end
        
        bound = z_critical * std_error
        upper_bound = bound
        lower_bound = -bound
    else
        upper_bound = T(Inf)
        lower_bound = T(-Inf)
    end
    
    return xcf, length_offsets, upper_bound, lower_bound
end

function get_time_series(filehead,filesave;fileparam = "xband1sublevels",use_gpu = true,savedata=true,plotdata=true)

    x=rand(Float32,10240,20480,1);
    if use_gpu && CUDA.functional()
        p = CUDA.CUFFT.plan_rfft(CuArray(x), (1,2));
        ip = CUDA.CUFFT.plan_irfft(p * CuArray(x), 10240, (1,2));
    else
        p = FFTW.plan_rfft(x, (1,2); flags=FFTW.MEASURE, timelimit=Inf);
        ip = FFTW.plan_irfft(p * x, 10240, (1,2); flags=FFTW.MEASURE, timelimit=Inf);
    end

    t = [11; 16; 21; 24:80];
    if isfile(filesave * "ts_d_h_Cub_Ewmax_" * fileparam * "_cg3hm.jld2")
        file = jldopen(filesave * "ts_d_h_Cub_Ewmax_" * fileparam * "_cg3hm.jld2", "r");
        width3 = [ones(1,3); file["widths"]];
        MLD3_front = [60*ones(1,3); file["MLDs_front"]];
        Cu4 = [zeros(1,4); file["Cus"]];
        Ewmax = [0; file["Ewmax"]];
        close(file)
        idx = [true; (.!(iszero.(Ewmax[2:end])))];
        t = [0; t][idx]
    else
        iterations = get_iterations_regex(filehead, fileparam)
        #iterations = iterations[t .>= 57]
        niter = length(iterations)
        width3 = zeros(niter,3)
        MLD3_front = zeros(niter,3)
        Ew3_front = zeros(niter,3)
        rcu3 = zeros(niter,3)
        Cu4 = zeros(niter,4)
        Ewmax = zeros(niter)
        @info "Iterations: $iterations"
        for (i,iteration) in enumerate(iterations)
            @info "Iteration: $iteration"
            # 1. Define the filename of the saved snapshot
            output_filename = filehead * "subdomains/" * fileparam * "_snapshot_iter$(iteration).jld2"

            # 2. Load the snapshot using the new function
            snapshot = load_subdomain_snapshot(output_filename;variables = ("T","MLD3","Ew3"),level=224);

            # Call the function to get the analysis results
            clevel, levels_T, points_x, points_y, MLD_front, Ew_front = analyze_frontal_properties(snapshot,iteration;use_gpu,plans=(p,ip));

            if length(levels_T) < 1
                break
            end
            center_idx = argmin(abs.(levels_T.- clevel))
            widths, arclength, curvature = compute_contour_properties(points_x, points_y, center_idx)
            if !isnothing(widths)
                d̅ = mean(widths)/2e3 # normalized width
                d_10percent = quantile(widths, 0.1)/2e3 # normalized width
                d_90percent = quantile(widths, 0.9)/2e3 # normalized width
            else
                break
            end
            h̅ = mean(MLD_front) # mean mixed layer depth
            h_10percent = quantile(MLD_front, 0.1) # 10th percentile mixed layer depth
            h_90percent = quantile(MLD_front, 0.9) # 90th percentile mixed layer depth
            E̅ = mean(Ew_front) # mean Ew
            E_10percent = quantile(Ew_front, 0.1) # 10th percentile Ew
            E_90percent = quantile(Ew_front, 0.9) # 90th percentile Ew
            Ewmax[i] = maximum(interior(snapshot[:Ew3]))
            @info "T levels: $levels_T"
            @info "Total arc length: $(arclength[end]/1e3) km"
            @info "Normalized width: $d̅ (10th percentile: $d_10percent, 90th percentile: $d_90percent)"
            @info "Mixed layer depth: $h̅ (10th percentile: $h_10percent, 90th percentile: $h_90percent)"
            @info "Depth-averaged w²/2: $E̅ (10th percentile: $E_10percent, 90th percentile: $E_90percent)"
            itp = interpolate((nodes(snapshot[:T])[1],nodes(snapshot[:T])[2]), interior(snapshot[:MLD3], :, :, 1), Gridded(Constant()))
            h_cfront = itp.(points_x[center_idx], points_y[center_idx])
            Cu = 2α * g * (max(levels_T...) - min(levels_T...)) / f^2 * h_cfront ./ widths .* curvature
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
            @info "Memory usage: " * (Sys.free_memory() |> Base.format_bytes) * " GB free"
        end

        if savedata
            jldopen(filesave * "ts_d_h_Cub_Ewmax_" * fileparam * "_cg3hm.jld2", "w") do file
                file["widths"] = width3
                file["MLDs_front"] = MLD3_front
                file["Cus"] = Cu4
                file["Ewmax"] = Ewmax 
                file["t"] = t
            end
        end
    end

    if plotdata
        fig = Figure(size=(640, 640))
        gl = fig[1, 1] = GridLayout()
        axis_kwargs = (titlealign = :left, limits=((0,61),nothing))
        axw = Axis(gl[1, 1]; title=L"\text{(a)}", ylabel=L"d/d_0", axis_kwargs...)
        axh = Axis(gl[2, 1]; title=L"\text{(b)}", ylabel=L"h/h_0", axis_kwargs...)
        axc = Axis(gl[3, 1]; title=L"\text{(c)}", ylabel=L"Cu_b", axis_kwargs...)
        axE = Axis(gl[4, 1]; title=L"\text{(d)}", ylabel=L"\overline{(w^2/2)}^h_{\max}/w_*^2", xlabel=L"t~\text{(h)}", axis_kwargs...)
        hidexdecorations!(axw, ticks=false)
        hidexdecorations!(axc, ticks=false)
        hidexdecorations!(axh, ticks=false)
        lines!(axw, t, width3[idx,1], color=wcolors[1])
        fill_between!(axw, t, width3[idx,2], width3[idx,3]; color = wcolors[1], alpha = 0.25)
        lines!(axh, t, MLD3_front[idx,1]/h₀, color=wcolors[2])
        fill_between!(axh, t, MLD3_front[idx,2]/h₀, MLD3_front[idx,3]/h₀; color = wcolors[2], alpha = 0.25)
        # lines!(axc, t, Cu3[:,1], color=wcolors[3])
        # fill_between!(axc, t, Cu3[:,2], Cu3[:,3]; color = wcolors[3], alpha = 0.25)
        lines!(axc, t, Cu4[idx,4], color=wcolors[3])
        # lines!(axb, t, Cu4[:,1], color=wcolors[4])
        lines!(axE, t, Ewmax[idx]/wₛ^2, color=wcolors[4])
        #fill_between!(axE, t, Ew3_front[:,2]/wₛ^2, Ew3_front[:,3]/wₛ^2; color = wcolors[5], alpha = 0.25)
        kwargs = (; markersize = 10, strokewidth = 1, color = :white)
        scatter!(axc, 44, Cu4[argmin(abs.(t.- 44)),4]; kwargs..., strokecolor = wcolors[3])
        scatter!(axE, 30, Ewmax[argmin(abs.(t.- 30))]/wₛ^2; kwargs..., strokecolor = wcolors[4])
        scatter!(axE, 49, Ewmax[argmin(abs.(t.- 49))]/wₛ^2; kwargs..., strokecolor = wcolors[4])
        rowgap!(gl, 3)
        resize_to_layout!(fig)
        save(filesave * "ts_d_h_Cub_Ewmax_" * fileparam * "_cg3hm.pdf", fig; pt_per_unit = 1)
    end
end

function plot_front_properties(filehead,fileparam,iteration,Nresample;use_gpu=true)
    @info "Iteration: $iteration"

    x=rand(Float32,10240,20480,1);
    if use_gpu && CUDA.functional()
        p = CUDA.CUFFT.plan_rfft(CuArray(x), (1,2));
        ip = CUDA.CUFFT.plan_irfft(p * CuArray(x), 10240, (1,2));
    else
        p = FFTW.plan_rfft(x, (1,2); flags=FFTW.MEASURE, timelimit=Inf);
        ip = FFTW.plan_irfft(p * x, 10240, (1,2); flags=FFTW.MEASURE, timelimit=Inf);
    end
    println(CUDA.pool_status())
    # 1. Define the filename of the saved snapshot
    output_filename = filehead * "subdomains/" * fileparam * "_snapshot_iter$(iteration).jld2"

    # 2. Load the snapshot using the new function
    snapshot = load_subdomain_snapshot(output_filename;variables = ("u","v","T","MLD3","Ew3"),level=224);

    x, y, z = nodes(snapshot[:T])
    idxEwm = argmax(interior(snapshot[:Ew3]))
    grid = snapshot[:u].grid
    T̅  = CenterField(grid);
    u̅  = XFaceField(grid);
    v̅  = YFaceField(grid);
    coarse_graining!(snapshot[:T] , T̅; kernel=:gaussian, cutoff=300, border = :ycircular, method = :spectral, use_gpu,plans=(p,ip));
    coarse_graining!(snapshot[:v] , v̅; kernel=:gaussian, cutoff=300, border = :ycircular, method = :spectral, use_gpu,plans=(p,ip));
    coarse_graining!(snapshot[:u] , u̅; kernel=:gaussian, cutoff=300, border = :ycircular, method = :spectral, use_gpu,plans=(p,ip));
    u̅g, v̅g = c2g(u̅), c2g(v̅)
    ζ̅g = (compute!(Field(∂x(v̅g) - ∂y(u̅g))));
    δ̅g = (compute!(Field(∂x(u̅g) + ∂y(v̅g))));
    ζ̅, δ̅ = g2c(ζ̅g), g2c(δ̅g)
    b̅g = compute!(Field(α * g * c2g(T̅)));
    B̅h = g2c(compute!(Field(-(∂x(b̅g)^2 * ∂x(u̅g) + ∂y(b̅g)^2 * ∂y(v̅g))-∂x(b̅g)*∂y(b̅g)*(∂x(v̅g) + ∂y(u̅g)))));
    Db̅² = g2c(compute!(Field(∂x(b̅g)^2 + ∂y(b̅g)^2)));

    # Call the function to get the analysis results
    clevel, levels_T, points_x, points_y, _ = analyze_frontal_properties(snapshot,iteration;use_gpu,plans=(p,ip));

    println(CUDA.pool_status())
    @info "T levels: $(levels_T)"
    center_idx = argmin(abs.(levels_T.- clevel))
    widths, arclength, curvature = compute_contour_properties(points_x, points_y, center_idx)
    @info "curvature head raw: $(curvature[1:10])"
    @info "Memory usage: " * (Sys.free_memory() |> Base.format_bytes) * " GB free"
    cps_x, cps_y, cd_Cu = resample_contour_periodic(points_x[center_idx], points_y[center_idx], Float32(1e5), Nresample;functions=(widths,curvature))
    @info "curvature resampled extrema: $(extrema(cd_Cu[2]))"
    @info "x head: $(cps_x[1:10]/1e3); y head: $(cps_y[1:10]/1e3)"
    @info "$(points_x[center_idx][1:10]/1e3), $(points_y[center_idx][1:10]/1e3)"

    Sn = Float32.(normal_strain_rate(cps_x, cps_y, etpU, etpV));
    @info "cps_x length: $(length(cps_x)), cps_y length: $(length(cps_y)), Sn length: $(length(Sn))"
    σn_cfront = Float32.(etpn.(cps_x, cps_y));
    itph = interpolate((Float32.(x),Float32.(y)), Float32.(interior(snapshot[:MLD3], :, :, 1)), Gridded(Constant(Interpolations.Periodic())))
    etph = extrapolate(itph, Interpolations.Periodic())
    h_cfront = Float32.(etph.(cps_x, cps_y));
    Cu = 2α * g * (max(levels_T...) - min(levels_T...)) / f^2 * h_cfront ./ cd_Cu[1] .* cd_Cu[2];
    @info "Cu extrema: $(extrema(Cu))"
    xζ, yζ, _ = nodes(ζ̅);
    itpζ = interpolate((Float32.(xζ),Float32.(yζ)), Float32.(interior(ζ̅, :, :, 1)), Gridded(Constant(Interpolations.Periodic())))
    etpζ = extrapolate(itpζ, Interpolations.Periodic())
    ζ_cfront = Float32.(etpζ.(cps_x, cps_y));
    xδ, yδ, _ = nodes(δ̅);
    itpδ = interpolate((Float32.(xδ),Float32.(yδ)), Float32.(interior(δ̅, :, :, 1)), Gridded(Constant(Interpolations.Periodic())))
    etpδ = extrapolate(itpδ, Interpolations.Periodic())
    δ_cfront = Float32.(etpδ.(cps_x, cps_y));
    xB, yB, _ = nodes(B̅h);
    itpB = interpolate((Float32.(xB),Float32.(yB)), Float32.(interior(B̅h, :, :, 1)), Gridded(Constant(Interpolations.Periodic())))
    etpB = extrapolate(itpB, Interpolations.Periodic())
    B_cfront = Float32.(etpB.(cps_x, cps_y));
    xDb, yDb, _ = nodes(Db̅²);
    itpDb = interpolate((Float32.(xDb),Float32.(yDb)), Float32.(interior(Db̅², :, :, 1)), Gridded(Constant(Interpolations.Periodic())))
    etpDb = extrapolate(itpDb, Interpolations.Periodic())
    D_cfront = Float32.(etpDb.(cps_x, cps_y));

    A=rand(Float32,10240,20480,4);
    if use_gpu && CUDA.functional()
        p = CUDA.CUFFT.plan_rfft(CuArray(A), (1,2));
        pk = CUDA.CUFFT.plan_rfft(CuArray(A[:,:,1:1]), (1,2));
        ip = CUDA.CUFFT.plan_irfft(p * CuArray(A), 10240, (1,2));
    else
        @error "CUDA not available"
    end
    println(CUDA.pool_status())
    fileparam = "xband1surf4"
    output_filename = filehead * "subdomains/" * fileparam * "_snapshot_iter$(iteration).jld2"
    snapshot2 = load_subdomain_snapshot(output_filename; variables = ("u", "v", "w"));
    u̅, v̅, w̅, τuu, τvv, τww = TKE(snapshot2; cutoff=300, border=:ycircular, Lx = snapshot2[:grid].Lx, Ly = snapshot2[:grid].Ly,
                                        method=:spectral,use_gpu,plans=(p,pk,ip));
    t0 = time()
    TKEₛg = (compute!(Field((τvv + τuu + τww)/2/wₛ^2)));
    @info "TKE compute time: $(time() - t0)"
    @info "Memory usage: " * (Sys.free_memory() |> Base.format_bytes) * " GB free"
    SKEₛ = ((u̅ - mean(u̅;dims=2))^2 + (v̅ - mean(v̅;dims=2))^2 + (w̅ - mean(w̅;dims=2))^2)/2/wₛ^2;
    SKEₛg = (compute!(Field(SKEₛ)));
    @info "SKE compute time: $(time() - t0)"
    @info "Memory usage: " * (Sys.free_memory() |> Base.format_bytes) * " GB free"
    TKEₛ, SKEₛ = g2c(TKEₛg), g2c(SKEₛg)
    idxTKEm = argmax(interior(TKEₛ))
    idxSKEm = argmax(interior(SKEₛ))
    xc, yc, zc = nodes(TKEₛ);
    xd, yd, zd = nodes(SKEₛ);
    @info "TKE maximum at x=$(xc[idxTKEm[1]]/1e3)km, y=$(yc[idxTKEm[2]]/1e3)km, z=$(zc[idxTKEm[3]])m"
    @info "SKE maximum at x=$(xd[idxSKEm[1]]/1e3)km, y=$(yd[idxSKEm[2]]/1e3)km, z=$(zd[idxSKEm[3]])m"
    k = length(zT)
    println(CUDA.pool_status())

    fig = Figure(size = (640, 450))
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
    hmc = heatmap!(axc, 1e-3xc, 1e-3yc, interior(TKEₛ, :, :, idxTKEm[3]); rasterize = true, colormap = :amp, colorrange = (0, max(interior(TKEₛ, :, :, idxTKEm[3])...)))
    hmd = heatmap!(axd, 1e-3xd, 1e-3yd, interior(SKEₛ, :, :, idxSKEm[3]); rasterize = true, colormap = :amp, colorrange = (0, max(interior(SKEₛ, :, :, idxSKEm[3])...)))
    scatter!(axb, 1e-3x[idxEwm[1]], 1e-3y[idxEwm[2]]; marker = :star4, markersize = 10, color = :white)
    scatter!(axc, 1e-3xc[idxTKEm[1]], 1e-3yc[idxTKEm[2]]; marker = :star4, markersize = 10, color = :black)
    scatter!(axd, 1e-3xd[idxSKEm[1]], 1e-3yd[idxSKEm[2]]; marker = :star4, markersize = 10, color = :black)
    Colorbar(fig[1, 2], hmT)
    Colorbar(fig[1, 4], hmb)
    Colorbar(fig[1, 6], hmc)
    Colorbar(fig[1, 8], hmd)
    for j = 1:3
        xs, ys = points_x[j], points_y[j]
        lvl = levels_T[j]
        CairoMakie.lines!(axT, 1e-3xs, 1e-3ys; label=@sprintf("%.2f°C", lvl), linewidth=0.8) 
    end
    #CairoMakie.lines!(axT, 1e-3cps_x, 1e-3cps_y; linewidth=0.8, linestyle=:dash, color=:black) 
    #axislegend(axT, labelsize=9, padding = (0f0, 0f0, 0f0, 0f0), patchlabelgap = 3)
    for j = 1:2:7
        colgap!(fig.layout, j, Relative(0)) 
    end
    for j = 2:2:6
        colgap!(fig.layout, j, 3) 
    end
    iarrow = 2:32:640
    arrows!(axb,xT[iarrow]/1e3.-50, yT[iarrow]/1e3, shift(interior(Ub,iarrow,iarrow,k)), shift(interior(Vb,iarrow,iarrow,k)), arrowsize = 3, lengthscale = 1e2,linecolor = :white, arrowcolor = :white, linewidth = 0.2,alpha=0.5)
    arrows!(axc,xT[iarrow]/1e3.-50, yT[iarrow]/1e3, shift(interior(Ub,iarrow,iarrow,k)), shift(interior(Vb,iarrow,iarrow,k)), arrowsize = 3, lengthscale = 1e2,linecolor = :black, arrowcolor = :black, linewidth = 0.2,alpha=0.5)
    arrows!(axd,xT[iarrow]/1e3.-50, yT[iarrow]/1e3, shift(interior(Ub,iarrow,iarrow,k)), shift(interior(Vb,iarrow,iarrow,k)), arrowsize = 3, lengthscale = 1e2,linecolor = :black, arrowcolor = :black, linewidth = 0.2,alpha=0.5)
    hlines!(axc, [18.75, 43.75, 68.75, 93.75], linestyle = :dash, color = :black, linewidth = 0.8)
    hlines!(axd, [18.75, 43.75, 68.75, 93.75], linestyle = :dash, color = :black, linewidth = 0.8)
    resize_to_layout!(fig)
    save(filesave * "T0_MLD3_TKE_SKE_" * fileparam * "_$(iteration).pdf", fig; pt_per_unit = 1)

    itpTKE = interpolate((Float32.(xc),Float32.(yc)), Float32.(interior(TKEₛ, :, :, length(zc))), Gridded(Constant(Interpolations.Periodic())))
    etpTKE = extrapolate(itpTKE, Interpolations.Periodic())
    TKE_cfront = Float32.(etpTKE.(cps_x, cps_y));
    itpSKE = interpolate((Float32.(xd),Float32.(yd)), Float32.(interior(SKEₛ, :, :, length(zd))), Gridded(Constant(Interpolations.Periodic())))
    etpSKE = extrapolate(itpSKE, Interpolations.Periodic())
    SKE_cfront = Float32.(etpSKE.(cps_x, cps_y));
    # itpEw = interpolate((x,y), interior(snapshot[:Ew3], :, :, 1), Gridded(Constant()))
    # Ew_cfront = itpEw.(cps_x, cps_y);

    # Calculate total contour length
    dx = diff(vcat(cps_x, cps_x[1]))
    dy = diff(vcat(cps_y, cps_y[1] + 1e5))  # Account for periodicity
    total_length = cumsum(sqrt.(dx.^2 + dy.^2))
    arclength = Float32.([0;total_length[1:end-1]])

    fig = Figure(size=(640, 720))
    limits = ((0, arclength[end]/1e5),nothing)
    axis_kwargs = (limits, xgridvisible = false,ygridvisible = false)
    axa = Axis(fig[1, 1]; titlealign = :left, title=L"\text{(a)}", ylabel=L"S_n/f", yticklabelcolor = wcolors[1], axis_kwargs...)
    axa2 = Axis(fig[1, 1]; ylabel=L"\sigma_n/f", yticklabelcolor = wcolors[2], yaxisposition = :right, axis_kwargs...)
    axb = Axis(fig[2, 1]; titlealign = :left, title=L"\text{(b)}", ylabel=L"h/h_0", yticklabelcolor = wcolors[1], axis_kwargs...)
    axb2 = Axis(fig[2, 1]; ylabel=L"d/d_0", yticklabelcolor = wcolors[2], yaxisposition = :right, axis_kwargs...)
    hidespines!(axb2);hidexdecorations!(axb2);
    axc = Axis(fig[3, 1]; titlealign = :left, title=L"\text{(c)}", ylabel=L"Cu", yticklabelcolor = wcolors[1], axis_kwargs...)
    axc2 = Axis(fig[3, 1]; ylabel=L"\overline{\zeta}/f", yticklabelcolor = wcolors[2], yaxisposition = :right, axis_kwargs...)
    axd = Axis(fig[4, 1]; titlealign = :left, title=L"\text{(d)}", ylabel=L"\overline{\delta}/f", yticklabelcolor = wcolors[1], axis_kwargs...)
    axd2 = Axis(fig[4, 1]; ylabel=L"\mathcal{F}_s/(f M_0^4)", yticklabelcolor = wcolors[2], yaxisposition = :right, axis_kwargs...)
    hidespines!(axd2);hidexdecorations!(axd2);
    axe = Axis(fig[5, 1]; titlealign = :left, title=L"\text{(e)}", xlabel=L"s/L_y", ylabel=L"\text{TKE}/w_*^2", yticklabelcolor = wcolors[1], axis_kwargs...)
    axe2 = Axis(fig[5, 1]; ylabel=L"\text{SKE}/w_*^2", yticklabelcolor = wcolors[2], yaxisposition = :right,axis_kwargs...)
    hidexdecorations!(axa, ticks=false)
    hidexdecorations!(axb, ticks=false)
    hidexdecorations!(axc, ticks=false)
    hidexdecorations!(axd, ticks=false)
    Sn_coarse, s_coarse = coarsen_binned_vectorized(Sn, arclength, 256)
    lines!(axa, arclength/1e5, Sn/f; linewidth=1, color = wcolors[1])
    lines!(axa, s_coarse/1e5, Sn_coarse/f; linewidth=1, color = wcolors[1], linestyle = :dash)
    lines!(axa2, arclength/1e5, σn_cfront/f; linewidth=1, color = wcolors[2])
    lines!(axb, arclength/1e5, h_cfront/60; linewidth=1, color = wcolors[1])
    lines!(axb2, arclength/1e5, cd_Cu[1]/2e3; linewidth=1, color = wcolors[2])
    lines!(axc, arclength/1e5, Cu; linewidth=1, color = wcolors[1])
    lines!(axc2, arclength/1e5, ζ_cfront/f; linewidth=1, color = wcolors[2])
    lines!(axd, arclength/1e5, δ_cfront/f; linewidth=1, color = wcolors[1])
    lines!(axd2, arclength/1e5, B_cfront/(f*M²₀^2); linewidth=1, color = wcolors[2])
    lines!(axe, arclength/1e5, TKE_cfront; linewidth=1, color = wcolors[1])
    #lines!(axe, arclength/1e5, Ew_cfront/wₛ^2; linewidth=1, color = wcolors[1], linestyle = :dash)
    lines!(axe2, arclength/1e5, SKE_cfront; linewidth=1, color = wcolors[2])
    rowgap!(fig.layout, 3)
    resize_to_layout!(fig)
    save(filesave * "along_front_" * fileparam * "_$(iteration).pdf", fig; pt_per_unit = 1)

    Ls = Float32.(total_length[end])
    cc_S_d, l_offsets, upper_bound, lower_bound = cross_correlation_periodic(σn_cfront, Float32.(cd_Cu[1]), Ls;confidence_level=Float32(0.99))
    cc_S_h, _, _, _ = cross_correlation_periodic(σn_cfront, h_cfront, Ls)
    cc_h_d, _, _, _ = cross_correlation_periodic(h_cfront,Float32.(cd_Cu[1]), Ls)
    cc_S_Cu, _, _, _ = cross_correlation_periodic(abs.(σn_cfront), Float32.(Cu), Ls)
    cc_S_ζ, _, _, _ = cross_correlation_periodic(abs.(σn_cfront), ζ_cfront, Ls)
    cc_Cu_ζ, _, _, _ = cross_correlation_periodic(Float32.(Cu), ζ_cfront, Ls)
    cc_S_δ, _, _, _ = cross_correlation_periodic(abs.(σn_cfront), δ_cfront, Ls)
    cc_S_B, _, _, _ = cross_correlation_periodic(abs.(σn_cfront), B_cfront, Ls)
    cc_δ_B, _, _, _ = cross_correlation_periodic(δ_cfront, B_cfront, Ls)
    cc_S_TKE, _, _, _ = cross_correlation_periodic(abs.(σn_cfront), TKE_cfront, Ls)
    cc_S_SKE, _, _, _ = cross_correlation_periodic(abs.(σn_cfront), SKE_cfront, Ls)
    cc_TKE_SKE, _, _, _ = cross_correlation_periodic(TKE_cfront, SKE_cfront, Ls)

    fig = Figure(size=(640, 640))
    xlimit = (l_offsets[1]/1e5, l_offsets[end]/1e5)
    ylimita = map(x -> x * 1.2, extrema([-cc_S_d; cc_S_h; -cc_h_d]))
    ylimitb = map(x -> x * 1.2, extrema([cc_S_ζ; cc_S_Cu; cc_Cu_ζ]))
    ylimitc = map(x -> x * 1.2, extrema([-cc_S_δ; cc_S_B; -cc_δ_B]))
    ylimitd = map(x -> x * 1.2, extrema([cc_S_TKE; cc_S_SKE; cc_TKE_SKE]))
    axis_kwargs = (titlealign = :left, xgridvisible = false)
    axa = Axis(fig[1, 1]; title=L"\text{(a)}", limits = (xlimit, ylimita), ylabel=L"CC(\cdot,\cdot)", axis_kwargs...)
    axb = Axis(fig[2, 1]; title=L"\text{(b)}", limits = (xlimit, ylimitb), ylabel=L"CC(\cdot,\cdot)", axis_kwargs...)
    axc = Axis(fig[3, 1]; title=L"\text{(c)}", limits = (xlimit, ylimitc), ylabel=L"CC(\cdot,\cdot)", axis_kwargs...)
    axd = Axis(fig[4, 1]; title=L"\text{(d)}", limits = (xlimit, ylimitd), xlabel=L"\Delta s/L_y", ylabel=L"CC(\cdot,\cdot)", axis_kwargs...)
    hidexdecorations!(axa, ticks=false)
    hidexdecorations!(axb, ticks=false)
    lines!(axa, l_offsets/1e5, cc_S_h; linewidth=1, label=L"(\sigma_n,h)", color = wcolors[1])
    lines!(axa, l_offsets/1e5, -cc_S_d; linewidth=1, label=L"(\sigma_n,-d)", color = wcolors[2])
    lines!(axa, l_offsets/1e5, -cc_h_d; linewidth=1, label=L"(h,-d)", color = wcolors[3])
    lines!(axb, l_offsets/1e5, cc_S_Cu; linewidth=1, label=L"(|\sigma_n|,Cu)", color = wcolors[1])
    lines!(axb, l_offsets/1e5, cc_S_ζ; linewidth=1, label=L"(|\sigma_n|,\overline{\zeta})", color = wcolors[2])
    lines!(axb, l_offsets/1e5, cc_Cu_ζ; linewidth=1, label=L"(Cu,\overline{\zeta})", color = wcolors[3])
    lines!(axc, l_offsets/1e5, -cc_S_δ; linewidth=1, label=L"(|\sigma_n|,-\overline{\delta})", color = wcolors[1])
    lines!(axc, l_offsets/1e5, cc_S_B; linewidth=1, label=L"(|\sigma_n|,\mathcal{F}_s)", color = wcolors[2])
    lines!(axc, l_offsets/1e5, -cc_δ_B; linewidth=1, label=L"(-\overline{\delta},\mathcal{F}_s)", color = wcolors[3])
    lines!(axd, l_offsets/1e5, cc_S_TKE; linewidth=1, label=L"(|\sigma_n|,\text{TKE})", color = wcolors[1])
    lines!(axd, l_offsets/1e5, cc_S_SKE; linewidth=1, label=L"(|\sigma_n|,\text{SKE})", color = wcolors[2])
    lines!(axd, l_offsets/1e5, cc_TKE_SKE; linewidth=1, label=L"(\text{TKE},\text{SKE})", color = wcolors[3])

    for ax in [axa, axb, axc, axd]
        fill_between!(ax, l_offsets/1e5, lower_bound, upper_bound; color=:black, alpha=0.5)
        axislegend(ax, backgroundcolor=nothing,framevisible=false, labelsize=10, position = :cb, patchsize = (15, 3), patchlabelgap = 3,orientation = :horizontal)
    end
    rowgap!(fig.layout, 3)
    resize_to_layout!(fig)
    save(filesave * "CCalong_front_" * fileparam * "_$(iteration).pdf", fig; pt_per_unit = 1)

    fig = Figure(size=(640, 540))
    axis_kwargs = (titlealign = :left, xgridvisible = false, ygridvisible = false)
    axa = Axis3(fig[1, 1]; azimuth = -0.275 * pi, title=L"\text{(a)}", xlabel = L"-| \nabla \overline{b}|^2 \overline{\delta}/(M_0^4 f)", ylabel=L"\mathcal{F}_s/(M_0^4 f)", zlabel=L"\sigma_n/f", axis_kwargs...)
    axb = Axis3(fig[1, 2]; title=L"\text{(b)}", xlabel=L"\text{TKE}/w_*^2", ylabel=L"\mathcal{F}_s/(M_0^4 f)", axis_kwargs...)
    axc = Axis3(fig[2, 1]; title=L"\text{(c)}", xlabel=L"\text{SKE}/w_*^2", ylabel=L"\mathcal{F}_s/(M_0^4 f)", zlabel=L"\sigma_n/f", axis_kwargs...)
    axd = Axis3(fig[2, 2]; title=L"\text{(d)}", xlabel=L"\text{SKE}/w_*^2", ylabel=L"\text{TKE}/w_*^2", axis_kwargs...)
    isigδ = δ_cfront .< minimum(δ_cfront)/3
    scatter!(axa, -D_cfront[.!isigδ].*δ_cfront[.!isigδ]/(f*M²₀^2), B_cfront[.!isigδ]/(f*M²₀^2), σn_cfront[.!isigδ]/f;  rasterize = true,colormap = :diff, color = σn_cfront[.!isigδ]/f, markersize=3, alpha = 0.1)
    scatter!(axa, -D_cfront[isigδ].*δ_cfront[isigδ]/(f*M²₀^2), B_cfront[isigδ]/(f*M²₀^2), σn_cfront[isigδ]/f;  rasterize = true,colormap = :diff, color = σn_cfront[isigδ]/f, markersize=10, alpha = 0.5)
    # Get colors from the colormap based on c values
    colors_from(c) = get(colorschemes[:diff], c, extrema(c))
    color_w_α(c,d) = [RGBAf(col.r, col.g, col.b, d[i]) for (i, col) in enumerate(colors_from(c))]
    isigTKE = TKE_cfront .> maximum(TKE_cfront)/3
    rTKE = Float64.(TKE_cfront/maximum(TKE_cfront))
    # scatter!(axb, TKE_cfront[.!isigTKE], B_cfront[.!isigTKE]/(f*M²₀^2), σn_cfront[.!isigTKE]/f;  rasterize = true,colormap = :diff, color = σn_cfront[.!isigTKE]/f, markersize=3, alpha = 0.1)
    scatter!(axb, TKE_cfront, B_cfront/(f*M²₀^2), σn_cfront/f;  rasterize = true, color = color_w_α(σn_cfront/f, rTKE), markersize=10*rTKE)
    isigSKE = SKE_cfront .> maximum(SKE_cfront)/3
    rSKE = Float64.(SKE_cfront/maximum(SKE_cfront))
    # scatter!(axc, SKE_cfront[.!isigSKE], B_cfront[.!isigSKE]/(f*M²₀^2), σn_cfront[.!isigSKE]/f;  rasterize = true,colormap = :diff, color = σn_cfront[.!isigSKE]/f, markersize=3, alpha = 0.1)
    scatter!(axc, SKE_cfront, B_cfront/(f*M²₀^2), σn_cfront/f;  rasterize = true, color = color_w_α(σn_cfront/f, rSKE), markersize=10*rSKE)
    isigKE = isigTKE .& isigSKE
    # scatter!(axd, SKE_cfront[.!isigKE], TKE_cfront[.!isigKE], σn_cfront[.!isigKE]/f;  rasterize = true,colormap = :diff, color = σn_cfront[.!isigKE]/f, markersize=3, alpha = 0.1)
    scatter!(axd, SKE_cfront, TKE_cfront, σn_cfront/f; rasterize = true, color = color_w_α(σn_cfront/f, rSKE.*rTKE), markersize=10*rTKE.*rSKE)
    # st = scatter!(axb, TKE_cfront, SKE_cfront; rasterize = true,colormap = :diff, color = σn_cfront, markersize=10, alpha = 0.3)
    # Colorbar(fig[1,3],st)
    colgap!(fig.layout, 1, 0)
    rowgap!(fig.layout, 1, 3)
    # colgap!(fig.layout, 2, 0)
    hidezdecorations!(axb, ticks=false)
    hidezdecorations!(axd, ticks=false)
    resize_to_layout!(fig)
    save(filesave * "cfront_scatters_" * fileparam * "_$(iteration).pdf", fig; pt_per_unit = 1)

    return arclength/1e5, Sn, σn_cfront, h_cfront, cd_Cu[1], cd_Cu[2], ζ_cfront, δ_cfront, TKE_cfront, SKE_cfront
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
xU,yU,zU = nodes(Ub);
itpU = interpolate((xU,yU), interior(Ub, :, :, length(zU)), Gridded(Linear(Interpolations.Periodic())))
etpU = extrapolate(itpU, Interpolations.Periodic())
xV,yV,zV = nodes(Vb);
itpV = interpolate((xV,yV), interior(Vb, :, :, length(zV)), Gridded(Linear(Interpolations.Periodic())))
etpV = extrapolate(itpV, Interpolations.Periodic())
σn = compute!(Field(-(∂x(Ub)-∂y(Vb))/2));
xn,yn,zn = nodes(σn);
itpn = interpolate((xn,yn), interior(σn, :, :, length(zn)), Gridded(Linear(Interpolations.Periodic())))
etpn = extrapolate(itpn, Interpolations.Periodic())

Nresample = 10240
s3 = zeros(Float32, 3, Nresample)
S3 = zeros(Float32, 3, Nresample)
σ3 = zeros(Float32, 3, Nresample)
d3 = zeros(Float32, 3, Nresample)
h3 = zeros(Float32, 3, Nresample)
C3 = zeros(Float32, 3, Nresample)
ζ3 = zeros(Float32, 3, Nresample)
δ3 = zeros(Float32, 3, Nresample)
TKE3 = zeros(Float32, 3, Nresample)
SKE3 = zeros(Float32, 3, Nresample)

fileparam = "xband1sublevels"
for (i, iteration) in enumerate([37003])
    s3[i,:],S3[i,:],σ3[i,:],h3[i,:],d3[i,:],C3[i,:],ζ3[i,:],δ3[i,:],TKE3[i,:],SKE3[i,:] = plot_front_properties(filehead,fileparam,iteration,Nresample)
    set_value!(; Lx = 1e5)
end

# jldopen(filesave * "along_front_" * fileparam * "_cg3hm_29h30h31h.jld2", "w") do file
#     file["s"] = s3
#     file["S"] = S3 
#     file["σ"] = σ3
#     file["h"] = h3
#     file["d"] = d3
#     file["C"] = C3
#     file["ζ"] = ζ3
#     file["δ"] = δ3
#     file["TKE"] = TKE3
#     file["SKE"] = SKE3
# end

# filename = filesave * "along_front_" * fileparam * "_cg3hm_29h30h31h.jld2";
# file = jldopen(filename, "r")
# s3 = file["s"];
# S3 = file["S"];
# σ3 = file["σ"];
# TKE3 = file["TKE"];
# SKE3 = file["SKE"];
# ζ3 = file["ζ"];
# δ3 = file["δ"];
# close(file)