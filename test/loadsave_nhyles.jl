using Oceananigans
using Oceananigans: location
using JLD2
using LESStudySetup
using LESStudySetup.Diagnostics
using LESStudySetup.Diagnostics: load_distributed_checkpoint_subdomain

# --- Simulation and Subdomain Parameters ---
filehead = "/orcd/data/abodner/002/nhyles_output/"
filename = filehead * "iteration5x/nonhydrostatic_checkpoint_"
iterations = [32207,32993,33772,34542]#72635 #52543#32207 #21955 at 11hr

# Define the desired subdomain in global coordinates.
# To create a 10km x 10km domain centered at (0,0) from the periodic
# full domain, we specify these limits. The function will correctly
# fetch data from the end of the full domain (e.g., 95km to 100km)
# for the negative coordinates.
xlims = (-12500, 12500)#(-3125, 3125)
ylims = (0.0, 1e5)#(5e4-3125, 53125)
zlims = nothing #(-252, 0.0) # nz vertical cells * 1.125m/cell
levels = [221, 224] #[162, 198, 221, 224]

# Load the subdomain data using the revised function
for iteration in iterations
    @info "Loading subdomain at iteration $iteration..."

    snapshot = load_distributed_checkpoint_subdomain(filename, iteration; 
                                                    xlims = xlims,
                                                    ylims = ylims, 
                                                    #zlims = zlims,
                                                    levels = levels,
                                                    getEw = true,
                                                    getMLD = 3)

    @show snapshot[:T].grid

    # The 'snapshot' Dict now contains the fields (u, v, w, T) and the grid
    # for your specified subdomain, running on the CPU.
    #u_sub = snapshot[:u]
    #v_sub = snapshot[:v]
    #w_sub = snapshot[:w]
    #T_sub = snapshot[:T]
    sub_grid = snapshot[:grid]

    @show sub_grid

    # --- Save the Subdomain to a New JLD2 File ---

    output_filename = filehead * "subdomains/xband1sublevels_snapshot_iter$(iteration).jld2"

    @info "Saving subdomain to $output_filename..."

    # We'll save the fields and grid in a format that's easy to reload.
    # It's good practice to save the raw data (interiors) and grid separately.
    jldopen(output_filename, "w") do file
        file["grid"] = sub_grid
        
        # Save each field's interior data and location metadata
        for field_name in keys(snapshot)
            if field_name != :grid
                field = snapshot[field_name]
                # Convert to standard Array on the CPU before saving
                field_data = Array(interior(field)) 
                
                file["fields/$field_name/data"] = field_data
                file["fields/$field_name/location"] = location(field)
                #file["fields/$field_name/metadata/grid"] = sub_grid # For context
            end
        end
        
        file["metadata/iteration"] = iteration
        file["metadata/xlims"] = xlims
        file["metadata/ylims"] = ylims
        if !isnothing(zlims)
            file["metadata/zlims"] = zlims
        elseif !isnothing(levels)
            file["metadata/levels"] = levels
        end
    end

    @info "Successfully saved subdomain data."
end