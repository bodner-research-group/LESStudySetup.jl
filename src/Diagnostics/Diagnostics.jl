module Diagnostics

export write_pointwise_diagnostics
export load_snapshots, propagate_function,
       ζ, ub, vb, wb, uw, vw, KE, MLD, BLD1D, PV

using Oceananigans
using Oceananigans

using Oceananigans
using Oceananigans.BuoyancyModels: buoyancy
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid
using Oceananigans.Architectures: device, architecture
using Oceananigans.Utils: launch!
using Oceananigans.Grids: Center, Face, inactive_node, znode
using Oceananigans.Operators: Δzᶜᶜᶠ, ζ₃ᶠᶠᶜ

using KernelAbstractions: @index, @kernel
using KernelAbstractions.Extras.LoopInfo: @unroll
using JLD2

import Oceananigans.Fields: compute!

using Oceananigans.Fields: OneField, condition_operand
using Oceananigans.AbstractOperations: materialize_condition!
using Oceananigans.Utils

using LESStudySetup

function write_pointwise_diagnostics(file_prefix; architecture = CPU())

    # Path to the snapshots
    filename = "hydrostatic_snapshots" * file_prefix * ".jld2"
    metadata = "experiment" * file_prefix * "_metadata.jld2"

    # Load in the snapshots
    snapshots = load_snapshots(filename; metadata, architecture)

    # Make sure parameters are correctly loaded
    @info parameters

    output_filename = "diagnostic" * file_prefix * ".jld2"

    UB = propagate_function(ub,  snapshots; filename = output_filename)
    VB = propagate_function(vb,  snapshots; filename = output_filename)
    WB = propagate_function(wb,  snapshots; filename = output_filename)
    UW = propagate_function(uw,  snapshots; filename = output_filename)
    VW = propagate_function(vw,  snapshots; filename = output_filename)
    Z  = propagate_function(ζ,   snapshots; filename = output_filename)
    D  = propagate_function(δ,   snapshots; filename = output_filename)
    Q  = propagate_function(PV,  snapshots; filename = output_filename)
    MX = propagate_function(MLD, snapshots; filename = output_filename)
    BD = propagate_function(BLD1D, snapshots; filename = output_filename)

    return (; UB, VB, WB, UW, VW, Z, D, Q, MX, BD)
end

const F = Face
const C = Center


function rewrite_variable_with_halos(filename, new_filename, variable_name; architecture = CPU())
    old_fts = FieldTimeSeries(filename, variable_name; architecture, backend = OnDisk())
    grid    = old_fts.grid
    times   = old_fts.times
    loc     = location(old_fts)
    new_fts = FieldTimeSeries{loc...}(grid, times; backend = OnDisk(), name = variable_name, path = new_filename)
    tmp     = Field{loc...}(grid)

    for t in eachindex(times)
        set!(tmp, old_fts[t])
        fill_halo_regions!(tmp)
        set!(new_fts, tmp, t)
    end

    return new_fts
end

function load_snapshots(filename; 
                        add_halos = false,
                        new_filename = nothing,
                        architecture = CPU(),
                        metadata = nothing,
                        variables = (:u, :v, :w, :T))

    snapshots = Dict()

    if add_halos
        for var in variables
            snapshots[var] = rewrite_variable_with_halos(filename, new_filename, string(var); architecture)
        end
    else
        for var in variables
            snapshots[var] = FieldTimeSeries(filename, string(var); architecture, backend = OnDisk())
        end
    end

    if !isnothing(metadata)
        params = jldopen(metadata)["parameters"]
        set_value!(params)
    end
    
    return snapshots
end

times(snapshots::Dict) = snapshots[first(keys(snapshots))].times

include("load_distributed_snapshot.jl")
include("mixed_layer.jl")
include("boundary_layer.jl")
include("pointwise_diagnostics.jl")
include("spectra.jl")
include("filtering.jl")

end