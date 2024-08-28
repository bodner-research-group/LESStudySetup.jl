module Diagnostics

export write_pointwise_diagnostics
export load_snapshots, propagate_function,
       ζ, ub, vb, wb, uw, vw, KE, MLD, PV, BLD1D

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

function load_snapshots(filename; 
                        architecture = CPU(),
                        metadata = nothing)

    snapshots = Dict()

    u = FieldTimeSeries(filename, "u"; architecture, backend = OnDisk())
    v = FieldTimeSeries(filename, "v"; architecture, backend = OnDisk())
    w = FieldTimeSeries(filename, "w"; architecture, backend = OnDisk())
    T = FieldTimeSeries(filename, "T"; architecture, backend = OnDisk())

    snapshots[:u] = u
    snapshots[:v] = v
    snapshots[:w] = w
    snapshots[:T] = T

    if !isnothing(metadata)
        params = jldopen(metadata)["parameters"]
        set_value!(params)
    end
    
    return snapshots
end

times(snapshots::Dict) = snapshots[first(keys(snapshots))].times

include("mixed_layer.jl")
include("boundary_layer.jl")
include("pointwise_diagnostics.jl")
include("spectra.jl")
include("filtering.jl")

end