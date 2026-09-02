#######################################################################
# hydrostatic_background_velocity_comparison.jl
#
# 150 days of HYDROSTATIC evolution on a coarse grid (Δh = 500 m, i.e.
# twice as fine as the 1 km default horizontal spacing), run three times.
# All three carry the eddy as a background velocity; they differ in how the
# baroclinic part of the eddy is represented:
#
#   with_background_forcing        full eddy velocity uᵢ in the background, eddy
#                                  temperature anomaly in the initial condition (Tᵢ)
#   without_baroclinic_temperature eddy temperature anomaly dropped from the initial
#                                  condition (Tᶠᶻ, the front alone) and the thermal
#                                  wind it would have supported carried in the
#                                  background twice over, u²ᶜ = uᵢ + (uᵢ - uᴮ)
#   with_barotropic_forcing        only the barotropic eddy velocity uᴮ in the
#                                  background, eddy temperature anomaly in the
#                                  initial condition (Tᵢ)
#
# The prognostic velocity is initialized with the front jet vᶠ in all three cases.
# The eddy never enters the initial conditions.
#######################################################################

using LESStudySetup
using LESStudySetup.Oceananigans.Units
using LESStudySetup.Oceananigans.OutputWriters: Checkpointer, JLD2Writer
using JLD2

# ── Architecture ─────────────────────────────────────────────────────
# Single GPU. For a multi-GPU run, prepend the MPI/GTL preamble used in
# `nonhydrostatic_experiment.jl` and replace `arch` with, e.g.
#     using MPI; MPI.Init()
#     arch = Distributed(GPU(), partition = Partition(proc, proc))
# making sure `proc` divides Nx = Ny = 200.
arch = GPU()

# ── Resolution ───────────────────────────────────────────────────────
Δh = 500.0    # horizontal spacing [m] — 2× finer than the 1 km default
Δz = 4.0      # vertical spacing [m] — the default

LESStudySetup.default_experimental_setup!(; Δh, Δz)
@info "Problem parameters" parameters

# ── Time stepping / output ───────────────────────────────────────────
stop_time            = 150days
output_frequency     = 3hours
average_frequency    = 1day
checkpoint_frequency = 5days

pickup = false   # set to `true` to restart each case from its last checkpoint

# ── Run the three configurations ─────────────────────────────────────
cases = (; with_background_forcing        = (; background_velocity = LESStudySetup.uᵢ,
                                              initial_temperature = LESStudySetup.Tᵢ),
           without_baroclinic_temperature = (; background_velocity = LESStudySetup.u²ᶜ,
                                              initial_temperature = LESStudySetup.Tᶠᶻ),
           with_barotropic_forcing        = (; background_velocity = LESStudySetup.uᴮ,
                                              initial_temperature = LESStudySetup.Tᵢ))

for (name, setup) in pairs(cases)
    case = string(name)
    @info "Hydrostatic run: $case" setup

    simulation = idealized_setup(arch; stop_time,
                                       hydrostatic_approximation = true,
                                       background_forcing = true,
                                       setup...)

    jldsave("hydrostatic_$(case)_metadata.jld2"; parameters)

    @info simulation

    model = simulation.model

    # Diffusivities are written only in the daily averages: at 3-hourly cadence over
    # 150 days they would dominate the output volume.
    snapshot_fields = merge(model.velocities, model.tracers)
    average_fields  = merge(snapshot_fields,
                            (; κu = model.closure_fields.κu,
                               κc = model.closure_fields.κc,
                               κe = model.closure_fields.κe))

    simulation.output_writers[:snapshots] = JLD2Writer(model, snapshot_fields;
        schedule           = TimeInterval(output_frequency),
        filename           = "hydrostatic_$(case)_snapshots",
        array_type         = Array{Float32},
        with_halos         = true,
        overwrite_existing = true)

    simulation.output_writers[:dailyaverages] = JLD2Writer(model, average_fields;
        schedule           = AveragedTimeInterval(average_frequency),
        filename           = "hydrostatic_$(case)_1Daverages",
        array_type         = Array{Float32},
        with_halos         = true,
        overwrite_existing = true)

    simulation.output_writers[:free_surface] = JLD2Writer(model, (; η = model.free_surface.displacement);
        schedule           = TimeInterval(output_frequency),
        filename           = "hydrostatic_$(case)_free_surface",
        array_type         = Array{Float32},
        with_halos         = true,
        overwrite_existing = true)

    simulation.output_writers[:checkpoint] = Checkpointer(model;
        schedule           = TimeInterval(checkpoint_frequency),
        prefix             = "hydrostatic_$(case)_checkpoint",
        overwrite_existing = true,
        cleanup            = true)

    run!(simulation; pickup)
end
