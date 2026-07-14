#######################################################################
# hydrostatic_bg_advection_comparison.jl
#
# 5 days of HYDROSTATIC evolution at 150 m horizontal spacing
# (same vertical spacing as the main nonhydrostatic experiment),
# run twice: WITHOUT and WITH the background advection term u′⋅∇U.
#######################################################################

using LESStudySetup
using LESStudySetup.Oceananigans.Units
using LESStudySetup.Oceananigans.Utils: ConsecutiveIterations
using LESStudySetup.Oceananigans.OutputWriters: Checkpointer, JLD2Writer
using JLD2

# ── Architecture ─────────────────────────────────────────────────────
# Single GPU by default. For a multi-GPU run, prepend the MPI/GTL preamble
# used in `nonhydrostatic_experiment.jl` and replace `arch` with, e.g.
#     using MPI; MPI.Init()
#     arch = Distributed(GPU(), partition = Partition(proc, proc))
# making sure `proc` divides Nx = ceil(Lx / Δh) (= 667 here, i.e. 23×29).
arch = GPU()

# ── Resolution ───────────────────────────────────────────────────────
Δh = 150.0    # horizontal spacing [m]
Δz = 1.125    # vertical spacing [m] — same as the main nonhydrostatic experiment

LESStudySetup.default_experimental_setup!(; Δh, Δz)
@info "Problem parameters" parameters

# ── Time stepping / output ───────────────────────────────────────────
stop_time            = 5days
output_frequency     = 30minutes
checkpoint_frequency = 6hours

background_forcing = true   # include the eddy field as a background velocity

# ── Run both configurations ──────────────────────────────────────────
for advect_background in (false, true)
    case = advect_background ? "with_bg_advection" : "without_bg_advection"
    @info "Hydrostatic run: $case (advect_background = $advect_background)"

    simulation = idealized_setup(arch; stop_time,
                                       hydrostatic_approximation = true,
                                       background_forcing,
                                       advect_background)

    jldsave("hydrostatic_$(case)_metadata.jld2"; parameters)

    @info simulation

    model         = simulation.model
    output_fields = merge(model.velocities, model.tracers)

    simulation.output_writers[:snapshots] = JLD2Writer(model, output_fields;
        schedule           = ConsecutiveIterations(TimeInterval(output_frequency)),
        filename           = "hydrostatic_$(case)_snapshots",
        array_type         = Array{Float32},
        overwrite_existing = true)

    simulation.output_writers[:checkpoint] = Checkpointer(model;
        schedule           = TimeInterval(checkpoint_frequency),
        prefix             = "hydrostatic_$(case)_checkpoint",
        overwrite_existing = true)

    run!(simulation)
end
