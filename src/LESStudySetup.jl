module LESStudySetup

export parameters
export idealized_setup, turbulence_generator_setup
export set_value!, set!

using Reexport
@reexport using Oceananigans
using Oceananigans.Units
using Printf
using Adapt

include("parameters.jl")
include("initial_conditions.jl")
include("background_field_forcing.jl")
include("model_setup.jl")
include("idealized_setup.jl")
include("Diagnostics/Diagnostics.jl")

using .Diagnostics

end