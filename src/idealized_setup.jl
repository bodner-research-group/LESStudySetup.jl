"""
    idealized_setup(arch; stop_time = 100days, hydrostatic_approximation = true)

Create and configure a simulation for an idealized LES setup.

## Arguments
- `arch`: The architecture to use for the simulation.
- `stop_time`: The duration of the simulation in days. Default is 100 days.
- `hydrostatic_approximation`: Whether to use the hydrostatic approximation. Default is `true`.

## Returns
A `Simulation` object configured for the idealized setup.

The function retrieves the problem constants and calculates the grid size based on the problem dimensions. 
It then creates a `RectilinearGrid` object with the specified architecture and grid parameters. 
The model type and specific keyword arguments are determined based on the hydrostatic approximation setting. 
The function sets up the boundary conditions for velocity and temperature fields, 
and creates a `Model` object with the specified grid, coriolis, buoyancy, boundary conditions, and keyword arguments.
The initial conditions for velocity and temperature fields are set, and the maximum velocity magnitude is calculated to determine the time step size. 
A `TimeStepWizard` object is created to control the time step size during the simulation. 
Finally, a `Simulation` object is created with the model, time step, and stop time, and the progress and wizard callbacks are added.
"""
function idealized_setup(arch; 
                         stop_time = 100days,
			             stop_iteration = Inf,
                         hydrostatic_approximation = false,
                         background_forcing = true) # by default we include the eddies as a background forcing 
    
    # Retrieving the problem constants
    Δh = parameters.Δh 
    Δz = parameters.Δz 
    Lx = parameters.Lx 
    Ly = parameters.Ly 
    Lz = parameters.Lz 
     α = parameters.α
     f = parameters.f  
    ρ₀ = parameters.ρ₀
    cₚ = parameters.cp
    τw = parameters.τw 
     θ = parameters.θ
     Q = parameters.Q

    # Calculating the grid-size
    Nx = ceil(Int, Lx / Δh)
    Ny = ceil(Int, Ly / Δh)
    Nz = ceil(Int, Lz / Δz)

    # Constructing the grid
    grid = RectilinearGrid(arch, 
                           size = (Nx, Ny, Nz), 
                           x = (0, Lx), 
                           y = (0, Ly), 
                           z = (-Lz, 0),
                           halo = (6, 6, 6))

    @info "Running on a grid with $Nx, $Ny, and $Nz cells"

    # ModelType can be either a `HydrostaticFreeSurfaceModel` or a `NonhydrostaticModel`
    ModelType = model_type(Val(hydrostatic_approximation))
    settings  = model_settings(ModelType, grid; background_forcing)

    coriolis = FPlane(; f)
    buoyancy = SeawaterBuoyancy(; equation_of_state = LinearEquationOfState(thermal_expansion = α), 
                                  constant_salinity = 35)
    
    # # Cooling in the middle of the domain and heating outside?
    # @inline Qtop(x, y, t, p) = - p.Q / p.ρ₀ / p.cₚ * cos(2π * x / p.Lx)

    u_top = FluxBoundaryCondition(τw * cosd(θ) / ρ₀)
    v_top = FluxBoundaryCondition(τw * sind(θ) / ρ₀)
    T_top = FluxBoundaryCondition(Q / ρ₀ / cₚ) # Positive fluxes at the top are cooling in Oceananigans

    u_bcs = FieldBoundaryConditions(top = u_top)
    v_bcs = FieldBoundaryConditions(top = v_top)
    T_bcs = FieldBoundaryConditions(top = T_top)

    boundary_conditions = (u = u_bcs, v = v_bcs, T = T_bcs)
    
    model = ModelType(; grid, 
                        coriolis,
                        buoyancy,
                        boundary_conditions,
                        settings...)

    if isforced(model)
        set!(model, v = vᶠ, T = Tᵢ) 
    else
        set!(model, u = uᵢ, v = vᵢᶠ, T = Tᵢ)
    end
    
    u, v, w = model.velocities
    
    um = maximum(abs, interior(u))
    vm = maximum(abs, interior(v))

    u_max = max(um, vm) 

    Δt = min(0.2 * Δh / u_max, 10)
    
    # HydrostaticFreeSurfaceModel uses quasi Adams-Bashforth-2 which requires
    # a very small timestep, the NonhydrostaticModel uses a Runge-Kutta-3 scheme,
    # which is heavier but more stable and can use larger timesteps.
    cfl = hydrostatic_approximation ? 0.25 : 0.75
    max_Δt = hydrostatic_approximation ? 1minutes : 3minutes

    wizard = TimeStepWizard(; cfl, max_change = 1.1, max_Δt)
    simulation = Simulation(model; Δt, stop_time, stop_iteration)

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

    return simulation
end

function default_experimental_setup!(; Δh=parameters.Δh, Δz=parameters.Δz)
    set_value!(; 
               # Grid,
              Δh = Δh,
              Δz = Δz,
              Lz = 252,
               # Forcing
               Q = 40.0,   # Cooling heat flux in W/m²
              τw = 0.1,    # Wind stress in N/m²
               θ = 30.0, # Wind stress angle in degrees (0 correspond to zonal wind stress)
              # Initial condition - frontal values
             M²₀ = 5e-7,
              T₀ = 20,
              m₀ = 60,
             Δmᶠ = 10,
             N²s = 5e-7,
             N²T = 2e-4,
             # Eddy forcing
             ΔTᵉ = 0.1,  # Eddy temperature difference
             Φ   = 0.01, # Barotropic eddy strength
             a   = 1.0,  # Eddy temperature magnitude
             Lf  = 0.9,  # Size of temperature front (large numbers correspond to steeper fronts)
             σ²  = 1.0)
    
    return nothing
end

function turbulence_generator_setup(arch; 
                                    stop_time = 10hours,
                                    background_forcing = false)

    # Retrieving the problem constants
    Δh = parameters.Δh 
    Δz = parameters.Δz 
    Lz = parameters.Lz 
     α = parameters.α
     f = parameters.f  
    ρ₀ = parameters.ρ₀
    cₚ = parameters.cp
    τw = parameters.τw 
     θ = parameters.θ
     
    # Reduced domain size (250 by 250 meters)
    Lx = Ly = 250

    # Remember to set the value!
    set_value!(; Lx, Ly)

    # Calculating the grid-size
    Nx = ceil(Int, Lx / Δh)
    Ny = ceil(Int, Ly / Δh)
    Nz = ceil(Int, Lz / Δz)

    # Constructing the grid
    grid = RectilinearGrid(arch, 
                           size = (Nx, Ny, Nz), 
                           x = (0, Lx), 
                           y = (0, Ly), 
                           z = (-Lz, 0),
                           halo = (6, 6, 6))

    @info "Running on a grid with $Nx, $Ny, and $Nz cells"

    coriolis = FPlane(; f)
    buoyancy = SeawaterBuoyancy(; equation_of_state = LinearEquationOfState(thermal_expansion = α), 
                                  constant_salinity = 35)
    
    # # Cooling in the middle of the domain and heating outside?
    # @inline Qtop(x, y, t, p) = - p.Q / p.ρ₀ / p.cₚ * cos(2π * x / p.Lx)

    u_top = FluxBoundaryCondition(τw * cosd(θ) / ρ₀)
    v_top = FluxBoundaryCondition(τw * sind(θ) / ρ₀)

    u_bcs = FieldBoundaryConditions(top = u_top)
    v_bcs = FieldBoundaryConditions(top = v_top)

    # We force only velocity!
    boundary_conditions = (u = u_bcs, v = v_bcs)
    
    model = NonhydrostaticModel(; grid, 
                                  coriolis,
                                  buoyancy,
                                  boundary_conditions,
                                  advection = WENO(; order = 9),
                                  tracers = :T)

    # We initialize with a fictitious
    # vertical profile that only depends on z 
    set!(model, T = Tᶻ) 
     
    # 10 seconds as an initial step does 
    # not seem preposterous
    Δt = 10
    
    # But let's always add a wizard to be sure!
    wizard = TimeStepWizard(cfl = 0.25, max_change = 1.1)

    simulation = Simulation(model; Δt, stop_time)

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

    return simulation
end