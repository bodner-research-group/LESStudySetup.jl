using LESStudySetup
using LESStudySetup.Diagnostics
using LESStudySetup.Diagnostics: extract_subdomain, save_subdomain_snapshot, load_subdomain_snapshot
using Oceananigans
using Oceananigans: location
using Test
using JLD2

@testset "LESStudySetup.jl" begin
    @testset "Subdomain extraction and save/load round-trip" begin
        # 1. Configure test parameters using set_value!
        #    Create a 10km x 10km x 100m domain with 1km horizontal, 10m vertical spacing
        #    This gives a 10x10x10 grid - small but testable
        set_value!(;
            Δh = 1000.0,     # 1km horizontal spacing
            Δz = 10.0,       # 10m vertical spacing
            Lx = 10000.0,    # 10km domain
            Ly = 10000.0,
            Lz = 100.0       # 100m depth
        )

        # 2. Create simulation using idealized_setup
        simulation = idealized_setup(CPU();
                                     stop_time = 1,
                                     hydrostatic_approximation = false)
        model = simulation.model
        grid = model.grid

        # 3. Build full-domain snapshot Dict from model fields
        full_snapshot = Dict{Symbol, Any}()
        full_snapshot[:u] = model.velocities.u
        full_snapshot[:v] = model.velocities.v
        full_snapshot[:w] = model.velocities.w
        full_snapshot[:T] = model.tracers.T
        full_snapshot[:grid] = grid

        # 4. Extract a SUBDOMAIN (e.g., middle 5km x 5km x 50m portion)
        #    Full domain: x ∈ [0, 10km], y ∈ [0, 10km], z ∈ [-100m, 0]
        #    Subdomain:   x ∈ [2km, 7km], y ∈ [2km, 7km], z ∈ [-80m, -30m]
        sub_xlims = (2000.0, 7000.0)
        sub_ylims = (2000.0, 7000.0)
        sub_zlims = (-80.0, -30.0)

        subdomain_snapshot = extract_subdomain(full_snapshot;
                                               xlims = sub_xlims,
                                               ylims = sub_ylims,
                                               zlims = sub_zlims)

        # 5. Verify subdomain dimensions
        sub_grid = subdomain_snapshot[:grid]
        @test sub_grid.Nx == 5   # 5km / 1km spacing = 5 cells
        @test sub_grid.Ny == 5
        @test sub_grid.Nz == 5   # 50m / 10m spacing = 5 cells

        # 6. Verify subdomain data matches expected slice from full domain
        #    Index ranges in full domain (1-indexed):
        #    x: indices 3:7 (cells centered at 2.5km, 3.5km, 4.5km, 5.5km, 6.5km)
        #    y: indices 3:7
        #    z: indices 3:7 (cells centered at -85m, -75m, -65m, -55m, -45m)
        full_T = interior(model.tracers.T)
        sub_T = interior(subdomain_snapshot[:T])
        @test sub_T ≈ full_T[3:7, 3:7, 3:7]

        # Verify u, v fields as well (x-face and y-face centered)
        # Note: u-field has extra x-point (Bounded topology), v-field has extra y-point
        full_u = interior(model.velocities.u)
        sub_u = interior(subdomain_snapshot[:u])
        @test sub_u ≈ full_u[3:8, 3:7, 3:7]  # x: 3:8 for XFaceField (6 points)

        full_v = interior(model.velocities.v)
        sub_v = interior(subdomain_snapshot[:v])
        @test sub_v ≈ full_v[3:7, 3:8, 3:7]  # y: 3:8 for YFaceField (6 points)

        # 7. Save subdomain to file
        test_file = tempname() * ".jld2"
        save_subdomain_snapshot(test_file, subdomain_snapshot;
                                iteration = 0,
                                xlims = sub_xlims,
                                ylims = sub_ylims,
                                zlims = sub_zlims)

        # 8. Load subdomain back
        loaded = load_subdomain_snapshot(test_file)

        # 9. Verify loaded data matches subdomain
        @test loaded[:grid].Nx == sub_grid.Nx
        @test loaded[:grid].Ny == sub_grid.Ny
        @test loaded[:grid].Nz == sub_grid.Nz
        @test interior(loaded[:T]) ≈ sub_T
        @test interior(loaded[:u]) ≈ interior(subdomain_snapshot[:u])
        @test interior(loaded[:v]) ≈ interior(subdomain_snapshot[:v])
        @test interior(loaded[:w]) ≈ interior(subdomain_snapshot[:w])

        # 10. Cleanup
        rm(test_file, force=true)
    end
end
