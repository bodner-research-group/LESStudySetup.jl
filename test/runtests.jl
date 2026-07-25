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

    @testset "Distributed subdomain extraction and save/load round-trip" begin
        # Use MPI.jl's mpiexec wrapper which is compatible with JLL libraries
        using MPI
        
        # Define the distributed test script as a string
        distributed_test_script = """
        using MPI
        MPI.Init()

        using Test
        using JLD2
        using LESStudySetup
        using LESStudySetup.Diagnostics: load_distributed_checkpoint,
                                          load_distributed_checkpoint_subdomain,
                                          save_subdomain_snapshot,
                                          load_subdomain_snapshot
        using Oceananigans
        using Oceananigans.OutputWriters: Checkpointer

        # Get temp directory from environment (passed from parent test)
        test_dir = ENV["TEST_TMPDIR"]

        # 1. Configure test parameters
        set_value!(;
            Δh = 1000.0,     # 1km horizontal spacing
            Δz = 10.0,       # 10m vertical spacing
            Lx = 8000.0,     # 8km domain (divisible by 2 ranks)
            Ly = 8000.0,
            Lz = 100.0       # 100m depth
        )

        # 2. Create distributed architecture (2x2 = 4 ranks)
        arch = Distributed(CPU(), partition = Partition(2, 2))

        # 3. Create simulation
        simulation = idealized_setup(arch;
                                     stop_time = 1,
                                     hydrostatic_approximation = false)
        model = simulation.model

        # 4. Save checkpoint after 1 iteration
        checkpoint_prefix = joinpath(test_dir, "test_checkpoint_\$(arch.local_rank)")
        simulation.output_writers[:checkpoint] = Checkpointer(model;
            schedule = IterationInterval(1),
            prefix = checkpoint_prefix,
            overwrite_existing = true)

        # Run for 1 iteration to trigger checkpoint
        run!(simulation)

        # 5. Synchronize all ranks before loading
        MPI.Barrier(MPI.COMM_WORLD)

        # 6. Only rank 0 performs verification
        if arch.local_rank == 0
            checkpoint_path = joinpath(test_dir, "test_checkpoint_")

            # Load full domain checkpoint at iteration 1
            full = load_distributed_checkpoint(checkpoint_path, 1)

            # Verify full domain dimensions
            @test size(interior(full[:T])) == (8, 8, 10)

            # Define subdomain limits that cross rank boundaries
            # With 8x8 grid, 2x2 partition, each rank has 4x4 cells
            # Rank boundaries: x=0,4000,8000 and y=0,4000,8000
            # This subdomain spans parts of all 4 ranks
            xlims = (2000.0, 6000.0)   # Crosses x rank boundary at 4000
            ylims = (2000.0, 6000.0)   # Crosses y rank boundary at 4000
            zlims = (-80.0, -30.0)

            # Load subdomain
            subdomain = load_distributed_checkpoint_subdomain(checkpoint_path, 1;
                xlims = xlims,
                ylims = ylims,
                zlims = zlims)

            # Verify subdomain dimensions
            @test subdomain[:grid].Nx == 4   # 4km / 1km = 4 cells
            @test subdomain[:grid].Ny == 4
            @test subdomain[:grid].Nz == 5   # 50m / 10m = 5 cells

            # Verify subdomain data matches expected slice from full domain
            # Index mapping: xlims (2km, 6km) -> indices 3:6
            #                ylims (2km, 6km) -> indices 3:6
            #                zlims (-80m, -30m) -> indices 3:7
            @test interior(subdomain[:T]) ≈ interior(full[:T])[3:6, 3:6, 3:7]

            # Load new subdomain with larger zlims (including original range)
            new_zlims = (-80.0, 0.0)  # 80m depth
            new_subdomain = load_distributed_checkpoint_subdomain(checkpoint_path, 1;
                xlims = xlims,
                ylims = ylims,
                zlims = new_zlims)

            # Verify new subdomain dimensions
            @test new_subdomain[:grid].Nx == 4
            @test new_subdomain[:grid].Ny == 4
            @test new_subdomain[:grid].Nz == 8   # 80m / 10m = 8 cells

            # Verify old subdomain matches slice of new subdomain
            # Old zlims (-80, -30) -> first 5 z-levels of new subdomain
            @test interior(subdomain[:T]) ≈ interior(new_subdomain[:T])[:, :, 1:5]

            # Save the new subdomain to file
            subdomain_file = joinpath(test_dir, "test_subdomain.jld2")
            save_subdomain_snapshot(subdomain_file, new_subdomain;
                iteration = 1,
                xlims = xlims,
                ylims = ylims,
                zlims = new_zlims)

            # Load the subdomain back
            loaded = load_subdomain_snapshot(subdomain_file)

            # Verify loaded data matches new subdomain
            @test loaded[:grid].Nx == new_subdomain[:grid].Nx
            @test loaded[:grid].Ny == new_subdomain[:grid].Ny
            @test loaded[:grid].Nz == new_subdomain[:grid].Nz
            @test interior(loaded[:T]) ≈ interior(new_subdomain[:T])
            @test interior(loaded[:u]) ≈ interior(new_subdomain[:u])
            @test interior(loaded[:v]) ≈ interior(new_subdomain[:v])
            @test interior(loaded[:w]) ≈ interior(new_subdomain[:w])

            @info "All distributed subdomain tests passed!"
        end

        # Final barrier before cleanup
        MPI.Barrier(MPI.COMM_WORLD)
        """

        # Create temporary directory for test files
        test_dir = mktempdir()

        # Write script to temp file
        script_file = joinpath(test_dir, "distributed_subdomain_test.jl")
        write(script_file, distributed_test_script)

        try
            # Set environment variable to pass temp dir path to MPI script
            withenv("TEST_TMPDIR" => test_dir) do
                # Execute with MPI.jl's mpiexec wrapper (compatible with JLL libraries)
                mpiexec() do mpiexec_cmd
                    run(`$mpiexec_cmd -n 4 julia --project -O0 $script_file`)
                end
            end
        finally
            # Cleanup
            rm(test_dir, recursive=true, force=true)
        end
    end

    @testset "Coarse-graining on subdomain grid" begin
        using LESStudySetup.Diagnostics: coarse_graining!
        using Statistics: var, mean
        using Oceananigans: fill_halo_regions!
        
        # Set up parameters for a 100km domain with 1km spacing
        set_value!(;
            Δh = 1000.0,
            Δz = 10.0,
            Lx = 100000.0,
            Ly = 100000.0,
            Lz = 100.0
        )
        
        # Create a SUBDOMAIN grid (10km x 10km) - different from parameters
        subdomain_grid = RectilinearGrid(CPU();
            size = (10, 10, 10),
            x = (0.0, 10000.0),
            y = (0.0, 10000.0),
            z = (-100.0, 0.0),
            topology = (Bounded, Bounded, Bounded)
        )
        
        # Verify subdomain grid differs from parameters
        @test subdomain_grid.Lx == 10000.0  # 10km, not 100km
        @test subdomain_grid.Lx != parameters.Lx
        
        # Create test fields on subdomain grid
        T_field = CenterField(subdomain_grid)
        T_filtered = CenterField(subdomain_grid, Float32)
        
        # Initialize with a pattern (warm blob in center)
        set!(T_field, (x, y, z) -> 10.0 + 5.0 * exp(-((x-5000)^2 + (y-5000)^2) / 2000^2))
        fill_halo_regions!(T_field)
        
        # Apply coarse-graining with 2km cutoff (should work on 10km subdomain)
        # This should NOT error despite parameters.Lx = 100km
        coarse_graining!(T_field, T_filtered; 
            kernel = :gaussian, 
            cutoff = 2000.0,  # 2km
            border = :reflect
        )
        
        # Verify output is valid (not NaN, smoothed)
        T_interior = interior(T_field)
        T_filt_interior = interior(T_filtered)
        
        @test !any(isnan, T_filt_interior)
        @test !any(isinf, T_filt_interior)
        
        # Filtered field should be smoother (lower variance)
        @test var(T_filt_interior) < var(T_interior)
        
        # Mean should be approximately preserved
        @test isapprox(mean(T_filt_interior), mean(T_interior), rtol=0.1)
        
        @info "Coarse-graining on subdomain test passed!"
    end

    @testset "Discrete levels save/load round-trip with save_subdomain_with_halo" begin
        using LESStudySetup.Diagnostics: save_subdomain_with_halo
        using Oceananigans: fill_halo_regions!
        
        # 1. Configure test parameters
        set_value!(;
            Δh = 1000.0,     # 1km horizontal spacing
            Δz = 10.0,       # 10m vertical spacing
            Lx = 10000.0,    # 10km domain
            Ly = 10000.0,
            Lz = 100.0       # 100m depth (10 z-levels)
        )
        
        # 2. Define discrete levels metadata (as would be used in quadrant_analysis_publication.jl)
        z_indices = [2, 5, 8]
        Nz_compact = length(z_indices)
        
        # 3. Create a "compact" grid with only these 3 levels
        # This simulates what create_compact_snapshot_from_levels() produces
        Δz_val = parameters.Δz
        Lz_val = parameters.Lz
        z_centers = [-Lz_val + (k - 1) * Δz_val + Δz_val/2 for k in z_indices]
        z_min = minimum(z_centers) - Δz_val/2
        z_max = maximum(z_centers) + Δz_val/2
        
        compact_grid = RectilinearGrid(CPU();
            size = (10, 10, Nz_compact),
            x = (0.0, 10000.0),
            y = (0.0, 10000.0),
            z = (z_min, z_max),
            topology = (Bounded, Bounded, Bounded))
        
        # 4. Create fields on compact grid with test data
        # This is what the compacted snapshot would contain after extraction
        u_compact = XFaceField(compact_grid)
        v_compact = YFaceField(compact_grid)
        w_compact = ZFaceField(compact_grid)
        T_compact = CenterField(compact_grid)
        
        # Fill with simple constant values for fast testing
        interior(u_compact) .= 0.1
        interior(v_compact) .= 0.2
        interior(w_compact) .= 0.01
        interior(T_compact) .= 10.5
        
        fill_halo_regions!(u_compact)
        fill_halo_regions!(v_compact)
        fill_halo_regions!(w_compact)
        fill_halo_regions!(T_compact)
        
        compact_snapshot = Dict{Symbol, Any}(
            :grid => compact_grid,
            :u => u_compact,
            :v => v_compact,
            :w => w_compact,
            :T => T_compact
        )
        
        # 5. Save with save_subdomain_with_halo using `levels` (not `zlims`)
        test_file = tempname() * ".jld2"
        save_subdomain_with_halo(test_file, compact_snapshot;
            core_xlims = (1000.0, 9000.0),
            core_ylims = (1000.0, 9000.0),
            halo_width = 1000.0,
            levels = z_indices,
            iteration = 0)
        
        # 6. Verify file contains levels metadata, not zlims
        jldopen(test_file, "r") do file
            @test haskey(file, "metadata/levels")
            @test file["metadata/levels"] == z_indices
            @test !haskey(file, "metadata/zlims")
        end
        
        # 7. Load back and verify dimensions match
        loaded = load_subdomain_snapshot(test_file)
        
        @test loaded[:grid].Nx == compact_grid.Nx
        @test loaded[:grid].Ny == compact_grid.Ny
        @test loaded[:grid].Nz == Nz_compact  # Critical: Nz should be 3, not 10
        
        @test interior(loaded[:T]) ≈ interior(T_compact)
        @test interior(loaded[:u]) ≈ interior(u_compact)
        @test interior(loaded[:v]) ≈ interior(v_compact)
        @test interior(loaded[:w]) ≈ interior(w_compact)
        
        # 8. Cleanup
        rm(test_file, force=true)
        
        @info "Discrete levels save/load round-trip test passed!"
    end
end
