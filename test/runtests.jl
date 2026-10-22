using Test
using LinearAlgebra
using TOML
using HubbardRenyiDQMC

include("square_lattice.jl")
include("test_ed.jl")
if get(ENV, "HUBBARD_TEST_DMRG", "false") == "true"
    include("test_dmrg_package.jl")
end

module SummarizeGridScript
include(joinpath(@__DIR__, "..", "analysis", "summarize_grid.jl"))
end

module GridLauncherScript
include(joinpath(@__DIR__, "..", "run", "run_grid_local.jl"))
end

module NodeRunnerScript
include(joinpath(@__DIR__, "..", "run", "run_lambda_node.jl"))
end

module EDReferenceScript
include(joinpath(@__DIR__, "..", "analysis", "compute_s2_ed.jl"))
end

module U4L4ThermalEDAuditScript
include(joinpath(@__DIR__, "..", "analysis", "audit_u4_l4_thermal_ed.jl"))
end

module DMRGReferenceScript
include(joinpath(@__DIR__, "..", "analysis", "compute_s2_itensor_dmrg.jl"))
end

module U4FSSAnalysisScript
include(joinpath(@__DIR__, "..", "analysis", "analyze_u4_fss.jl"))
end

const TEST_METADATA = (
    Ninc = 1,
    S2_ED = nothing,
    root_seeds = [7],
    N_bins = 2,
    run_name = "fixture",
)

function write_bins(path, rows; header = SummarizeGridScript.BINS_HEADER)
    open(path, "w") do io
        println(io, header)
        for row in rows
            println(io, row)
        end
    end
end

const K0_ROWS = [
    "fixture,0,0.0,7,1,1.0,1.0,0.1,0.80,",
    "fixture,0,0.0,7,2,1.0,1.0,0.1,0.82,",
]
const K1_ROWS = [
    "fixture,1,1.0,7,1,1.0,2.0,1.0,,1.20",
    "fixture,1,1.0,7,2,1.0,2.0,1.0,,1.22",
]
const K0_SUMMARY_ROW = join(
    [
        "fixture", "2", "1.0", "2.0", "0.0", "4.0", "0.1", "2", "1",
        "1", "1", "1.0", "0", "0.0", "7", "11", "12", "plus", "0.0",
        "1.0", "1.0", "0.1", "0.0", "0.1", "1.0", "0.1", "1.0", "0.1",
        "0.5", "0.0", "1.0",
    ],
    ",",
)

function write_valid_grid(dir)
    write_bins(joinpath(dir, "fixture_k000_seed7_bins.csv"), K0_ROWS)
    write_bins(joinpath(dir, "fixture_k001_seed7_bins.csv"), K1_ROWS)
end

function write_parameter_fixture(path; boundary_condition = nothing)
    open(path, "w") do io
        println(io, "[model]")
        println(io, "L = 8")
        println(io, "t = 1.0")
        println(io, "U = 2.0")
        println(io, "mu = 0.0")
        println(io, "beta = 32.0")
        println(io, "dtau = 0.05")
        println(io, "n_stab = 10")
        println(io, "N_sweeps = 1")
        boundary_condition === nothing ||
            println(io, "boundary_condition = \"$boundary_condition\"")
        println(io, "\n[subsystem]")
        println(io, "nA = 4")
        println(io, "\n[lambda_grid]")
        println(io, "Ninc = 1")
        println(io, "S2_ED = nan")
        println(io, "\n[sampling]")
        println(io, "N_burnin = 1")
        println(io, "N_measurements = 2")
        println(io, "N_bins = 2")
        println(io, "root_seeds = [7]")
        println(io, "\n[output]")
        println(io, "run_name = \"fixture\"")
        println(io, "raw_dir = \"raw\"")
        println(io, "processed_dir = \"processed\"")
        println(io, "figures_dir = \"figures\"")
    end
end

function write_u4_comparison_fixture(path; overrides = Dict{Tuple{Int,String},String}())
    header = [
        "L", "nA", "beta", "boundary", "U", "DQMC_S2", "DQMC_error",
        "reference_method", "reference_S2", "DQMC_minus_reference", "z_DQMC_only",
    ]
    open(path, "w") do io
        println(io, join(header, ','))
        for L in 4:2:22
            reference = 1.1 + 0.18 * log(L / pi) + 0.04 * (-1.0)^(L ÷ 2) / L
            dqmc = reference + 0.002 * sin(L)
            error = 0.01 + 0.0002 * L
            values = Dict(
                "L" => string(L),
                "nA" => string(L ÷ 2),
                "beta" => string(4L),
                "boundary" => "pbc",
                "U" => "4",
                "DQMC_S2" => string(dqmc),
                "DQMC_error" => string(error),
                "reference_method" => (L <= 12 ? "ED" : "DMRG"),
                "reference_S2" => string(reference),
                "DQMC_minus_reference" => string(dqmc - reference),
                "z_DQMC_only" => string((dqmc - reference) / error),
            )
            for column in header
                values[column] = get(overrides, (L, column), values[column])
            end
            println(io, join((values[column] for column in header), ','))
        end
    end
end

@testset "simulation-free workflow checks" begin
    @testset "PBC and APBC kinetic matrices" begin
        @test normalize_boundary_condition(:pbc) == :pbc
        @test normalize_boundary_condition("apbc") == :apbc
        @test_throws ArgumentError normalize_boundary_condition("periodic")
        @test_throws ArgumentError normalize_boundary_condition(:APBC)

        L = 8
        t = 1.0
        pbc = hubbard_chain_kinetic_matrix(t, L)
        apbc = hubbard_chain_kinetic_matrix(t, L, :apbc)

        expected_pbc = zeros(L, L)
        for site in 1:(L - 1)
            expected_pbc[site, site + 1] = -t
            expected_pbc[site + 1, site] = -t
        end
        expected_pbc[L, 1] = -t
        expected_pbc[1, L] = -t

        @test pbc == expected_pbc
        @test ishermitian(apbc)
        @test apbc[1, L] == t
        @test apbc[L, 1] == t
        for site in 1:(L - 1)
            @test apbc[site, site + 1] == pbc[site, site + 1]
            @test apbc[site + 1, site] == pbc[site + 1, site]
        end
        @test count(x -> !iszero(x), apbc - pbc) == 2

        sublattice = Diagonal([isodd(site) ? 1.0 : -1.0 for site in 1:L])
        @test sublattice * apbc * sublattice == -apbc

        expected_apbc_spectrum = sort([
            -2t * cos(2π * (n + 0.5) / L) for n in 0:(L - 1)
        ])
        @test sort(eigvals(Symmetric(apbc))) ≈ expected_apbc_spectrum atol = 1e-12

        one_particle_masks = EDReferenceScript.fixed_particle_masks(L, 1)
        one_particle_index = EDReferenceScript.mask_index(one_particle_masks)
        for boundary in (:pbc, :apbc)
            spin_hops = EDReferenceScript.build_spin_hops(
                one_particle_masks,
                one_particle_index,
                L,
                t;
                boundary = string(boundary),
            )
            ed_kinetic = zeros(L, L)
            for source in eachindex(spin_hops)
                for (destination, value) in spin_hops[source]
                    ed_kinetic[destination, source] += value
                end
            end
            @test ed_kinetic == hubbard_chain_kinetic_matrix(t, L, boundary)
        end
    end

    @testset "ED and DMRG reference metadata" begin
        ed_cfg = EDReferenceScript.parse_cli([
            "--L", "2", "--nA", "1", "--lanczos-maxiter", "17",
            "--lanczos-tol", "1e-9",
        ])
        @test ed_cfg.lanczos_maxiter == 17
        @test ed_cfg.lanczos_tol == 1e-9
        ed_space = EDReferenceScript.build_space(2, 1.0, 2.0, 0.0, 1, 1)
        @test_throws ErrorException EDReferenceScript.ground_state_lanczos(
            ed_space;
            maxiter = 1,
            tol = 1e-300,
        )

        ed_psi = ones(EDReferenceScript.dimension(ed_space))
        ed_psi ./= norm(ed_psi)
        ed_io = IOBuffer()
        EDReferenceScript.write_results(
            ed_io,
            ed_cfg,
            ed_space,
            ed_psi;
            energy = -1.25,
            iterations = 7,
            residual = 2e-12,
            lanczos_tol = 1e-9,
        )
        ed_lines = split(chomp(String(take!(ed_io))), '\n')
        @test split(ed_lines[1], ',') == [
            "L", "t", "U", "mu", "Nup", "Ndn", "nA", "boundary", "S2_ED",
            "energy", "iterations", "residual", "lanczos_tol",
        ]
        @test length(split(ed_lines[2], ',')) == 13

        schedule_cfg = DMRGReferenceScript.parse_cli([
            "--nsweeps", "4", "--maxdim-schedule", "20,50,50,100",
        ])
        @test schedule_cfg.maxdim_schedule == [20, 50, 50, 100]
        @test schedule_cfg.maxdim_cap == 100
        @test DMRGReferenceScript.effective_maxdim_schedule(schedule_cfg) ==
              [20, 50, 50, 100]
        @test_throws ErrorException DMRGReferenceScript.parse_cli([
            "--nsweeps", "3", "--maxdim-schedule", "20,10,30",
        ])
        @test_throws ErrorException DMRGReferenceScript.parse_cli([
            "--nsweeps", "3", "--maxdim-schedule", "20,40",
        ])
        @test_throws ErrorException DMRGReferenceScript.parse_cli([
            "--nsweeps", "3", "--maxdim-schedule", "20,40,80",
            "--maxdim-cap", "4800",
        ])

        dmrg_io = IOBuffer()
        DMRGReferenceScript.write_result(
            dmrg_io,
            schedule_cfg,
            -8.0,
            1.5,
            1.0,
            90,
            4e-10,
        )
        dmrg_lines = split(chomp(String(take!(dmrg_io))), '\n')
        @test split(dmrg_lines[1], ',')[13:end] == [
            "nsweeps", "cutoff", "maxdim_cap", "maxdim_schedule",
            "energy_density", "final_sweep_max_truncerr",
        ]
        @test length(split(dmrg_lines[2], ',')) == 18
        @test split(dmrg_lines[2], ',')[16] == "20;50;50;100"
        @test parse(Float64, split(dmrg_lines[2], ',')[17]) == -1.0
        @test parse(Float64, split(dmrg_lines[2], ',')[18]) == 4e-10

        mktempdir() do dir
            ed_path = joinpath(dir, "ed.csv")
            EDReferenceScript.write_results_atomic(
                ed_path,
                ed_cfg,
                ed_space,
                ed_psi;
                energy = -1.25,
                iterations = 7,
                residual = 2e-12,
                lanczos_tol = 1e-9,
            )
            @test isfile(ed_path)
            @test_throws ErrorException EDReferenceScript.write_results_atomic(
                ed_path,
                ed_cfg,
                ed_space,
                ed_psi;
                energy = -1.25,
                iterations = 7,
                residual = 2e-12,
                lanczos_tol = 1e-9,
            )

            dmrg_path = joinpath(dir, "dmrg.csv")
            DMRGReferenceScript.write_result_atomic(
                dmrg_path,
                schedule_cfg,
                -8.0,
                1.5,
                1.0,
                90,
                4e-10,
            )
            @test isfile(dmrg_path)
            @test_throws ErrorException DMRGReferenceScript.write_result_atomic(
                dmrg_path,
                schedule_cfg,
                -8.0,
                1.5,
                1.0,
                90,
                4e-10,
            )
        end
    end

    @testset "U4 L4 thermal ED audit" begin
        audit = U4L4ThermalEDAuditScript.U4L4ThermalEDAudit

        kinetic0, potential0, masks = audit.full_fock_matrices(
            4,
            1.0,
            0.0,
            0.0,
        )
        dimension = length(masks)^2
        mixed_rho = Matrix{Float64}(I, dimension, dimension) ./ dimension
        mixed_S2, mixed_purity, mixed_rhoA = audit.renyi2(
            mixed_rho,
            masks,
            4,
            2,
        )
        @test mixed_S2 ≈ 4 * log(2) atol = 1e-12
        @test mixed_purity ≈ 1 / 16 atol = 1e-14
        @test mixed_rhoA ≈ Matrix{Float64}(I, 16, 16) ./ 16

        exact_rho0, _ = audit.exact_thermal_density(
            kinetic0 + Diagonal(potential0),
            1.0,
        )
        trotter_rho0, slices0 = audit.asymmetric_trotter_density(
            kinetic0,
            potential0,
            1.0,
            0.1,
        )
        @test slices0 == 10
        @test trotter_rho0 ≈ exact_rho0 atol = 2e-12

        kinetic, potential, _ = audit.full_fock_matrices(4, 1.0, 4.0, 0.0)
        one_slice_rho, one_slice = audit.asymmetric_trotter_density(
            kinetic,
            potential,
            0.05,
            0.05,
        )
        direct_slice = Diagonal(exp.(-0.05 .* potential)) * exp(-0.05 .* kinetic)
        direct_slice ./= tr(direct_slice)
        @test one_slice == 1
        @test one_slice_rho ≈ direct_slice atol = 2e-12

        result = audit.run_audit()
        @test result.slices == 320
        @test result.fixed_ground_S2 ≈ 0.7521220360017 atol = 5e-11
        @test result.thermal_exact_S2 ≈ 0.7712028368091 atol = 5e-11
        @test result.trotter_asymmetric_S2 ≈ 0.7689951876109 atol = 5e-11
    end

    @testset "boundary condition parameter loading" begin
        mktempdir() do dir
            default_path = joinpath(dir, "default.toml")
            apbc_path = joinpath(dir, "apbc.toml")
            invalid_path = joinpath(dir, "invalid.toml")
            write_parameter_fixture(default_path)
            write_parameter_fixture(apbc_path; boundary_condition = "apbc")
            write_parameter_fixture(invalid_path; boundary_condition = "periodic")

            default_config = NodeRunnerScript.load_config(default_path)
            apbc_config = NodeRunnerScript.load_config(apbc_path)
            @test default_config.boundary_condition == :pbc
            @test apbc_config.boundary_condition == :apbc
            @test NodeRunnerScript.raw_run_name(default_config) == "fixture"
            @test NodeRunnerScript.raw_run_name(apbc_config) == "fixture__bc_apbc"
            grid_apbc_config = GridLauncherScript.load_config(apbc_path)
            @test GridLauncherScript.raw_run_name(grid_apbc_config) ==
                  "fixture__bc_apbc"
            @test SummarizeGridScript.load_metadata(apbc_path).run_name ==
                  "fixture__bc_apbc"
            @test_throws ArgumentError NodeRunnerScript.load_config(invalid_path)
            @test_throws ErrorException GridLauncherScript.load_config(invalid_path)
            @test_throws ErrorException SummarizeGridScript.load_metadata(invalid_path)
        end
    end

    @testset "Grover identity and seeds" begin
        G1 = [0.8 0.1; -0.2 0.6]
        G2 = [0.7 -0.1; 0.05 0.4]
        identity2 = Matrix{Float64}(I, 2, 2)
        @test grover_matrix(G1, G2, [1, 2]) ≈
              G1 * G2 + (identity2 - G1) * (identity2 - G2)

        seed1 = NodeRunnerScript.derived_seed(1234, 4, 1)
        seed2 = NodeRunnerScript.derived_seed(1234, 4, 2)
        @test seed1 == NodeRunnerScript.derived_seed(1234, 4, 1)
        @test seed1 != seed2
    end

    @testset "scripts are importable" begin
        @test isdefined(SummarizeGridScript, :main)
        @test isdefined(GridLauncherScript, :main)
        @test isdefined(NodeRunnerScript, :main)
        @test NodeRunnerScript.parse_cli([
            "params.toml",
            "--lambda-index",
            "1",
            "--force",
        ]).force
    end

    @testset "U4 FSS production parameter matrix" begin
        for L in 4:2:22
            beta = 4L
            nA = L ÷ 2
            path = joinpath(
                @__DIR__,
                "..",
                "parameters",
                "L$(L)_U4_beta$(beta)_nA$(nA)_Ninc16_burnin8000.toml",
            )
            @test isfile(path)
            config = TOML.parsefile(path)
            model = config["model"]
            subsystem = config["subsystem"]
            grid = config["lambda_grid"]
            sampling = config["sampling"]

            @test model["L"] == L
            @test model["U"] == 4.0
            @test model["beta"] == beta
            @test model["dtau"] == 0.05
            @test model["n_stab"] == 10
            @test model["N_sweeps"] == 1
            @test model["boundary_condition"] == "pbc"
            @test subsystem["nA"] == nA
            @test grid["Ninc"] == 16
            @test sampling["N_burnin"] == 8000
            @test sampling["N_measurements"] == 4000
            @test sampling["N_bins"] == 40
            @test sampling["root_seeds"] == [1234, 11234, 21234, 31234]
        end
    end

    @testset "U4 FSS post-campaign analysis" begin
        mktempdir() do dir
            input = joinpath(dir, "comparison.csv")
            fits_path = joinpath(dir, "fits.csv")
            plot_path = joinpath(dir, "s2.svg")
            write_u4_comparison_fixture(input)

            result = U4FSSAnalysisScript.U4FSSAnalysis.analyze(
                input,
                fits_path,
                plot_path,
            )
            @test length(result.observations) == 10
            @test length(result.fits) == 10
            @test count(fit -> fit.dataset == "DQMC" && fit.method == "WLS", result.fits) == 5
            @test count(fit -> fit.dataset == "reference" && fit.method == "OLS", result.fits) == 5
            @test count(fit -> fit.model == "total", result.fits) == 2
            @test count(fit -> fit.model == "branch", result.fits) == 4
            @test count(fit -> fit.model == "oscillation", result.fits) == 4
            reference_p1 = only(fit for fit in result.fits if
                fit.dataset == "reference" && fit.model == "oscillation" && fit.p == 1)
            @test reference_p1.c_eff_diagnostic ≈ 0.72 atol = 1e-11
            @test reference_p1.intercept ≈ 1.1 atol = 1e-11
            @test reference_p1.oscillation_amplitude ≈ 0.04 atol = 1e-11
            @test length(readlines(fits_path)) == 11
            @test occursin("c_eff is a finite-window diagnostic only", read(fits_path, String))
            svg = read(plot_path, String)
            @test startswith(svg, "<svg")
            @test occursin("DQMC +/- 1 sigma", svg)
            @test occursin("ED reference", svg)
            @test occursin("DMRG reference", svg)

            original_fits = read(fits_path)
            original_plot = read(plot_path)
            @test_throws ErrorException U4FSSAnalysisScript.U4FSSAnalysis.analyze(
                input,
                fits_path,
                plot_path,
            )
            @test read(fits_path) == original_fits
            @test read(plot_path) == original_plot
        end

        for (column, invalid_value) in (
            ("U", "2"),
            ("boundary", "apbc"),
            ("nA", "3"),
            ("beta", "15"),
            ("DQMC_error", "0"),
        )
            mktempdir() do dir
                input = joinpath(dir, "invalid.csv")
                write_u4_comparison_fixture(
                    input;
                    overrides = Dict((4, column) => invalid_value),
                )
                @test_throws ErrorException U4FSSAnalysisScript.U4FSSAnalysis.read_comparison(input)
            end
        end
    end

    @testset "strict raw CSV validation" begin
        mktempdir() do dir
            write_valid_grid(dir)
            nodes = SummarizeGridScript.read_bins(dir, TEST_METADATA)
            @test sort!(collect(keys(nodes))) == [0, 1]
            @test length(nodes[0].plus) == 2

            result = redirect_stdout(devnull) do
                SummarizeGridScript.summarize(dir, TEST_METADATA)
            end
            @test isfinite(result.S2)
            @test result.Ninc == 1
            @test_throws ErrorException SummarizeGridScript.ratio_from_pairs([
                (weight = 0.0, phase = 1.0),
            ])
        end

        mktempdir() do dir
            write_bins(joinpath(dir, "bad_bins.csv"), K0_ROWS; header = "wrong,header")
            @test_throws ErrorException SummarizeGridScript.read_bins(dir, TEST_METADATA)
        end

        mktempdir() do dir
            write_bins(
                joinpath(dir, "duplicate_bins.csv"),
                [K0_ROWS[1], K0_ROWS[1]],
            )
            @test_throws ErrorException SummarizeGridScript.read_bins(dir, TEST_METADATA)
        end

        mktempdir() do dir
            bad_rows = copy(K0_ROWS)
            bad_rows[1] = "fixture,0,0.0,7,1,1.0,1.0,0.1,NaN,"
            write_bins(joinpath(dir, "nonfinite_bins.csv"), bad_rows)
            @test_throws ErrorException SummarizeGridScript.read_bins(dir, TEST_METADATA)
        end

        mktempdir() do dir
            write_bins(joinpath(dir, "fixture_k000_seed7_bins.csv"), K0_ROWS)
            @test_throws ErrorException SummarizeGridScript.read_bins(dir, TEST_METADATA)
        end

        mktempdir() do dir
            write_bins(joinpath(dir, "fixture_k000_seed7_bins.csv"), K0_ROWS[1:1])
            write_bins(joinpath(dir, "fixture_k001_seed7_bins.csv"), K1_ROWS)
            @test_throws ErrorException SummarizeGridScript.read_bins(dir, TEST_METADATA)
        end
    end

    @testset "launcher job selection" begin
        config = (
            Ninc = 1,
            N_bins = 2,
            root_seeds = [7],
            run_name = "fixture",
            boundary_condition = "pbc",
        )
        cli = GridLauncherScript.parse_cli(["params.toml", "--all-seeds", "--max-jobs", "2"])

        mktempdir() do dir
            jobs = GridLauncherScript.build_jobs(config, cli, dir)
            @test jobs == [
                (lambda_index = 0, root_seed = 7),
                (lambda_index = 1, root_seed = 7),
            ]

            summary = GridLauncherScript.summary_path(config, dir, 0, 7)
            write(summary, "fixture\n")
            @test_throws ErrorException GridLauncherScript.build_jobs(config, cli, dir)

            bins = GridLauncherScript.bins_path(config, dir, 0, 7)
            write(bins, "fixture\n")
            @test_throws ErrorException GridLauncherScript.build_jobs(config, cli, dir)

            write(summary, GridLauncherScript.SUMMARY_HEADER * "\n" * K0_SUMMARY_ROW * "\n")
            write_bins(bins, ["fixture", "fixture"])
            @test_throws ErrorException GridLauncherScript.build_jobs(config, cli, dir)

            write_bins(bins, K0_ROWS)
            remaining = GridLauncherScript.build_jobs(config, cli, dir)
            @test remaining == [(lambda_index = 1, root_seed = 7)]
        end

        duplicate_cli = GridLauncherScript.parse_cli([
            "params.toml",
            "--lambda-indices",
            "0,0",
        ])
        mktempdir() do dir
            @test_throws ErrorException GridLauncherScript.build_jobs(
                config,
                duplicate_cli,
                dir,
            )
        end
    end
end
