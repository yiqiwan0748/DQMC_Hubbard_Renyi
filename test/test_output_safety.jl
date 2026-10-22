using Test
using TOML

module SquareOutputScript
include(joinpath(@__DIR__, "..", "run", "run_square_hubbard.jl"))
end

module DMRGOutputScript
include(joinpath(@__DIR__, "..", "examples", "dmrg_square.jl"))
end

const MOCK_SQUARE_RESULT = (
    energy=-4.0, energy_error=0.1, energy_per_site=-1.0,
    energy_per_site_error=0.025, density=1.0, density_error=0.0,
    double_occupancy=0.2, double_occupancy_error=0.01,
    average_sign=1.0, acceptance=0.5, deltaG_max=1e-14, bins=2,
    bin_data=(sign=[1.0, 1.0], energy=[-4.1, -3.9],
        density=[1.0, 1.0], double_occupancy=[0.19, 0.21]),
)

const MOCK_DMRG_RESULT = (
    Nup=2, Ndn=2, energy=-6.1, S2=0.9, cut=2, maxlinkdim=16,
    nsweeps=10, sweep_energy_change=1e-13, energy_variance=0.0,
    final_truncation_error=1e-15, root_seed=1234,
)

@testset "Square DQMC reserves a new output directory before sampling" begin
    mktempdir() do directory
        settings = joinpath(@__DIR__, "..", "parameters", "square_4x4.toml")
        output = joinpath(directory, "nested", "run")
        args = [settings, "--output-dir", output]
        calls = Ref(0)
        function sample(; kwargs...)
            calls[] += 1
            @test kwargs[:root_seed] == 1234
            @test isfile(joinpath(output, "parameters.toml"))
            # A second run attempts to start while the first still samples.
            competitor_called = Ref(false)
            @test_throws ErrorException SquareOutputScript.square_hubbard_main(args;
                run_simulation=(; kwargs...) -> begin
                    competitor_called[] = true
                    error("must not sample")
                end)
            @test !competitor_called[]
            return MOCK_SQUARE_RESULT
        end
        SquareOutputScript.square_hubbard_main(args; run_simulation=sample)
        @test calls[] == 1
        @test TOML.parsefile(joinpath(output, "parameters.toml")) == TOML.parsefile(settings)
        summary = read(joinpath(output, "summary.csv"), String)
        @test startswith(summary, "observable,mean,standard_error\n")
        @test countlines(joinpath(output, "bins.csv")) == 3
        @test_throws ErrorException SquareOutputScript.square_hubbard_main(args;
            run_simulation=sample)
        @test calls[] == 1
        @test read(joinpath(output, "summary.csv"), String) == summary

        failed_output = joinpath(directory, "failed")
        @test_throws ErrorException SquareOutputScript.square_hubbard_main(
            [settings, "--output-dir", failed_output];
            run_simulation=(; kwargs...) -> error("sampling failed"))
        @test isfile(joinpath(failed_output, "parameters.toml"))
        @test !ispath(joinpath(failed_output, "summary.csv"))
    end
end

@testset "DMRG reserves its CSV and records the complete model and solver settings" begin
    mktempdir() do directory
        settings = joinpath(@__DIR__, "..", "examples", "dmrg_square.toml")
        output = joinpath(directory, "nested", "dmrg.csv")
        args = [settings, "--output", output]
        calls = Ref(0)
        function solve(K; kwargs...)
            calls[] += 1
            @test kwargs[:root_seed] == 1234
            @test isfile(output)
            @test isempty(read(output, String))
            competitor_called = Ref(false)
            @test_throws ErrorException DMRGOutputScript.main(args;
                solve_ground_state=(K; kwargs...) -> begin
                    competitor_called[] = true
                    error("must not solve")
                end)
            @test !competitor_called[]
            return MOCK_DMRG_RESULT
        end
        DMRGOutputScript.main(args; solve_ground_state=solve)
        @test calls[] == 1
        lines = readlines(output)
        @test length(lines) == 2
        header, values = split.(lines, ',')
        @test length(header) == length(values)
        row = Dict(zip(header, values))
        @test row["t"] == "1.0"
        @test row["periodic_x"] == "false"
        @test row["periodic_y"] == "false"
        @test parse(Float64, row["cutoff"]) == 1e-12
        @test row["maxdim_schedule"] == "16;32;64;128"
        @test_throws ErrorException DMRGOutputScript.main(args; solve_ground_state=solve)
        @test calls[] == 1
        @test readlines(output) == lines

        failed_output = joinpath(directory, "failed.csv")
        @test_throws ErrorException DMRGOutputScript.main([settings, "--output", failed_output];
            solve_ground_state=(K; kwargs...) -> error("solver failed"))
        @test isfile(failed_output)
        @test isempty(read(failed_output, String))
        @test_throws ErrorException DMRGOutputScript.reserve_dmrg_output(failed_output)

        # O_EXCL must reject even a dangling symlink instead of following it.
        if !Sys.iswindows()
            link_output = joinpath(directory, "link.csv")
            absent_target = joinpath(directory, "absent.csv")
            symlink(absent_target, link_output)
            @test_throws ErrorException DMRGOutputScript.reserve_dmrg_output(link_output)
            @test !ispath(absent_target)
        end
    end
end
