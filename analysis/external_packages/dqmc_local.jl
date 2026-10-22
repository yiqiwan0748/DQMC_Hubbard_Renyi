#!/usr/bin/env julia
# Run with the public package project, in a separate Julia process from SmoQy.
using HubbardRenyiDQMC
include(joinpath(@__DIR__, "dqmc_common.jl"))

function main(args = ARGS)
    length(args) == 1 || error("Usage: julia --project=PUBLIC_PACKAGE dqmc_local.jl OUTPUT")
    BLAS.set_num_threads(1)
    paths = prepare_output(args[1], "local")
    K = square_lattice_kinetic_matrix(2, 2; periodic = (false, false))
    expected = [0. -1 -1 0; -1 0 0 -1; -1 0 0 -1; 0 -1 -1 0]
    K == expected || error("Unexpected four-bond square hopping matrix")
    pkgroot = dirname(dirname(pathof(HubbardRenyiDQMC)))
    extras = Dict{String,Any}(
        "package_source" => pkgroot,
        "package_git_commit" => readchomp(`git -C $pkgroot rev-parse HEAD`),
        "package_working_tree_dirty" => !isempty(readchomp(`git -C $pkgroot status --porcelain`)),
        "K_rows" => [collect(K[i,:]) for i in axes(K,1)],
        "source_sha256" => Dict(name => bytes2hex(sha256(read(joinpath(pkgroot,name))))
            for name in ("src/square_lattice.jl", "src/dqmc_kernel.jl")))
    open(paths[3], "w") do io
        TOML.print(io, metadata("HubbardRenyiDQMC", copy(DQMC_SETTINGS), extras))
    end
    open(paths[1], "w") do bins
        println(bins, "engine,U,root_seed,bin,sign,energy_per_site,density,double_occupancy")
        open(paths[2], "w") do summary
            println(summary, "engine,U,root_seed,observable,mean,standard_error,average_sign,deltaG_max")
            for U in (0.0, 2.0), seed in DQMC_ROOT_SEEDS
                result = run_square_dqmc(; Lx = 2, Ly = 2, periodic = (false,false),
                    U, mu = 0.0, beta = 2.0, dtau = 0.05, n_stab = 10,
                    root_seed = seed, burnin = 1000, measurements = 20000, bins = 40)
                data = (; sign = result.bin_data.sign,
                    energy_per_site = result.bin_data.energy ./ 4,
                    density = result.bin_data.density,
                    double_occupancy = result.bin_data.double_occupancy)
                write_results(summary, bins, "HubbardRenyiDQMC", U, seed, data, result.deltaG_max)
            end
        end
    end
end
abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
