#!/usr/bin/env julia
# Explicit output paths prevent replacement of previous benchmark evidence.
using HubbardRenyiDQMC
using LinearAlgebra
using Random
using Statistics
using TOML
using Printf
using SHA

function save_metadata(path, settings)
    ispath(path) && error("Refusing to overwrite $path")
    metadata = Dict{String, Any}(
        "julia_version" => string(VERSION), "kernel" => string(Sys.KERNEL),
        "architecture" => string(Sys.ARCH), "cpu" => Sys.cpu_info()[1].model,
        "julia_threads" => Threads.nthreads(), "blas_threads" => BLAS.get_num_threads(),
        "blas" => string(BLAS.get_config()),
        "git_commit" => readchomp(`git rev-parse HEAD`),
        "working_tree_dirty" => !isempty(readchomp(`git status --porcelain`)),
        "source_sha256" => Dict(name => bytes2hex(sha256(read(joinpath(@__DIR__, "..", name))))
            for name in ("src/dqmc_kernel.jl", "src/grover_rank1_update.jl",
                "src/square_lattice.jl", "analysis/cv_benchmark_dqmc.jl")),
        "hamiltonian" => "sum(K_ij cdag_i c_j) + U sum((nup-1/2)(ndn-1/2)) - mu N",
        "settings" => settings)
    open(path, "w") do io
        TOML.print(io, metadata)
    end
end

function thermal_validation(output)
    paths = [joinpath(output, name) for name in
        ("dqmc_ed_thermal.csv", "dqmc_ed_bins.csv", "dqmc_ed_metadata.toml")]
    any(ispath, paths) && error("Thermal benchmark evidence already exists; choose a new directory")
    mkpath(output)
    seeds = [1103, 2207, 3313, 4409]
    reference = ed_thermal(square_lattice_kinetic_matrix(2, 2); U = 2.0, mu = 0.0, beta = 2.0)
    save_metadata(paths[3], Dict("Lx" => 2, "Ly" => 2, "U" => 2.0, "mu" => 0.0,
        "beta" => 2.0, "periodic" => [true, true], "dtau" => [0.1, 0.05],
        "burnin" => 1000, "measurements_per_seed" => 20000, "bins_per_seed" => 40,
        "root_seeds" => seeds, "ensemble" => "grand canonical", "reference" => "full-spectrum exact diagonalization"))
    open(paths[1], "w") do summary
        println(summary, "dtau,root_seed,observable,dqmc,standard_error,ed,difference,z_sampling_only,average_sign,deltaG_max")
        open(paths[2], "w") do binio
            println(binio, "dtau,root_seed,bin,sign,energy,density,double_occupancy")
            for dtau in (0.1, 0.05)
                combined = Dict(key => Float64[] for key in (:sign, :energy, :density, :double_occupancy))
                for root_seed in seeds
                    result = run_square_dqmc(; Lx = 2, Ly = 2, U = 2.0, mu = 0.0,
                        beta = 2.0, dtau, periodic = (true, true), root_seed,
                        burnin = 1000, measurements = 20000, bins = 40)
                    for name in (:energy, :density, :double_occupancy)
                        value, error = getproperty(result, name), getproperty(result, Symbol(name, "_error"))
                        exact = getproperty(reference, name)
                        println(summary, join((dtau, root_seed, name, value, error, exact,
                            value - exact, error > 0 ? (value - exact) / error : NaN,
                            result.average_sign, result.deltaG_max), ','))
                    end
                    for bin in 1:result.bins
                        println(binio, join((dtau, root_seed, bin,
                            (getproperty(result.bin_data, name)[bin] for name in
                                (:sign, :energy, :density, :double_occupancy))...), ','))
                    end
                    for name in keys(combined)
                        append!(combined[name], getproperty(result.bin_data, name))
                    end
                    @printf("Thermal ED check: dtau=%.2f seed=%d energy=%.8f +/- %.8f\n",
                        dtau, root_seed, result.energy, result.energy_error)
                    flush(stdout)
                end
                for name in (:energy, :density, :double_occupancy)
                    estimate, error = HubbardRenyiDQMC.jdqmcm.jackknife((x, s) -> x / s,
                        combined[name], combined[:sign])
                    exact = getproperty(reference, name)
                    println(summary, join((dtau, "pooled", name, estimate, error, exact,
                        estimate - exact, error > 0 ? (estimate - exact) / error : NaN,
                        mean(combined[:sign]), ""), ','))
                    @printf("Pooled dtau=%.2f %s=%.10f +/- %.10f; ED=%.10f\n", dtau, name, estimate, error, exact)
                end
            end
        end
    end
end

function update_rounds(f, r1, r2, rounds, C, lambda)
    deltaG = 0.0
    for _ in 1:rounds
        r1, deltaG, _ = f(r1, r2, deltaG, C, lambda)
        r2, deltaG, _ = f(r2, r1, deltaG, C, lambda)
    end
    return r1, r2, deltaG
end

function time_kernel(f, base1, base2, rounds, C, lambda)
    r1, r2 = deepcopy(base1), deepcopy(base2)
    # Compile and initialize the selected path before timing.
    r1, r2, _ = update_rounds(f, r1, r2, 2, C, lambda)
    GC.gc()
    timing = @timed update_rounds(f, r1, r2, rounds, C, lambda)
    return (; seconds = timing.time, bytes = timing.bytes, gc_seconds = timing.gctime,
        deltaG_max = timing.value[3])
end

function check_matched_kernels(base1, base2, C, lambda)
    direct = update_rounds(update_replica_fused!, deepcopy(base1), deepcopy(base2), 3, C, lambda)
    rank1 = update_rounds(update_replica_rank1_fused!, deepcopy(base1), deepcopy(base2), 3, C, lambda)
    max_error = 0.0
    for replica in 1:2
        direct[replica].s == rank1[replica].s || error("Matched-kernel auxiliary fields disagree")
        rand(copy(direct[replica].rng), UInt64) == rand(copy(rank1[replica].rng), UInt64) ||
            error("Matched-kernel RNG streams disagree")
        for name in (:logdetGup, :logdetGdn)
            isapprox(getproperty(direct[replica], name), getproperty(rank1[replica], name);
                rtol = 2e-8, atol = 2e-9) || error("Matched-kernel $name disagrees")
        end
        for name in (:sgndetGup, :sgndetGdn)
            getproperty(direct[replica], name) == getproperty(rank1[replica], name) ||
                error("Matched-kernel $name disagrees")
        end
        for name in (:Gup, :Gdn, :Gup_00, :Gdn_00, :Gup_τ0, :Gup_0τ, :Gdn_τ0, :Gdn_0τ)
            A, B = getproperty(direct[replica], name), getproperty(rank1[replica], name)
            isapprox(A, B; rtol = 2e-8, atol = 2e-9) || error("Matched-kernel $name disagrees")
            max_error = max(max_error, maximum(abs.(A - B)))
        end
    end
    return max_error
end

function rank1_timing(output)
    paths = [joinpath(output, name) for name in ("rank1_timing.csv", "rank1_metadata.toml")]
    any(ispath, paths) && error("Timing evidence already exists; choose a new directory")
    mkpath(output)
    rounds, repeats = 100, 7
    save_metadata(paths[2], Dict("rounds" => rounds, "repeats" => repeats,
        "warmup_rounds_per_repeat" => 2, "lambda" => 0.5,
        "scope" => "two-replica stabilized updates, excluding initialization and measurements",
        "baseline" => "direct Grover determinant rebuilds; Green functions already use rank-one updates",
        "optimized" => "rank-one Grover determinant/inverse tracker",
        "square" => Dict("Lx" => 4, "Ly" => 4, "U" => 4.0, "beta" => 4.0,
            "dtau" => 0.1, "nA" => 8, "n_stab" => 10),
        "chain" => Dict("L" => 8, "U" => 2.0, "beta" => 32.0,
            "dtau" => 0.05, "nA" => 4, "n_stab" => 10)))
    open(paths[1], "w") do io
        println(io, "geometry,repeat,order,kernel,rounds,seconds,allocated_bytes,gc_seconds,deltaG_max,matched_G_max_error")
        for geometry in (:square, :chain)
            C = collect(1:(geometry == :square ? 8 : 4))
            direct_times, rank1_times = Float64[], Float64[]
            for repeat in 1:repeats
                seed1, seed2 = 62001 + 100repeat, 62002 + 100repeat
                if geometry == :square
                    K = square_lattice_kinetic_matrix(4, 4)
                    base1 = initializer_from_kinetic(K, 4.0, 0.0, 4.0, 0.1, seed1)
                    base2 = initializer_from_kinetic(K, 4.0, 0.0, 4.0, 0.1, seed2)
                else
                    base1 = initializer(1.0, 2.0, 0.0, 32.0, 8, 0.05, seed1)
                    base2 = initializer(1.0, 2.0, 0.0, 32.0, 8, 0.05, seed2)
                end
                error = check_matched_kernels(base1, base2, C, 0.5)
                kernels = isodd(repeat) ? (:direct, :rank1) : (:rank1, :direct)
                for (order, kernel) in enumerate(kernels)
                    f = kernel == :direct ? update_replica_fused! : update_replica_rank1_fused!
                    result = time_kernel(f, base1, base2, rounds, C, 0.5)
                    push!(kernel == :direct ? direct_times : rank1_times, result.seconds)
                    println(io, join((geometry, repeat, order, kernel, rounds,
                        result.seconds, result.bytes, result.gc_seconds, result.deltaG_max, error), ','))
                    flush(io)
                end
                @printf("%s matched repeat %d/7 complete\n", geometry, repeat)
                flush(stdout)
            end
            @printf("%s: ratio of mean direct/rank-one time = %.6fx; median paired ratio = %.6fx\n",
                geometry, mean(direct_times) / mean(rank1_times), median(direct_times ./ rank1_times))
        end
    end
end

function main(args = ARGS)
    length(args) == 2 || error("Usage: julia --project=. .\\analysis\\cv_benchmark_dqmc.jl thermal|timing OUTPUT_DIRECTORY")
    BLAS.set_num_threads(1)
    args[1] == "thermal" ? thermal_validation(args[2]) :
        args[1] == "timing" ? rank1_timing(args[2]) : error("Choose thermal or timing")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
