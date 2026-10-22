using LinearAlgebra, Statistics, TOML, SHA, Pkg, Dates

const DQMC_ROOT_SEEDS = [17011, 27011, 37013, 47017]
const DQMC_SETTINGS = Dict{String,Any}(
    "Lx" => 2, "Ly" => 2, "sites" => 4, "t" => 1.0,
    "periodic" => [false, false], "beta" => 2.0, "dtau" => 0.05,
    "U_values" => [0.0, 2.0], "mu" => 0.0, "n_stab" => 10,
    "burnin" => 1000, "measurements" => 20000, "bins" => 40,
    "sweeps_per_measurement" => 1, "root_seeds" => DQMC_ROOT_SEEDS,
    "symmetric" => false, "checkerboard" => false,
    "hamiltonian" => "sum(K_ij cdag_i c_j) + U sum((nup-1/2)(ndn-1/2)) - mu N",
    "ensemble" => "grand canonical", "energy_normalization" => "per site",
    "propagator" => "exp(-dtau V) exp(-dtau K)")

function ratio_stats(x, signs)
    n = length(x)
    estimate = sum(x) / sum(signs)
    deleted = [(sum(x) - x[i]) / (sum(signs) - signs[i]) for i in 1:n]
    error = sqrt((n-1) / n * sum(abs2, deleted .- mean(deleted)))
    return estimate, error
end

function prepare_output(output, engine)
    mkpath(output)
    paths = [joinpath(output, engine * suffix) for suffix in
             ("_bins.csv", "_summary.csv", "_metadata.toml")]
    any(ispath, paths) && error("Evidence already exists; choose a new output directory")
    return paths
end

function metadata(engine, settings, extras = Dict{String,Any}())
    packages = Dict{String,Any}()
    for (_, dependency) in Pkg.dependencies()
        isnothing(dependency.version) && continue
        packages[dependency.name] = Dict("version" => string(dependency.version),
            "tree_hash" => isnothing(dependency.tree_hash) ? "" : string(dependency.tree_hash))
    end
    return merge(Dict{String,Any}(
        "engine" => engine, "timestamp_utc" => string(now(UTC)),
        "julia_version" => string(VERSION), "kernel" => string(Sys.KERNEL),
        "architecture" => string(Sys.ARCH), "cpu" => Sys.cpu_info()[1].model,
        "julia_threads" => Threads.nthreads(), "blas_threads" => BLAS.get_num_threads(),
        "blas" => string(BLAS.get_config()), "packages" => packages,
        "settings" => settings,
        "script_sha256" => Dict(basename(path) => bytes2hex(sha256(read(path)))
            for path in (abspath(PROGRAM_FILE), joinpath(@__DIR__, "dqmc_common.jl")))), extras)
end

function write_results(summary_io, bin_io, engine, U, seed, data, deltaG)
    names = (:energy_per_site, :density, :double_occupancy)
    for name in names
        value, error = ratio_stats(getproperty(data, name), data.sign)
        println(summary_io, join((engine, U, seed, name, value, error, mean(data.sign), deltaG), ','))
        println("$engine U=$U seed=$seed $name=$value +/- $error")
    end
    for i in eachindex(data.sign)
        println(bin_io, join((engine, U, seed, i, data.sign[i],
            (getproperty(data, name)[i] for name in names)...), ','))
    end
    flush(summary_io); flush(bin_io); flush(stdout)
end
