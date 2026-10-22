#!/usr/bin/env julia

using Printf
using TOML
using LinearAlgebra
import JDQMCMeasurements as jdqmcm

const PROJECT_ROOT = dirname(dirname(@__FILE__))
const MODULE_PATH = joinpath(PROJECT_ROOT, "src", "HubbardRenyiDQMC.jl")

function usage()
    println("""
    Run one Hubbard Renyi lambda node.

    Usage:
      julia run/run_lambda_node.jl parameters/L10_U2_beta40_Ninc8.toml --lambda-index 4 --root-seed 1234

    Required:
      parameter file
      --lambda-index K       integer lambda node, 0 <= K <= Ninc

    Optional:
      --root-seed SEED       one user-facing seed; defaults to the first seed in the parameter file
      --output-dir DIR       override output.raw_dir from the parameter file
      --force                replace an existing complete node intentionally
      --help                 print this message

    Output:
      One summary CSV and one bins CSV under the raw result directory.
    """)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__) &&
   (isempty(ARGS) || "--help" in ARGS || "-h" in ARGS)
    usage()
    exit(0)
end

include(MODULE_PATH)
using .HubbardRenyiDQMC

function parse_cli(args)
    if isempty(args) || "--help" in args || "-h" in args
        usage()
        exit(0)
    end

    param_file = args[1]
    lambda_index = nothing
    root_seed = nothing
    output_dir = nothing
    force = false

    i = 2
    while i <= length(args)
        arg = args[i]
        if arg == "--lambda-index"
            i += 1
            i <= length(args) || error("--lambda-index requires a value")
            lambda_index = parse(Int, args[i])
        elseif arg == "--root-seed"
            i += 1
            i <= length(args) || error("--root-seed requires a value")
            root_seed = parse(Int, args[i])
        elseif arg == "--output-dir"
            i += 1
            i <= length(args) || error("--output-dir requires a value")
            output_dir = args[i]
        elseif arg == "--force"
            force = true
        else
            error("Unknown argument: $arg")
        end
        i += 1
    end

    lambda_index === nothing && error("Missing --lambda-index")

    return (
        param_file = param_file,
        lambda_index = lambda_index,
        root_seed = root_seed,
        output_dir = output_dir,
        force = force,
    )
end

function get_number(table, key, ::Type{T}) where {T}
    haskey(table, key) || error("Missing parameter: $key")
    return convert(T, table[key])
end

function resolve_project_path(path)
    isfile(path) && return path
    project_path = joinpath(PROJECT_ROOT, path)
    isfile(project_path) && return project_path
    return path
end

function load_config(path)
    config = TOML.parsefile(resolve_project_path(path))

    model = config["model"]
    subsystem = config["subsystem"]
    lambda_grid = config["lambda_grid"]
    sampling = config["sampling"]
    output = config["output"]
    boundary_condition = normalize_boundary_condition(
        get(model, "boundary_condition", "pbc"),
    )

    return (
        L = get_number(model, "L", Int),
        t = get_number(model, "t", Float64),
        U = get_number(model, "U", Float64),
        mu = get_number(model, "mu", Float64),
        beta = get_number(model, "beta", Float64),
        dtau = get_number(model, "dtau", Float64),
        n_stab = get_number(model, "n_stab", Int),
        N_sweeps = get_number(model, "N_sweeps", Int),
        boundary_condition = boundary_condition,
        nA = get_number(subsystem, "nA", Int),
        Ninc = get_number(lambda_grid, "Ninc", Int),
        S2_ED = get_number(lambda_grid, "S2_ED", Float64),
        N_burnin = get_number(sampling, "N_burnin", Int),
        N_measurements = get_number(sampling, "N_measurements", Int),
        N_bins = get_number(sampling, "N_bins", Int),
        root_seeds = [convert(Int, seed) for seed in sampling["root_seeds"]],
        run_name = String(output["run_name"]),
        raw_dir = String(output["raw_dir"]),
        processed_dir = String(output["processed_dir"]),
        figures_dir = String(output["figures_dir"])
    )
end

function raw_run_name(config)
    config.boundary_condition == :pbc && return config.run_name
    return "$(config.run_name)__bc_$(config.boundary_condition)"
end

function splitmix64(x::UInt64)
    z = x + 0x9e3779b97f4a7c15
    z = xor(z, z >> 30) * 0xbf58476d1ce4e5b9
    z = xor(z, z >> 27) * 0x94d049bb133111eb
    return xor(z, z >> 31)
end

function derived_seed(root_seed::Integer, lambda_index::Integer, replica_id::Integer)
    x = UInt64(root_seed)
    x = xor(x, UInt64(lambda_index + 1) * 0x9e3779b97f4a7c15)
    x = xor(x, UInt64(replica_id) * 0xbf58476d1ce4e5b9)
    return Int(splitmix64(x) % UInt64(typemax(Int)))
end

function ratio_statistics(avg_weight, avg_phase)
    ratio, ratio_err =
        jdqmcm.jackknife(
            (weight_bin, phase_bin) -> real(weight_bin / phase_bin),
            avg_weight,
            avg_phase
        )

    entropy, entropy_err =
        jdqmcm.jackknife(
            (weight_bin, phase_bin) -> -log(abs(weight_bin / phase_bin)),
            avg_weight,
            avg_phase
        )

    return (
        ratio = real(ratio),
        ratio_err = ratio_err,
        entropy = real(entropy),
        entropy_err = entropy_err
    )
end

function make_node_measurements!(
    sgndetGup,
    sgndetGdn,
    bin,
    avg_phase,
    avg_den,
    avg_num,
    avg_step_plus,
    avg_step_minus,
    gCup,
    gCdn,
    lambda_value,
    delta_lambda,
    plus_defined,
    minus_defined
)
    sgn = sign(sgndetGup * sgndetGdn)
    X_raw = det(gCup) * det(gCdn)

    if !isfinite(real(X_raw)) || !isfinite(imag(X_raw))
        error("Encountered nonfinite Grover estimator X = $X_raw")
    end

    imaginary_tolerance = 1e-10 * max(1.0, abs(real(X_raw)))
    if abs(imag(X_raw)) > imaginary_tolerance
        error("Grover estimator has a significant imaginary part: X = $X_raw")
    end

    X = real(X_raw)
    if X <= 0.0
        error("Power-increment measurement requires X > 0, but X = $X")
    end

    logX = log(X)
    reweight = exp(-lambda_value * logX)

    avg_phase[bin] += sgn
    avg_den[bin] += sgn * reweight
    avg_num[bin] += sgn * X * reweight

    if plus_defined
        avg_step_plus[bin] += sgn * exp(delta_lambda * logX)
    end

    if minus_defined
        avg_step_minus[bin] += sgn * exp(-delta_lambda * logX)
    end

    return nothing
end

function run_lambda_node_simulation(config, lambda_index, root_seed)
    @assert config.Ninc >= 1
    @assert 0 <= lambda_index <= config.Ninc
    @assert config.N_measurements % config.N_bins == 0

    delta_lambda = 1.0 / config.Ninc
    lambda_value = lambda_index * delta_lambda
    plus_defined = lambda_index < config.Ninc
    minus_defined = lambda_index > 0

    seed1 = derived_seed(root_seed, lambda_index, 1)
    seed2 = derived_seed(root_seed, lambda_index, 2)

    t_start = time()
    deltaG = 0.0
    acceptance_rate = 0.0
    N_binsize = div(config.N_measurements, config.N_bins)
    C = collect(1:config.nA)

    replica1 = initializer(
        config.t,
        config.U,
        config.mu,
        config.beta,
        config.L,
        config.dtau,
        seed1;
        n_stab = config.n_stab,
        boundary_condition = config.boundary_condition,
    )

    replica2 = initializer(
        config.t,
        config.U,
        config.mu,
        config.beta,
        config.L,
        config.dtau,
        seed2;
        n_stab = config.n_stab,
        boundary_condition = config.boundary_condition,
    )

    avg_phase = zeros(Float64, config.N_bins)
    avg_den = zeros(Float64, config.N_bins)
    avg_num = zeros(Float64, config.N_bins)
    avg_step_plus = zeros(Float64, config.N_bins)
    avg_step_minus = zeros(Float64, config.N_bins)

    for update in 1:config.N_burnin
        replica1, deltaG, acc1 =
            update_replica_rank1_fused!(replica1, replica2, deltaG, C, lambda_value)
        replica2, deltaG, acc2 =
            update_replica_rank1_fused!(replica2, replica1, deltaG, C, lambda_value)
        acceptance_rate += (acc1 + acc2) / 2
    end

    for bin in 1:config.N_bins
        for measurement in 1:N_binsize
            for sweep in 1:config.N_sweeps
                replica1, deltaG, acc1 =
                    update_replica_rank1_fused!(replica1, replica2, deltaG, C, lambda_value)
                replica2, deltaG, acc2 =
                    update_replica_rank1_fused!(replica2, replica1, deltaG, C, lambda_value)
                acceptance_rate += (acc1 + acc2) / 2
            end

            gCup = grover_matrix(replica1.Gup_00, replica2.Gup_00, C)
            gCdn = grover_matrix(replica1.Gdn_00, replica2.Gdn_00, C)

            make_node_measurements!(
                replica1.sgndetGup * replica2.sgndetGup,
                replica1.sgndetGdn * replica2.sgndetGdn,
                bin,
                avg_phase,
                avg_den,
                avg_num,
                avg_step_plus,
                avg_step_minus,
                gCup,
                gCdn,
                lambda_value,
                delta_lambda,
                plus_defined,
                minus_defined
            )
        end

        avg_phase[bin] /= N_binsize
        avg_den[bin] /= N_binsize
        avg_num[bin] /= N_binsize
        if plus_defined
            avg_step_plus[bin] /= N_binsize
        end
        if minus_defined
            avg_step_minus[bin] /= N_binsize
        end
    end

    N_update_rounds =
        config.N_burnin +
        config.N_bins * N_binsize * config.N_sweeps

    acceptance_rate /= N_update_rounds
    runtime = time() - t_start

    endpoint_stats = power_increment_statistics(avg_phase, avg_den, avg_num)
    plus_stats = plus_defined ? ratio_statistics(avg_step_plus, avg_phase) : nothing
    minus_stats = minus_defined ? ratio_statistics(avg_step_minus, avg_phase) : nothing

    return (
        lambda_index = lambda_index,
        lambda_value = lambda_value,
        delta_lambda = delta_lambda,
        plus_defined = plus_defined,
        minus_defined = minus_defined,
        root_seed = root_seed,
        seed1 = seed1,
        seed2 = seed2,
        C = C,
        acceptance = acceptance_rate,
        deltaG = deltaG,
        runtime = runtime,
        avg_phase = avg_phase,
        avg_den = avg_den,
        avg_num = avg_num,
        avg_step_plus = avg_step_plus,
        avg_step_minus = avg_step_minus,
        endpoint_stats = endpoint_stats,
        plus_stats = plus_stats,
        minus_stats = minus_stats
    )
end

function write_summary_csv(path, config, result)
    open(path, "w") do io
        println(
            io,
            join(
                [
                    "run_name",
                    "L",
                    "t",
                    "U",
                    "mu",
                    "beta",
                    "dtau",
                    "n_stab",
                    "N_sweeps",
                    "nA",
                    "Ninc",
                    "delta_lambda",
                    "lambda_index",
                    "lambda",
                    "root_seed",
                    "seed_replica_1",
                    "seed_replica_2",
                    "direction",
                    "lambda_from",
                    "lambda_to",
                    "ratio",
                    "ratio_err",
                    "entropy_increment",
                    "entropy_increment_err",
                    "Z1_over_Zlambda",
                    "Z1_over_Zlambda_err",
                    "Zlambda_over_Z0",
                    "Zlambda_over_Z0_err",
                    "acceptance",
                    "deltaG",
                    "runtime_seconds"
                ],
                ","
            )
        )

        function write_direction(direction, lambda_from, lambda_to, stats)
            println(
                io,
                join(
                    [
                        raw_run_name(config),
                        config.L,
                        @sprintf("%.8f", config.t),
                        @sprintf("%.8f", config.U),
                        @sprintf("%.8f", config.mu),
                        @sprintf("%.8f", config.beta),
                        @sprintf("%.8f", config.dtau),
                        config.n_stab,
                        config.N_sweeps,
                        config.nA,
                        config.Ninc,
                        @sprintf("%.12f", result.delta_lambda),
                        result.lambda_index,
                        @sprintf("%.12f", result.lambda_value),
                        result.root_seed,
                        result.seed1,
                        result.seed2,
                        direction,
                        @sprintf("%.12f", lambda_from),
                        @sprintf("%.12f", lambda_to),
                        @sprintf("%.12e", stats.ratio),
                        @sprintf("%.12e", stats.ratio_err),
                        @sprintf("%.12e", stats.entropy),
                        @sprintf("%.12e", stats.entropy_err),
                        @sprintf("%.12e", result.endpoint_stats.Z1_over_Zlambda),
                        @sprintf("%.12e", result.endpoint_stats.Z1_over_Zlambda_err),
                        @sprintf("%.12e", result.endpoint_stats.Zlambda_over_Z0),
                        @sprintf("%.12e", result.endpoint_stats.Zlambda_over_Z0_err),
                        @sprintf("%.8f", result.acceptance),
                        @sprintf("%.12e", result.deltaG),
                        @sprintf("%.6f", result.runtime)
                    ],
                    ","
                )
            )
        end

        if result.plus_defined
            write_direction(
                "plus",
                result.lambda_value,
                result.lambda_value + result.delta_lambda,
                result.plus_stats
            )
        end

        if result.minus_defined
            write_direction(
                "minus",
                result.lambda_value,
                result.lambda_value - result.delta_lambda,
                result.minus_stats
            )
        end
    end

    return nothing
end

function write_bins_csv(path, config, result)
    open(path, "w") do io
        println(
            io,
            "run_name,lambda_index,lambda,root_seed,bin,phase,den,num,step_plus,step_minus"
        )

        for bin in eachindex(result.avg_phase)
            println(
                io,
                join(
                    [
                        raw_run_name(config),
                        result.lambda_index,
                        @sprintf("%.12f", result.lambda_value),
                        result.root_seed,
                        bin,
                        @sprintf("%.12e", result.avg_phase[bin]),
                        @sprintf("%.12e", result.avg_den[bin]),
                        @sprintf("%.12e", result.avg_num[bin]),
                        result.plus_defined ?
                            @sprintf("%.12e", result.avg_step_plus[bin]) :
                            "",
                        result.minus_defined ?
                            @sprintf("%.12e", result.avg_step_minus[bin]) :
                            ""
                    ],
                    ","
                )
            )
        end
    end

    return nothing
end

function write_node_outputs(summary_csv, bins_csv, config, result; force = false)
    if !force && (ispath(summary_csv) || ispath(bins_csv))
        error("Node output already exists; use --force only for an intentional rerun")
    end

    summary_tmp, summary_io = mktemp(dirname(summary_csv))
    bins_tmp, bins_io = mktemp(dirname(bins_csv))
    close(summary_io)
    close(bins_io)

    try
        write_summary_csv(summary_tmp, config, result)
        write_bins_csv(bins_tmp, config, result)
        mv(summary_tmp, summary_csv; force = force)
        mv(bins_tmp, bins_csv; force = force)
    finally
        isfile(summary_tmp) && rm(summary_tmp; force = true)
        isfile(bins_tmp) && rm(bins_tmp; force = true)
    end

    return nothing
end

function main(args = ARGS)
    cli = parse_cli(args)
    BLAS.set_num_threads(1)
    config = load_config(cli.param_file)

    if cli.lambda_index < 0 || cli.lambda_index > config.Ninc
        error("--lambda-index must satisfy 0 <= K <= Ninc, with Ninc = $(config.Ninc)")
    end

    root_seed =
        cli.root_seed === nothing ?
        first(config.root_seeds) :
        cli.root_seed

    output_dir =
        cli.output_dir === nothing ?
        joinpath(PROJECT_ROOT, config.raw_dir) :
        cli.output_dir

    mkpath(output_dir)

    prefix = joinpath(
        output_dir,
        @sprintf("%s_k%03d_seed%d", config.run_name, cli.lambda_index, root_seed),
    )
    summary_csv = prefix * "_summary.csv"
    bins_csv = prefix * "_bins.csv"

    if !cli.force && (ispath(summary_csv) || ispath(bins_csv))
        error("Node output already exists; refusing to run without --force")
    end

    println()
    println("========================================")
    println("Hubbard Renyi lambda-node run")
    println("========================================")
    println("parameter file = ", cli.param_file)
    println("run name       = ", config.run_name)
    println("boundary       = ", config.boundary_condition)
    println("lambda index   = ", cli.lambda_index, " / ", config.Ninc)
    println("root seed      = ", root_seed)
    println("BLAS threads   = ", BLAS.get_num_threads())
    println("output dir     = ", output_dir)
    println("========================================")
    flush(stdout)

    result = run_lambda_node_simulation(config, cli.lambda_index, root_seed)
    write_node_outputs(summary_csv, bins_csv, config, result; force = cli.force)

    println()
    println("Finished lambda-node run")
    println("summary CSV = ", summary_csv)
    println("bins CSV    = ", bins_csv)
    if result.plus_defined
        println(
            "R_plus  = ",
            result.plus_stats.ratio,
            " +/- ",
            result.plus_stats.ratio_err
        )
    end
    if result.minus_defined
        println(
            "R_minus = ",
            result.minus_stats.ratio,
            " +/- ",
            result.minus_stats.ratio_err
        )
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
