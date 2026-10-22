#!/usr/bin/env julia

using Printf
using TOML

const PROJECT_ROOT = dirname(dirname(@__FILE__))
const NODE_SCRIPT = joinpath(PROJECT_ROOT, "run", "run_lambda_node.jl")
const BINS_HEADER =
    "run_name,lambda_index,lambda,root_seed,bin,phase,den,num,step_plus,step_minus"
const SUMMARY_HEADER = join(
    [
        "run_name", "L", "t", "U", "mu", "beta", "dtau", "n_stab",
        "N_sweeps", "nA", "Ninc", "delta_lambda", "lambda_index",
        "lambda", "root_seed", "seed_replica_1", "seed_replica_2",
        "direction", "lambda_from", "lambda_to", "ratio", "ratio_err",
        "entropy_increment", "entropy_increment_err", "Z1_over_Zlambda",
        "Z1_over_Zlambda_err", "Zlambda_over_Z0", "Zlambda_over_Z0_err",
        "acceptance", "deltaG", "runtime_seconds",
    ],
    ",",
)

function usage()
    println("""
    Run Hubbard Renyi lambda nodes in a rolling local process pool.

    Usage:
      julia run/run_grid_local.jl parameters/L10_U2_beta40_Ninc8.toml --root-seed 1234 --max-jobs 4

    Defaults:
      - runs all lambda indices 0:Ninc
      - uses the first root seed in the parameter file
      - skips nodes whose complete raw outputs already exist
      - launches up to 4 Julia processes at a time

    Options:
      --root-seed SEED          run one root seed
      --all-seeds               run every root seed in the parameter file
      --lambda-indices LIST     comma list, for example 0,1,2,6,7,8
      --max-jobs N              number of simultaneous Julia processes
      --output-dir DIR          write raw outputs outside the configured directory
      --force                   rerun even if raw node outputs exist
      --help                    print this message
    """)
end

function resolve_project_path(path)
    isfile(path) && return path
    project_path = joinpath(PROJECT_ROOT, path)
    isfile(project_path) && return project_path
    return path
end

function get_number(table, key, ::Type{T}) where {T}
    haskey(table, key) || error("Missing parameter: $key")
    return convert(T, table[key])
end

function load_config(path)
    config = TOML.parsefile(resolve_project_path(path))
    model = config["model"]
    lambda_grid = config["lambda_grid"]
    sampling = config["sampling"]
    output = config["output"]
    boundary_condition = String(get(model, "boundary_condition", "pbc"))
    boundary_condition in ("pbc", "apbc") ||
        error("boundary_condition must be \"pbc\" or \"apbc\"")

    return (
        Ninc = get_number(lambda_grid, "Ninc", Int),
        N_bins = get_number(sampling, "N_bins", Int),
        root_seeds = [convert(Int, seed) for seed in sampling["root_seeds"]],
        run_name = String(output["run_name"]),
        raw_dir = String(output["raw_dir"]),
        boundary_condition = boundary_condition,
    )
end

function raw_run_name(config)
    config.boundary_condition == "pbc" && return config.run_name
    return "$(config.run_name)__bc_$(config.boundary_condition)"
end

function parse_csv_int(text, path, row, field)
    value = tryparse(Int, text)
    value === nothing && error("Invalid $field in $path:$row")
    return value
end

function parse_csv_finite(text, path, row, field)
    value = tryparse(Float64, text)
    value === nothing && error("Invalid $field in $path:$row")
    isfinite(value) || error("Nonfinite $field in $path:$row")
    return value
end

function validate_summary_csv(path, config, lambda_index, root_seed)
    lines = readlines(path)
    expected_rows = lambda_index in (0, config.Ninc) ? 2 : 3
    length(lines) == expected_rows || error("Wrong row count in $path")
    first(lines) == SUMMARY_HEADER || error("Wrong header in $path")

    lambda_value = lambda_index / config.Ninc
    delta_lambda = 1 / config.Ninc
    directions = Set{String}()

    for (data_row, line) in enumerate(Iterators.drop(lines, 1))
        row = data_row + 1
        fields = split(line, ","; keepempty = true)
        length(fields) == 31 || error("Malformed summary row in $path:$row")
        fields[1] == raw_run_name(config) || error("Wrong run name in $path:$row")
        parse_csv_int(fields[11], path, row, "Ninc") == config.Ninc ||
            error("Wrong Ninc in $path:$row")
        parse_csv_int(fields[13], path, row, "lambda index") == lambda_index ||
            error("Wrong lambda index in $path:$row")
        parse_csv_int(fields[15], path, row, "root seed") == root_seed ||
            error("Wrong root seed in $path:$row")

        for column in 2:17
            parse_csv_finite(fields[column], path, row, "numeric field $column")
        end
        for column in 19:31
            parse_csv_finite(fields[column], path, row, "numeric field $column")
        end

        parsed_lambda = parse_csv_finite(fields[14], path, row, "lambda")
        isapprox(parsed_lambda, lambda_value; rtol = 0.0, atol = 1e-12) ||
            error("Wrong lambda value in $path:$row")

        direction = fields[18]
        direction in ("plus", "minus") || error("Wrong direction in $path:$row")
        direction in directions && error("Duplicate direction in $path:$row")
        push!(directions, direction)

        lambda_from = parse_csv_finite(fields[19], path, row, "lambda_from")
        lambda_to = parse_csv_finite(fields[20], path, row, "lambda_to")
        expected_to = direction == "plus" ?
            lambda_value + delta_lambda :
            lambda_value - delta_lambda
        isapprox(lambda_from, lambda_value; rtol = 0.0, atol = 1e-12) ||
            error("Wrong lambda_from in $path:$row")
        isapprox(lambda_to, expected_to; rtol = 0.0, atol = 1e-12) ||
            error("Wrong lambda_to in $path:$row")
    end

    expected_directions =
        lambda_index == 0 ? Set(["plus"]) :
        lambda_index == config.Ninc ? Set(["minus"]) :
        Set(["plus", "minus"])
    directions == expected_directions || error("Missing direction in $path")
    return nothing
end

function validate_bins_csv(path, config, lambda_index, root_seed)
    lines = readlines(path)
    length(lines) == config.N_bins + 1 || error("Wrong row count in $path")
    first(lines) == BINS_HEADER || error("Wrong header in $path")

    lambda_value = lambda_index / config.Ninc
    seen_bins = Set{Int}()
    for (data_row, line) in enumerate(Iterators.drop(lines, 1))
        row = data_row + 1
        fields = split(line, ","; keepempty = true)
        length(fields) == 10 || error("Malformed bins row in $path:$row")
        fields[1] == raw_run_name(config) || error("Wrong run name in $path:$row")
        parse_csv_int(fields[2], path, row, "lambda index") == lambda_index ||
            error("Wrong lambda index in $path:$row")
        parse_csv_int(fields[4], path, row, "root seed") == root_seed ||
            error("Wrong root seed in $path:$row")

        parsed_lambda = parse_csv_finite(fields[3], path, row, "lambda")
        isapprox(parsed_lambda, lambda_value; rtol = 0.0, atol = 1e-12) ||
            error("Wrong lambda value in $path:$row")
        bin = parse_csv_int(fields[5], path, row, "bin")
        bin in seen_bins && error("Duplicate bin in $path:$row")
        push!(seen_bins, bin)

        for column in 6:8
            parse_csv_finite(fields[column], path, row, "numeric field $column")
        end
        if lambda_index < config.Ninc
            parse_csv_finite(fields[9], path, row, "step_plus")
        else
            isempty(fields[9]) || error("Unexpected step_plus in $path:$row")
        end
        if lambda_index > 0
            parse_csv_finite(fields[10], path, row, "step_minus")
        else
            isempty(fields[10]) || error("Unexpected step_minus in $path:$row")
        end
    end

    seen_bins == Set(1:config.N_bins) || error("Incomplete bins in $path")
    return nothing
end

function validate_node_outputs(config, summary, bins, lambda_index, root_seed)
    validate_summary_csv(summary, config, lambda_index, root_seed)
    validate_bins_csv(bins, config, lambda_index, root_seed)
    return nothing
end

function parse_int_list(text)
    isempty(strip(text)) && return Int[]
    return [parse(Int, strip(piece)) for piece in split(text, ",")]
end

function parse_cli(args)
    param_file = args[1]
    root_seed = nothing
    all_seeds = false
    lambda_indices = nothing
    max_jobs = 4
    output_dir = nothing
    force = false

    i = 2
    while i <= length(args)
        arg = args[i]
        if arg == "--root-seed"
            i += 1
            i <= length(args) || error("--root-seed requires a value")
            root_seed = parse(Int, args[i])
        elseif arg == "--all-seeds"
            all_seeds = true
        elseif arg == "--lambda-indices"
            i += 1
            i <= length(args) || error("--lambda-indices requires a value")
            lambda_indices = parse_int_list(args[i])
        elseif arg == "--max-jobs"
            i += 1
            i <= length(args) || error("--max-jobs requires a value")
            max_jobs = parse(Int, args[i])
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

    max_jobs >= 1 || error("--max-jobs must be at least 1")

    return (
        param_file = param_file,
        root_seed = root_seed,
        all_seeds = all_seeds,
        lambda_indices = lambda_indices,
        max_jobs = max_jobs,
        output_dir = output_dir,
        force = force
    )
end

function summary_path(config, output_dir, lambda_index, root_seed)
    filename = @sprintf(
        "%s_k%03d_seed%d_summary.csv",
        config.run_name,
        lambda_index,
        root_seed
    )
    return joinpath(output_dir, filename)
end

function bins_path(config, output_dir, lambda_index, root_seed)
    filename = @sprintf(
        "%s_k%03d_seed%d_bins.csv",
        config.run_name,
        lambda_index,
        root_seed
    )
    return joinpath(output_dir, filename)
end

function log_path(log_dir, config, lambda_index, root_seed)
    filename = @sprintf(
        "%s_k%03d_seed%d.log",
        config.run_name,
        lambda_index,
        root_seed
    )
    return joinpath(log_dir, filename)
end

function build_jobs(config, cli, output_dir)
    seeds =
        cli.all_seeds ?
        config.root_seeds :
        [cli.root_seed === nothing ? first(config.root_seeds) : cli.root_seed]

    lambda_indices =
        cli.lambda_indices === nothing ?
        collect(0:config.Ninc) :
        cli.lambda_indices

    length(unique(seeds)) == length(seeds) || error("Duplicate root seeds")
    length(unique(lambda_indices)) == length(lambda_indices) ||
        error("Duplicate lambda indices")

    for k in lambda_indices
        0 <= k <= config.Ninc ||
            error("lambda index $k is outside 0:$(config.Ninc)")
    end

    jobs = []
    for seed in seeds
        for k in lambda_indices
            summary = summary_path(config, output_dir, k, seed)
            bins = bins_path(config, output_dir, k, seed)
            summary_exists = isfile(summary)
            bins_exists = isfile(bins)

            if summary_exists && bins_exists && !cli.force
                validate_node_outputs(config, summary, bins, k, seed)
                println("Skipping existing node: k=$k seed=$seed")
            elseif summary_exists != bins_exists && !cli.force
                error(
                    "Partial raw-output node detected for k=$k seed=$seed.\n" *
                    "Expected summary CSV: $summary\n" *
                    "Expected bins CSV: $bins\n" *
                    "Inspect these files manually, or rerun with --force only " *
                    "if you intentionally want to overwrite raw outputs."
                )
            else
                push!(jobs, (lambda_index = k, root_seed = seed))
            end
        end
    end

    return jobs
end

function launch_job(param_file, output_dir, log_dir, config, job; force = false)
    active_project = Base.active_project()
    active_project === nothing && error("No active Julia project")

    log_file = log_path(log_dir, config, job.lambda_index, job.root_seed)
    io = open(log_file, "w")

    cmd = `$(Base.julia_cmd()) --project=$(dirname(active_project)) $NODE_SCRIPT $param_file --lambda-index $(job.lambda_index) --root-seed $(job.root_seed) --output-dir $output_dir`
    force && (cmd = `$cmd --force`)

    println("Starting k=$(job.lambda_index) seed=$(job.root_seed)")
    println(io, "Command: ", cmd)
    flush(io)

    proc = run(pipeline(cmd, stdout = io, stderr = io), wait = false)

    return (
        job = job,
        proc = proc,
        io = io,
        log_file = log_file
    )
end

function wait_for_job!(completed, item)
    @async begin
        wait(item.proc)
        put!(completed, item)
    end
    return nothing
end

function finish_job(item)
    close(item.io)

    if success(item.proc)
        println("Finished k=$(item.job.lambda_index) seed=$(item.job.root_seed)")
        return 0
    end

    println(
        "FAILED k=$(item.job.lambda_index) seed=$(item.job.root_seed); log = ",
        item.log_file
    )
    return 1
end

function run_pool(param_file, output_dir, log_dir, config, jobs, max_jobs; force = false)
    completed = Channel{Any}(max_jobs)
    next_job = 1
    active_jobs = 0
    failures = 0

    while next_job <= length(jobs) || active_jobs > 0
        while next_job <= length(jobs) && active_jobs < max_jobs
            item = launch_job(
                param_file,
                output_dir,
                log_dir,
                config,
                jobs[next_job],
                force = force,
            )
            wait_for_job!(completed, item)
            next_job += 1
            active_jobs += 1
        end

        item = take!(completed)
        failures += finish_job(item)
        active_jobs -= 1
    end

    return failures
end

function main(args = ARGS)
    if isempty(args) || "--help" in args || "-h" in args
        usage()
        return
    end

    cli = parse_cli(args)
    config = load_config(cli.param_file)

    param_file = resolve_project_path(cli.param_file)
    output_dir =
        cli.output_dir === nothing ?
        joinpath(PROJECT_ROOT, config.raw_dir) :
        abspath(cli.output_dir)
    log_dir = joinpath(output_dir, "logs")
    mkpath(output_dir)
    mkpath(log_dir)

    jobs = build_jobs(config, cli, output_dir)

    println()
    println("========================================")
    println("Local lambda-grid launcher")
    println("========================================")
    println("parameter file = ", cli.param_file)
    println("run name       = ", config.run_name)
    println("nodes to run   = ", length(jobs))
    println("max jobs       = ", cli.max_jobs)
    println("output dir     = ", output_dir)
    println("log dir        = ", log_dir)
    println("========================================")

    isempty(jobs) && begin
        println("Nothing to run.")
        return
    end

    failures = run_pool(
        param_file,
        output_dir,
        log_dir,
        config,
        jobs,
        cli.max_jobs,
        force = cli.force,
    )

    println()
    println("========================================")
    println("Local launcher finished")
    println("failures = ", failures)
    println("========================================")

    failures == 0 || exit(1)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
