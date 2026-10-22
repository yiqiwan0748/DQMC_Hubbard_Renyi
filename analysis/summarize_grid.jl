#!/usr/bin/env julia

using Printf
using TOML

const PROJECT_ROOT = dirname(dirname(@__FILE__))
const BINS_HEADER =
    "run_name,lambda_index,lambda,root_seed,bin,phase,den,num,step_plus,step_minus"
const EMPTY_METADATA = (
    Ninc = nothing,
    S2_ED = nothing,
    root_seeds = nothing,
    N_bins = nothing,
    run_name = nothing,
)

mutable struct NodeBins
    lambda::Float64
    phase::Vector{Float64}
    plus::Vector{Float64}
    minus::Vector{Float64}
end

function usage()
    println("""
    Summarize Hubbard Renyi lambda-node CSV output.

    Usage:
      julia analysis/summarize_grid.jl results/raw/L10_U2_beta40_Ninc8 --parameter-file parameters/L10_U2_beta40_Ninc8.toml

    Output:
      - per-step forward ratios
      - R_plus(lambda_k) * R_minus(lambda_{k+1}) consistency checks
      - final Z(1)/Z(0)
      - final S2_DQMC with propagated uncertainty

    Options:
      --parameter-file FILE   read Ninc and S2_ED from a TOML file
      --output FILE           write processed CSV summary
      --force                 allow overwriting an existing output file
      --help                  print this message
    """)
end

function resolve_project_path(path)
    isfile(path) || isdir(path) ? path : joinpath(PROJECT_ROOT, path)
end

function parse_cli(args)
    raw_dir = args[1]
    parameter_file = nothing
    output_file = nothing
    force = false

    i = 2
    while i <= length(args)
        arg = args[i]
        if arg == "--parameter-file"
            i += 1
            i <= length(args) || error("--parameter-file requires a value")
            parameter_file = args[i]
        elseif arg == "--output"
            i += 1
            i <= length(args) || error("--output requires a value")
            output_file = args[i]
        elseif arg == "--force"
            force = true
        else
            error("Unknown argument: $arg")
        end
        i += 1
    end

    return (
        raw_dir = resolve_project_path(raw_dir),
        parameter_file = parameter_file === nothing ? nothing : resolve_project_path(parameter_file),
        output_file = output_file,
        force = force
    )
end

function load_metadata(parameter_file)
    parameter_file === nothing && return EMPTY_METADATA

    config = TOML.parsefile(parameter_file)
    model = config["model"]
    grid = config["lambda_grid"]
    sampling = config["sampling"]
    output = config["output"]
    boundary_condition = String(get(model, "boundary_condition", "pbc"))
    boundary_condition in ("pbc", "apbc") ||
        error("boundary_condition must be \"pbc\" or \"apbc\"")
    run_name = String(output["run_name"])
    raw_run_name = boundary_condition == "pbc" ?
        run_name : "$(run_name)__bc_$(boundary_condition)"
    S2_ED = convert(Float64, grid["S2_ED"])
    return (
        Ninc = convert(Int, grid["Ninc"]),
        S2_ED = isfinite(S2_ED) ? S2_ED : nothing,
        root_seeds = [convert(Int, seed) for seed in sampling["root_seeds"]],
        N_bins = convert(Int, sampling["N_bins"]),
        run_name = raw_run_name,
    )
end

function parse_finite(text, field, file, line_number)
    isempty(text) && error("Empty $field in $file:$line_number")
    value = parse(Float64, text)
    isfinite(value) || error("Nonfinite $field in $file:$line_number")
    return value
end

function parse_optional_finite(text, field, file, line_number)
    isempty(text) && return NaN
    return parse_finite(text, field, file, line_number)
end

function get_node!(nodes, k, lambda_value)
    if !haskey(nodes, k)
        nodes[k] = NodeBins(lambda_value, Float64[], Float64[], Float64[])
    elseif !isapprox(nodes[k].lambda, lambda_value; rtol = 0.0, atol = 1e-12)
        error("Inconsistent lambda values for index $k")
    end
    return nodes[k]
end

function read_bins(raw_dir, metadata = EMPTY_METADATA)
    isdir(raw_dir) || error("Raw result directory not found: $raw_dir")

    files = sort([
        joinpath(raw_dir, name)
        for name in readdir(raw_dir)
        if endswith(name, "_bins.csv") && isfile(joinpath(raw_dir, name))
    ])

    isempty(files) && error("No *_bins.csv files found in $raw_dir")

    nodes = Dict{Int, NodeBins}()
    seen_rows = Set{Tuple{Int, Int, Int}}()
    node_files = Dict{Tuple{Int, Int}, String}()
    bins_by_node = Dict{Tuple{Int, Int}, Set{Int}}()
    observed_run_name = nothing

    for file in files
        open(file, "r") do io
            header = readline(io)
            header == BINS_HEADER || error("Unexpected CSV header in $file")

            file_node = nothing
            row_count = 0
            for (data_line, line) in enumerate(eachline(io))
                line_number = data_line + 1
                row_count += 1
                fields = split(line, ","; keepempty = true)
                length(fields) == 10 || error("Malformed row in $file:$line_number")

                run_name = fields[1]
                isempty(run_name) && error("Empty run_name in $file:$line_number")
                k = parse(Int, fields[2])
                lambda_value = parse_finite(fields[3], "lambda", file, line_number)
                root_seed = parse(Int, fields[4])
                bin = parse(Int, fields[5])
                phase = parse_finite(fields[6], "phase", file, line_number)
                parse_finite(fields[7], "den", file, line_number)
                parse_finite(fields[8], "num", file, line_number)
                plus = parse_optional_finite(
                    fields[9],
                    "step_plus",
                    file,
                    line_number,
                )
                minus = parse_optional_finite(
                    fields[10],
                    "step_minus",
                    file,
                    line_number,
                )

                k >= 0 || error("Negative lambda index in $file:$line_number")
                bin >= 1 || error("Invalid bin number in $file:$line_number")
                if metadata.Ninc !== nothing
                    k <= metadata.Ninc ||
                        error("Lambda index $k exceeds Ninc=$(metadata.Ninc)")
                end
                if metadata.N_bins !== nothing
                    bin <= metadata.N_bins ||
                        error("Bin $bin exceeds N_bins=$(metadata.N_bins)")
                end
                if metadata.run_name !== nothing
                    run_name == metadata.run_name ||
                        error("Unexpected run_name '$run_name' in $file:$line_number")
                elseif observed_run_name !== nothing
                    run_name == observed_run_name ||
                        error("Mixed run names in $raw_dir")
                end
                observed_run_name = run_name

                node_key = (k, root_seed)
                file_node === nothing || file_node == node_key ||
                    error("Mixed lambda nodes in $file")
                file_node = node_key

                if haskey(node_files, node_key) && node_files[node_key] != file
                    error("Multiple bin files for k=$k seed=$root_seed")
                end
                node_files[node_key] = file

                row_key = (k, root_seed, bin)
                row_key in seen_rows &&
                    error("Duplicate row for k=$k seed=$root_seed bin=$bin")
                push!(seen_rows, row_key)
                push!(get!(bins_by_node, node_key, Set{Int}()), bin)

                node = get_node!(nodes, k, lambda_value)
                push!(node.phase, phase)
                push!(node.plus, plus)
                push!(node.minus, minus)
            end
            row_count > 0 || error("No data rows in $file")
        end
    end

    Ninc = metadata.Ninc === nothing ? maximum(keys(nodes)) : metadata.Ninc
    for (k, node) in nodes
        if k < Ninc
            all(isfinite, node.plus) || error("Missing step_plus values for lambda index $k")
        else
            all(isnan, node.plus) || error("Unexpected step_plus values for lambda index $k")
        end
        if k > 0
            all(isfinite, node.minus) || error("Missing step_minus values for lambda index $k")
        else
            all(isnan, node.minus) || error("Unexpected step_minus values for lambda index $k")
        end
    end

    if metadata.root_seeds !== nothing
        expected_nodes = Set(
            (k, seed) for k in 0:Ninc for seed in metadata.root_seeds
        )
        actual_nodes = Set(keys(node_files))
        missing_nodes = sort!(collect(setdiff(expected_nodes, actual_nodes)))
        unexpected_nodes = sort!(collect(setdiff(actual_nodes, expected_nodes)))
        isempty(missing_nodes) || error("Missing bin files for nodes: $missing_nodes")
        isempty(unexpected_nodes) ||
            error("Unexpected bin files for nodes: $unexpected_nodes")
    end

    if metadata.N_bins !== nothing
        expected_bins = Set(1:metadata.N_bins)
        for node_key in sort!(collect(keys(bins_by_node)))
            bins_by_node[node_key] == expected_bins ||
                error("Incomplete bin coverage for node $node_key")
        end
    end

    return nodes
end

function jackknife(values, estimator)
    n = length(values)
    n >= 2 || error("Need at least two bins for jackknife")

    estimates = Vector{Float64}(undef, n)
    for leave_out in 1:n
        sample = [values[i] for i in 1:n if i != leave_out]
        estimates[leave_out] = estimator(sample)
    end

    avg = sum(estimates) / n
    err = sqrt((n - 1) / n * sum((x - avg)^2 for x in estimates))

    return avg, err
end

function valid_pairs(weights, phases)
    pairs = [
        (weight = weights[i], phase = phases[i])
        for i in eachindex(weights)
        if isfinite(weights[i])
    ]
    return pairs
end

function ratio_from_pairs(pairs)
    numerator = sum(pair.weight for pair in pairs) / length(pairs)
    denominator = sum(pair.phase for pair in pairs) / length(pairs)
    denominator == 0.0 && error("Mean phase is zero")
    ratio = real(numerator / denominator)
    isfinite(ratio) || error("Nonfinite ratio")
    ratio == 0.0 && error("Mean ratio is zero")
    return ratio
end

function entropy_from_pairs(pairs)
    return -log(abs(ratio_from_pairs(pairs)))
end

function ratio_statistics(weights, phases)
    pairs = valid_pairs(weights, phases)
    isempty(pairs) && return nothing

    ratio, ratio_err = jackknife(pairs, ratio_from_pairs)
    entropy, entropy_err = jackknife(pairs, entropy_from_pairs)

    return (
        ratio = ratio,
        ratio_err = ratio_err,
        entropy = entropy,
        entropy_err = entropy_err,
        n_bins = length(pairs)
    )
end

function product_with_error(a, da, b, db)
    product = a * b
    err = abs(product) * sqrt((da / a)^2 + (db / b)^2)
    return product, err
end

function print_forward_table(plus_stats)
    println()
    println("Forward adjacent ratios")
    println("k    R_plus +/- err              entropy +/- err")
    for k in sort(collect(keys(plus_stats)))
        stats = plus_stats[k]
        println(
            @sprintf(
                "%-4d %.12f +/- %.12f   %.12f +/- %.12f",
                k,
                stats.ratio,
                stats.ratio_err,
                stats.entropy,
                stats.entropy_err
            )
        )
    end
end

function print_consistency_table(plus_stats, minus_stats)
    println()
    println("Bidirectional consistency checks")
    println("check                         product +/- err        z")

    for k in sort(collect(keys(plus_stats)))
        haskey(minus_stats, k + 1) || continue

        plus = plus_stats[k]
        minus = minus_stats[k + 1]
        product, err = product_with_error(
            plus.ratio,
            plus.ratio_err,
            minus.ratio,
            minus.ratio_err
        )
        z = (product - 1.0) / err

        println(
            @sprintf(
                "R_plus(%d) * R_minus(%d)    %.12f +/- %.12f   %+.3f",
                k,
                k + 1,
                product,
                err,
                z
            )
        )
    end
end

function summarize(raw_dir, metadata)
    nodes = read_bins(raw_dir, metadata)
    inferred_Ninc = maximum(keys(nodes))
    Ninc = metadata.Ninc === nothing ? inferred_Ninc : metadata.Ninc

    plus_stats = Dict{Int, Any}()
    minus_stats = Dict{Int, Any}()

    for k in sort(collect(keys(nodes)))
        node = nodes[k]
        plus = ratio_statistics(node.plus, node.phase)
        minus = ratio_statistics(node.minus, node.phase)

        plus === nothing || (plus_stats[k] = plus)
        minus === nothing || (minus_stats[k] = minus)
    end

    required_plus = collect(0:(Ninc - 1))
    missing_plus = [k for k in required_plus if !haskey(plus_stats, k)]
    isempty(missing_plus) ||
        error("Missing forward ratios for lambda indices: $(missing_plus)")

    Z1_over_Z0 = prod(plus_stats[k].ratio for k in required_plus)
    rel_var = sum(
        (plus_stats[k].ratio_err / plus_stats[k].ratio)^2
        for k in required_plus
    )
    Z1_over_Z0_err = abs(Z1_over_Z0) * sqrt(rel_var)

    S2 = sum(plus_stats[k].entropy for k in required_plus)
    S2_err = sqrt(sum(plus_stats[k].entropy_err^2 for k in required_plus))

    print_forward_table(plus_stats)
    print_consistency_table(plus_stats, minus_stats)

    println()
    println("Complete lambda-grid result")
    println(@sprintf("Z1/Z0   = %.12f +/- %.12f", Z1_over_Z0, Z1_over_Z0_err))
    println(@sprintf("S2_DQMC = %.12f +/- %.12f", S2, S2_err))

    if metadata.S2_ED !== nothing
        delta = S2 - metadata.S2_ED
        z = delta / S2_err
        println(@sprintf("S2_ED   = %.12f", metadata.S2_ED))
        println(@sprintf("DeltaS2 = %+.12f", delta))
        println(@sprintf("z       = %+.3f", z))
    end

    return (
        plus_stats = plus_stats,
        minus_stats = minus_stats,
        Z1_over_Z0 = Z1_over_Z0,
        Z1_over_Z0_err = Z1_over_Z0_err,
        S2 = S2,
        S2_err = S2_err,
        Ninc = Ninc
    )
end

function write_processed_summary(path, result, metadata; force = false)
    if isfile(path) && !force
        error("Output file already exists: $path. Use --force to overwrite.")
    end

    mkpath(dirname(path))

    open(path, "w") do io
        println(io, "quantity,value,error")
        println(io, @sprintf("Z1_over_Z0,%.12e,%.12e", result.Z1_over_Z0, result.Z1_over_Z0_err))
        println(io, @sprintf("S2_DQMC,%.12e,%.12e", result.S2, result.S2_err))
        if metadata.S2_ED !== nothing
            println(io, @sprintf("S2_ED,%.12e,", metadata.S2_ED))
            println(io, @sprintf("DeltaS2,%.12e,", result.S2 - metadata.S2_ED))
        end
    end

    return nothing
end

function main(args = ARGS)
    if isempty(args) || "--help" in args || "-h" in args
        usage()
        return
    end

    cli = parse_cli(args)
    metadata = load_metadata(cli.parameter_file)
    result = summarize(cli.raw_dir, metadata)

    output_file =
        cli.output_file === nothing ?
        joinpath(PROJECT_ROOT, "results", "processed", "lambda_grid_summary.csv") :
        resolve_project_path(cli.output_file)

    write_processed_summary(output_file, result, metadata; force = cli.force)
    println()
    println("Wrote processed summary:")
    println(output_file)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
