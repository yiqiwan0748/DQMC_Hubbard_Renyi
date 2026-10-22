#!/usr/bin/env julia

module U4FSSAnalysis

using LinearAlgebra
using Printf

const REQUIRED_COLUMNS = [
    "L",
    "nA",
    "beta",
    "boundary",
    "U",
    "DQMC_S2",
    "DQMC_error",
    "reference_method",
    "reference_S2",
]
const EXPECTED_SIZES = collect(4:2:22)
const BRANCH_0 = "L mod 4 = 0"
const BRANCH_2 = "L mod 4 = 2"
const DIAGNOSTIC_NOTE =
    "c_eff is a finite-window diagnostic only, not an asymptotic central-charge determination"

struct Observation
    L::Int
    dqmc_s2::Float64
    dqmc_error::Float64
    reference_method::String
    reference_s2::Float64
end

struct FitResult
    dataset::String
    method::String
    model::String
    branch::String
    p::Union{Nothing,Int}
    n_points::Int
    n_parameters::Int
    dof::Int
    L_min::Int
    L_max::Int
    condition_number::Float64
    c_eff_diagnostic::Float64
    c_eff_nominal_se::Float64
    intercept::Float64
    intercept_nominal_se::Float64
    oscillation_amplitude::Union{Nothing,Float64}
    oscillation_amplitude_nominal_se::Union{Nothing,Float64}
    statistic_name::String
    statistic_value::Float64
end

struct Options
    input::String
    fit_output::String
    plot_output::String
end

function usage(io::IO = stdout)
    println(io, "Usage:")
    println(io, "  julia analysis/analyze_u4_fss.jl --input COMPARISON.csv \\")
    println(io, "    --fit-output FITS.csv --plot-output S2_VS_L.svg")
    println(io)
    println(io, "Both outputs must be new paths; existing files are never overwritten.")
end

function parse_options(args::Vector{String})
    args == ["--help"] && (usage(); return nothing)
    values = Dict{String,String}()
    i = 1
    while i <= length(args)
        option = args[i]
        option in ("--input", "--fit-output", "--plot-output") ||
            error("Unknown argument: $option")
        haskey(values, option) && error("Duplicate option: $option")
        i == length(args) && error("Missing value for $option")
        startswith(args[i + 1], "--") && error("Missing value for $option")
        isempty(args[i + 1]) && error("Empty value for $option")
        values[option] = args[i + 1]
        i += 2
    end
    for option in ("--input", "--fit-output", "--plot-output")
        haskey(values, option) || error("Missing required option: $option")
    end
    paths = abspath.([values["--input"], values["--fit-output"], values["--plot-output"]])
    length(unique(paths)) == 3 || error("Input and output paths must be distinct")
    return Options(paths...)
end

# Small RFC-4180 reader: quoted commas, doubled quotes, CRLF, and embedded
# newlines are accepted. The campaign output itself uses only plain fields.
function parse_csv(text::String)
    rows = Vector{Vector{String}}()
    row = String[]
    field = IOBuffer()
    quoted = false
    chars = collect(text)
    i = 1
    while i <= length(chars)
        char = chars[i]
        if quoted
            if char == '"'
                if i < length(chars) && chars[i + 1] == '"'
                    write(field, '"')
                    i += 1
                else
                    quoted = false
                end
            else
                write(field, char)
            end
        elseif char == '"'
            position(field) == 0 || error("Malformed CSV: quote inside unquoted field")
            quoted = true
        elseif char == ','
            push!(row, String(take!(field)))
        elseif char == '\n'
            push!(row, String(take!(field)))
            push!(rows, row)
            row = String[]
        elseif char != '\r'
            write(field, char)
        end
        i += 1
    end
    quoted && error("Malformed CSV: unterminated quoted field")
    if position(field) > 0 || !isempty(row)
        push!(row, String(take!(field)))
        push!(rows, row)
    end
    isempty(rows) && error("Comparison CSV is empty")
    return rows
end

function finite_float(value::String, column::String, row_number::Int)
    parsed = tryparse(Float64, strip(value))
    parsed === nothing && error("Row $row_number has invalid $column: $(repr(value))")
    isfinite(parsed) || error("Row $row_number has non-finite $column")
    return parsed
end

function integer_value(value::String, column::String, row_number::Int)
    parsed = tryparse(Int, strip(value))
    parsed === nothing && error("Row $row_number has invalid integer $column: $(repr(value))")
    return parsed
end

function read_comparison(path::String)
    isfile(path) || error("Comparison CSV does not exist: $path")
    rows = parse_csv(read(path, String))
    length(rows) == 11 ||
        error("Expected one header plus 10 data rows for L=4:2:22; found $(length(rows) - 1) data rows")
    header = strip.(rows[1])
    length(unique(header)) == length(header) || error("Comparison CSV has duplicate column names")
    missing_columns = setdiff(REQUIRED_COLUMNS, header)
    isempty(missing_columns) || error("Comparison CSV is missing columns: $(join(missing_columns, ", "))")
    column = Dict(name => something(findfirst(==(name), header)) for name in REQUIRED_COLUMNS)

    observations = Observation[]
    for (offset, fields) in enumerate(rows[2:end])
        row_number = offset + 1
        length(fields) == length(header) ||
            error("Row $row_number has $(length(fields)) fields; expected $(length(header))")
        L = integer_value(fields[column["L"]], "L", row_number)
        nA = integer_value(fields[column["nA"]], "nA", row_number)
        beta = finite_float(fields[column["beta"]], "beta", row_number)
        interaction = finite_float(fields[column["U"]], "U", row_number)
        boundary = lowercase(strip(fields[column["boundary"]]))
        dqmc_s2 = finite_float(fields[column["DQMC_S2"]], "DQMC_S2", row_number)
        dqmc_error = finite_float(fields[column["DQMC_error"]], "DQMC_error", row_number)
        reference_method = uppercase(strip(fields[column["reference_method"]]))
        reference_s2 = finite_float(fields[column["reference_S2"]], "reference_S2", row_number)

        L in EXPECTED_SIZES || error("Row $row_number has unexpected L=$L")
        nA == L ÷ 2 || error("Row $row_number violates nA=L/2 for L=$L")
        isapprox(beta, 4L; atol = 1e-12, rtol = 0.0) ||
            error("Row $row_number violates beta=4L for L=$L")
        boundary == "pbc" || error("Row $row_number is not PBC")
        isapprox(interaction, 4.0; atol = 1e-12, rtol = 0.0) ||
            error("Row $row_number is not U=4")
        dqmc_error > 0.0 || error("Row $row_number has nonpositive DQMC_error")
        expected_method = L <= 12 ? "ED" : "DMRG"
        reference_method == expected_method ||
            error("Row $row_number expected $expected_method reference at L=$L, found $reference_method")

        push!(observations, Observation(L, dqmc_s2, dqmc_error, reference_method, reference_s2))
    end
    sort!(observations; by = observation -> observation.L)
    [observation.L for observation in observations] == EXPECTED_SIZES ||
        error("Comparison CSV must contain each size L=4:2:22 exactly once")
    return observations
end

function design_matrix(points::Vector{Observation}, model::String, p::Union{Nothing,Int})
    X = hcat([0.25 * log(point.L / pi) for point in points], ones(length(points)))
    if model == "oscillation"
        p in (1, 2) || error("Oscillation fit requires p=1 or p=2")
        oscillator = [(-1.0)^(point.L ÷ 2) / point.L^p for point in points]
        X = hcat(X, oscillator)
    elseif model != "total" && model != "branch"
        error("Unknown model: $model")
    end
    return Matrix{Float64}(X)
end

function fit_one(
    observations::Vector{Observation},
    dataset::String,
    model::String;
    branch::String = "all",
    p::Union{Nothing,Int} = nothing,
)
    branch_residue = branch == BRANCH_0 ? 0 : branch == BRANCH_2 ? 2 : nothing
    branch == "all" || branch_residue !== nothing || error("Unknown branch: $branch")
    points = branch == "all" ? observations :
             [point for point in observations if mod(point.L, 4) == branch_residue]
    X = design_matrix(points, model, p)
    y = dataset == "DQMC" ? [point.dqmc_s2 for point in points] :
        dataset == "reference" ? [point.reference_s2 for point in points] :
        error("Unknown dataset: $dataset")
    sigma = [point.dqmc_error for point in points]
    weighted = dataset == "DQMC"
    solve_X = weighted ? X ./ reshape(sigma, :, 1) : X
    solve_y = weighted ? y ./ sigma : y
    n_points, n_parameters = size(X)
    dof = n_points - n_parameters
    dof > 0 || error("Fit $dataset/$model/$branch has nonpositive degrees of freedom")

    decomposition = svd(solve_X)
    tolerance = maximum(size(solve_X)) * eps(Float64) * maximum(decomposition.S)
    count(>(tolerance), decomposition.S) == n_parameters ||
        error("Fit $dataset/$model/$branch is rank deficient")
    coefficients = decomposition.V *
                   ((transpose(decomposition.U) * solve_y) ./ decomposition.S)
    inverse_normal = decomposition.V * Diagonal(1.0 ./ decomposition.S .^ 2) *
                     transpose(decomposition.V)
    residual = y - X * coefficients
    if weighted
        statistic_name = "chi2_per_dof"
        statistic_value = sum(abs2, residual ./ sigma) / dof
        covariance = inverse_normal
        method = "WLS"
    else
        statistic_name = "rss_per_dof"
        statistic_value = sum(abs2, residual) / dof
        covariance = statistic_value * inverse_normal
        method = "OLS"
    end
    standard_errors = sqrt.(max.(0.0, diag(covariance)))
    all(isfinite, coefficients) && all(isfinite, standard_errors) && isfinite(statistic_value) ||
        error("Fit $dataset/$model/$branch produced a non-finite result")

    return FitResult(
        dataset,
        method,
        model,
        branch,
        p,
        n_points,
        n_parameters,
        dof,
        minimum(point.L for point in points),
        maximum(point.L for point in points),
        decomposition.S[1] / decomposition.S[end],
        coefficients[1],
        standard_errors[1],
        coefficients[2],
        standard_errors[2],
        model == "oscillation" ? coefficients[3] : nothing,
        model == "oscillation" ? standard_errors[3] : nothing,
        statistic_name,
        statistic_value,
    )
end

function run_fits(observations::Vector{Observation})
    fits = FitResult[]
    for dataset in ("DQMC", "reference")
        push!(fits, fit_one(observations, dataset, "total"))
        for branch in (BRANCH_0, BRANCH_2)
            push!(fits, fit_one(observations, dataset, "branch"; branch))
        end
        for p in (1, 2)
            push!(fits, fit_one(observations, dataset, "oscillation"; p))
        end
    end
    length(fits) == 10 || error("Internal fit-coverage error")
    return fits
end

function csv_field(value)
    text = value === nothing ? "" :
           value isa AbstractFloat ? @sprintf("%.17g", value) : string(value)
    if any(character -> character in (',', '"', '\n', '\r'), text)
        return "\"" * replace(text, "\"" => "\"\"") * "\""
    end
    return text
end

const FIT_HEADER = [
    "dataset", "fit_method", "model", "branch", "p", "n_points", "n_parameters",
    "dof", "L_min", "L_max", "condition_number", "c_eff_diagnostic",
    "c_eff_nominal_se", "intercept", "intercept_nominal_se",
    "oscillation_amplitude", "oscillation_amplitude_nominal_se", "statistic_name",
    "statistic_value", "notes",
]

function render_fit_summary(fits::Vector{FitResult})
    io = IOBuffer()
    println(io, join(FIT_HEADER, ','))
    for fit in fits
        uncertainty_note = fit.dataset == "DQMC" ?
            "DQMC WLS nominal SE uses reported MC errors only; beta/dtau systematics are excluded" :
            "reference OLS nominal SE is regression scatter, not ED precision or DMRG uncertainty"
        row = Any[
            fit.dataset, fit.method, fit.model, fit.branch, fit.p, fit.n_points,
            fit.n_parameters, fit.dof, fit.L_min, fit.L_max, fit.condition_number,
            fit.c_eff_diagnostic, fit.c_eff_nominal_se, fit.intercept,
            fit.intercept_nominal_se, fit.oscillation_amplitude,
            fit.oscillation_amplitude_nominal_se, fit.statistic_name,
            fit.statistic_value, "$DIAGNOSTIC_NOTE; $uncertainty_note",
        ]
        println(io, join(csv_field.(row), ','))
    end
    return String(take!(io))
end

function render_svg(observations::Vector{Observation})
    width, height = 960, 650
    left, right, top, bottom = 85.0, 35.0, 75.0, 85.0
    plot_width = width - left - right
    plot_height = height - top - bottom
    x_values = log.([point.L for point in observations])
    y_values = vcat(
        [point.dqmc_s2 - point.dqmc_error for point in observations],
        [point.dqmc_s2 + point.dqmc_error for point in observations],
        [point.reference_s2 for point in observations],
    )
    x_min, x_max = extrema(x_values)
    y_min, y_max = extrema(y_values)
    y_padding = max(0.05 * (y_max - y_min), 0.02)
    y_min -= y_padding
    y_max += y_padding
    xmap(value) = left + (value - x_min) / (x_max - x_min) * plot_width
    ymap(value) = top + (y_max - value) / (y_max - y_min) * plot_height

    io = IOBuffer()
    println(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"$width\" height=\"$height\" viewBox=\"0 0 $width $height\">")
    println(io, "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>")
    println(io, "<text x=\"$(width / 2)\" y=\"30\" text-anchor=\"middle\" font-family=\"Arial,sans-serif\" font-size=\"20\">U=4 PBC half-chain Renyi S2</text>")
    println(io, "<text x=\"$(width / 2)\" y=\"53\" text-anchor=\"middle\" font-family=\"Arial,sans-serif\" font-size=\"13\" fill=\"#555\">L=4:2:22, beta=4L; finite-window diagnostics only</text>")

    for tick in range(y_min, y_max; length = 6)
        y = ymap(tick)
        @printf(io, "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"#e5e5e5\"/>\n", left, y, left + plot_width, y)
        @printf(io, "<text x=\"%.2f\" y=\"%.2f\" text-anchor=\"end\" dominant-baseline=\"middle\" font-family=\"Arial,sans-serif\" font-size=\"12\">%.3f</text>\n", left - 9, y, tick)
    end
    for point in observations
        x = xmap(log(point.L))
        @printf(io, "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"#b8b8b8\"/>\n", x, top, x, top + plot_height)
        @printf(io, "<text x=\"%.2f\" y=\"%.2f\" text-anchor=\"middle\" font-family=\"Arial,sans-serif\" font-size=\"12\">%d</text>\n", x, top + plot_height + 22, point.L)
    end
    @printf(io, "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"black\" stroke-width=\"1.5\"/>\n", left, top + plot_height, left + plot_width, top + plot_height)
    @printf(io, "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"black\" stroke-width=\"1.5\"/>\n", left, top, left, top + plot_height)

    dqmc_polyline = join((@sprintf("%.2f,%.2f", xmap(log(point.L)), ymap(point.dqmc_s2)) for point in observations), " ")
    println(io, "<polyline points=\"$dqmc_polyline\" fill=\"none\" stroke=\"#1f77b4\" stroke-width=\"1.2\" opacity=\"0.55\"/>")
    for point in observations
        x = xmap(log(point.L))
        y = ymap(point.dqmc_s2)
        y_low = ymap(point.dqmc_s2 - point.dqmc_error)
        y_high = ymap(point.dqmc_s2 + point.dqmc_error)
        @printf(io, "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"#1f77b4\" stroke-width=\"1.8\"/>\n", x, y_low, x, y_high)
        @printf(io, "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"#1f77b4\" stroke-width=\"1.8\"/>\n", x - 5, y_low, x + 5, y_low)
        @printf(io, "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"#1f77b4\" stroke-width=\"1.8\"/>\n", x - 5, y_high, x + 5, y_high)
        @printf(io, "<circle cx=\"%.2f\" cy=\"%.2f\" r=\"4.5\" fill=\"white\" stroke=\"#1f77b4\" stroke-width=\"2\"/>\n", x, y)

        reference_y = ymap(point.reference_s2)
        if point.reference_method == "ED"
            @printf(io, "<rect x=\"%.2f\" y=\"%.2f\" width=\"8\" height=\"8\" fill=\"#e68613\" stroke=\"#8f4d00\"/>\n", x - 4, reference_y - 4)
        else
            coordinates = @sprintf("%.2f,%.2f %.2f,%.2f %.2f,%.2f %.2f,%.2f", x, reference_y - 5, x + 5, reference_y, x, reference_y + 5, x - 5, reference_y)
            println(io, "<polygon points=\"$coordinates\" fill=\"#d62728\" stroke=\"#7f1010\"/>")
        end
    end

    @printf(io, "<text x=\"%.2f\" y=\"%.2f\" text-anchor=\"middle\" font-family=\"Arial,sans-serif\" font-size=\"14\">L (logarithmic axis)</text>\n", left + plot_width / 2, height - 27)
    @printf(io, "<text x=\"23\" y=\"%.2f\" text-anchor=\"middle\" transform=\"rotate(-90 23 %.2f)\" font-family=\"Arial,sans-serif\" font-size=\"14\">S2</text>\n", top + plot_height / 2, top + plot_height / 2)
    legend_x, legend_y = left + 15, top + 18
    @printf(io, "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"#1f77b4\" stroke-width=\"1.8\"/><circle cx=\"%.2f\" cy=\"%.2f\" r=\"4\" fill=\"white\" stroke=\"#1f77b4\" stroke-width=\"2\"/><text x=\"%.2f\" y=\"%.2f\" dominant-baseline=\"middle\" font-family=\"Arial,sans-serif\" font-size=\"12\">DQMC +/- 1 sigma</text>\n", legend_x, legend_y, legend_x + 22, legend_y, legend_x + 11, legend_y, legend_x + 30, legend_y)
    @printf(io, "<rect x=\"%.2f\" y=\"%.2f\" width=\"8\" height=\"8\" fill=\"#e68613\"/><text x=\"%.2f\" y=\"%.2f\" dominant-baseline=\"middle\" font-family=\"Arial,sans-serif\" font-size=\"12\">ED reference</text>\n", legend_x + 7, legend_y + 14, legend_x + 30, legend_y + 18)
    coordinates = @sprintf("%.2f,%.2f %.2f,%.2f %.2f,%.2f %.2f,%.2f", legend_x + 11, legend_y + 31, legend_x + 16, legend_y + 36, legend_x + 11, legend_y + 41, legend_x + 6, legend_y + 36)
    println(io, "<polygon points=\"$coordinates\" fill=\"#d62728\"/><text x=\"$(legend_x + 30)\" y=\"$(legend_y + 36)\" dominant-baseline=\"middle\" font-family=\"Arial,sans-serif\" font-size=\"12\">DMRG reference</text>")
    println(io, "</svg>")
    return String(take!(io))
end

function stage_text(target::String, contents::String)
    directory = dirname(target)
    mkpath(directory)
    for attempt in 1:100
        temporary = joinpath(directory, ".$(basename(target)).tmp.$(getpid()).$(time_ns()).$attempt")
        try
            io = Base.Filesystem.open(
                temporary,
                Base.JL_O_WRONLY | Base.JL_O_CREAT | Base.JL_O_EXCL,
                0o644,
            )
            try
                write(io, contents)
                flush(io)
            finally
                close(io)
            end
            return temporary
        catch exception
            isfile(temporary) && rm(temporary)
            exception isa Base.IOError || rethrow()
        end
    end
    error("Could not allocate a temporary output beside $target")
end

function publish_new_outputs(outputs::Vector{Pair{String,String}})
    targets = first.(outputs)
    length(unique(targets)) == length(targets) || error("Output paths must be distinct")
    for target in targets
        ispath(target) && error("Refusing to overwrite existing output: $target")
    end
    staged = Pair{String,String}[]
    try
        for (target, contents) in outputs
            push!(staged, stage_text(target, contents) => target)
        end
        for (temporary, target) in staged
            Base.Filesystem.hardlink(temporary, target)
        end
    finally
        for (temporary, _) in staged
            isfile(temporary) && rm(temporary)
        end
    end
    return nothing
end

function analyze(input::String, fit_output::String, plot_output::String)
    input_path, fit_path, plot_path = abspath.((input, fit_output, plot_output))
    length(unique((input_path, fit_path, plot_path))) == 3 ||
        error("Input and output paths must be distinct")
    observations = read_comparison(input_path)
    fits = run_fits(observations)
    summary = render_fit_summary(fits)
    plot = render_svg(observations)
    publish_new_outputs([fit_path => summary, plot_path => plot])
    return (observations = observations, fits = fits, fit_output = fit_path, plot_output = plot_path)
end

function main(args::Vector{String} = ARGS)
    options = parse_options(args)
    options === nothing && return nothing
    result = analyze(options.input, options.fit_output, options.plot_output)
    println("PASS: validated U=4 PBC L=4:2:22 comparison")
    println("Fit summary: $(result.fit_output)")
    println("SVG plot: $(result.plot_output)")
    println(DIAGNOSTIC_NOTE)
    return nothing
end

end # module U4FSSAnalysis

if abspath(PROGRAM_FILE) == @__FILE__
    U4FSSAnalysis.main()
end
