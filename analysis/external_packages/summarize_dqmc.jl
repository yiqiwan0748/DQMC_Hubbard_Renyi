#!/usr/bin/env julia
# Read only the newly generated CSV evidence; refuse to replace a prior report.
include(joinpath(@__DIR__, "dqmc_common.jl"))

function read_bins(path)
    rows = NamedTuple[]
    for line in Iterators.drop(eachline(path), 1)
        fields = split(line, ',')
        push!(rows, (; engine = fields[1], U = parse(Float64,fields[2]),
            seed = parse(Int,fields[3]), bin = parse(Int,fields[4]),
            sign = parse(Float64,fields[5]), energy_per_site = parse(Float64,fields[6]),
            density = parse(Float64,fields[7]), double_occupancy = parse(Float64,fields[8])))
    end
    return rows
end

function pool(rows, U, name; merge_bins = 1)
    selected = filter(row -> row.U == U, rows)
    signs, values = Float64[], Float64[]
    seed_means = Float64[]
    for seed in sort(unique(row.seed for row in selected))
        seedrows = filter(row -> row.seed == seed, selected)
        length(seedrows) % merge_bins == 0 || error("Rebinning must divide the original bins")
        for start in 1:merge_bins:length(seedrows)
            section = seedrows[start:start+merge_bins-1]
            push!(signs, mean(row.sign for row in section))
            push!(values, mean(getproperty(row,name) for row in section))
        end
        push!(seed_means, sum(getproperty(row,name) for row in seedrows) / sum(row.sign for row in seedrows))
    end
    estimate, error = ratio_stats(values, signs)
    seed_error = std(seed_means) / sqrt(length(seed_means))
    return (; estimate, error, seed_error, average_sign = mean(signs), bins = length(signs))
end

function main(args = ARGS)
    length(args) == 1 || error("Usage: julia summarize_dqmc.jl EVIDENCE_DIRECTORY")
    output = joinpath(args[1], "comparison.csv")
    ispath(output) && error("Refusing to replace previous comparison")
    a = read_bins(joinpath(args[1], "local_bins.csv"))
    b = read_bins(joinpath(args[1], "smoqy_bins.csv"))
    all_pass = true
    open(output, "w") do io
        println(io, "U,observable,local_mean,local_standard_error,smoqy_mean,smoqy_standard_error,difference,combined_error,z,absolute_tolerance,passed,local_seed_standard_error,smoqy_seed_standard_error,local_rebin2_error,smoqy_rebin2_error,average_sign_local,average_sign_smoqy,exact_U0")
        for U in (0.0,2.0), name in (:energy_per_site,:density,:double_occupancy)
            local_stats, smoqy_stats = pool(a,U,name), pool(b,U,name)
            local_rebin, smoqy_rebin = pool(a,U,name;merge_bins=2), pool(b,U,name;merge_bins=2)
            difference = local_stats.estimate - smoqy_stats.estimate
            combined = hypot(local_stats.error,smoqy_stats.error)
            tolerance = max(5combined, 1e-10)
            passed = abs(difference) <= tolerance
            exact = U == 0 ? (name == :energy_per_site ? -tanh(2.0) :
                name == :density ? 1.0 : 0.25) : NaN
            U == 0 && (passed &= max(abs(local_stats.estimate-exact),abs(smoqy_stats.estimate-exact)) <= 1e-10)
            all_pass &= passed
            z = combined > 1e-12 ? difference/combined : NaN
            println(io, join((U,name,local_stats.estimate,local_stats.error,smoqy_stats.estimate,
                smoqy_stats.error,difference,combined,z,tolerance,passed,
                local_stats.seed_error,smoqy_stats.seed_error,local_rebin.error,smoqy_rebin.error,
                local_stats.average_sign,smoqy_stats.average_sign,exact), ','))
            println("U=$U $name: local=$(local_stats.estimate) +/- $(local_stats.error); SmoQy=$(smoqy_stats.estimate) +/- $(smoqy_stats.error); z=$z; pass=$passed")
        end
    end
    all_pass || error("One or more package comparisons failed")
end
abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
