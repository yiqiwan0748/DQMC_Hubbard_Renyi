#!/usr/bin/env julia

using ITensors
using ITensorMPS
using Printf
using Random

const DEFAULT_BOUNDARY = "pbc"
const DEFAULT_SEED = 1234

function usage()
    println("""
    Compute a DMRG/MPS ground-state S2 benchmark for the 1D Hubbard chain.

    Usage:
      julia --project=. analysis/compute_s2_itensor_dmrg.jl --L 8 --U 2.0 --nA 4 --boundary pbc --seed 1234

    Notes:
      - This is a DMRG/MPS benchmark, not ED.
      - Use --boundary pbc or --boundary apbc to match ED/DQMC.
      - The legacy --periodic true|false option remains supported.
      - The subsystem is C = 1:nA, matching the current ED/DQMC convention.

    Options:
      --seed INT          initial-MPS seed; default $DEFAULT_SEED
      --nsweeps INT       number of DMRG sweeps; default 14
      --cutoff FLOAT      truncation cutoff; default 1e-12
      --maxdim-cap INT    final bond-dimension cap; default 4800
      --maxdim-schedule LIST
                          explicit comma-separated schedule, one positive
                          nondecreasing value per sweep; overrides the default
                          ramp (cannot be combined with --maxdim-cap)
      --output PATH       atomically write the final CSV row to PATH
    """)
end

function parse_maxdim_schedule(text)
    stripped = strip(text)
    isempty(stripped) && error("--maxdim-schedule cannot be empty")
    pieces = split(stripped, ","; keepempty = true)
    any(piece -> isempty(strip(piece)), pieces) &&
        error("--maxdim-schedule contains an empty value: $text")
    schedule = [parse(Int, strip(piece)) for piece in pieces]
    all(>(0), schedule) ||
        error("--maxdim-schedule values must be positive")
    all(schedule[i] <= schedule[i + 1] for i in 1:(length(schedule) - 1)) ||
        error("--maxdim-schedule must be nondecreasing")
    return schedule
end

function parse_boundary(text; allow_obc = false)
    boundary = lowercase(strip(text))
    allowed = allow_obc ? ("pbc", "apbc", "obc") : ("pbc", "apbc")
    choices = allow_obc ? "pbc, apbc, or obc" : "pbc or apbc"
    boundary in allowed || error("--boundary must be $choices, got: $text")
    return boundary
end

function option_value(args, index, option)
    index <= length(args) || error("$option requires a value")
    return args[index]
end

function parse_cli(args)
    L = 8
    t = 1.0
    U = 2.0
    mu = 0.0
    nA = 4
    boundary = DEFAULT_BOUNDARY
    boundary_was_set = false
    periodic_was_set = false
    seed = DEFAULT_SEED
    nsweeps = 14
    cutoff = 1e-12
    maxdim_cap = 4800
    maxdim_cap_was_set = false
    maxdim_schedule = nothing
    output = nothing

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--help" || arg == "-h"
            usage()
            exit(0)
        elseif arg == "--L"
            i += 1
            L = parse(Int, option_value(args, i, arg))
        elseif arg == "--t"
            i += 1
            t = parse(Float64, option_value(args, i, arg))
        elseif arg == "--U"
            i += 1
            U = parse(Float64, option_value(args, i, arg))
        elseif arg == "--mu"
            i += 1
            mu = parse(Float64, option_value(args, i, arg))
        elseif arg == "--nA"
            i += 1
            nA = parse(Int, option_value(args, i, arg))
        elseif arg == "--boundary"
            i += 1
            boundary = parse_boundary(option_value(args, i, arg))
            boundary_was_set = true
        elseif arg == "--periodic"
            i += 1
            periodic = parse(Bool, option_value(args, i, arg))
            boundary = periodic ? "pbc" : "obc"
            periodic_was_set = true
        elseif arg == "--seed"
            i += 1
            seed = parse(Int, option_value(args, i, arg))
        elseif arg == "--nsweeps"
            i += 1
            nsweeps = parse(Int, option_value(args, i, arg))
        elseif arg == "--cutoff"
            i += 1
            cutoff = parse(Float64, option_value(args, i, arg))
        elseif arg == "--maxdim-cap"
            i += 1
            maxdim_cap = parse(Int, option_value(args, i, arg))
            maxdim_cap_was_set = true
        elseif arg == "--maxdim-schedule"
            i += 1
            maxdim_schedule = parse_maxdim_schedule(option_value(args, i, arg))
        elseif arg == "--output"
            i += 1
            output = option_value(args, i, arg)
        else
            error("Unknown argument: $arg")
        end
        i += 1
    end

    boundary_was_set && periodic_was_set &&
        error("Use either --boundary or --periodic, not both")
    nsweeps >= 1 || error("--nsweeps must be positive")
    isfinite(cutoff) && cutoff > 0 || error("--cutoff must be finite and positive")
    maxdim_cap_was_set && maxdim_schedule !== nothing &&
        error("Use either --maxdim-cap or --maxdim-schedule, not both")
    if maxdim_schedule === nothing
        maxdim_cap >= 4800 || error("--maxdim-cap must be at least 4800")
    else
        length(maxdim_schedule) == nsweeps || error(
            "--maxdim-schedule must contain exactly --nsweeps=$nsweeps values, " *
            "got $(length(maxdim_schedule))",
        )
        maxdim_cap = last(maxdim_schedule)
    end

    return (
        L = L,
        t = t,
        U = U,
        mu = mu,
        nA = nA,
        boundary = boundary,
        periodic = boundary != "obc",
        seed = seed,
        nsweeps = nsweeps,
        cutoff = cutoff,
        maxdim_cap = maxdim_cap,
        maxdim_schedule = maxdim_schedule,
        output = output
    )
end

function effective_maxdim_schedule(cfg)
    if cfg.maxdim_schedule !== nothing
        return cfg.maxdim_schedule
    end

    ramp = [
        20, 50, 100, 200, 400, 800, 1200, 1600, 2400, 3200,
        cfg.maxdim_cap,
    ]
    if cfg.nsweeps <= length(ramp)
        return ramp[1:cfg.nsweeps]
    end
    return vcat(ramp, fill(cfg.maxdim_cap, cfg.nsweeps - length(ramp)))
end

function hubbard_mpo(
    sites;
    L::Int,
    t::Float64,
    U::Float64,
    mu::Float64,
    boundary::AbstractString = DEFAULT_BOUNDARY
)
    boundary = parse_boundary(boundary; allow_obc = true)

    os = OpSum()

    for j in 1:(L - 1)
        os += -t, "Cdagup", j, "Cup", j + 1
        os += -t, "Cdagup", j + 1, "Cup", j

        os += -t, "Cdagdn", j, "Cdn", j + 1
        os += -t, "Cdagdn", j + 1, "Cdn", j
    end

    if boundary != "obc"
        seam_hopping = boundary == "pbc" ? -t : t
        os += seam_hopping, "Cdagup", L, "Cup", 1
        os += seam_hopping, "Cdagup", 1, "Cup", L

        os += seam_hopping, "Cdagdn", L, "Cdn", 1
        os += seam_hopping, "Cdagdn", 1, "Cdn", L
    end

    for j in 1:L
        os += U, "Nupdn", j
        if mu != 0.0
            os += -mu, "Nup", j
            os += -mu, "Ndn", j
        end
    end

    return MPO(os, sites)
end

function half_filled_sz0_state(L::Int)
    return [isodd(j) ? "Up" : "Dn" for j in 1:L]
end

function renyi2_bond_entropy(psi::MPS, b::Int)
    orthogonalize!(psi, b)

    sites = siteinds(psi)
    left_inds = b == 1 ? (sites[b],) : (linkind(psi, b - 1), sites[b])

    _, S, _ = svd(psi[b], left_inds)

    trrho2 = 0.0
    normcheck = 0.0

    for n in 1:dim(S, 1)
        s = real(S[n, n])
        p = s^2
        normcheck += p
        trrho2 += p^2
    end

    trrho2 /= normcheck^2

    return -log(trrho2), normcheck
end

function run_case(cfg)
    sites = siteinds("Electron", cfg.L; conserve_qns = true)

    H = hubbard_mpo(
        sites;
        L = cfg.L,
        t = cfg.t,
        U = cfg.U,
        mu = cfg.mu,
        boundary = cfg.boundary
    )

    state = half_filled_sz0_state(cfg.L)
    Random.seed!(cfg.seed)
    psi0 = random_mps(sites, state; linkdims = 20)

    maxdim = effective_maxdim_schedule(cfg)
    observer = DMRGObserver()

    energy, psi = dmrg(
        H,
        psi0;
        nsweeps = cfg.nsweeps,
        maxdim = maxdim,
        cutoff = cfg.cutoff,
        outputlevel = 1,
        observer = observer,
    )
    sweep_truncerrors = truncerrors(observer)
    isempty(sweep_truncerrors) &&
        error("DMRG observer returned no sweep truncation errors")
    final_sweep_max_truncerr = last(sweep_truncerrors)

    S2, normcheck = renyi2_bond_entropy(psi, cfg.nA)

    return energy, S2, normcheck, maxlinkdim(psi), final_sweep_max_truncerr
end

function write_result(
    io,
    cfg,
    energy,
    S2,
    normcheck,
    maxlink,
    final_sweep_max_truncerr,
)
    println(
        io,
        "L,t,U,mu,nA,periodic,boundary,seed,energy,S2_MPS,normcheck,maxlinkdim," *
        "nsweeps,cutoff,maxdim_cap,maxdim_schedule,energy_density," *
        "final_sweep_max_truncerr",
    )
    schedule_text = join(effective_maxdim_schedule(cfg), ";")
    @printf(
        io,
        "%d,%.8f,%.8f,%.8f,%d,%s,%s,%d,%.16e,%.16e,%.16e,%d,%d,%.16e,%d,%s,%.16e,%.16e\n",
        cfg.L,
        cfg.t,
        cfg.U,
        cfg.mu,
        cfg.nA,
        string(cfg.periodic),
        cfg.boundary,
        cfg.seed,
        energy,
        S2,
        normcheck,
        maxlink,
        cfg.nsweeps,
        cfg.cutoff,
        cfg.maxdim_cap,
        schedule_text,
        energy / cfg.L,
        final_sweep_max_truncerr,
    )
end

function write_result_atomic(
    path,
    cfg,
    energy,
    S2,
    normcheck,
    maxlink,
    final_sweep_max_truncerr,
)
    output_path = abspath(path)
    (ispath(output_path) || islink(output_path)) &&
        error("Output path already exists: $output_path")
    mkpath(dirname(output_path))
    temp_path, temp_io = mktemp(dirname(output_path))
    try
        write_result(
            temp_io,
            cfg,
            energy,
            S2,
            normcheck,
            maxlink,
            final_sweep_max_truncerr,
        )
        close(temp_io)
        mv(temp_path, output_path)
    finally
        isopen(temp_io) && close(temp_io)
        isfile(temp_path) && rm(temp_path; force = true)
    end
end

function main(args = ARGS)
    cfg = parse_cli(args)

    @assert iseven(cfg.L)
    @assert 1 <= cfg.nA < cfg.L

    energy, S2, normcheck, maxlink, final_sweep_max_truncerr = run_case(cfg)

    if cfg.output === nothing
        write_result(
            stdout,
            cfg,
            energy,
            S2,
            normcheck,
            maxlink,
            final_sweep_max_truncerr,
        )
    else
        write_result_atomic(
            cfg.output,
            cfg,
            energy,
            S2,
            normcheck,
            maxlink,
            final_sweep_max_truncerr,
        )
        println("Wrote DMRG S2 reference CSV: ", abspath(cfg.output))
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
