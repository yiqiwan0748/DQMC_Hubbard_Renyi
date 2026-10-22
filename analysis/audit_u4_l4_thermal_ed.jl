#!/usr/bin/env julia

module U4L4ThermalEDAudit

using LinearAlgebra
using Printf

module GroundStateED
include(joinpath(@__DIR__, "compute_s2_ed.jl"))
end

const DEFAULT_L = 4
const DEFAULT_T = 1.0
const DEFAULT_U = 4.0
const DEFAULT_MU = 0.0
const DEFAULT_BETA = 16.0
const DEFAULT_DTAU = 0.05
const DEFAULT_NA = 2
const DEFAULT_BOUNDARY = "pbc"
const DEFAULT_DQMC_S2 = 0.9550849355244
const DEFAULT_DQMC_ERROR = 0.01192389684662

function usage()
    println("""
    Exact L=4 audit of the U=4 DQMC Renyi target.

    Usage:
      julia analysis/audit_u4_l4_thermal_ed.jl [--output PATH]

    Options:
      --output PATH    write a new CSV; otherwise print it
      --help           print this message

    Fixed production point:
      L=4, t=1, U=4, mu=0, beta=16, dtau=0.05, nA=2, PBC

    The thermal calculation is grand canonical and uses the particle-hole
    symmetric interaction U(n_up-1/2)(n_down-1/2). The Trotter result uses
    the production ordering exp(-dtau*V) exp(-dtau*K).
    """)
end

function parse_cli(args)
    output = nothing
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--help" || arg == "-h"
            usage()
            return nothing
        elseif arg == "--output"
            i += 1
            i <= length(args) || error("--output requires a path")
            output = args[i]
        else
            error("Unknown argument: $arg")
        end
        i += 1
    end
    return (; output)
end

function full_fock_matrices(L, t, U, mu; boundary = DEFAULT_BOUNDARY)
    spin_dimension = 1 << L
    masks = collect(UInt(0):UInt(spin_dimension - 1))
    index = GroundStateED.mask_index(masks)
    hops = GroundStateED.build_spin_hops(masks, index, L, t; boundary)
    dimension = spin_dimension^2

    kinetic = zeros(Float64, dimension, dimension)
    potential = zeros(Float64, dimension)

    for iu in eachindex(masks), idn in eachindex(masks)
        source = (iu - 1) * spin_dimension + idn
        up = masks[iu]
        down = masks[idn]

        for (iu_to, amplitude) in hops[iu]
            destination = (iu_to - 1) * spin_dimension + idn
            kinetic[destination, source] += amplitude
        end
        for (idn_to, amplitude) in hops[idn]
            destination = (iu - 1) * spin_dimension + idn_to
            kinetic[destination, source] += amplitude
        end

        particles = count_ones(up) + count_ones(down)
        double_occupancy = count_ones(up & down)
        potential[source] = U * (
            double_occupancy - particles / 2 + L / 4
        ) - mu * particles
    end

    ishermitian(kinetic) || error("Many-body kinetic matrix is not Hermitian")
    return kinetic, potential, masks
end

function subsystem_labels(masks, L, nA)
    sitesA = collect(0:(nA - 1))
    sitesB = collect(nA:(L - 1))
    spin_dimension = length(masks)
    dimension = spin_dimension^2
    a_label = Vector{Int}(undef, dimension)
    b_label = Vector{Int}(undef, dimension)
    signs = Vector{Float64}(undef, dimension)

    for iu in eachindex(masks), idn in eachindex(masks)
        state = (iu - 1) * spin_dimension + idn
        upA = GroundStateED.compressed_bits(masks[iu], sitesA)
        upB = GroundStateED.compressed_bits(masks[iu], sitesB)
        dnA = GroundStateED.compressed_bits(masks[idn], sitesA)
        dnB = GroundStateED.compressed_bits(masks[idn], sitesB)
        a_label[state] = Int(upA | (dnA << nA)) + 1
        b_label[state] = Int(upB | (dnB << (L - nA))) + 1
        signs[state] = GroundStateED.subsystem_sign(upA, upB, dnA)
    end

    return a_label, b_label, signs
end

function reduced_operator(rho, masks, L, nA)
    size(rho, 1) == size(rho, 2) || error("rho must be square")
    size(rho, 1) == length(masks)^2 || error("rho has the wrong dimension")
    1 <= nA <= L || error("nA must satisfy 1 <= nA <= L")

    a_label, b_label, signs = subsystem_labels(masks, L, nA)
    dimension_A = 1 << (2 * nA)
    dimension_B = 1 << (2 * (L - nA))
    states_by_B = [Int[] for _ in 1:dimension_B]
    for state in eachindex(b_label)
        push!(states_by_B[b_label[state]], state)
    end

    rhoA = zeros(eltype(rho), dimension_A, dimension_A)
    for states in states_by_B
        for source in states, destination in states
            rhoA[a_label[source], a_label[destination]] +=
                signs[source] * signs[destination] * rho[source, destination]
        end
    end
    return rhoA
end

function renyi2(rho, masks, L, nA)
    rhoA = reduced_operator(rho, masks, L, nA)
    trace_error = abs(tr(rhoA) - 1)
    trace_error <= 5e-11 || error("Reduced density trace error: $trace_error")
    purity = real(tr(rhoA * rhoA))
    isfinite(purity) && purity > 0 || error("Invalid reduced purity: $purity")
    return -log(purity), purity, rhoA
end

function exact_thermal_density(hamiltonian, beta)
    decomposition = eigen(Symmetric(hamiltonian))
    ground_energy = first(decomposition.values)
    weights = exp.(-beta .* (decomposition.values .- ground_energy))
    probabilities = weights ./ sum(weights)
    rho = decomposition.vectors * Diagonal(probabilities) * decomposition.vectors'
    return rho, decomposition
end

function asymmetric_trotter_density(kinetic, potential, beta, dtau)
    slices = round(Int, beta / dtau)
    isapprox(slices * dtau, beta; atol = 1e-12, rtol = 0.0) ||
        error("beta/dtau must be an integer")

    expm_half_kinetic = exp(-dtau .* kinetic ./ 2)
    expp_half_kinetic = exp(+dtau .* kinetic ./ 2)
    expm_potential = Diagonal(exp.(-dtau .* potential))

    # The production slice is T = exp(-dtau*V) exp(-dtau*K). It is similar to
    # the positive symmetric matrix S below:
    # T = exp(+dtau*K/2) S exp(-dtau*K/2).
    symmetric_slice = Symmetric(
        expm_half_kinetic * expm_potential * expm_half_kinetic,
    )
    decomposition = eigen(symmetric_slice)
    minimum(decomposition.values) > 0 || error("Trotter slice is not positive")
    log_weights = slices .* log.(decomposition.values)
    log_weights .-= maximum(log_weights)
    probabilities = exp.(log_weights)
    probabilities ./= sum(probabilities)
    symmetric_density =
        decomposition.vectors * Diagonal(probabilities) * decomposition.vectors'
    rho = expp_half_kinetic * symmetric_density * expm_half_kinetic
    rho ./= tr(rho)
    return rho, slices
end

function fixed_sector_ground_state(L, t, U, mu, nA; boundary = DEFAULT_BOUNDARY)
    space = GroundStateED.build_space(
        L,
        t,
        U,
        mu,
        L ÷ 2,
        L ÷ 2;
        boundary,
    )
    energy, psi, iterations, residual = GroundStateED.ground_state_lanczos(
        space;
        maxiter = 300,
        tol = 1e-12,
    )
    S2 = GroundStateED.renyi2_entropy(psi, space, nA)
    return (; energy, S2, iterations, residual)
end

function run_audit()
    L = DEFAULT_L
    t = DEFAULT_T
    U = DEFAULT_U
    mu = DEFAULT_MU
    beta = DEFAULT_BETA
    dtau = DEFAULT_DTAU
    nA = DEFAULT_NA
    boundary = DEFAULT_BOUNDARY

    kinetic, potential, masks = full_fock_matrices(L, t, U, mu; boundary)
    hamiltonian = kinetic + Diagonal(potential)
    thermal_rho, decomposition = exact_thermal_density(hamiltonian, beta)
    ground_vector = decomposition.vectors[:, 1]
    ground_rho = ground_vector * ground_vector'
    trotter_rho, slices = asymmetric_trotter_density(
        kinetic,
        potential,
        beta,
        dtau,
    )

    ground_S2, _, _ = renyi2(ground_rho, masks, L, nA)
    thermal_S2, _, _ = renyi2(thermal_rho, masks, L, nA)
    trotter_S2, _, trotter_rhoA = renyi2(trotter_rho, masks, L, nA)
    fixed = fixed_sector_ground_state(L, t, U, mu, nA; boundary)

    abs(ground_S2 - fixed.S2) <= 5e-11 || error(
        "Full-Fock and fixed-sector ground-state S2 disagree: " *
        "$ground_S2 versus $(fixed.S2)",
    )
    expected_shifted_energy = fixed.energy - U * L / 4
    abs(first(decomposition.values) - expected_shifted_energy) <= 5e-10 || error(
        "Particle-hole energy-shift check failed",
    )

    thermal_shift = thermal_S2 - ground_S2
    trotter_shift = trotter_S2 - thermal_S2
    dqmc_minus_trotter = DEFAULT_DQMC_S2 - trotter_S2
    dqmc_z = dqmc_minus_trotter / DEFAULT_DQMC_ERROR
    trotter_asymmetry = norm(trotter_rhoA - trotter_rhoA') / norm(trotter_rhoA)

    return (
        L,
        t,
        U,
        mu,
        beta,
        dtau,
        slices,
        nA,
        boundary,
        fixed_ground_S2 = fixed.S2,
        full_ground_S2 = ground_S2,
        thermal_exact_S2 = thermal_S2,
        trotter_asymmetric_S2 = trotter_S2,
        thermal_minus_ground = thermal_shift,
        trotter_minus_thermal = trotter_shift,
        dqmc_S2 = DEFAULT_DQMC_S2,
        dqmc_error = DEFAULT_DQMC_ERROR,
        dqmc_minus_trotter,
        dqmc_z_vs_trotter = dqmc_z,
        trotter_reduced_asymmetry = trotter_asymmetry,
    )
end

const OUTPUT_HEADER = join((
    "L", "t", "U", "mu", "beta", "dtau", "Ltau", "nA", "boundary",
    "fixed_ground_S2", "full_ground_S2", "thermal_exact_S2",
    "trotter_asymmetric_S2", "thermal_minus_ground", "trotter_minus_thermal",
    "dqmc_S2", "dqmc_error", "dqmc_minus_trotter", "dqmc_z_vs_trotter",
    "trotter_reduced_asymmetry",
), ',')

function write_result(io, result)
    println(io, OUTPUT_HEADER)
    println(io, join((
        result.L,
        @sprintf("%.8f", result.t),
        @sprintf("%.8f", result.U),
        @sprintf("%.8f", result.mu),
        @sprintf("%.8f", result.beta),
        @sprintf("%.8f", result.dtau),
        result.slices,
        result.nA,
        result.boundary,
        @sprintf("%.15e", result.fixed_ground_S2),
        @sprintf("%.15e", result.full_ground_S2),
        @sprintf("%.15e", result.thermal_exact_S2),
        @sprintf("%.15e", result.trotter_asymmetric_S2),
        @sprintf("%.15e", result.thermal_minus_ground),
        @sprintf("%.15e", result.trotter_minus_thermal),
        @sprintf("%.15e", result.dqmc_S2),
        @sprintf("%.15e", result.dqmc_error),
        @sprintf("%.15e", result.dqmc_minus_trotter),
        @sprintf("%.8f", result.dqmc_z_vs_trotter),
        @sprintf("%.15e", result.trotter_reduced_asymmetry),
    ), ','))
end

function write_result_atomic(path, result)
    output = abspath(path)
    (ispath(output) || islink(output)) && error("Output already exists: $output")
    mkpath(dirname(output))
    temporary, io = mktemp(dirname(output))
    try
        write_result(io, result)
        close(io)
        mv(temporary, output)
    finally
        isopen(io) && close(io)
        isfile(temporary) && rm(temporary; force = true)
    end
    return output
end

function main(args = ARGS)
    cli = parse_cli(args)
    cli === nothing && return
    result = run_audit()
    if cli.output === nothing
        write_result(stdout, result)
    else
        output = write_result_atomic(cli.output, result)
        println("Wrote U=4 L=4 thermal ED audit: $output")
    end
end

end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    U4L4ThermalEDAudit.main()
end
