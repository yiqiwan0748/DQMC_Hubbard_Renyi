#!/usr/bin/env julia

using LinearAlgebra
using Printf
using Random

const DEFAULT_L = 10
const DEFAULT_T = 1.0
const DEFAULT_U = 2.0
const DEFAULT_MU = 0.0
const DEFAULT_NA = "1:5"
const DEFAULT_BOUNDARY = "pbc"
const DEFAULT_LANCZOS_MAXITER = 300
const DEFAULT_LANCZOS_TOL = 1e-11

function usage()
    println("""
    Compute exact-diagonalization S2 references for a 1D Hubbard ring.

    Usage:
      julia analysis/compute_s2_ed.jl --L 10 --t 1.0 --U 2.0 --mu 0.0 --nA 1:5 --boundary pbc --output results/processed/s2_ed.csv

    Options:
      --L INT          chain length, must be even; default $DEFAULT_L
      --t FLOAT        nearest-neighbor hopping amplitude; default $DEFAULT_T
      --U FLOAT        onsite Hubbard interaction; default $DEFAULT_U
      --mu FLOAT       chemical potential; default $DEFAULT_MU
      --nA LIST_OR_RANGE
                       subsystem sizes, for example 1:5 or 1,2,3,4,5; default $DEFAULT_NA
      --boundary pbc|apbc
                       boundary condition; default $DEFAULT_BOUNDARY
      --lanczos-maxiter INT
                       maximum Lanczos iterations; default $DEFAULT_LANCZOS_MAXITER
      --lanczos-tol FLOAT
                       required Ritz residual tolerance; default $DEFAULT_LANCZOS_TOL
      --output PATH    write CSV to PATH; if omitted, write CSV to stdout
      --help           print this message

    The calculation uses the half-filled sector Nup = L/2, Ndn = L/2.
    """)
end

function parse_boundary(text)
    boundary = lowercase(strip(text))
    boundary in ("pbc", "apbc") ||
        error("--boundary must be pbc or apbc, got: $text")
    return boundary
end

boundary_phase(boundary) = parse_boundary(boundary) == "pbc" ? 1.0 : -1.0

function parse_int_list_or_range(text)
    text = strip(text)
    isempty(text) && error("--nA cannot be empty")

    if occursin(":", text)
        pieces = [parse(Int, strip(piece)) for piece in split(text, ":")]
        if length(pieces) == 2
            return collect(pieces[1]:pieces[2])
        elseif length(pieces) == 3
            return collect(pieces[1]:pieces[2]:pieces[3])
        else
            error("Invalid range for --nA: $text")
        end
    end

    return [parse(Int, strip(piece)) for piece in split(text, ",")]
end

function parse_cli(args)
    L = DEFAULT_L
    t = DEFAULT_T
    U = DEFAULT_U
    mu = DEFAULT_MU
    nA_values = parse_int_list_or_range(DEFAULT_NA)
    boundary = DEFAULT_BOUNDARY
    lanczos_maxiter = DEFAULT_LANCZOS_MAXITER
    lanczos_tol = DEFAULT_LANCZOS_TOL
    output = nothing

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--help" || arg == "-h"
            usage()
            exit(0)
        elseif arg == "--L"
            i += 1
            i <= length(args) || error("--L requires a value")
            L = parse(Int, args[i])
        elseif arg == "--t"
            i += 1
            i <= length(args) || error("--t requires a value")
            t = parse(Float64, args[i])
        elseif arg == "--U"
            i += 1
            i <= length(args) || error("--U requires a value")
            U = parse(Float64, args[i])
        elseif arg == "--mu"
            i += 1
            i <= length(args) || error("--mu requires a value")
            mu = parse(Float64, args[i])
        elseif arg == "--nA"
            i += 1
            i <= length(args) || error("--nA requires a value")
            nA_values = parse_int_list_or_range(args[i])
        elseif arg == "--boundary"
            i += 1
            i <= length(args) || error("--boundary requires a value")
            boundary = parse_boundary(args[i])
        elseif arg == "--lanczos-maxiter"
            i += 1
            i <= length(args) || error("--lanczos-maxiter requires a value")
            lanczos_maxiter = parse(Int, args[i])
        elseif arg == "--lanczos-tol"
            i += 1
            i <= length(args) || error("--lanczos-tol requires a value")
            lanczos_tol = parse(Float64, args[i])
        elseif arg == "--output"
            i += 1
            i <= length(args) || error("--output requires a value")
            output = args[i]
        else
            error("Unknown argument: $arg")
        end
        i += 1
    end

    iseven(L) || error("--L must be even so Nup = Ndn = L/2")
    L <= Sys.WORD_SIZE - 2 || error("--L is too large for UInt bit masks on this Julia build")
    all(nA -> 1 <= nA <= L, nA_values) ||
        error("Every nA must satisfy 1 <= nA <= L")
    lanczos_maxiter >= 1 || error("--lanczos-maxiter must be positive")
    isfinite(lanczos_tol) && lanczos_tol > 0 ||
        error("--lanczos-tol must be finite and positive")

    return (
        L = L,
        t = t,
        U = U,
        mu = mu,
        Nup = div(L, 2),
        Ndn = div(L, 2),
        nA_values = nA_values,
        boundary = boundary,
        lanczos_maxiter = lanczos_maxiter,
        lanczos_tol = lanczos_tol,
        output = output
    )
end

function fixed_particle_masks(L, N)
    limit = UInt(1) << L
    return [mask for mask in UInt(0):(limit - UInt(1)) if count_ones(mask) == N]
end

function mask_index(masks)
    return Dict(mask => i for (i, mask) in enumerate(masks))
end

function occupied(mask, site)
    return !iszero(mask & (UInt(1) << site))
end

function fermion_hop_sign(mask, from_site, to_site)
    lo = min(from_site, to_site) + 1
    hi = max(from_site, to_site) - 1
    if lo > hi
        return 1.0
    end

    between = UInt(0)
    for site in lo:hi
        between |= UInt(1) << site
    end

    return isodd(count_ones(mask & between)) ? -1.0 : 1.0
end

function build_spin_hops(masks, index, L, t; boundary = DEFAULT_BOUNDARY)
    hops = [Tuple{Int, Float64}[] for _ in eachindex(masks)]
    seam_phase = boundary_phase(boundary)

    # Each bond connects i to i+1 modulo L. APBC multiplies only the wrap
    # bond by -1, changing its hopping coefficient from -t to +t.
    for (basis_index, mask) in enumerate(masks)
        for site in 0:(L - 1)
            next_site = mod(site + 1, L)
            bond_phase = site == L - 1 ? seam_phase : 1.0
            for (from_site, to_site) in ((site, next_site), (next_site, site))
                if occupied(mask, from_site) && !occupied(mask, to_site)
                    newmask = xor(mask, UInt(1) << from_site)
                    newmask = xor(newmask, UInt(1) << to_site)
                    sign = fermion_hop_sign(mask, from_site, to_site)
                    push!(hops[basis_index], (index[newmask], -t * bond_phase * sign))
                end
            end
        end
    end

    return hops
end

struct HubbardEDSpace
    L::Int
    t::Float64
    U::Float64
    mu::Float64
    Nup::Int
    Ndn::Int
    up_masks::Vector{UInt}
    dn_masks::Vector{UInt}
    up_hops::Vector{Vector{Tuple{Int, Float64}}}
    dn_hops::Vector{Vector{Tuple{Int, Float64}}}
    diagonal::Vector{Float64}
end

function build_space(L, t, U, mu, Nup, Ndn; boundary = DEFAULT_BOUNDARY)
    up_masks = fixed_particle_masks(L, Nup)
    dn_masks = fixed_particle_masks(L, Ndn)
    up_index = mask_index(up_masks)
    dn_index = mask_index(dn_masks)
    up_hops = build_spin_hops(up_masks, up_index, L, t; boundary = boundary)
    dn_hops = build_spin_hops(dn_masks, dn_index, L, t; boundary = boundary)

    nup = length(up_masks)
    ndn = length(dn_masks)
    diagonal = Vector{Float64}(undef, nup * ndn)

    # In this fixed-number half-filled sector, -mu * (Nup + Ndn) is a constant.
    # The particle-hole symmetric Hubbard constant, if used elsewhere, also only
    # shifts energies in this sector and does not change the ground-state vector
    # or the Renyi entropy.
    for iu in eachindex(up_masks)
        for idn in eachindex(dn_masks)
            index = (iu - 1) * ndn + idn
            double_occupancy = count_ones(up_masks[iu] & dn_masks[idn])
            diagonal[index] = U * double_occupancy - mu * (Nup + Ndn)
        end
    end

    return HubbardEDSpace(
        L,
        t,
        U,
        mu,
        Nup,
        Ndn,
        up_masks,
        dn_masks,
        up_hops,
        dn_hops,
        diagonal
    )
end

function dimension(space::HubbardEDSpace)
    return length(space.diagonal)
end

function apply_hamiltonian!(y, x, space::HubbardEDSpace)
    fill!(y, 0.0)

    ndn = length(space.dn_masks)
    for iu in eachindex(space.up_masks)
        for idn in eachindex(space.dn_masks)
            source = (iu - 1) * ndn + idn
            amplitude = x[source]
            y[source] += space.diagonal[source] * amplitude

            for (iu_to, hop) in space.up_hops[iu]
                destination = (iu_to - 1) * ndn + idn
                y[destination] += hop * amplitude
            end

            for (idn_to, hop) in space.dn_hops[idn]
                destination = (iu - 1) * ndn + idn_to
                y[destination] += hop * amplitude
            end
        end
    end

    return y
end

function ground_state_lanczos(
    space::HubbardEDSpace;
    maxiter = DEFAULT_LANCZOS_MAXITER,
    tol = DEFAULT_LANCZOS_TOL,
)
    maxiter >= 1 || error("Lanczos maxiter must be positive")
    isfinite(tol) && tol > 0 ||
        error("Lanczos tolerance must be finite and positive")
    n = dimension(space)
    Random.seed!(1234)

    q = randn(n)
    q ./= norm(q)
    qprev = zeros(Float64, n)
    z = similar(q)

    basis = Vector{Vector{Float64}}()
    alpha = Float64[]
    beta = Float64[]
    beta_prev = 0.0

    for iter in 1:maxiter
        push!(basis, copy(q))
        apply_hamiltonian!(z, q, space)

        if iter > 1
            z .-= beta_prev .* qprev
        end

        a = dot(q, z)
        z .-= a .* q

        # Full reorthogonalization keeps the small ED reference stable and
        # deterministic for the benchmark-sized Hilbert spaces used here.
        for previous in basis
            z .-= dot(previous, z) .* previous
        end

        b = norm(z)
        push!(alpha, a)

        tridiagonal = SymTridiagonal(alpha, beta)
        eigen_result = eigen(tridiagonal)
        ground_position = argmin(eigen_result.values)
        energy = eigen_result.values[ground_position]
        ritz = eigen_result.vectors[:, ground_position]
        residual = b * abs(ritz[end])

        if residual < tol
            psi = zeros(Float64, n)
            for i in eachindex(basis)
                psi .+= ritz[i] .* basis[i]
            end
            psi ./= norm(psi)
            return energy, psi, iter, residual
        end

        if iter == maxiter
            error(
                @sprintf(
                    "Lanczos did not converge in %d iterations: residual %.16e >= tolerance %.16e",
                    maxiter,
                    residual,
                    tol,
                ),
            )
        end

        push!(beta, b)
        qprev = q
        q = z ./ b
        beta_prev = b
    end

    error("Lanczos failed to return a ground state")
end

function compressed_bits(mask, sites)
    compressed = UInt(0)
    for (position, site) in enumerate(sites)
        if occupied(mask, site)
            compressed |= UInt(1) << (position - 1)
        end
    end
    return compressed
end

function subsystem_sign(upA, upB, dnA)
    # The fixed-sector basis is ordered as up(A), up(B), down(A), down(B).
    # For the spatial bipartition we group subsystem modes as
    # A = up(A), down(A) and B = up(B), down(B). Moving down(A) past up(B)
    # gives (-1)^(N_down_A * N_up_B).
    crossings = count_ones(dnA) * count_ones(upB)
    return isodd(crossings) ? -1.0 : 1.0
end

function renyi2_entropy_sites(psi, space::HubbardEDSpace, sitesA)
    sitesB = [site for site in 0:(space.L - 1) if !(site in sitesA)]
    nA = length(sitesA)
    nB = length(sitesB)
    dimA = 1 << (2 * nA)
    ndn = length(space.dn_masks)
    groups = Dict{Int, Vector{Tuple{Int, Float64}}}()

    # Compress arbitrary site sets into dense subsystem/environment bit labels.
    # This lets large contiguous regions use the smaller pure-state complement
    # without allocating rhoA for the larger side.
    for iu in eachindex(space.up_masks)
        upA = compressed_bits(space.up_masks[iu], sitesA)
        upB = compressed_bits(space.up_masks[iu], sitesB)
        for idn in eachindex(space.dn_masks)
            coefficient = psi[(iu - 1) * ndn + idn]
            iszero(coefficient) && continue

            dnA = compressed_bits(space.dn_masks[idn], sitesA)
            dnB = compressed_bits(space.dn_masks[idn], sitesB)
            a_index = Int(upA | (dnA << nA)) + 1
            b_index = Int(upB | (dnB << nB)) + 1
            signed_coefficient =
                subsystem_sign(upA, upB, dnA) * coefficient

            entries = get!(groups, b_index, Tuple{Int, Float64}[])
            push!(entries, (a_index, signed_coefficient))
        end
    end

    rhoA = zeros(Float64, dimA, dimA)
    for entries in values(groups)
        for (a1, c1) in entries
            for (a2, c2) in entries
                rhoA[a1, a2] += c1 * c2
            end
        end
    end

    tr_rhoA_squared = sum(abs2, rhoA)
    return -log(tr_rhoA_squared)
end

function renyi2_entropy(psi, space::HubbardEDSpace, nA)
    if nA <= div(space.L, 2)
        # DQMC convention: A is sites 1:nA, represented by zero-based bits
        # 0:(nA - 1).
        return renyi2_entropy_sites(psi, space, collect(0:(nA - 1)))
    end

    # The ground state is pure, so S2(A) = S2(B). For A = sites 1:nA, the
    # complement is sites nA+1:L, represented by zero-based bits nA:(L - 1).
    # Computing rho for this smaller side avoids the dense 2^(2*nA) allocation
    # that would otherwise fail for large nA.
    return renyi2_entropy_sites(psi, space, collect(nA:(space.L - 1)))
end

function write_results(
    io,
    config,
    space,
    psi;
    energy,
    iterations,
    residual,
    lanczos_tol,
)
    println(
        io,
        "L,t,U,mu,Nup,Ndn,nA,boundary,S2_ED,energy,iterations,residual,lanczos_tol",
    )
    for nA in config.nA_values
        S2 = renyi2_entropy(psi, space, nA)
        println(
            io,
            join(
                [
                    config.L,
                    @sprintf("%.12g", config.t),
                    @sprintf("%.12g", config.U),
                    @sprintf("%.12g", config.mu),
                    config.Nup,
                    config.Ndn,
                    nA,
                    config.boundary,
                    @sprintf("%.15f", S2),
                    @sprintf("%.16e", energy),
                    iterations,
                    @sprintf("%.16e", residual),
                    @sprintf("%.16e", lanczos_tol),
                ],
                ","
            )
        )
    end
end

function write_results_atomic(
    path,
    config,
    space,
    psi;
    energy,
    iterations,
    residual,
    lanczos_tol,
)
    output_path = abspath(path)
    (ispath(output_path) || islink(output_path)) &&
        error("Output path already exists: $output_path")
    mkpath(dirname(output_path))
    temp_path, temp_io = mktemp(dirname(output_path))
    try
        write_results(
            temp_io,
            config,
            space,
            psi;
            energy = energy,
            iterations = iterations,
            residual = residual,
            lanczos_tol = lanczos_tol,
        )
        close(temp_io)
        mv(temp_path, output_path)
    finally
        isopen(temp_io) && close(temp_io)
        isfile(temp_path) && rm(temp_path; force = true)
    end
end

function main(args = ARGS)
    config = parse_cli(args)
    space = build_space(
        config.L,
        config.t,
        config.U,
        config.mu,
        config.Nup,
        config.Ndn;
        boundary = config.boundary
    )
    energy, psi, iterations, residual = ground_state_lanczos(
        space;
        maxiter = config.lanczos_maxiter,
        tol = config.lanczos_tol,
    )

    @info "Computed ground state" energy iterations residual dimension = dimension(space)

    if config.output === nothing
        write_results(
            stdout,
            config,
            space,
            psi;
            energy = energy,
            iterations = iterations,
            residual = residual,
            lanczos_tol = config.lanczos_tol,
        )
    else
        write_results_atomic(
            config.output,
            config,
            space,
            psi;
            energy = energy,
            iterations = iterations,
            residual = residual,
            lanczos_tol = config.lanczos_tol,
        )
        println("Wrote ED S2 reference CSV: ", abspath(config.output))
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
