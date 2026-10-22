using Test
using Random
using LinearAlgebra
using Printf

const PROJECT_ROOT = dirname(@__DIR__)

include(joinpath(PROJECT_ROOT, "src", "HubbardRenyiDQMC.jl"))
using .HubbardRenyiDQMC

function explicit_grover_from_blocks(G1C, G2C)
    nC = size(G1C, 1)
    identity_C = Matrix{eltype(G1C)}(I, nC, nC)
    return G1C * G2C + (identity_C - G1C) * (identity_C - G2C)
end

function test_tracker_matches_explicit(tracker, G1C, G2C; atol = 2e-9, rtol = 2e-8)
    expected_grover = explicit_grover_from_blocks(G1C, G2C)
    expected_logabsdet, expected_sign = logabsdet(expected_grover)

    @test isapprox(tracker.G1C, G1C; atol, rtol)
    @test isapprox(tracker.G2C, G2C; atol, rtol)
    @test isapprox(tracker.grover, expected_grover; atol, rtol)
    @test isapprox(tracker.logabsdet_value, expected_logabsdet; atol, rtol)

    if tracker.usable_inverse && !iszero(expected_sign)
        @test isapprox(
            tracker.inverse_grover * expected_grover,
            Matrix{Float64}(I, size(G1C, 1), size(G1C, 1));
            atol,
            rtol,
        )
    end

    return nothing
end

function validate_determinant_lemma()
    rng = Xoshiro(0x5eed)

    @testset "Grover rank-one determinant identity" begin
        for nC in 1:8
            for _ in 1:50
                G1C = randn(rng, nC, nC)
                G2C = randn(rng, nC, nC)
                p = randn(rng, nC)
                q = randn(rng, nC)

                grover_old = explicit_grover_from_blocks(G1C, G2C)
                abs(det(grover_old)) > 1e-8 || continue

                grover_new = explicit_grover_from_blocks(G1C + p * transpose(q), G2C)
                r = transpose(2 * G2C - I) * q

                direct_ratio = det(grover_new) / det(grover_old)
                lemma_ratio = 1 + transpose(r) * (grover_old \ p)

                @test isapprox(
                    direct_ratio,
                    lemma_ratio;
                    rtol = 5e-10,
                    atol = 5e-11,
                )
            end
        end
    end

    return nothing
end

function validate_tracker()
    rng = Xoshiro(0x71ac)
    fast_path_cases = 0

    @testset "Grover rank-one tracker single proposals" begin
        for N in 4:8
            C = collect(1:2:N)
            nC = length(C)

            for _ in 1:40
                G1 = randn(rng, N, N)
                G2 = randn(rng, N, N)
                G0tau = randn(rng, N, N)
                Gtau0 = randn(rng, N, N)
                site = rand(rng, 1:N)
                ratio = 0.5 + rand(rng)
                delta = 0.2 * randn(rng)

                tracker = GroverRank1Tracker(Float64, nC)
                set_passive_grover_block!(tracker, G2, C)
                rebuild_active_grover!(tracker, G1, C)
                tracker.usable_inverse || continue

                old_logabsdet = tracker.logabsdet_value
                proposed_logabsdet = prepare_grover_rank1_proposal!(
                    tracker,
                    G0tau,
                    Gtau0,
                    ratio,
                    delta,
                    site,
                    C,
                )
                tracker.pending_explicit_rebuild && continue
                fast_path_cases += 1

                p = (delta / ratio) .* G0tau[C, site]
                q = vec(Gtau0[site, C])
                G1C_new = G1[C, C] + p * transpose(q)
                grover_new = explicit_grover_from_blocks(G1C_new, G2[C, C])
                direct_logabsdet, direct_sign = logabsdet(grover_new)

                @test isapprox(
                    proposed_logabsdet,
                    direct_logabsdet;
                    rtol = 5e-10,
                    atol = 5e-11,
                )
                @test isapprox(
                    proposed_logabsdet - old_logabsdet,
                    log(abs(det(grover_new) / det(explicit_grover_from_blocks(G1[C, C], G2[C, C]))));
                    rtol = 5e-10,
                    atol = 5e-11,
                )

                accept_grover_rank1_proposal!(tracker)
                test_tracker_matches_explicit(tracker, G1C_new, G2[C, C])
            end
        end

        @test fast_path_cases >= 100
    end

    return nothing
end

function validate_tracker_sequences()
    rng = Xoshiro(0x5e91)
    N = 8
    C = [1, 3, 6, 8]
    nC = length(C)
    G1 = 0.4 .* Matrix{Float64}(I, N, N) + 0.02 .* randn(rng, N, N)
    G2 = 0.6 .* Matrix{Float64}(I, N, N) + 0.02 .* randn(rng, N, N)
    G1C = copy(G1[C, C])
    G2C = copy(G2[C, C])
    tracker = GroverRank1Tracker(Float64, nC)
    set_passive_grover_block!(tracker, G2, C)
    rebuild_active_grover!(tracker, G1, C)

    @testset "Sequential accepted rank-one updates" begin
        for update_index in 1:4nC
            G0tau = randn(rng, N, N)
            Gtau0 = randn(rng, N, N)
            site = rand(rng, 1:N)
            ratio = 0.5 + rand(rng)
            delta = 0.05 * randn(rng)

            proposed_logabsdet = prepare_grover_rank1_proposal!(
                tracker, G0tau, Gtau0, ratio, delta, site, C,
            )
            @test !tracker.pending_explicit_rebuild
            p = (delta / ratio) .* G0tau[C, site]
            q = vec(Gtau0[site, C])
            G1C += p * transpose(q)
            expected_logabsdet, _ =
                logabsdet(explicit_grover_from_blocks(G1C, G2C))
            @test isapprox(
                proposed_logabsdet,
                expected_logabsdet;
                rtol = 2e-8,
                atol = 2e-9,
            )

            accept_grover_rank1_proposal!(tracker)
            test_tracker_matches_explicit(tracker, G1C, G2C)
            if update_index % nC == 0
                @test tracker.accepted_since_rebuild == 0
            end
        end
    end

    @testset "Rejected proposal leaves tracker unchanged" begin
        G1C_before = copy(tracker.G1C)
        grover_before = copy(tracker.grover)
        inverse_before = copy(tracker.inverse_grover)
        logabsdet_before = tracker.logabsdet_value

        prepare_grover_rank1_proposal!(
            tracker,
            randn(rng, N, N),
            randn(rng, N, N),
            1.1,
            0.07,
            2,
            C;
            denominator_tolerance = Inf,
        )
        @test tracker.pending_explicit_rebuild
        @test tracker.G1C == G1C_before
        @test tracker.grover == grover_before
        @test tracker.inverse_grover == inverse_before
        @test tracker.logabsdet_value == logabsdet_before

        G0tau = randn(rng, N, N)
        Gtau0 = randn(rng, N, N)
        ratio = 0.9
        delta = -0.04
        site = 7
        prepare_grover_rank1_proposal!(
            tracker, G0tau, Gtau0, ratio, delta, site, C,
        )
        G1C += (delta / ratio) .* G0tau[C, site] * transpose(vec(Gtau0[site, C]))
        accept_grover_rank1_proposal!(tracker)
        test_tracker_matches_explicit(tracker, G1C, G2C)
    end

    @testset "Forced explicit fallback accepts exact state" begin
        G0tau = randn(rng, N, N)
        Gtau0 = randn(rng, N, N)
        ratio = 1.2
        delta = 0.03
        site = 4
        prepare_grover_rank1_proposal!(
            tracker,
            G0tau,
            Gtau0,
            ratio,
            delta,
            site,
            C;
            denominator_tolerance = Inf,
        )
        @test tracker.pending_explicit_rebuild
        G1C += (delta / ratio) .* G0tau[C, site] * transpose(vec(Gtau0[site, C]))
        accept_grover_rank1_proposal!(tracker)
        test_tracker_matches_explicit(tracker, G1C, G2C)
    end

    @testset "Ill-conditioned Grover matrix disables inverse shortcut" begin
        ill_tracker = GroverRank1Tracker(Float64, 3)
        identity3 = Matrix{Float64}(I, 3, 3)
        ill_G2 = Diagonal([1.0, 1e-10, 1e-12]) |> Matrix
        set_passive_grover_block!(ill_tracker, ill_G2, 1:3)
        rebuild_active_grover!(ill_tracker, identity3, 1:3)
        @test !ill_tracker.usable_inverse

        G0tau = Matrix{Float64}(I, 3, 3)
        Gtau0 = Matrix{Float64}(I, 3, 3)
        proposed_logabsdet = prepare_grover_rank1_proposal!(
            ill_tracker, G0tau, Gtau0, 1.0, 0.01, 1, 1:3,
        )
        @test ill_tracker.pending_explicit_rebuild
        expected_G1C = identity3 + 0.01 .* G0tau[:, 1] * transpose(Gtau0[1, :])
        expected_logabsdet, _ =
            logabsdet(explicit_grover_from_blocks(expected_G1C, ill_G2))
        @test isapprox(proposed_logabsdet, expected_logabsdet; rtol = 1e-10, atol = 1e-11)
        accept_grover_rank1_proposal!(ill_tracker)
        @test !ill_tracker.usable_inverse
        test_tracker_matches_explicit(ill_tracker, expected_G1C, ill_G2)

        second_logabsdet = prepare_grover_rank1_proposal!(
            ill_tracker, G0tau, Gtau0, 1.0, -0.02, 2, 1:3,
        )
        @test ill_tracker.pending_explicit_rebuild
        second_G1C = expected_G1C - 0.02 .* G0tau[:, 2] * transpose(Gtau0[2, :])
        second_expected_logabsdet, _ =
            logabsdet(explicit_grover_from_blocks(second_G1C, ill_G2))
        @test isapprox(
            second_logabsdet,
            second_expected_logabsdet;
            rtol = 1e-10,
            atol = 1e-11,
        )
    end

    return nothing
end

function compare_replica_states(direct, rank1; atol = 2e-9, rtol = 2e-8)
    @test direct.s == rank1.s

    for field in (
        :Gup, :Gdn,
        :Gup_00, :Gdn_00,
        :Gup_τ0, :Gdn_τ0,
        :Gup_0τ, :Gdn_0τ,
    )
        @test isapprox(
            getproperty(direct, field),
            getproperty(rank1, field);
            atol,
            rtol,
        )
    end

    @test isapprox(direct.logdetGup, rank1.logdetGup; atol, rtol)
    @test isapprox(direct.logdetGdn, rank1.logdetGdn; atol, rtol)
    @test direct.sgndetGup == rank1.sgndetGup
    @test direct.sgndetGdn == rank1.sgndetGdn
    return nothing
end

function validate_same_seed_update_rounds(; rounds = 12)
    t = 1.0
    U = 2.0
    mu = 0.0
    beta = 4.0
    L = 4
    dtau = 0.1
    n_stab = 10
    C = collect(1:2)
    lambda = 0.5

    base1 = initializer(t, U, mu, beta, L, dtau, 61001; n_stab)
    base2 = initializer(t, U, mu, beta, L, dtau, 61002; n_stab)

    direct1, direct2 = deepcopy(base1), deepcopy(base2)
    rank1_1, rank1_2 = deepcopy(base1), deepcopy(base2)
    direct_deltaG = 0.0
    rank1_deltaG = 0.0

    @testset "Same-seed direct versus rank-one update" begin
        for _ in 1:rounds
            direct1, direct_deltaG, direct_acceptance1 =
                update_replica_fused!(direct1, direct2, direct_deltaG, C, lambda)
            rank1_1, rank1_deltaG, rank1_acceptance1 =
                update_replica_rank1_fused!(rank1_1, rank1_2, rank1_deltaG, C, lambda)

            @test direct_acceptance1 == rank1_acceptance1
            compare_replica_states(direct1, rank1_1)

            direct2, direct_deltaG, direct_acceptance2 =
                update_replica_fused!(direct2, direct1, direct_deltaG, C, lambda)
            rank1_2, rank1_deltaG, rank1_acceptance2 =
                update_replica_rank1_fused!(rank1_2, rank1_1, rank1_deltaG, C, lambda)

            @test direct_acceptance2 == rank1_acceptance2
            compare_replica_states(direct2, rank1_2)

            direct_grover_up = grover_matrix(direct1.Gup_00, direct2.Gup_00, C)
            rank1_grover_up = grover_matrix(rank1_1.Gup_00, rank1_2.Gup_00, C)
            direct_grover_dn = grover_matrix(direct1.Gdn_00, direct2.Gdn_00, C)
            rank1_grover_dn = grover_matrix(rank1_1.Gdn_00, rank1_2.Gdn_00, C)
            @test isapprox(direct_grover_up, rank1_grover_up; rtol = 2e-8, atol = 2e-9)
            @test isapprox(direct_grover_dn, rank1_grover_dn; rtol = 2e-8, atol = 2e-9)
            @test isapprox(
                det(direct_grover_up) * det(direct_grover_dn),
                det(rank1_grover_up) * det(rank1_grover_dn);
                rtol = 2e-8,
                atol = 2e-9,
            )
        end

        @test isapprox(direct_deltaG, rank1_deltaG; rtol = 2e-8, atol = 2e-10)
        @test rand(copy(direct1.rng), UInt64) == rand(copy(rank1_1.rng), UInt64)
        @test rand(copy(direct2.rng), UInt64) == rand(copy(rank1_2.rng), UInt64)
    end

    return nothing
end

function time_update_rounds(update_function; rounds, seed_offset)
    t = 1.0
    U = 2.0
    mu = 0.0
    beta = 32.0
    L = 8
    dtau = 0.05
    n_stab = 10
    C = collect(1:4)
    lambda = 0.5

    replica1 = initializer(t, U, mu, beta, L, dtau, 62001 + seed_offset; n_stab)
    replica2 = initializer(t, U, mu, beta, L, dtau, 62002 + seed_offset; n_stab)
    deltaG = 0.0

    replica1, deltaG, _ = update_function(replica1, replica2, deltaG, C, lambda)
    replica2, deltaG, _ = update_function(replica2, replica1, deltaG, C, lambda)

    GC.gc()
    timing = @timed begin
        for _ in 1:rounds
            replica1, deltaG, _ = update_function(replica1, replica2, deltaG, C, lambda)
            replica2, deltaG, _ = update_function(replica2, replica1, deltaG, C, lambda)
        end
    end

    return timing.time, deltaG, timing.bytes, timing.gctime
end

function benchmark_update_rounds(; rounds = 20, repeats = 5)
    BLAS.set_num_threads(1)
    direct_times = Float64[]
    rank1_times = Float64[]
    direct_bytes = Int[]
    rank1_bytes = Int[]
    direct_gc_times = Float64[]
    rank1_gc_times = Float64[]

    for repeat in 1:repeats
        if isodd(repeat)
            direct_seconds, _, direct_allocated, direct_gc = time_update_rounds(
                update_replica_fused!;
                rounds,
                seed_offset = 100 * repeat,
            )
            rank1_seconds, _, rank1_allocated, rank1_gc = time_update_rounds(
                update_replica_rank1_fused!;
                rounds,
                seed_offset = 100 * repeat,
            )
        else
            rank1_seconds, _, rank1_allocated, rank1_gc = time_update_rounds(
                update_replica_rank1_fused!;
                rounds,
                seed_offset = 100 * repeat,
            )
            direct_seconds, _, direct_allocated, direct_gc = time_update_rounds(
                update_replica_fused!;
                rounds,
                seed_offset = 100 * repeat,
            )
        end
        push!(direct_times, direct_seconds)
        push!(rank1_times, rank1_seconds)
        push!(direct_bytes, direct_allocated)
        push!(rank1_bytes, rank1_allocated)
        push!(direct_gc_times, direct_gc)
        push!(rank1_gc_times, rank1_gc)
    end

    direct_mean = sum(direct_times) / length(direct_times)
    rank1_mean = sum(rank1_times) / length(rank1_times)

    println("Direct seconds: ", join((@sprintf("%.6f", x) for x in direct_times), ", "))
    println("Rank-one seconds: ", join((@sprintf("%.6f", x) for x in rank1_times), ", "))
    @printf(
        "Mean allocation: direct %.3f MiB, rank-one %.3f MiB\n",
        sum(direct_bytes) / length(direct_bytes) / 2.0^20,
        sum(rank1_bytes) / length(rank1_bytes) / 2.0^20,
    )
    @printf(
        "Mean GC time: direct %.6f s, rank-one %.6f s\n",
        sum(direct_gc_times) / length(direct_gc_times),
        sum(rank1_gc_times) / length(rank1_gc_times),
    )
    @printf("Mean speedup: %.6fx\n", direct_mean / rank1_mean)
    return nothing
end

function main()
    validate_determinant_lemma()
    validate_tracker()
    validate_tracker_sequences()
    validate_same_seed_update_rounds()

    if "--benchmark" in ARGS
        benchmark_update_rounds()
    end

    println("Grover rank-one validation completed.")
end

main()
