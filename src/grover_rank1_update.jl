"""
Exact rank-one Grover-determinant updates for the active replica.

For

    g(G1, G2) = G1 * G2 + (I - G1) * (I - G2),

a local update `G1' = G1 + p*qᵀ` gives

    g' = g + p * (qᵀ * (2G2 - I)).

The matrix determinant lemma therefore evaluates `det(g') / det(g)` in
`O(nA^2)`. The tracker is rebuilt initially and at unequal-time stabilization
boundaries, and falls back to an explicit rebuild when a Sherman--Morrison
denominator is small.

The proposal ratio needs only `1 + rᵀ*g⁻¹*p`. The separate `g⁻ᵀ*r` product is
used only by the Sherman--Morrison inverse update, so it is evaluated only
after the proposal is accepted.

Between stabilization boundaries the tracker, rather than the full `G00`, is
the authoritative subsystem state. The full `G00` is reconstructed at every
sweep endpoint, so callers must measure only after a complete replica update.
"""

mutable struct GroverRank1Tracker{T<:Real}
    G1C::Matrix{T}
    G2C::Matrix{T}
    identity_minus_G2::Matrix{T}
    twice_G2_minus_I::Matrix{T}
    grover::Matrix{T}
    inverse_grover::Matrix{T}
    identity_minus_G1_work::Matrix{T}
    product_work::Matrix{T}
    logabsdet_value::Float64
    usable_inverse::Bool
    condition_estimate::Float64
    accepted_since_rebuild::Int

    p::Vector{T}
    q::Vector{T}
    r::Vector{T}
    inverse_times_p::Vector{T}
    transpose_inverse_times_r::Vector{T}

    trial_G1C::Matrix{T}
    trial_grover::Matrix{T}
    trial_inverse_grover::Matrix{T}
    pending_ratio::T
    pending_logabsdet::Float64
    pending_usable_inverse::Bool
    pending_condition_estimate::Float64
    pending_explicit_rebuild::Bool
end

function GroverRank1Tracker(::Type{T}, nC::Integer) where {T<:Real}
    nC >= 0 || throw(ArgumentError("subsystem size must be nonnegative"))
    matrix() = zeros(T, nC, nC)
    vector() = zeros(T, nC)

    return GroverRank1Tracker(
        matrix(), matrix(), matrix(), matrix(), matrix(), matrix(), matrix(), matrix(),
        0.0, false, Inf, 0,
        vector(), vector(), vector(), vector(), vector(),
        matrix(), matrix(), matrix(),
        one(T), 0.0, false, Inf, false,
    )
end

function _copy_subsystem_block!(dest, source, C)
    nC = length(C)
    size(dest) == (nC, nC) || throw(DimensionMismatch("subsystem block size changed"))

    @inbounds for j in 1:nC
        source_j = C[j]
        for i in 1:nC
            dest[i, j] = source[C[i], source_j]
        end
    end

    return dest
end

function _form_grover_from_blocks!(
    dest,
    G1C,
    G2C,
    identity_minus_G2,
    identity_minus_G1_work,
    product_work,
)
    mul!(dest, G1C, G2C)
    nC = size(dest, 1)

    @inbounds for j in 1:nC
        for i in 1:nC
            identity_minus_G1_work[i, j] =
                (i == j ? one(eltype(dest)) : zero(eltype(dest))) - G1C[i, j]
        end
    end

    mul!(product_work, identity_minus_G1_work, identity_minus_G2)
    @inbounds @simd for index in eachindex(dest, product_work)
        dest[index] += product_work[index]
    end

    return dest
end

function _matrix_inf_norm(matrix)
    largest = 0.0
    @inbounds for i in axes(matrix, 1)
        row_sum = 0.0
        @simd for j in axes(matrix, 2)
            row_sum += abs(matrix[i, j])
        end
        largest = max(largest, row_sum)
    end
    return largest
end

function _inverse_residual_inf_norm(matrix, inverse_matrix)
    largest = 0.0
    n = size(matrix, 1)

    @inbounds for i in 1:n
        row_sum = 0.0
        for j in 1:n
            value = zero(eltype(matrix))
            @simd for k in 1:n
                value += matrix[i, k] * inverse_matrix[k, j]
            end
            value -= i == j ? one(value) : zero(value)
            row_sum += abs(value)
        end
        largest = max(largest, row_sum)
    end

    return largest
end

function _assess_inverse(matrix, inverse_matrix)
    condition_estimate =
        _matrix_inf_norm(matrix) * _matrix_inf_norm(inverse_matrix)
    residual = max(
        _inverse_residual_inf_norm(matrix, inverse_matrix),
        _inverse_residual_inf_norm(inverse_matrix, matrix),
    )
    tolerance = 256 * eps(Float64) * max(1.0, condition_estimate)
    condition_limit = inv(sqrt(eps(Float64)))
    usable =
        isfinite(condition_estimate) &&
        isfinite(residual) &&
        condition_estimate <= condition_limit &&
        residual <= max(1e-12, tolerance)
    return condition_estimate, usable
end

function _inverse_and_logabsdet!(inverse_dest, matrix)
    copyto!(inverse_dest, matrix)
    factorization = lu!(inverse_dest; check = false)
    logabs, sign_value = logabsdet(factorization)
    usable = isfinite(logabs) && !iszero(sign_value)

    if usable
        try
            LinearAlgebra.LAPACK.getri!(
                factorization.factors,
                factorization.ipiv,
            )
        catch error
            if error isa LinearAlgebra.SingularException ||
               error isa LinearAlgebra.LAPACKException
                usable = false
            else
                rethrow()
            end
        end
    end

    condition_estimate = Inf
    if usable
        condition_estimate, usable = _assess_inverse(matrix, inverse_dest)
    end

    return Float64(logabs), usable, condition_estimate
end

function set_passive_grover_block!(tracker::GroverRank1Tracker, G2, C)
    _copy_subsystem_block!(tracker.G2C, G2, C)
    nC = length(C)

    @inbounds for j in 1:nC
        for i in 1:nC
            identity_value =
                i == j ? one(eltype(tracker.G2C)) : zero(eltype(tracker.G2C))
            tracker.identity_minus_G2[i, j] =
                identity_value - tracker.G2C[i, j]
            tracker.twice_G2_minus_I[i, j] =
                2 * tracker.G2C[i, j] - identity_value
        end
    end

    return tracker
end

function rebuild_active_grover!(tracker::GroverRank1Tracker, G1, C)
    _copy_subsystem_block!(tracker.G1C, G1, C)
    _form_grover_from_blocks!(
        tracker.grover,
        tracker.G1C,
        tracker.G2C,
        tracker.identity_minus_G2,
        tracker.identity_minus_G1_work,
        tracker.product_work,
    )
    tracker.logabsdet_value, tracker.usable_inverse, tracker.condition_estimate =
        _inverse_and_logabsdet!(tracker.inverse_grover, tracker.grover)
    tracker.accepted_since_rebuild = 0
    return tracker
end

@inline function _bilinear_no_conjugation(left, right)
    value = zero(promote_type(eltype(left), eltype(right)))
    @inbounds @simd for i in eachindex(left, right)
        value += left[i] * right[i]
    end
    return value
end

function _rank1_add_no_conjugation!(matrix, left, right, scale = one(eltype(matrix)))
    @inbounds for j in eachindex(right)
        scaled_right = scale * right[j]
        @simd for i in eachindex(left)
            matrix[i, j] += left[i] * scaled_right
        end
    end
    return matrix
end

function _prepare_explicit_trial!(tracker::GroverRank1Tracker)
    copyto!(tracker.trial_G1C, tracker.G1C)
    _rank1_add_no_conjugation!(tracker.trial_G1C, tracker.p, tracker.q)
    _form_grover_from_blocks!(
        tracker.trial_grover,
        tracker.trial_G1C,
        tracker.G2C,
        tracker.identity_minus_G2,
        tracker.identity_minus_G1_work,
        tracker.product_work,
    )

    tracker.pending_logabsdet,
    tracker.pending_usable_inverse,
    tracker.pending_condition_estimate =
        _inverse_and_logabsdet!(tracker.trial_inverse_grover, tracker.trial_grover)
    tracker.pending_explicit_rebuild = true
    return tracker.pending_logabsdet
end

function prepare_grover_rank1_proposal!(
    tracker::GroverRank1Tracker,
    G0tau,
    Gtau0,
    ratio,
    delta,
    site::Int,
    C;
    denominator_tolerance::Float64 = sqrt(eps(Float64)),
)
    coefficient = delta / ratio

    @inbounds for a in eachindex(C)
        subsystem_site = C[a]
        tracker.p[a] = coefficient * G0tau[subsystem_site, site]
        tracker.q[a] = Gtau0[site, subsystem_site]
    end

    if !tracker.usable_inverse
        return _prepare_explicit_trial!(tracker)
    end

    mul!(tracker.r, transpose(tracker.twice_G2_minus_I), tracker.q)
    mul!(tracker.inverse_times_p, tracker.inverse_grover, tracker.p)

    determinant_ratio =
        one(eltype(tracker.G1C)) +
        _bilinear_no_conjugation(tracker.r, tracker.inverse_times_p)

    summation_scale = 1.0
    @inbounds @simd for a in eachindex(tracker.r, tracker.inverse_times_p)
        summation_scale += abs(tracker.r[a] * tracker.inverse_times_p[a])
    end
    scaled_tolerance =
        denominator_tolerance *
        sqrt(max(1.0, tracker.condition_estimate)) *
        summation_scale

    if !isfinite(determinant_ratio) || abs(determinant_ratio) <= scaled_tolerance
        return _prepare_explicit_trial!(tracker)
    end

    tracker.pending_ratio = determinant_ratio
    tracker.pending_logabsdet =
        tracker.logabsdet_value + log(abs(determinant_ratio))
    tracker.pending_usable_inverse = true
    tracker.pending_condition_estimate = tracker.condition_estimate
    tracker.pending_explicit_rebuild = false
    return tracker.pending_logabsdet
end

function accept_grover_rank1_proposal!(tracker::GroverRank1Tracker)
    if tracker.pending_explicit_rebuild
        copyto!(tracker.G1C, tracker.trial_G1C)
        copyto!(tracker.grover, tracker.trial_grover)
        if tracker.pending_usable_inverse
            copyto!(tracker.inverse_grover, tracker.trial_inverse_grover)
        end
        tracker.condition_estimate = tracker.pending_condition_estimate
        tracker.accepted_since_rebuild = 0
    else
        mul!(
            tracker.transpose_inverse_times_r,
            transpose(tracker.inverse_grover),
            tracker.r,
        )
        _rank1_add_no_conjugation!(tracker.G1C, tracker.p, tracker.q)
        _rank1_add_no_conjugation!(tracker.grover, tracker.p, tracker.r)
        _rank1_add_no_conjugation!(
            tracker.inverse_grover,
            tracker.inverse_times_p,
            tracker.transpose_inverse_times_r,
            -inv(tracker.pending_ratio),
        )
        tracker.accepted_since_rebuild += 1
    end

    tracker.logabsdet_value = tracker.pending_logabsdet
    tracker.usable_inverse = tracker.pending_usable_inverse

    condition_limit = inv(sqrt(eps(Float64)))
    if tracker.usable_inverse && !tracker.pending_explicit_rebuild
        tracker.condition_estimate =
            _matrix_inf_norm(tracker.grover) *
            _matrix_inf_norm(tracker.inverse_grover)

        if !isfinite(tracker.condition_estimate) ||
           tracker.condition_estimate > condition_limit
            tracker.logabsdet_value,
            tracker.usable_inverse,
            tracker.condition_estimate =
                _inverse_and_logabsdet!(tracker.inverse_grover, tracker.grover)
            tracker.accepted_since_rebuild = 0
        end
    end

    # The condition-norm check above already rebuilds on nonfinite entries.
    # The more expensive two-sided residual check is needed only periodically.
    check_interval = max(1, size(tracker.G1C, 1))
    if tracker.usable_inverse && tracker.accepted_since_rebuild >= check_interval
        tracker.condition_estimate, tracker.usable_inverse =
            _assess_inverse(tracker.grover, tracker.inverse_grover)

        if !tracker.usable_inverse
            tracker.logabsdet_value,
            tracker.usable_inverse,
            tracker.condition_estimate =
                _inverse_and_logabsdet!(tracker.inverse_grover, tracker.grover)
        end
        tracker.accepted_since_rebuild = 0
    end

    return tracker
end

function _stabilization_rebuilds_G00(fermion_greens_calculator)
    (; l, Lτ, forward, n_stab, N_stab) = fermion_greens_calculator
    _, interval_slice =
        jdqmcf.stabilization_interval(fermion_greens_calculator)

    if forward
        return l == Lτ || (interval_slice == n_stab && N_stab > 1)
    else
        return l == 1 || (interval_slice == 1 && N_stab > 1)
    end
end

function _local_update_replica_rank1_fused!(
    Gup, Gup_tau0, Gup_0tau, logdetGup, sgndetGup, Bup,
    Gdn, Gdn_tau0, Gdn_0tau, logdetGdn, sgndetGdn, Bdn,
    s, mu, alpha, dtau, rng, perm, u, v,
    C, lambda,
    tracker_up, tracker_dn,
)
    shuffle!(rng, perm)

    accepted = 0
    N = length(perm)

    for i in perm
        delta_V_up = -alpha / dtau * (-2 * s[i])
        delta_V_dn = +alpha / dtau * (-2 * s[i])

        ratio_up, delta_up = jdqmcf.local_update_det_ratio(Gup, delta_V_up, i, dtau)
        ratio_dn, delta_dn = jdqmcf.local_update_det_ratio(Gdn, delta_V_dn, i, dtau)

        logabsdet_up_new = prepare_grover_rank1_proposal!(
            tracker_up, Gup_0tau, Gup_tau0, ratio_up, delta_up, i, C,
        )
        logabsdet_dn_new = prepare_grover_rank1_proposal!(
            tracker_dn, Gdn_0tau, Gdn_tau0, ratio_dn, delta_dn, i, C,
        )

        logX_old = tracker_up.logabsdet_value + tracker_dn.logabsdet_value
        logX_new = logabsdet_up_new + logabsdet_dn_new
        ratio_lambda = fused_ratio_from_logs(logX_new, logX_old, lambda)
        probability = min(1.0, abs(ratio_up * ratio_dn) * ratio_lambda)

        if rand(rng) < probability
            accepted += 1
            s[i] = -s[i]

            local_update_Gτ0!(Gup_tau0, Gup, ratio_up, delta_up, i, u, v)
            local_update_G0τ!(Gup_0tau, Gup, ratio_up, delta_up, i, u, v)
            local_update_Gτ0!(Gdn_tau0, Gdn, ratio_dn, delta_dn, i, u, v)
            local_update_G0τ!(Gdn_0tau, Gdn, ratio_dn, delta_dn, i, u, v)

            logdetGup, sgndetGup = jdqmcf.local_update_greens!(
                Gup, logdetGup, sgndetGup, Bup, ratio_up, delta_up, i, u, v,
            )
            logdetGdn, sgndetGdn = jdqmcf.local_update_greens!(
                Gdn, logdetGdn, sgndetGdn, Bdn, ratio_dn, delta_dn, i, u, v,
            )

            accept_grover_rank1_proposal!(tracker_up)
            accept_grover_rank1_proposal!(tracker_dn)
        end
    end

    return logdetGup, sgndetGup, logdetGdn, sgndetGdn, accepted / N
end

function local_update_replica_rank1_fused!(
    Gup::Matrix{T}, logdetGup, sgndetGup, Bup, fermion_greens_calc_up,
    Gdn::Matrix{T}, logdetGdn, sgndetGdn, Bdn, fermion_greens_calc_dn,
    Gup_00, Gup_tau0, Gup_0tau,
    Gdn_00, Gdn_tau0, Gdn_0tau,
    G2up_00, G2dn_00,
    s, mu, alpha, dtau, deltaG, rng,
    C, lambda,
) where {T<:Real}
    Ltau = length(Bup)
    N = size(Gup, 1)

    u = zeros(T, N)
    v = zeros(T, N)
    perm = collect(1:N)

    tracker_up = GroverRank1Tracker(T, length(C))
    tracker_dn = GroverRank1Tracker(T, length(C))
    set_passive_grover_block!(tracker_up, G2up_00, C)
    set_passive_grover_block!(tracker_dn, G2dn_00, C)

    acceptance_rate = 0.0

    jdqmcf.initialize_unequaltime_greens!(Gup_tau0, Gup_0tau, Gup, Gup_00)
    jdqmcf.initialize_unequaltime_greens!(Gdn_tau0, Gdn_0tau, Gdn, Gdn_00)
    rebuild_active_grover!(tracker_up, Gup_00, C)
    rebuild_active_grover!(tracker_dn, Gdn_00, C)

    for l in fermion_greens_calc_up
        jdqmcf.propagate_unequaltime_greens!(
            Gup_tau0, Gup_0tau, Gup,
            fermion_greens_calc_up, Bup,
        )
        jdqmcf.propagate_unequaltime_greens!(
            Gdn_tau0, Gdn_0tau, Gdn,
            fermion_greens_calc_dn, Bdn,
        )

        s_l = @view s[:, l]
        logdetGup, sgndetGup,
        logdetGdn, sgndetGdn,
        acceptance_rate_l = _local_update_replica_rank1_fused!(
            Gup, Gup_tau0, Gup_0tau, logdetGup, sgndetGup, Bup[l],
            Gdn, Gdn_tau0, Gdn_0tau, logdetGdn, sgndetGdn, Bdn[l],
            s_l, mu, alpha, dtau, rng, perm, u, v,
            C, lambda,
            tracker_up, tracker_dn,
        )

        acceptance_rate += acceptance_rate_l / Ltau

        logdetGup, sgndetGup, deltaGup, _ = jdqmcf.stabilize_unequaltime_greens!(
            Gup_tau0, Gup_0tau, Gup, Gup_00,
            logdetGup, sgndetGup,
            fermion_greens_calc_up, Bup,
            update_B̄ = true,
        )
        logdetGdn, sgndetGdn, deltaGdn, _ = jdqmcf.stabilize_unequaltime_greens!(
            Gdn_tau0, Gdn_0tau, Gdn, Gdn_00,
            logdetGdn, sgndetGdn,
            fermion_greens_calc_dn, Bdn,
            update_B̄ = true,
        )

        # Between stabilization boundaries, the Grover weight consumes only
        # G00[C,C], which the trackers update exactly. The framework rebuilds
        # the full G00 matrices from the propagators at these boundaries and
        # at sweep completion, so resynchronize the tracked blocks only then.
        if _stabilization_rebuilds_G00(fermion_greens_calc_up)
            rebuild_active_grover!(tracker_up, Gup_00, C)
            rebuild_active_grover!(tracker_dn, Gdn_00, C)
        end

        deltaG = max(deltaG, deltaGup, deltaGdn)
        iterate(fermion_greens_calc_dn, fermion_greens_calc_up)
    end

    return logdetGup, sgndetGup, logdetGdn, sgndetGdn, deltaG, acceptance_rate
end

function update_replica_rank1_fused!(active, passive, deltaG, C, lambda)
    if lambda == 0
        return update_replica!(active, deltaG)
    end

    logdetGup, sgndetGup,
    logdetGdn, sgndetGdn,
    deltaG, acceptance_rate = local_update_replica_rank1_fused!(
        active.Gup, active.logdetGup, active.sgndetGup,
        active.Bup, active.fermion_greens_calc_up,
        active.Gdn, active.logdetGdn, active.sgndetGdn,
        active.Bdn, active.fermion_greens_calc_dn,
        active.Gup_00, active.Gup_τ0, active.Gup_0τ,
        active.Gdn_00, active.Gdn_τ0, active.Gdn_0τ,
        passive.Gup_00, passive.Gdn_00,
        active.s, active.μ, active.α, active.Δτ,
        deltaG, active.rng,
        C, lambda,
    )

    active = merge(active, (
        logdetGup = logdetGup,
        sgndetGup = sgndetGup,
        logdetGdn = logdetGdn,
        sgndetGdn = sgndetGdn,
    ))

    return active, deltaG, acceptance_rate
end
