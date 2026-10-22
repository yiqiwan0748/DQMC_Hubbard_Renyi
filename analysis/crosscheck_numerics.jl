using HubbardRenyiDQMC
using LinearAlgebra
using Test

BLAS.set_num_threads(1)

# Independent full-Fock-space construction with explicit creation/annihilation
# matrices, rather than the package's separate spin-mask hopping routine.
function independent_hubbard(K, U, mu)
    n = size(K, 1)
    dim = 1 << (2n)
    destroy = [zeros(dim, dim) for _ in 1:(2n)]
    for mode in 1:(2n), bits in 0:(dim - 1)
        bit = 1 << (mode - 1)
        iszero(bits & bit) && continue
        parity = isodd(count_ones(bits & (bit - 1))) ? -1.0 : 1.0
        destroy[mode][xor(bits, bit) + 1, bits + 1] = parity
    end
    number = [c' * c for c in destroy]
    H = zeros(dim, dim)
    doubles = zeros(dim, dim)
    particles = sum(number)
    for i in 1:n, j in 1:n
        H .+= K[i, j] .* (destroy[i]' * destroy[j] +
                          destroy[n + i]' * destroy[n + j])
    end
    for i in 1:n
        doubles .+= number[i] * number[n + i]
        H .+= U .* ((number[i] - 0.5I) * (number[n + i] - 0.5I))
    end
    H .-= mu .* particles
    return H, particles, doubles
end

@testset "Independent full-Fock Hubbard reference" begin
    K = [0.1 -0.7 0.2; -0.7 -0.2 -1.1; 0.2 -1.1 0.3]
    n = size(K, 1)
    beta, mu = 0.8, 0.23
    for U in (-2.0, 0.0, 3.1)
        H, particles, doubles = independent_hubbard(K, U, mu)
        for Nup in 0:n, Ndn in 0:n
            indices = [bits + 1 for bits in 0:((1 << (2n)) - 1)
                       if count_ones(bits & ((1 << n) - 1)) == Nup &&
                          count_ones(bits >> n) == Ndn]
            exact = ed_ground_state(K; U, mu, Nup, Ndn)
            @test exact.spectrum ≈ eigvals(Symmetric(H[indices, indices])) atol=2e-12 rtol=0
        end
        eig = eigen(Symmetric(H))
        weights = exp.(-beta .* (eig.values .- first(eig.values)))
        probabilities = weights ./ sum(weights)
        expected_energy = dot(probabilities, eig.values)
        expected_density = dot(probabilities, diag(eig.vectors' * particles * eig.vectors)) / n
        expected_double = dot(probabilities, diag(eig.vectors' * doubles * eig.vectors)) / n
        expected_logZ = log(sum(weights)) - beta * first(eig.values)
        thermal = ed_thermal(K; U, mu, beta)
        @test thermal.energy ≈ expected_energy atol=2e-12 rtol=0
        @test thermal.density ≈ expected_density atol=2e-12 rtol=0
        @test thermal.double_occupancy ≈ expected_double atol=2e-12 rtol=0
        @test thermal.logZ ≈ expected_logZ atol=2e-12 rtol=0
    end
end

@testset "Custom interacting ED-DMRG agreement" begin
    K = [0.1 -0.7 0.2; -0.7 -0.2 -1.1; 0.2 -1.1 0.3]
    for (Nup, Ndn) in ((2, 1), (1, 1), (0, 0), (3, 3), (2, 0))
        exact = ed_ground_state(K; U=3.1, mu=0.23, Nup, Ndn)
        mps = dmrg_ground_state(K; U=3.1, mu=0.23, Nup, Ndn,
            nsweeps=5, maxdim=[8, 16, 32])
        @test mps.energy ≈ exact.energy atol=1e-10 rtol=0
        @test mps.S2 ≈ ed_renyi2(exact, [1]) atol=1e-8 rtol=0
        @test sum(mps.density_up) ≈ Nup atol=1e-12 rtol=0
        @test sum(mps.density_dn) ≈ Ndn atol=1e-12 rtol=0
    end
end

include(joinpath(@__DIR__, "..", "test", "test_dmrg_package.jl"))
