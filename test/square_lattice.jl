using Test
using LinearAlgebra
using HubbardRenyiDQMC

@testset "Editable square lattice" begin
    K = square_lattice_kinetic_matrix(2, 2)
    @test K == [0 -1 -1 0; -1 0 0 -1; -1 0 0 -1; 0 -1 -1 0]
    @test square_lattice_kinetic_matrix(2, 2; periodic = (false, false)) == K
    @test count(!iszero, square_lattice_kinetic_matrix(3, 2; periodic = (false, false))) == 14
    K4 = square_lattice_kinetic_matrix(4, 4)
    @test all(sum(!iszero, K4; dims = 2) .== 4)
    @test sort(eigvals(Symmetric(K4))) ≈ sort([-2cos(2pi*x/4)-2cos(2pi*y/4) for x in 0:3 for y in 0:3])
    @test_throws ArgumentError square_lattice_kinetic_matrix(0, 2)
    @test_throws ArgumentError square_lattice_kinetic_matrix(2, 2; periodic = [true])
    @test_throws ArgumentError initializer_from_kinetic(K, -1, 0, 2, 0.1, 1234)
    @test_throws ArgumentError initializer_from_kinetic(K, 2, 0, 2, 0.13, 1234)
    @test_throws ArgumentError initializer_from_kinetic([0.0 -1.0; 0.0 0.0], 2, 0, 2, 0.1, 1234)
end

@testset "Free square-lattice DQMC equals Fermi-Dirac result" begin
    beta, mu = 2.0, 0.3
    K = square_lattice_kinetic_matrix(2, 2)
    levels = eigvals(Symmetric(K))
    occupations = 1 ./ (1 .+ exp.(beta .* (levels .- mu)))
    expected_density = 2sum(occupations) / 4
    expected_energy = 2sum((levels .- mu) .* occupations)
    result = run_square_dqmc(Lx = 2, Ly = 2, U = 0, mu = mu, beta = beta,
        dtau = 0.1, burnin = 0, measurements = 4, bins = 2)
    @test result.energy ≈ expected_energy atol = 1e-10
    @test result.density ≈ expected_density atol = 1e-11
    @test result.double_occupancy ≈ (expected_density / 2)^2 atol = 1e-11
    @test result.average_sign == 1
end
