using Test
using LinearAlgebra

@testset "ED: analytic Hubbard dimer" begin
    for t in (0.5, 1.0, 1.7), U in (0.0, 2.0, 4.0, 8.0), mu in (-0.3, 0.0)
        K = [0.0 -t; -t 0.0]
        result = ed_ground_state(K; U, mu, Nup=1, Ndn=1)
        gap = sqrt(U^2 + 16t^2)
        # Centering shifts the conventional half-filled energy by -U/2.
        exact_spectrum = sort([-gap/2, -U/2, U/2, gap/2] .- 2mu)
        @test result.spectrum ≈ exact_spectrum atol=2e-12 rtol=0
        @test result.residual < 2e-12
        @test result.density == 1.0
        @test result.double_occupancy ≈ (1 - U/gap)/4 atol=2e-13 rtol=0
        p = (1 - U/gap)/2
        @test ed_renyi2(result, [1]) ≈ -log((p^2 + (1-p)^2)/2) atol=2e-12 rtol=0
        @test abs(ed_renyi2(result, Int[])) < 2e-12
        @test abs(ed_renyi2(result, [1, 2])) < 2e-12
    end
end

@testset "ED: free square and thermal limits" begin
    K = [0.0 -1 0 -1; -1 0 -1 0; 0 -1 0 -1; -1 0 -1 0]
    beta, mu = 1.3, 0.2
    eigenvalues = eigvals(Symmetric(K)) .- mu
    occupations = 1 ./ (1 .+ exp.(beta .* eigenvalues))
    thermal = ed_thermal(K; U=0.0, beta, mu)
    @test thermal.energy ≈ 2dot(eigenvalues, occupations) atol=2e-12 rtol=0
    @test thermal.density ≈ 2sum(occupations)/4 atol=2e-13 rtol=0
    @test thermal.double_occupancy ≈ (sum(occupations)/4)^2 atol=2e-13 rtol=0
    @test thermal.logZ ≈ 2sum(log1p.(exp.(-beta .* eigenvalues))) atol=2e-12 rtol=0
    infinite_temperature = ed_thermal(K; U=4.0, beta=0.0)
    @test abs(infinite_temperature.energy) < 2e-12
    @test infinite_temperature.density ≈ 1.0 atol=2e-13 rtol=0
    @test infinite_temperature.double_occupancy ≈ 0.25 atol=2e-13 rtol=0
    @test infinite_temperature.logZ ≈ 4log(4) atol=2e-12 rtol=0
    ground = ed_ground_state(K; U=4.0)
    @test ed_renyi2(ground, [1, 3]) ≈ ed_renyi2(ground, [2, 4]) atol=2e-12 rtol=0
    @test ed_renyi2(ground, [1]) ≈ ed_renyi2(ground, [2, 3, 4]) atol=2e-12 rtol=0
end

@testset "ED: input limits" begin
    @test_throws ArgumentError ed_ground_state(zeros(2, 3))
    @test_throws ArgumentError ed_ground_state([0.0 1; 0 0])
    @test_throws ArgumentError ed_ground_state(zeros(2, 2); Nup=3)
    @test_throws ArgumentError ed_ground_state(zeros(4, 4); max_dimension=2)
    @test_throws ArgumentError ed_thermal(zeros(7, 7))
    @test_throws ArgumentError ed_thermal(zeros(2, 2); beta=-1)
    ground = ed_ground_state([0.0 -1; -1 0])
    @test_throws ArgumentError ed_renyi2(ground, [1, 1])
    @test_throws ArgumentError ed_renyi2(ground, [3])
end

@testset "ED: noncontiguous fermionic spatial cut" begin
    K = zeros(6, 6)
    for (i, j) in ((1,2),(2,3),(4,5),(5,6),(1,4),(2,5),(3,6))
        K[i,j] = K[j,i] = -1.0
    end
    one_body = eigen(Symmetric(K))
    occupied = one_body.vectors[:, 1:3]
    sites = [1, 3, 5]
    correlation = (occupied * occupied')[sites, sites]
    probabilities = eigvals(Symmetric(correlation))
    exact_entropy = -2sum(log.(probabilities.^2 .+ (1 .- probabilities).^2))
    ground = ed_ground_state(K; U=0.0)
    @test ed_renyi2(ground, sites) ≈ exact_entropy atol=2e-12 rtol=0
end

@testset "ED: editable hoppings and site energies" begin
    K = [0.1 -0.7 0.2; -0.7 -0.2 -1.1; 0.2 -1.1 0.3]
    beta, mu = 0.8, -0.25
    one_body = eigen(Symmetric(K))
    levels = one_body.values .- mu
    occupations = 1 ./ (1 .+ exp.(beta .* levels))
    local_density_one_spin = abs2.(one_body.vectors) * occupations
    thermal = ed_thermal(K; U=0.0, beta, mu)
    @test thermal.energy ≈ 2dot(levels, occupations) atol=2e-12 rtol=0
    @test thermal.density ≈ 2sum(occupations)/3 atol=2e-12 rtol=0
    @test thermal.double_occupancy ≈ sum(abs2, local_density_one_spin)/3 atol=2e-12 rtol=0
end
