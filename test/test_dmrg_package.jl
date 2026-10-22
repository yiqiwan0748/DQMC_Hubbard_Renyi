using Test
using LinearAlgebra
using HubbardRenyiDQMC

@testset "ITensorMPS Hubbard ground states" begin
    K = square_lattice_kinetic_matrix(2, 1; periodic=(false, false))
    result = dmrg_ground_state(K; U=4.0, nsweeps=4, maxdim=[4, 8])
    @test result.energy ≈ -sqrt(32.0) / 2 atol=1e-12
    @test result.S2 ≈ -log(3 / 8) atol=1e-12
    @test result.density_up ≈ [0.5, 0.5] atol=1e-12
    @test result.density_dn ≈ [0.5, 0.5] atol=1e-12
    @test result.energy_variance < 1e-12

    # A nonzero chemical potential and a one-particle sector test the centered
    # interaction and particle counts independently of the half-filled case.
    one_particle = dmrg_ground_state(K; U=4.0, mu=0.3, Nup=1, Ndn=0,
        nsweeps=4, maxdim=[4, 8])
    @test one_particle.energy ≈ -1.3 atol=1e-12
    @test sum(one_particle.density_up) ≈ 1 atol=1e-12
    @test sum(one_particle.density_dn) ≈ 0 atol=1e-12

    square = square_lattice_kinetic_matrix(2, 2; periodic=(false, false))
    exact = ed_ground_state(square; U=4.0)
    mps = dmrg_ground_state(square; U=4.0, nsweeps=8, maxdim=[8, 16, 32])
    @test mps.energy ≈ exact.energy atol=1e-10
    @test mps.S2 ≈ ed_renyi2(exact, [1, 2]) atol=1e-8

    # Odd lattices must default to half filling, as in the ED interface.
    odd = square_lattice_kinetic_matrix(3, 1; periodic=(false, false))
    odd_exact = ed_ground_state(odd; U=4.0)
    odd_mps = dmrg_ground_state(odd; U=4.0, nsweeps=4, maxdim=[8, 16])
    @test (odd_mps.Nup, odd_mps.Ndn) == (odd_exact.basis.Nup, odd_exact.basis.Ndn) == (2, 1)
    @test sum(odd_mps.density_up + odd_mps.density_dn) ≈ 3 atol=1e-12
    @test odd_mps.energy ≈ odd_exact.energy atol=1e-10
    @test odd_mps.S2 ≈ ed_renyi2(odd_exact, [1]) atol=1e-8

    for particles in (0, 2)
        edge = dmrg_ground_state(K; U=4.0, Nup=particles, Ndn=particles,
            nsweeps=2, maxdim=[4, 8])
        @test edge.energy ≈ 2.0 atol=1e-12
        @test abs(edge.S2) < 1e-12
        @test sum(edge.density_up) ≈ particles atol=1e-12
        @test sum(edge.density_dn) ≈ particles atol=1e-12
    end

    @test_throws ArgumentError dmrg_ground_state([0.0 -1; 0 0])
    @test_throws ArgumentError dmrg_ground_state(K; Nup=3)
    @test_throws ArgumentError dmrg_ground_state(K; cut=0)
    @test_throws ArgumentError dmrg_ground_state(K; nsweeps=1)
end
