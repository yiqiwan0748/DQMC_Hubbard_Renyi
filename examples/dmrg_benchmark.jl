#!/usr/bin/env julia
using HubbardRenyiDQMC
using LinearAlgebra
using Printf

function benchmark(output)
    (ispath(output) || islink(output)) && error("Output already exists: $output")
    U, t = 4.0, 1.0
    cases = [
        (name="dimer", Lx=2, Ly=1, sweeps=8, maxdim=[4, 8], seed=1234),
        (name="square_2x2", Lx=2, Ly=2, sweeps=8, maxdim=[8, 16, 32, 64], seed=1234),
        (name="square_2x2_convergence", Lx=2, Ly=2, sweeps=12, maxdim=[16, 32, 64, 128], seed=5678),
        (name="rectangle_3x2", Lx=3, Ly=2, sweeps=8, maxdim=[16, 32, 64, 128], seed=1234),
        (name="rectangle_3x2_convergence", Lx=3, Ly=2, sweeps=12, maxdim=[32, 64, 128, 256], seed=5678),
    ]
    mkpath(dirname(abspath(output)))
    open(output, "w") do io
        println(io, "case,Lx,Ly,U,t,Nup,Ndn,standard,reference_energy,dmrg_energy,abs_energy_error,reference_S2,dmrg_S2,abs_S2_error,cut,nsweeps,maxdim_cap,maxlinkdim,final_sweep_energy_change,energy_variance,final_truncation_error,root_seed,dmrg_elapsed_seconds")
        for c in cases
            N = c.Lx * c.Ly
            K = square_lattice_kinetic_matrix(c.Lx, c.Ly; t, periodic=(false, false))
            exact = ed_ground_state(K; U, Nup=N÷2, Ndn=N÷2)
            exact_S2 = ed_renyi2(exact, collect(1:(N÷2)))
            standard = "LAPACK_exact_diagonalization"
            reference_energy = exact.energy
            if c.name == "dimer"
                standard = "Carrascal_2015_exact_Hubbard_dimer"
                reference_energy = -sqrt(U^2 + 16t^2) / 2
                ionic_weight = (1 - U / sqrt(U^2 + 16t^2)) / 2
                exact_S2 = -log((ionic_weight^2 + (1 - ionic_weight)^2) / 2)
                @assert abs(exact.energy - reference_energy) < 1e-12
            end
            result = dmrg_ground_state(K; U, cut=N÷2, nsweeps=c.sweeps,
                maxdim=c.maxdim, root_seed=c.seed, outputlevel=1)
            row = (c.name, c.Lx, c.Ly, U, t, result.Nup, result.Ndn, standard,
                reference_energy, result.energy, abs(result.energy - reference_energy),
                exact_S2, result.S2, abs(result.S2 - exact_S2), result.cut,
                result.nsweeps, last(c.maxdim), result.maxlinkdim,
                result.sweep_energy_change, result.energy_variance,
                result.final_truncation_error, result.root_seed, result.elapsed_seconds)
            println(io, join(row, ','))
            flush(io)
            @printf("%s: energy error %.3e; S2 error %.3e\n", c.name, row[11], row[14])
            @assert row[11] < 1e-8 "Ground-state energy did not converge to ED"
            @assert row[14] < 1e-7 "Renyi entropy did not converge to ED"
        end
    end
    println("Wrote ", abspath(output))
end

function main(args=ARGS)
    length(args) == 1 || error("Usage: julia --project=. examples/dmrg_benchmark.jl NEW_OUTPUT.csv")
    benchmark(args[1])
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
