#!/usr/bin/env julia
# Run the published ITensorMPS wrapper, without changing its implementation.
using HubbardRenyiDQMC
using ITensorMPS
using ITensors
using LinearAlgebra
using Statistics
using TOML

length(ARGS) == 1 || error("Usage: julia --project=output/github_release analysis/external_packages/dmrg_public.jl NEW_OUTPUT.csv")
output = abspath(only(ARGS))
ispath(output) && error("Output already exists: $output")
mkpath(dirname(output))
BLAS.set_num_threads(1)
cfg = TOML.parsefile(joinpath(@__DIR__, "dmrg_cases.toml"))
revision = readchomp(`git -C $(dirname(dirname(pathof(HubbardRenyiDQMC)))) rev-parse HEAD`)
header = "case,setting,package,package_revision,Lx,Ly,t,U,mu,Nup,Ndn,cut,root_seed,nsweeps,maxdim_cap,cutoff,energy,S2,density,double_occupancy,maxlinkdim,energy_variance,sweep_energy_change,final_truncation_error,elapsed_seconds,ed_energy,ed_S2,ed_double_occupancy,ed_gap"
io = Base.Filesystem.open(output, Base.JL_O_WRONLY | Base.JL_O_CREAT | Base.JL_O_EXCL, 0o644)
try
    println(io, header)
    flush(io)
    for case in cfg["cases"]
        K = square_lattice_kinetic_matrix(case["Lx"], case["Ly"]; t=case["t"], periodic=(false, false))
        exact = ed_ground_state(K; U=case["U"], mu=case["mu"], Nup=case["Nup"], Ndn=case["Ndn"])
        exact_s2 = ed_renyi2(exact, 1:case["cut"])
        for settings in cfg["settings"]
            println("Public DMRG: ", case["name"], " / ", settings["name"])
            flush(stdout)
            result = dmrg_ground_state(K; U=case["U"], mu=case["mu"], Nup=case["Nup"], Ndn=case["Ndn"],
                cut=case["cut"], nsweeps=settings["nsweeps"], maxdim=settings["maxdim"],
                cutoff=settings["cutoff"], root_seed=settings["root_seed"], outputlevel=0)
            double_occupancy = mean(expect(result.state, "Nupdn"))
            density = mean(result.density_up + result.density_dn)
            row = (case["name"], settings["name"], "HubbardRenyiDQMC_ITensorMPS", revision,
                case["Lx"], case["Ly"], case["t"], case["U"], case["mu"], case["Nup"], case["Ndn"], case["cut"],
                settings["root_seed"], settings["nsweeps"], last(settings["maxdim"]), settings["cutoff"],
                result.energy, result.S2, density, double_occupancy, result.maxlinkdim,
                result.energy_variance, result.sweep_energy_change, result.final_truncation_error,
                result.elapsed_seconds, exact.energy, exact_s2, exact.double_occupancy,
                exact.spectrum[2] - exact.spectrum[1])
            println(io, join(row, ','))
            flush(io)
        end
    end
finally
    close(io)
end
println("Wrote ", output)
