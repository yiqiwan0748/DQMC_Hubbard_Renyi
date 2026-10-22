#!/usr/bin/env julia
using HubbardRenyiDQMC
using TOML
using Printf

function main(args=ARGS)
    if "--help" in args || "-h" in args
        println("Usage: julia --project=. examples/dmrg_square.jl [settings.toml] [--output new.csv]")
        println("Edit lattice dimensions, boundary conditions, particle counts, and DMRG settings in the TOML file.")
        return
    end
    input = joinpath(@__DIR__, "dmrg_square.toml")
    output = nothing
    i = 1
    while i <= length(args)
        if args[i] == "--output"
            i += 1
            i <= length(args) || error("--output requires a new CSV path")
            output = abspath(args[i])
        elseif startswith(args[i], "--")
            error("Unknown option: $(args[i])")
        elseif i == 1
            input = args[i]
        else
            error("Unexpected argument: $(args[i])")
        end
        i += 1
    end
    output === nothing || (!ispath(output) && !islink(output)) || error("Output already exists: $output")
    cfg = TOML.parsefile(input)
    lattice, model, sector, solver = cfg["lattice"], cfg["model"], cfg["sector"], cfg["dmrg"]
    K = square_lattice_kinetic_matrix(lattice["Lx"], lattice["Ly"];
        t=lattice["t"], periodic=Tuple(lattice["periodic"]))
    result = dmrg_ground_state(K; U=model["U"], mu=model["mu"],
        Nup=sector["Nup"], Ndn=sector["Ndn"], cut=solver["cut"],
        nsweeps=solver["nsweeps"], maxdim=solver["maxdim"], cutoff=solver["cutoff"],
        root_seed=solver["root_seed"], outputlevel=1)
    function write_result(io)
        println(io, "Lx,Ly,Nup,Ndn,U,mu,energy,S2,cut,maxlinkdim,nsweeps,final_sweep_energy_change,energy_variance,final_truncation_error,root_seed")
        println(io, join((lattice["Lx"], lattice["Ly"], result.Nup, result.Ndn,
            model["U"], model["mu"], result.energy, result.S2, result.cut,
            result.maxlinkdim, result.nsweeps, result.sweep_energy_change,
            result.energy_variance, result.final_truncation_error, result.root_seed), ','))
    end
    if output === nothing
        write_result(stdout)
    else
        mkpath(dirname(output))
        open(write_result, output; write=true, create=true, truncate=false)
        println("Wrote ", output)
    end
    @printf("Energy = %.12f; S2(1:%d) = %.12f\n", result.energy, result.cut, result.S2)
    return result
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
