#!/usr/bin/env julia
using HubbardRenyiDQMC
using TOML
using Printf

function reserve_dmrg_output(output)
    mkpath(dirname(output))
    try
        return Base.Filesystem.open(output,
            Base.JL_O_WRONLY | Base.JL_O_CREAT | Base.JL_O_EXCL, 0o644)
    catch err
        (ispath(output) || islink(output)) && error("Output already exists: $output")
        rethrow(err)
    end
end

function main(args=ARGS; solve_ground_state=dmrg_ground_state)
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
    output_io = output === nothing ? nothing : reserve_dmrg_output(output)
    try
        result = solve_ground_state(K; U=model["U"], mu=model["mu"],
            Nup=sector["Nup"], Ndn=sector["Ndn"], cut=solver["cut"],
            nsweeps=solver["nsweeps"], maxdim=solver["maxdim"], cutoff=solver["cutoff"],
            root_seed=solver["root_seed"], outputlevel=1)
        function write_result(io)
            println(io, "Lx,Ly,Nup,Ndn,U,mu,energy,S2,cut,maxlinkdim,nsweeps,final_sweep_energy_change,energy_variance,final_truncation_error,root_seed,t,periodic_x,periodic_y,cutoff,maxdim_schedule")
            println(io, join((lattice["Lx"], lattice["Ly"], result.Nup, result.Ndn,
                model["U"], model["mu"], result.energy, result.S2, result.cut,
                result.maxlinkdim, result.nsweeps, result.sweep_energy_change,
                result.energy_variance, result.final_truncation_error, result.root_seed,
                lattice["t"], lattice["periodic"]..., solver["cutoff"],
                join(solver["maxdim"], ';')), ','))
        end
        if output === nothing
            write_result(stdout)
        else
            buffer = IOBuffer()
            write_result(buffer)
            write(output_io, take!(buffer))
            println("Wrote ", output)
        end
        @printf("Energy = %.12f; S2(1:%d) = %.12f\n", result.energy, result.cut, result.S2)
        return result
    finally
        output_io === nothing || close(output_io)
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
