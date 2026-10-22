#!/usr/bin/env julia
using HubbardRenyiDQMC
using TOML
using Printf
using LinearAlgebra

function square_hubbard_main(args = ARGS; run_simulation = run_square_dqmc)
    if isempty(args) || "--help" in args
        println("Run square-lattice Hubbard DQMC: julia --project=. .\\run\\run_square_hubbard.jl .\\parameters\\square_4x4.toml [--output-dir DIR]")
        return
    end
    length(args) in (1, 3) || error("Provide a TOML file, optionally followed by --output-dir DIR")
    length(args) == 3 && args[2] != "--output-dir" && error("Unknown option: $(args[2])")
    cfg = TOML.parsefile(args[1])
    model, sampling = cfg["model"], cfg["sampling"]
    output = length(args) == 3 ? args[3] : cfg["output"]["directory"]
    # Reserve the directory before sampling so concurrent runs cannot share it.
    mkpath(dirname(abspath(output)))
    try
        mkdir(output)
    catch err
        (ispath(output) || islink(output)) &&
            error("Output already exists: $output. Choose a new --output-dir.")
        rethrow(err)
    end
    open(joinpath(output, "parameters.toml"), "w") do io
        TOML.print(io, cfg)
    end
    BLAS.set_num_threads(1)
    result = run_simulation(;
        Lx = model["Lx"], Ly = model["Ly"], t = model["t"], U = model["U"],
        mu = model["mu"], beta = model["beta"], dtau = model["dtau"],
        periodic = model["periodic"], n_stab = model["n_stab"],
        root_seed = sampling["root_seed"], burnin = sampling["burnin"],
        measurements = sampling["measurements"], bins = sampling["bins"],
        sweeps_per_measurement = sampling["sweeps_per_measurement"])
    open(joinpath(output, "summary.csv"), "w") do io
        println(io, "observable,mean,standard_error")
        for name in (:energy, :energy_per_site, :density, :double_occupancy)
            println(io, name, ",", getproperty(result, name), ",",
                getproperty(result, Symbol(name, "_error")))
        end
        println(io, "average_sign,", result.average_sign, ",")
        println(io, "acceptance,", result.acceptance, ",")
        println(io, "deltaG_max,", result.deltaG_max, ",")
    end
    open(joinpath(output, "bins.csv"), "w") do io
        println(io, "bin,sign,energy,density,double_occupancy")
        for bin in 1:result.bins
            values = [getproperty(result.bin_data, key)[bin] for key in
                (:sign, :energy, :density, :double_occupancy)]
            println(io, bin, ",", join(values, ','))
        end
    end
    @printf("Energy per site: %.8f +/- %.8f\n", result.energy_per_site, result.energy_per_site_error)
    @printf("Density: %.8f +/- %.8f\n", result.density, result.density_error)
    @printf("Double occupancy: %.8f +/- %.8f\n", result.double_occupancy, result.double_occupancy_error)
    @printf("Average sign: %.6f; maximum Green-function correction: %.3e\n", result.average_sign, result.deltaG_max)
    println("Saved bins, summary and parameters to ", abspath(output))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && square_hubbard_main()
