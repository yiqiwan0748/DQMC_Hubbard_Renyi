using HubbardRenyiDQMC
using LinearAlgebra
using Printf
using TOML
using Dates

# Independent closed-form levels of the two-site Hubbard model. The N=2
# singlet energy is Eq. (36) of Carrascal et al., JPCM 27, 393001 (2015).
# The remaining particle sectors follow from one-particle bonding/antibonding
# levels and particle-hole symmetry. No ED matrix enters these references.
function analytic_dimer_thermal(t, U, mu, beta)
    gap = sqrt(U^2 + 16t^2)
    energy, particles, doubles = Float64[], Float64[], Float64[]
    function level(e, N, d, degeneracy=1)
        append!(energy, fill(e - mu*N, degeneracy))
        append!(particles, fill(N, degeneracy))
        append!(doubles, fill(d, degeneracy))
    end
    level(U/2, 0, 0)
    level(-t, 1, 0, 2)
    level(t, 1, 0, 2)
    level(-gap/2, 2, (1 - U/gap)/2)
    level(-U/2, 2, 0, 3)
    level(U/2, 2, 1)
    level(gap/2, 2, (1 + U/gap)/2)
    level(-t, 3, 1, 2)
    level(t, 3, 1, 2)
    level(U/2, 4, 2)
    weights = exp.(-beta .* (energy .- minimum(energy)))
    probabilities = weights ./ sum(weights)
    return (; energy=dot(probabilities, energy),
            density=dot(probabilities, particles)/2,
            double_occupancy=dot(probabilities, doubles)/2,
            logZ=log(sum(weights)) - beta*minimum(energy))
end

function benchmark_ed(output_dir)
    files = ["ed_dimer_ground.csv", "ed_dimer_thermal.csv", "ed_free_square.csv",
             "ed_square_references.csv", "ed_metadata.toml"]
    for file in files
        ispath(joinpath(output_dir, file)) && error("Output already exists: $(joinpath(output_dir, file))")
    end
    mkpath(output_dir)
    ground_errors, thermal_errors, free_errors = Float64[], Float64[], Float64[]
    open(joinpath(output_dir, files[1]), "w") do io
        println(io, "t,U,mu,energy,analytic_energy,energy_abs_error,spectrum_max_abs_error,double_occupancy_abs_error,S2_abs_error,residual")
        for t in (0.5, 1.0, 1.7), U in (0.0, 2.0, 4.0, 8.0), mu in (-0.3, 0.0)
            result = ed_ground_state([0.0 -t; -t 0.0]; U, mu)
            gap = sqrt(U^2 + 16t^2)
            spectrum = sort([-gap/2, -U/2, U/2, gap/2] .- 2mu)
            energy_error = abs(result.energy - first(spectrum))
            spectrum_error = maximum(abs.(result.spectrum .- spectrum))
            double_error = abs(result.double_occupancy - (1 - U/gap)/4)
            p = (1 - U/gap)/2
            entropy_error = abs(ed_renyi2(result, [1]) + log((p^2+(1-p)^2)/2))
            push!(ground_errors, spectrum_error)
            println(io, join((t,U,mu,result.energy,first(spectrum),energy_error,
                             spectrum_error,double_error,entropy_error,result.residual), ','))
        end
    end
    open(joinpath(output_dir, files[2]), "w") do io
        println(io, "t,U,mu,beta,energy,analytic_energy,energy_abs_error,density_abs_error,double_occupancy_abs_error,logZ_abs_error")
        for t in (0.5, 1.0, 1.7), U in (0.0, 2.0, 4.0, 8.0), mu in (-0.3, 0.0), beta in (0.2, 2.0)
            result = ed_thermal([0.0 -t; -t 0.0]; U, mu, beta)
            exact = analytic_dimer_thermal(t, U, mu, beta)
            errors = abs.([result.energy-exact.energy, result.density-exact.density,
                           result.double_occupancy-exact.double_occupancy, result.logZ-exact.logZ])
            append!(thermal_errors, errors)
            println(io, join((t,U,mu,beta,result.energy,exact.energy,errors...), ','))
        end
    end
    K = square_lattice_kinetic_matrix(2, 2; periodic=(false, false))
    open(joinpath(output_dir, files[3]), "w") do io
        println(io, "Lx,Ly,U,mu,beta,energy,free_energy,energy_abs_error,density_abs_error,double_occupancy_abs_error,logZ_abs_error")
        for mu in (-0.3, 0.0, 0.2), beta in (0.2, 1.3, 4.0)
            result = ed_thermal(K; U=0.0, mu, beta)
            levels = eigvals(Symmetric(K)) .- mu
            occupations = 1 ./ (1 .+ exp.(beta .* levels))
            energy = 2dot(levels, occupations)
            errors = abs.([result.energy-energy, result.density-sum(occupations)/2,
                           result.double_occupancy-(sum(occupations)/4)^2,
                           result.logZ-2sum(log1p.(exp.(-beta .* levels)))])
            append!(free_errors, errors)
            println(io, join((2,2,0.0,mu,beta,result.energy,energy,errors...), ','))
        end
    end
    open(joinpath(output_dir, files[4]), "w") do io
        println(io, "Lx,Ly,t,U,mu,beta,ensemble,energy,energy_per_site,density,double_occupancy,S2_prefix_half")
        for (Lx, Ly) in ((2, 1), (2, 2), (3, 2))
            K = square_lattice_kinetic_matrix(Lx, Ly; periodic=(false, false))
            for U in (2.0, 4.0)
                ground = ed_ground_state(K; U)
                entropy = ed_renyi2(ground, 1:(Lx*Ly÷2))
                println(io, join((Lx,Ly,1.0,U,0.0,"Inf","canonical_ground",
                                 ground.energy,ground.energy_per_site,ground.density,
                                 ground.double_occupancy,entropy), ','))
                thermal = ed_thermal(K; U, beta=2.0)
                println(io, join((Lx,Ly,1.0,U,0.0,2.0,"grand_canonical_thermal",
                                 thermal.energy,thermal.energy_per_site,thermal.density,
                                 thermal.double_occupancy,""), ','))
            end
        end
    end
    tolerance = 2e-12
    all(errors -> maximum(errors) < tolerance, (ground_errors,thermal_errors,free_errors)) ||
        error("An ED benchmark exceeded the absolute tolerance $tolerance")
    metadata = Dict(
        "recorded_utc" => Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SSZ"),
        "julia_version" => string(VERSION), "os" => string(Sys.KERNEL),
        "cpu" => first(Sys.cpu_info()).model, "julia_threads" => Threads.nthreads(),
        "blas_threads" => BLAS.get_num_threads(), "blas" => string(BLAS.get_config()),
        "ground_models" => length(ground_errors), "thermal_models" => length(thermal_errors)÷4,
        "free_square_models" => length(free_errors)÷4,
        "ground_spectrum_max_abs_error" => maximum(ground_errors),
        "thermal_observable_max_abs_error" => maximum(thermal_errors),
        "free_square_observable_max_abs_error" => maximum(free_errors),
        "absolute_tolerance" => tolerance,
        "interaction" => "U(n_up-1/2)(n_down-1/2)-mu(n_up+n_down)",
        "reference" => "Carrascal et al., J. Phys.: Condens. Matter 27, 393001 (2015), Eq. 36",
        "reference_url" => "https://arxiv.org/html/1502.02194v1#S2.E36")
    open(joinpath(output_dir, files[5]), "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    @printf("24 dimer spectra: max absolute error %.3e\n", maximum(ground_errors))
    @printf("48 thermal dimer cases: max absolute observable error %.3e\n", maximum(thermal_errors))
    @printf("9 free square cases: max absolute observable error %.3e\n", maximum(free_errors))
    println("Wrote ED evidence to ", abspath(output_dir))
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) <= 1 || error("Usage: julia --project=. examples/benchmark_ed.jl [NEW_OUTPUT_DIRECTORY]")
    output_dir = isempty(ARGS) ? joinpath("results", "cv_benchmarks_20261021") : only(ARGS)
    benchmark_ed(output_dir)
end
