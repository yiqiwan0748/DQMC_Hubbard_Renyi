using HubbardRenyiDQMC
using LinearAlgebra
using Printf

"""
Evaluate the finite-time-step partition function by exact many-body matrix
products on a small square. This is a deterministic check of DQMC's
discretized target, separate from the continuous-time Hamiltonian reference.
"""
function exact_trotter_observables(K; U=2.0, mu=0.0, beta=2.0, dtau=0.05)
    n = size(K, 1)
    n <= 6 || error("Use at most six sites for this dense benchmark")
    slices = round(Int, beta/dtau)
    slices >= 1 && isapprox(slices*dtau, beta; atol=1e-12, rtol=0) ||
        error("beta/dtau must be a positive integer")
    partition, energy_trace, double_trace, particle_trace = 0.0, 0.0, 0.0, 0.0
    for Nup in 0:n, Ndn in 0:n
        kinetic, _, _ = HubbardRenyiDQMC._ed_sector(K, 0.0, 0.0, Nup, Ndn, 4096)
        H, double_counts, _ = HubbardRenyiDQMC._ed_sector(K, U, mu, Nup, Ndn, 4096)
        potential = diag(H - kinetic)
        # exp(-dtau*K) exp(-dtau*V) and the reverse ordering give equal
        # partition, K, and V expectation traces by cyclicity. D commutes
        # with V, so its trace is also invariant under reversing the order.
        slice = exp(-dtau*kinetic) * Diagonal(exp.(-dtau .* potential))
        density_unnormalized = slice^slices
        sector_partition = tr(density_unnormalized)
        partition += sector_partition
        energy_trace += tr(H * density_unnormalized)
        double_trace += dot(double_counts, diag(density_unnormalized))
        particle_trace += (Nup + Ndn)*sector_partition
    end
    return (; energy=energy_trace/partition, density=particle_trace/(partition*n),
            double_occupancy=double_trace/(partition*n), logZ=log(partition))
end

function benchmark_ed_trotter(output)
    ispath(output) && error("Output already exists: $output")
    mkpath(dirname(abspath(output)))
    K = square_lattice_kinetic_matrix(2, 2; periodic=(false, false))
    exact = ed_thermal(K; U=2.0, mu=0.0, beta=2.0)
    open(output, "w") do io
        println(io, "Lx,Ly,t,U,mu,beta,dtau,energy,energy_minus_continuous,density,double_occupancy,double_minus_continuous,logZ")
        println(io, join((2,2,1.0,2.0,0.0,2.0,0.0,exact.energy,0.0,exact.density,
                         exact.double_occupancy,0.0,exact.logZ), ','))
        for dtau in (0.1, 0.05, 0.025)
            result = exact_trotter_observables(K; dtau)
            @assert abs(result.density-1.0) < 1e-12
            println(io, join((2,2,1.0,2.0,0.0,2.0,dtau,result.energy,
                             result.energy-exact.energy,result.density,result.double_occupancy,
                             result.double_occupancy-exact.double_occupancy,result.logZ), ','))
            @printf("dtau=%.3f E=%.12f D=%.12f\n", dtau, result.energy, result.double_occupancy)
        end
    end
    println("Wrote exact finite-step reference to ", abspath(output))
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) <= 1 || error("Usage: julia --project=. examples/benchmark_ed_trotter.jl [NEW_OUTPUT_CSV]")
    output = isempty(ARGS) ? joinpath("results", "cv_benchmarks_20261021", "ed_trotter_square.csv") : only(ARGS)
    benchmark_ed_trotter(output)
end
