#!/usr/bin/env julia
# Official SmoQyDQMC model construction, high-level local updates and native measurements.
# See https://smoqysuite.github.io/SmoQyDQMC.jl/stable/tutorials/hubbard_square/
using SmoQyDQMC, Random
import SmoQyDQMC.LatticeUtilities as lu
import SmoQyDQMC.JDQMCFramework as dqmcf
import SmoQyDQMC.JDQMCMeasurements as dqmcm
include(joinpath(@__DIR__, "dqmc_common.jl"))

function run_smoqy(U, seed)
    beta, dtau, n_stab = 2.0, 0.05, 10
    burnin, measurements, nbins = 1000, 20000, 40
    # A single four-orbital cell represents the open square exactly. This avoids
    # duplicate wrapping bonds for a periodic axis of length two.
    cell = lu.UnitCell(lattice_vecs = [[2.,0.],[0.,2.]],
        basis_vecs = [[0.,0.],[1.,0.],[0.,1.],[1.,1.]])
    # ModelGeometry requires periodic cell coordinates. No bond connects to an
    # adjacent cell, so this one-cell representation still has the open square H.
    lattice = lu.Lattice(L = [1,1], periodic = [true,true])
    geometry = ModelGeometry(cell, lattice)
    bonds = [lu.Bond(orbitals = pair, displacement = [0,0])
        for pair in ((1,2),(1,3),(2,4),(3,4))]
    foreach(bond -> add_bond!(geometry, bond), bonds)
    tbmodel = TightBindingModel(model_geometry = geometry, t_bonds = bonds,
        t_mean = ones(4), t_std = zeros(4), ϵ_mean = zeros(4),
        ϵ_std = zeros(4), μ = 0.0)
    hubmodel = HubbardModel(ph_sym_form = true, U_orbital = collect(1:4),
        U_mean = fill(U,4), U_std = zeros(4))
    rng = Xoshiro(seed)
    tb = TightBindingParameters(tight_binding_model = tbmodel,
        model_geometry = geometry, rng = rng)
    hub = HubbardParameters(model_geometry = geometry, hubbard_model = hubmodel, rng = rng)
    hs = HubbardSpinHirschHST(β = beta, Δτ = dtau,
        hubbard_parameters = hub, rng = rng)
    up = FermionPathIntegral(tight_binding_parameters = tb, β = beta, Δτ = dtau)
    dn = FermionPathIntegral(tight_binding_parameters = tb, β = beta, Δτ = dtau)
    initialize!(up, dn, hub)
    initialize!(up, dn, hs)
    Bup = initialize_propagators(up; symmetric = false, checkerboard = false)
    Bdn = initialize_propagators(dn; symmetric = false, checkerboard = false)
    K = zeros(4,4)
    for b in axes(tb.neighbor_table,2)
        i,j = tb.neighbor_table[:,b]
        K[i,j] -= tb.t[b]
        K[j,i] -= tb.t[b]
    end
    expected = [0. -1 -1 0; -1 0 0 -1; -1 0 0 -1; 0 -1 -1 0]
    K == expected || error("SmoQy hopping matrix differs from open four-bond square")
    calcup = dqmcf.FermionGreensCalculator(Bup, beta, dtau, n_stab)
    calcdn = dqmcf.FermionGreensCalculator(Bdn, beta, dtau, n_stab)
    Gup, Gdn = zeros(4,4), zeros(4,4)
    logup, signup = dqmcf.calculate_equaltime_greens!(Gup, calcup)
    logdn, signdn = dqmcf.calculate_equaltime_greens!(Gdn, calcdn)
    deltaG, deltaTheta = 0.0, 0.0
    signs, energies, densities, doubles = (zeros(nbins) for _ in 1:4)
    binsize = measurements ÷ nbins
    for sweep in 1:(burnin + measurements)
        (_, logup, signup, logdn, signdn, deltaG, deltaTheta) = local_updates!(
            Gup, logup, signup, Gdn, logdn, signdn, hs;
            fermion_path_integral_up = up, fermion_path_integral_dn = dn,
            fermion_greens_calculator_up = calcup, fermion_greens_calculator_dn = calcdn,
            Bup, Bdn, δG_max = 1e-6, δG = deltaG, δθ = deltaTheta, rng,
            update_stabilization_frequency = false)
        sweep <= burnin && continue
        bin = (sweep - burnin - 1) ÷ binsize + 1
        sgn = sign(real(signup * signdn))
        density = dqmcm.measure_n(Gup) + dqmcm.measure_n(Gdn)
        double = dqmcm.measure_double_occ(Gup, Gdn)
        # SmoQy's hopping measurement is per bond of a given type; sum the
        # physical bond and orbital energies before dividing by four sites.
        kinetic = sum(length(tb.hopping_slices[id]) *
            (SmoQyDQMC.measure_hopping_energy(tb, up, Gup, id) +
             SmoQyDQMC.measure_hopping_energy(tb, dn, Gdn, id)) for id in 1:4)
        interaction = sum(SmoQyDQMC.measure_hubbard_energy(hub, Gup, Gdn, id) for id in 1:4)
        energy = real(kinetic + interaction) / 4
        signs[bin] += sgn / binsize
        energies[bin] += sgn * energy / binsize
        densities[bin] += sgn * real(density) / binsize
        doubles[bin] += sgn * real(double) / binsize
    end
    return (; sign = signs, energy_per_site = energies, density = densities,
        double_occupancy = doubles), deltaG, K, string(typeof(Bup[1]))
end

function main(args = ARGS)
    length(args) == 1 || error("Usage: julia --project=SMOQY_ENV dqmc_smoqy.jl OUTPUT")
    BLAS.set_num_threads(1)
    paths = prepare_output(args[1], "smoqy")
    settings = copy(DQMC_SETTINGS)
    settings["root_seeds"] = DQMC_ROOT_SEEDS .+ 500000
    extras = Dict{String,Any}(
        "upstream_source" => "https://github.com/SmoQySuite/SmoQyDQMC.jl",
        "upstream_tutorial" => "https://smoqysuite.github.io/SmoQyDQMC.jl/stable/tutorials/hubbard_square/",
        "unit_cell_orbitals" => [[0.,0.],[1.,0.],[0.,1.],[1.,1.]],
        "computational_cell_periodic" => [true,true],
        "unit_cell_count" => [1,1], "bond_orbitals" => [[1,2],[1,3],[2,4],[3,4]],
        "bond_displacements" => [[0,0] for _ in 1:4],
        "updates" => "official SmoQyDQMC.local_updates!; no reflection updates",
        "measurements" => "native SmoQyDQMC hopping/Hubbard energies and JDQMCMeasurements density/double occupancy")
    open(paths[1], "w") do bins
        println(bins, "engine,U,root_seed,bin,sign,energy_per_site,density,double_occupancy")
        open(paths[2], "w") do summary
            println(summary, "engine,U,root_seed,observable,mean,standard_error,average_sign,deltaG_max")
            for U in (0.0, 2.0), seed in settings["root_seeds"]
                data, deltaG, K, propagator = run_smoqy(U, seed)
                extras["K_rows"] = [collect(K[i,:]) for i in axes(K,1)]
                extras["propagator_type"] = propagator
                write_results(summary, bins, "SmoQyDQMC", U, seed, data, deltaG)
            end
        end
    end
    open(paths[3], "w") do io
        TOML.print(io, metadata("SmoQyDQMC", settings, extras))
    end
end
abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
