using HubbardRenyiDQMC
using LinearAlgebra
using TOML

BLAS.set_num_threads(1)
length(ARGS) == 2 || error("usage: ed_julia.jl cases.toml output.toml")
ispath(ARGS[2]) && error("Output already exists; choose a new evidence path: $(ARGS[2])")
cases = TOML.parsefile(ARGS[1])

function kinetic(case)
    if haskey(case, "K")
        return reduce(vcat, permutedims.(case["K"]))
    end
    Lx, Ly = case["shape"]
    return square_lattice_kinetic_matrix(Lx, Ly; t=case["t"], periodic=case["periodic"])
end

out = Dict{String,Any}("julia_version" => string(VERSION),
    "package_version" => string(pkgversion(HubbardRenyiDQMC)),
    "source_path" => pathof(HubbardRenyiDQMC),
    "source_commit" => readchomp(`git -C $(dirname(dirname(pathof(HubbardRenyiDQMC)))) rev-parse HEAD`),
    "ground" => Any[], "thermal" => Any[])
for case in cases["ground"]
    K = kinetic(case)
    exact = ed_ground_state(K; U=case["U"], mu=case["mu"],
        Nup=case["Nup"], Ndn=case["Ndn"])
    renyi = [ed_renyi2(exact, cut) for cut in case["cuts"]]
    push!(out["ground"], Dict("id" => case["id"], "K" => collect.(eachrow(K)),
        "energy" => exact.energy, "density" => exact.density,
        "double_occupancy" => exact.double_occupancy, "residual" => exact.residual,
        "spectrum" => exact.spectrum, "renyi2" => renyi))
    println("ED ground: ", case["id"], " energy=", exact.energy)
end
for case in cases["thermal"]
    K = kinetic(case)
    exact = ed_thermal(K; U=case["U"], mu=case["mu"], beta=case["beta"])
    push!(out["thermal"], Dict("id" => case["id"], "K" => collect.(eachrow(K)),
        "energy" => exact.energy, "density" => exact.density,
        "double_occupancy" => exact.double_occupancy, "logZ" => exact.logZ))
    println("ED thermal: ", case["id"], " energy=", exact.energy)
end
mkpath(dirname(abspath(ARGS[2])))
open(ARGS[2], "w") do io
    TOML.print(io, out; sorted=true)
end
