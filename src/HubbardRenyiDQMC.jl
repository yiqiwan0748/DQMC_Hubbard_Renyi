module HubbardRenyiDQMC

include("dqmc_kernel.jl")
include("grover_rank1_update.jl")
include("square_lattice.jl")
include("exact_diagonalization.jl")

export GroverRank1Tracker
export accept_grover_rank1_proposal!
export grover_matrix
export hubbard_chain_kinetic_matrix
export initializer
export initializer_from_kinetic
export run_square_dqmc
export square_lattice_kinetic_matrix
export dmrg_ground_state
export ed_ground_state
export ed_thermal
export ed_renyi2
export normalize_boundary_condition
export power_increment_statistics
export prepare_grover_rank1_proposal!
export rebuild_active_grover!
export set_passive_grover_block!
export update_replica!
export update_replica_fused!
export update_replica_rank1_fused!

const DMRG_LOAD_LOCK = ReentrantLock()

"""
    dmrg_ground_state(K; kwargs...)

Compute a Hubbard ground state with ITensorMPS. The optional backend loads on
first use; see `examples/dmrg_square.jl` for a complete small-lattice example.
"""
function dmrg_ground_state(K; kwargs...)
    lock(DMRG_LOAD_LOCK) do
        if !isdefined(@__MODULE__, :HubbardDMRG)
            include(joinpath(@__DIR__, "dmrg.jl"))
        end
    end
    return Base.invokelatest() do
        HubbardDMRG.dmrg_ground_state(K; kwargs...)
    end
end

end
