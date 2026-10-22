# Current Scope

The planned implementation work is complete:

- The Grover determinant is updated with the validated rank-one identity.
- The active code and exact Julia environment live in this repository.
- The local launcher resumes valid grids and stops on partial or invalid data.
- The raw-data summarizer validates complete coverage before fitting.
- The scientific result is consolidated in
  [`dqmc_renyi_ee.md`](dqmc_renyi_ee.md).

Larger-size U=2 FSS runs and L=16 subsystem scans are not part of the current
scope. The active production step is the detached U=4 PBC FSS calculation and
ED/DMRG comparison through L=22, specified in
[`u4_fss_L4_to_L22_campaign.md`](u4_fss_L4_to_L22_campaign.md). The earlier U=4
screening contrast remains documented in
[`u4_mott_crossover.md`](u4_mott_crossover.md). Future performance changes must
preserve the proposal-weight identity and pass same-seed
direct-versus-optimized checks.
