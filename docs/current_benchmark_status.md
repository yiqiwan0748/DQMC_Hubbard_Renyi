# Current Benchmark Status

This note summarizes the current high-stat L=10, U=2 lambda-grid benchmark results.

Reference exact diagonalization value:

```text
S2_ED = 1.716551661433
```

## Beta convergence at dtau = 0.05, Ninc = 16

| Check | S2_DQMC | Error | DeltaS2 | z |
|---|---:|---:|---:|---:|
| beta20 dtau0.05 Ninc16 | 1.716774793264 | 0.003078757854 | +0.000223131831 | +0.072 |
| beta30 dtau0.05 Ninc16 | 1.716525821965 | 0.003094484199 | -0.000025839469 | -0.008 |
| beta40 dtau0.05 Ninc16 | 1.716441955736 | 0.003040377094 | -0.000109705697 | -0.036 |
| beta60 dtau0.05 Ninc16 | 1.720125232459 | 0.002984680310 | +0.003573571025 | +1.197 |

The beta20, beta30, and beta40 results agree with ED essentially at the center of the error bars. The beta60 point is higher than ED, but only by about 1.2 sigma, so it is not statistically significant at the current precision.

## Trotter check at beta = 30, Ninc = 16

| Check | S2_DQMC | Error | DeltaS2 | z |
|---|---:|---:|---:|---:|
| beta30 dtau0.05 Ninc16 | 1.716525821965 | 0.003094484199 | -0.000025839469 | -0.008 |
| beta30 dtau0.10 Ninc16 | 1.712723681943 | 0.002910235085 | -0.003827979490 | -1.315 |

The dtau0.10 result is lower than the dtau0.05 result, but it remains within about 1.3 sigma of ED. At this precision, there is no clear evidence of a Trotter-step failure.

## Current conclusion

All current high-stat checks agree with exact diagonalization within statistical uncertainty. The strongest remaining hints are beta60 being mildly high and dtau0.10 being mildly low, but neither reaches 2 sigma.
