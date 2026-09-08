# mateR2 manuscript analysis

Analysis code, figures and summary outputs for:

> Sard, N.M. *mateR2: simulating pedigrees using a bipartite mating network
> framework.* (in review)

The simulation engine lives in a separate package,
[nicksard/mateR2](https://github.com/nicksard/mateR2). This repository holds
only the analyses that produce the manuscript's figures, tables and reported
numbers.

## Running it

```r
renv::restore()                  # installs every dependency, including mateR2
```

```bash
Rscript scripts/00_run_all.R     # steps 01-09
Rscript scripts/00_run_all.R 10  # the weight sweep, opt-in
```

`renv.lock` pins mateR2 to the public commit tagged `v0.1.0` and pins CRAN to a
dated snapshot, so `restore()` rebuilds the environment rather than whatever
happens to be current. It does not install R itself; the analyses were run under
R 4.4.3.

Set `MATER2_WORKERS` to control parallelism (default 4). More than 4 to 6 is
counterproductive: at N_P = 1600 each worker holds four million-row chain
histories in memory.

**Runtime**, on four workers: about 19 minutes from an empty state, or about 6
minutes once `data/outputs/chains/` exists, since step 02 is restartable and
skips chains it has already written. Delete that directory to force
regeneration. The weight sweep in step 10 is a separate 180-job run and is not
part of the default pipeline.

## The pipeline

| Step | Script | Produces |
|---|---|---|
| 01 | `01_setup_scenarios.R` | the 30-scenario x 3-seed-set job grid |
| 02 | `02_run_mcmc_chains.R` | 90 chain files (the long step; restartable) |
| 03 | `03_convergence_stats.R` | convergence summaries, median and envelope |
| 04 | `04_generate_mateR2_figures.R` | Figures 3 and S1 |
| 05 | `05_topology_validation.R` | Figure 4 |
| 06 | `06_sibling_dynamics.R` | Figures 5 and 6 |
| 07 | `07_figure_1_Bipartite_Network_visuals.R` | Figure 1 |
| 08 | `08_figure_2_worked_Example.R` | Figure 2 |
| 09 | `09_sibling_asymmetry_mechanism.R` | Table 1, sibling asymmetry results |
| 10 | `10_weight_sensitivity.R` | Figure S2 (opt-in) |

## What is and is not in the repository

Summary outputs under `data/outputs/` are tracked. Rendered figures and the raw
MCMC chains are not: `figures/` and `data/outputs/chains/` are in `.gitignore`
because both regenerate from the code above, and the chains run to hundreds of
megabytes. Cloning this repository and running the pipeline reproduces every
figure in the manuscript.

`scripts/audit_provenance.R` checks that every output traces to a single
coherent execution: that no output predates the script that produces it or the
installed package build, that every output has a producing script, and that all
90 chains carry the same mateR2 version. It exits non-zero on failure, so it can
gate a release step.

## Reproducibility

The 90 chains reproduce **bit-identically** from a clean state. Re-running the
full grid against a fresh install of mateR2 v0.1.0 returned every R-hat,
effective sample size and demographic error identical to full double precision;
the only difference in any output was recorded wall-clock runtime.

The pipeline has also been run under Linux with R 4.3.3 and Windows with R
4.4.3, giving identical diagnostics to every reported digit. This depends on
mateR2 v0.1.0 or later, where the C++ sampler draws from R's own RNG stream;
earlier versions seeded from the system clock and could not be reproduced from a
script at all.

## Licence

MIT. See `LICENSE`.
