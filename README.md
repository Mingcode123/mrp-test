# Comparing Small Area Estimation Methods Using the Household Pulse Survey

## Overview

This project uses Household Pulse Survey (HPS) Week 1 public-use microdata to compare small area estimation methods under unequal-probability sampling.

I treat the cleaned full dataset as a known finite population. In each replication, I draw approximately one-fifteenth of the records, estimate the proportion of respondents expecting employment income loss in each state, and compare the estimates with the known full-data proportions.

Using Horvitz–Thompson (HT), non-spatial Bayesian (NSB), and spatial Bayesian (SB) methods as benchmarks, I compare six multilevel regression and poststratification (MRP) specifications. The analysis focuses on three questions:

- Can different MRP interaction structures improve state-level estimation?
- Is there a clear relationship between state sample size and estimation error?
- How do the methods compare in accuracy and computational cost?

## My contributions

- Cleaned and prepared 71,213 survey records and constructed poststratification cells defined by state, age group, sex, and race/ethnicity.
- Implemented repeated unequal-probability sampling and calculated inclusion probabilities, inverse-probability weights, and scaled weights for Bayesian model fitting.
- Specified and fitted six MRP models using `rstanarm`, comparing demographic and state-by-demographic interactions.
- Implemented HT estimation and integrated published NSB and SB Stan models into a common sampling, fitting, and poststratification workflow.
- Evaluated mean squared error, bias, posterior interval coverage, and model-fitting time; interpreted the results through bias–variance decomposition, sample-size comparisons, and geographic visualizations.
- Saved intermediate fits and estimates so repeated computation could run in batches and resume from unfinished replications.

The NSB and SB Stan model definitions use Sun, Parker, and Holan’s public implementation, credited below.

## Data and sampling design

The cleaned dataset contains 71,213 respondents across the 48 contiguous states and Washington, D.C. The binary outcome indicates whether a respondent expects employment income loss.

Each replication draws approximately 4,748 respondents using the Midzuno method. The sampling size measure combines the original survey weight and the outcome, giving respondents unequal inclusion probabilities.

HT uses inverse inclusion probabilities directly. NSB, SB, and MRP use weights scaled to sum to the sample size during model fitting.

Poststratification cells cross-classify state, five age groups, sex, and five race/ethnicity categories. Cell-level predicted probabilities are aggregated to states using cell counts from the full dataset.

The evaluation target is the known state-level proportion within the cleaned survey records. Poststratification counts come from the same full dataset; no external population table is used.

## Methods

The three benchmarks are:

- **HT:** Direct estimation using inverse inclusion probabilities and known finite-population state totals.
- **NSB:** A weighted Bayesian logistic model with fixed demographic effects and state random effects.
- **SB:** An extension of NSB with spatial effects based on state adjacency.

Baseline MRP includes varying intercepts for state, age group, sex, and race/ethnicity. The six specifications are:

| Model  | Additional varying interaction effects                |
| ------ | ----------------------------------------------------- |
| MRP    | None                                                  |
| MRP2   | State × race/ethnicity                                |
| MRP2.1 | Age × race/ethnicity; state × race/ethnicity          |
| MRP3   | Age × race/ethnicity; sex × race/ethnicity            |
| MRP4   | Age × race/ethnicity; sex × race/ethnicity; age × sex |
| MRP5   | All interactions in MRP4, plus state × race/ethnicity |

This design compares lower-dimensional demographic interactions with higher-dimensional state-by-demographic interactions. For example, state × race/ethnicity contains 245 combinations, compared with 25 for age × race/ethnicity.

## Main findings

The final report compares the first 10 replications across all nine methods using MSE, average signed bias, Bayesian posterior interval coverage, and fitting time summarized from recorded Stan chain timings.

- **Model-based estimates reduced MSE.** NSB, SB, and the MRP specifications all had lower MSE than HT, with state-level estimates more concentrated around the full-data proportions.
- **NSB and MRP produced similar estimates.** They share the main demographic information and use poststratification to obtain state estimates. In this dataset, their different parameterizations did not produce substantial differences in estimates.
- **More complex interactions did not clearly improve MSE.** Some specifications reduced average signed bias without reducing overall MSE. Further bias–variance decomposition showed why average bias alone was insufficient to assess model performance.
- **State sample size showed no clear monotonic relationship with error.** States with larger samples did not necessarily have smaller errors, and sample size alone did not explain differences in estimation performance across states.
- **Computational costs differed substantially.** SB had the lowest MSE but required longer average fitting time. NSB was the fastest Bayesian method and achieved accuracy similar to the MRP specifications.

These findings are limited to this dataset and the 10-replication comparison. Further work could increase the number of replications, examine stronger demographic heterogeneity and spatial structure, and incorporate spatial effects into MRP. None of the six current MRP specifications includes spatial adjacency effects.

## Repository contents

| File or directory | Purpose |
| --- | --- |
| `dataCleanHPS.R` | Clean data and construct poststratification cells and state-level targets |
| `sampleIds.R` | Generate repeated samples and their inclusion-probability weights |
| `MRP.R` | Fit six MRP specifications and perform state-level poststratification |
| `sunMethods.R` | Run HT, NSB, and SB benchmarks |
| `plot.R` | Evaluate saved estimates and generate comparison figures and summary tables |
| `stan/` | Published Stan definitions used for NSB and SB |
| `data/` | Survey data, processed data, and state adjacency information |
| `results/` | Saved samples, estimates for all nine methods, and selected fitted model objects |
| `results/method_chain_runtime.csv` | Recorded Stan chain timings used for computational-cost comparisons |
| `results/plot/` | Comparison figures and summary tables for the final analysis |

## Final analysis and repository version

The final report compares the first 10 replications across all nine methods. The repository includes the saved estimates, plotting script, figures, and summary tables for this comparison.

Extended runs are also retained for selected methods, with the main sampling and fitting scripts configured for up to 100 replications. The saved estimate tables contain:

| Method | Saved replications |
| --- | --- |
| HT | 100 |
| NSB | 100 |
| SB | 29 |
| MRP | 25 |
| MRP2 | 25 |
| MRP2.1 | 50 |
| MRP3 | 16 |
| MRP4 | 10 |
| MRP5 | 10 |

The plotting script sets `B_eval = 10` and uses replications 1–10 for every method, matching the comparison in the final report. Additional saved replications are not included in the reported findings.

## Running the workflow

Run all scripts from the repository root.

To regenerate figures and summary tables from the saved results, install `tidyverse` and `maps`, then run `plot.R`. It reads the included processed data, sample IDs, estimate tables in `.rds` format, and recorded chain timings, and writes outputs to `results/plot/`. No model refitting is needed for this step.

To run the sampling and model-fitting workflow, the additional packages are `sampling`, `rstan`, and `rstanarm`. The Stan benchmarks require a working compilation toolchain. Run the scripts in this order:

1. `dataCleanHPS.R`
2. `sampleIds.R`
3. `sunMethods.R`
4. `MRP.R`
5. `plot.R`

Review `B`, `rep_start`, and `rep_end` before running. The fitting scripts skip replications already present in saved estimate files, so running the current checkout may reuse existing results. For a complete refit, use a separate working copy with fresh result directories.

Runtime comparisons use the supplied `results/method_chain_runtime.csv`. The fitting scripts do not automatically refresh this table; comparisons based on new fits require updated chain timings.

## Data and model sources

The data are from the U.S. Census Bureau’s Household Pulse Survey Week 1 public-use microdata.

NSB and SB use `Binomial.stan` and `Binomial_ICAR.stan`, respectively, from Sun, Parker, and Holan’s implementation accompanying *Analysis of Household Pulse Survey Public-Use Microdata via Unit-Level Models for Informative Sampling*.

[Original Stan model repository](https://github.com/QuarkofDorothy/Analysis-of-HPS-Public-Use-Microdata-via-Unit-Level-Models-for-Informative-Sampling)