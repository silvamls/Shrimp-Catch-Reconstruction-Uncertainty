# shrimp-catch-reconstruction-uncertainty

Data, code and results of a study on how the **catch-reconstruction (gap-filling) method** changes the reconstructed catch series of **Brazilian marine shrimp** and the resulting **ZBRT depletion indicator (B/K)**.

The unit of analysis is the **total shrimp catch of each of the 17 coastal states**. Regions (North, Northeast, Southeast, South) and Brazil are sums of states, used only for summaries and figures. The study does not assess stocks or species.

## Questions and steps

| Step | Question | Main output |
|---|---|---|
| 1. Compare published reconstructions | Do Freire et al. (2021), Silva and Andrade (2023) and MPA (2026) agree in the state-years they had to fill (years without an observed value)? | `Table_1`, `Table_2`, `Fig_1`–`Fig_3`, `Fig_9` |
| 2. Hide and refill | How large is the error of four gap-filling methods when known years are hidden, and is it larger at the start or end of a series than inside it (edge-error hypothesis)? | `Table_3`, `Table_S1`, `Fig_4` |
| 3. Rebuild series and apply ZBRT | Do the four methods change the B/K estimated by ZBRT by more than ZBRT's own 95% interval, in more than half of the years in more than half of the states (divergence hypothesis)? | `Table_4`, `Table_S2`, `Fig_5`–`Fig_8` |

The four gap-filling methods are applied identically to every state series, using only the years of the series itself:

- **M0 naive**: the nearest observed value is carried over.
- **M1 interpolation**: linear interpolation inside gaps and a LOESS curve at the edges (family of Freire et al. 2021).
- **M2 gamma GLM**: gamma GLM with log link, year² and sin(year), chosen by AIC (family of Silva and Andrade 2023).
- **M3 neural network**: small feedforward network with year as the only input.

Observed catch = national bulletins 1962–2007 plus MPA collected landings 2023–2024. The analysis plan is in `analysis_plan_OSF_EN.docx` (also as PDF).

## Repository structure

```
analysis_plan_OSF_EN.docx / .pdf        statistical analysis plan
Code and results/
  Catch_reconstruction_uncertainty.R    the single script for the whole study
  FINAL_RECONSTRUCTED_Brazil_1950_2015_CommercialEtapaII_04072021.xlsx   Freire et al. (2021) database
  Banco_reconstrucao_marinha_2026.xlsx  MPA (2026) reconstruction 1950-2022 and collected landings 2023-2024
  data_AL.csv ... data_SP.csv           Silva and Andrade (2023) GLM outputs by state (17 databases)
  Camaroes_final-3.csv                  observed bulletin landings by state and year, 1962-2007
  data/                                 unified databases built by the script (observed and reconstructed catch by state-year, coverage, rebuilt series)
  models/                               fitted models saved by the script (reused on later runs)
  results/                              tables (.xlsx), figures (.png at 300 dpi and .pdf) and session_info.txt
```

## How to run

1. Open R with the working directory set to `Code and results`.
2. Install the packages listed below.
3. Run `Catch_reconstruction_uncertainty.R`. The script builds the unified database, runs Steps 1 to 3, and writes `data/`, `models/` and `results/`. Every step, parameter and choice is explained in comments inside the script. Models already saved in `models/` are reused; delete them to refit.

## Computational environment

R 4.6.1 on Windows 11 (full details in `results/session_info.txt`).

Packages: `data.table`, `readxl`, `openxlsx`, `ggplot2`, `lme4`, `lmerTest`, `emmeans`, `callr`, `keras` and `tensorflow` (neural network, with a Python backend through `reticulate`), and `datalimited2` (ZBRT, Zhou et al. 2017).

## Data sources and limitations

- The observed bulletins were not standardised across periods: the institution responsible and the collection method changed over time, and the bulletins ended in 2007 (state totals only in 2008–2011). Years with all zeros in the bulletins are treated as no data.
- The 2023–2024 MPA values come from a collection that differs from the earlier ones.
- The reconstructions of Freire et al. (2021), Silva and Andrade (2023) and MPA (2026) are third-party or earlier products and are included here to make the analysis reproducible. Please cite the original sources when using them.
- A zero in the Freire et al. and MPA series means that the source has no record for that state-year.
- The four methods are deliberately simple and use no covariates, so that they are comparable; more elaborate specifications could change the ranking.

## License

Code, analysis plan and results: CC-BY 4.0, as stated in the script header.

## Citation

A citation will be added when the paper is published.

## Contact

Matheus L. S. Silva, PGRPAq, Universidade Federal Rural de Pernambuco (UFRPE), CAMMARADA project.
