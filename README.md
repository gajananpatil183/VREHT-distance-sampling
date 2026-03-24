# VREHT Distance Sampling Simulation

This repository contains R code for the simulation study presented in the paper:

"Extending the Vertically Reflected Exponentiated Hyperbolic Tangent Detection Function for Grouped Line-Transect Data"

## Description

The code simulates grouped distance sampling data under different detection scenarios and fits the following models:

- Half-Normal (HN)
- Hazard-rate (HR)
- VREHT (proposed model)

Performance is evaluated using:
- Probability of detection (Pa)
- Standard error of Pa (Delta method)
- Absolute Relative Bias (ARB)
- Normalized RMSE (NRMSE)
- Model selection criteria (AIC)

## How to Run

1. Open R / RStudio
2. Run the script:

```r
source("simulation_VREHT_group.R")
