# VREHT-CDS-grouped-distance-sampling

This repository contains the R scripts used in the simulation study and analyses for the paper:

**"Extending the Vertically Reflected Exponentiated Hyperbolic Tangent Detection Function for Grouped Line-Transect Data"**

## Description

The code implements a simulation study to evaluate the performance of the proposed **VREHT detection function** in comparison with the classical **Half-Normal (HN)** and **Hazard-Rate (HR)** detection models for grouped line-transect distance sampling data.

The simulation considers a mixture detection function:

0.5 × HN(20) + 0.5 × HR(10, 7)

The competing models fitted to the simulated data are:

- Half-Normal (HN)
- Hazard-Rate (HR)
- VREHT detection function

## How to run

Run the R script:

simulation_mixture_HN_HR.R

The script performs the simulation study and generates the corresponding summary results.

## Author

Gajanan Patil
