# MF-PFM-2D-Materials

MATLAB implementation of the multi-flake phase-field model described in:

**Phase-Field Modeling of Multi-Flake Growth and Grain-Boundary Formation in CVD-Grown 2D Materials**

Muhammad Hassaan Ali, Won-Kyu Lee, and Yongwoo Kwon

*RSC Advances* (manuscript under revision)

## Overview

This repository provides MATLAB code for simulating the growth and interaction of multiple two-dimensional (2D) material flakes using a coupled three-field phase-field model.

The model incorporates:

- Phase-field evolution of growing flakes
- Adatom supply, surface diffusion, desorption, and attachment
- Crystallographic orientation evolution
- Growth competition, flake impingement, and grain-boundary formation

The framework enables simulations of multiple flakes with independently prescribed crystallographic orientations.

## Requirements

MATLAB R2022a or later is recommended. The original simulations were performed using MATLAB R2025b.

The core simulation solver uses base MATLAB. Some analysis and figure-generation scripts require the Image Processing Toolbox.

## Repository Structure

| Directory or file | Description |
|---|---|
| `KWC_EulerCN_standalone_v2_Fig2a.m` | Standalone example for Fig. 2(a) |
| `KWC_EulerCN_standalone_v2_Fig6.m` | Standalone example for Fig. 6(a–c) |
| `solver/` | Main simulation solver and standalone implementations |
| `cases/` | Simulation parameters and execution scripts |
| `analysis/` | Morphology and grain-boundary analysis |
| `figures_code/` | Figure-generation scripts |
| `results/` | Available simulation figures and quantitative analysis results |
| `paths_config.m` | MATLAB path configuration |

## Representative MATLAB Examples

Two standalone scripts are provided in the repository root directory.

**Example 1 — Fig. 2(a)**

`KWC_EulerCN_standalone_v2_Fig2a.m`

Simulates the growth of a single star-like flake with threefold anisotropy.

Parameters: D = 5, F = 0.002, tau_phi = 2, and g = 0.1.

**Example 2 — Fig. 6(a–c)**

`KWC_EulerCN_standalone_v2_Fig6.m`

Simulates the collision of two triangular flakes with different crystallographic orientations and the formation of a grain boundary.

Parameters: D = 10, F = 0.002, tau_phi = 5, and g = 0.3.

Both scripts contain the required numerical procedures and local functions.

### Running the Standalone Examples

1. Download the desired MATLAB script.
2. Open it in MATLAB.
3. Run the script.

The simulations display the phase field, adatom concentration, and crystallographic orientation during growth.

## Additional Simulation Cases

The `solver/mfpfm_run.m` function provides the main solver used for the parameter studies and additional simulations.

The `cases/` directory includes:

- `case_params.m` — Model parameters for named simulation cases
- `case_list.m` — Simulation cases grouped by figure or analysis
- `run_case.m` — Execute a selected simulation
- `run_all.m` — Batch execution and postprocessing
- `nucleus_survival.m` — Seed survival analysis

To initialize the MATLAB paths, open MATLAB in the repository root directory and run:

`paths_config();`

To execute a representative test case:

`run_case('BX_20', 'test', true);`

Full simulations can require substantial computation time.

## Analysis and Visualization

The `analysis/` directory contains routines for computing morphology descriptors, grain-boundary statistics, and other quantities reported in the manuscript.

The `figures_code/` directory contains scripts for generating the main figures and supplementary figures.

The `results/` directory provides available output figures and numerical analysis results.

Some figure-generation and analysis routines require saved simulation data. When these data are not available, the corresponding simulation cases must first be executed.

## Numerical Implementation

The governing equations are solved using finite differences with periodic boundary conditions.

The phase-field and orientation-field equations are advanced explicitly. The adatom concentration equation uses a Crank–Nicolson formulation.

The function-based solver also supports FFT-based solution of the periodic concentration equation.

The model parameters and initial conditions can be modified to investigate different growth conditions and flake configurations.

## Reproducibility

This repository provides the numerical implementation, simulation configurations, and analysis scripts supporting the accompanying manuscript.

The two standalone examples offer direct demonstrations of single-flake growth and two-flake impingement. Additional simulations can be reproduced using the parameter definitions and execution scripts, subject to the required computational resources.

The accompanying manuscript and Supplementary Information provide further descriptions of the governing equations, numerical procedures, parameters, and simulation results.
