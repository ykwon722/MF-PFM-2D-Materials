# MF-PFM-2D-Materials

MATLAB implementation of the multi-flake phase-field model described in the manuscript:

**Phase-Field Modeling of Multi-Flake Growth and Grain-Boundary Formation in CVD-Grown 2D Materials**

Muhammad Hassaan Ali, Won-Kyu Lee, and Yongwoo Kwon

*RSC Advances* (under revision)

## Overview

This repository provides MATLAB code for simulating the growth and interaction of two-dimensional (2D) material flakes using a coupled three-field phase-field model.

The model incorporates:

- Phase-field evolution of growing flakes
- Adatom supply, surface diffusion, desorption, and attachment
- Crystallographic orientation evolution
- Growth competition, flake impingement, and grain-boundary formation

## MATLAB Examples

Two standalone MATLAB scripts are provided.

| File | Description |
|---|---|
| `KWC_EulerCN_standalone_v2_Fig2a.m` | Single-flake growth with threefold symmetry, corresponding to Fig. 2(a) |
| `KWC_EulerCN_standalone_v2_Fig6.m` | Collision of two triangular flakes with different crystallographic orientations, corresponding to Fig. 6(a–c) |

Each script contains the numerical solver and required local functions.

## How to Run

1. Download the desired `.m` file.
2. Open the file in MATLAB.
3. Run the script.

The scripts display the phase field, adatom concentration, and crystallographic orientation during the simulation.

Simulation parameters, initial conditions, and output settings can be modified directly in the scripts.

## Numerical Implementation

The model employs finite differences with periodic boundary conditions. The phase-field and orientation-field equations are advanced explicitly, while the adatom concentration equation is solved using a Crank–Nicolson scheme.

## Availability

The scripts are provided as representative simulation examples. Additional model parameters, numerical analyses, and simulation results are described in the accompanying manuscript and its Supplementary Information.
