# Upper critical dimension for dirty Weyl semimetal to metal quantum phase transitions

### Data and codes accompanying the manuscript "Upper critical dimension for dirty Weyl semimetal-to-metal quantum phase transitions" ([arXiv:2609.30265 (2026)](https://arxiv.org/abs/2609.30265))

## Overview

This repository contains extensive raw data and necessray code and scripts that does kernel polynomial method (KPM) [^1] calculations using `kwant` [^2] and data analysis for our numerical study on the upper-critical dimension for disordered Weyl semimetal-to-metal transitions. 

From our raw data, one can perform several scaling analyses and data collapses. Key scaling analyses include:
  1. BCS like scaling for average density of states (ADOS) at zero energy as a function of disorder strength for 2D Weyl fermions: $\rho(0) ~ \exp(-\lambda / W),$ where $\lambda$ is a non-universal fitting parameter.
  2. ADOS scaling at criticality for $d \geq 3$: $\rho(E) \tilde{}  |E|^{\frac{d}{z} - 1}$, where $z$ is the dynamic exponent.
  3. ADOS at zero energy scaling as a function of the reduced distance from quantum critical point (QCP) for $d \geq 3$: $\rho(0) ~ \delta ^{\beta}, where $\beta$ is the order parameter exponent.
With $z$ and $\beta$ obtained from this series of scaling analyses, one may evaluate the correlation length exponent $\nu$ as $\nu = \beta / (d - z)$ for $d \geq 3$.

With our data as well as the critical exponents obtained from scalings, one can perform the following data collapses to confirm the validity of the analyses:
  1. Data collpase of $\rho(E)\delta^{-\nu (d - z)}$ against $|E|\delta^{-\nu z}$ over full range of $W$ on the largest system size of one certain dimension, and
  2. Finite-size data collapse of $\rho(0)L^{d - z}$ against $\delta L ^{1 / \nu}$ on the metallic side (i.e., for $W > W_{c}$).

## Repository contents

The detailed layout and repository contents are as follows:

```
Upper-critical-dimension-for-dirty-Weyl-semimetal-to-metal-quantum-phase-transitions/
├── README.md                                      # this file
├── Dirac_Weyl_KPM/                                # the directory containing KPM code
|   |── Dirac_Weyl_KPM_ver0.6                      # the KPM code 
|   |   |── datafiles/                             # the directory containing data files produced by our code
|   |   |── errorfiles/                            # Any numerically erroneous data file is stored in this directory 
|   |   |── tempfiles/                             # the directory that stores temporarily data files from bins of disorder realizations
|   |   |── core_routines.py                       # Core KPM routines and other important numerical routines
|   |   |── Dirac_Weyl_main.py                     # the main python file
|   |   |── Hamiltonians.py                        # the Hamiltonians
|   |   |── input.yaml                             # the input file for a KPM calculation
|   |   └── postroutines.py                        # post-routines for writing ADOS and some other information
|   |── Dirac_Weyl_KPM_MPI_ver0.7.2                # the embarrassingly-parallelized version of KPM code with MPI
|   |   |── datafiles/
|   |   |── errorfiles/
|   |   |── tempfiles/
|   |   |── core_routines.py
|   |   |── Dirac_Weyl_main.py
|   |   |── Hamiltonians.py
|   |   └── input.yaml
|   └────── postroutines.py
└── Weyl_semimetal_upper_critical_dim.zip          # the raw data
```


## Filename conventions


## Data format


## Numerical methods


## Using the data

## Citations:
[^1]: A. Weiße, G. Wellein, A. Alvermann, and H. Fehske, The kernel polynomial method, Rev. Mod. Phys. 78, 275 (2006)
[^2]: C. W. Groth, M. Wimmer, A. R. Akhmerov, and X. Waintal, Kwant: a software package for quantum transport, New J. Phys. 16, 063065 (2014)
