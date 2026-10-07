# Upper critical dimension for dirty Weyl semimetal to metal quantum phase transitions

### Data and codes accompanying the manuscript "Upper critical dimension for dirty Weyl semimetal-to-metal quantum phase transitions" ([arXiv:2609.30265 (2026)](https://arxiv.org/abs/2609.30265))

## Overview

This repository contains extensive raw data and necessray code and scripts that does kernel polynomial method (KPM) [^1] calculations using `kwant` [^2] and data analysis for our numerical study on the upper-critical dimension for disordered Weyl semimetal-to-metal transitions. 

From our raw data, one can perform several scaling analyses and data collapses. Key scaling analyses include:
  1. BCS like scaling for average density of states (ADOS) at zero energy as a function of disorder strength for 2D Weyl fermions: $\rho(0) ~ \exp(-\lambda / W),$ where $\lambda$ is a non-universal fitting parameter.
  2. ADOS linear scaling for $W$ far from $W_{c}$ (for $d \geq 3$) and $W \ll 1$ for $d = 2$, in order to confirm $\rho(E) \sim |E|^{d - 1}$ in the clean limit. 
  3. ADOS scaling at criticality for $d \geq 3$: $\rho(E) \sim |E|^{\frac{d}{z} - 1}$, where $z$ is the dynamic exponent.
  4. ADOS at zero energy scaling as a function of the reduced distance from quantum critical point (QCP) for $d \geq 3$: $\rho(0) \sim \delta ^{\beta}, where $\beta$ is the order parameter exponent.
With $z$ and $\beta$ obtained from this series of scaling analyses, one may evaluate the correlation length exponent $\nu$ as $\nu = \beta / (d - z)$ for $d \geq 3$.

With our data as well as the critical exponents obtained from scalings, one can perform the following data collapses to confirm the validity of the analyses:
  1. Data collpase of $\rho(E)\delta^{-\nu (d - z)}$ against $|E|\delta^{-\nu z}$ over full range of $W$ on the largest system size of one certain dimension, and
  2. Finite-size data collapse of $\rho(0)L^{d - z}$ against $\delta L ^{1 / \nu}$ on the metallic side (i.e., for $W > W_{c}$).

## Repository contents and layouts

The detailed layout and repository contents are as follows:

```
Upper-critical-dimension-for-dirty-Weyl-semimetal-to-metal-quantum-phase-transitions/
├── README.md                                      # this file
|
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
|
└── Weyl_semimetal_upper_critical_dim.zip          # the raw data and scripts for post-precessing/scaling analysis 
```

All the raw data and scripts for post-precessing/scaling analysis are zipped in `Weyl_semimetal_upper_critical_dim.zip`. After unzipping, the layout is:

```
Weyl_semimetal_upper_critical_dim/
├── ADOS_2D/                                       # The master directory that stores all raw data for 2D Weyl fermions (here, L = 2500)
|   ├── L_2500/                                    # The sub-directory that stores all data for one particular system size (here, L = 2500)
|   |   ├── datafiles_x/                           # The directory that stores all data files for one particular job, where "x" (as a three-digit integer) stands for a job index. 
|   |   ├── ...
|   |   └── job_table_ADOS_2D_L2500.xlsx           # The spreadsheet that maps the job to a disorder strength, and some other key input parameters
|   └── preliminary_results/                       # The master directory that stores system size-dependent ADOS vs. W data. 
|       └── ADOS_vs_W_L2500/                       # The sub-directory that stores ADOS vs. W data for a specific system size (here, L = 2500)
|           └── ADOS_vs_W_L2500.txt                # the ADOS vs. W data
| 
├── ADOS_3D/                                       # The master directory that stores all raw data for 3D Weyl fermions 
|   ├── L_0100/                                    # The sub-directory that stores all data for one particular system size (here, L = 100)
|   |   ├── datafiles_x/                           # The directory that stores all data files for one particular job, where "x" (as a three-digit integer) stands for a job index. 
|   |   ├── ...
|   |   └── job_table_ADOS_3D_L0100.xlsx           # The spreadsheet that maps the job to a disorder strength, and some other key input parameters
|   ├── L_0120/                                    # The master directory that stores all raw data for 3D Weyl fermions (here, L = 120)
|   |   └── ...
|   ├── ...
|   └── preliminary_results/                       # The master directory that stores system size-dependent ADOS vs. W data. 
|       └── ...
|
├── ADOS_4D/                                       # The master directory that stores all raw data for 4D Weyl fermions 
|   └── ...
|
├── ADOS_5D/                                       # The master directory that stores all raw data for 5D Weyl fermions 
|   └── ...
|
├── ADOS_6D/                                       # The master directory that stores all raw data for 6D Weyl fermions 
|   └── ...
|
├── Comprehensive_plots/                           # The directory that stores all the comprehensive plots that are on our manuscript.
|
└── post_processing_routines                       # The master directory that contains all the scripts used for data post-processing, scaling analysis, and data collapsing
    ├── ADOS_processing/                           # The directory containing scripts for ADOS post-processing
    |   └── ADOS_extractor.py                      # If running in a manually parallel manner, the script that averages DOS from each disorder realization to form the ADOS
    ├── Plotters/                                  # The directory that contains the plotting scripts for figures in our manuscript
    |   ├── ADOS_plotters.py                       # The script that plots ADOS
    |   ├── Comprehensive_ADOS_scaling_plotter.py  # The script that plots the comprehensive ADOS scaling plot (i.e., Fig. 2 in our manuscript) 
    |   ├── Comprehensive_data_collapsing.py       # The script that plots the comprehensive data collapsing plot (i.e., Fig. 3 in our manuscript) 
    |   └── E_rescaling_and_BCS_scaling.py         # The script that plots the energy rescaling from all dimensions along with BCS scaling for 2D Weyl fermions (i.e., Fig. 1 in our manuscript) 
    └── Scaling_analysis/                          # The directory that contains the scripts for scaling analysis
        ├── Collasing_ADOS_0_vs_L.py               # The script for finite-size data collapsing
        ├── Collapsing_ADOS_vs_E.py                # The script for data collapsing for ADOS against E
        ├── Find_Wc.py                             # (deprecated) The script for finding Wc via fitting
        ├── Fit_ADOS_vs_delta.py                   # The script for scaling of ADOS(0) vs. delta
        ├── Fit_ADOS_vs_E.py                       # The script for scaling of ADOS vs. E
        ├── Fit_ADOS_vs_W_BCS.py                   # The script for BCS-like scaling for 2D Weyl fermions
        └── Rescaling_ADOS_vs_Epow.py              # The script for rescaling of ADOS for W << Wc for Weyl fermions in all dimensions
```

## Filename conventions
#### * For data files
#### * For folders

## Data format


## Numerical methods


## Using the data

## Citations:
[^1]: A. Weiße, G. Wellein, A. Alvermann, and H. Fehske, The kernel polynomial method, Rev. Mod. Phys. 78, 275 (2006)
[^2]: C. W. Groth, M. Wimmer, A. R. Akhmerov, and X. Waintal, Kwant: a software package for quantum transport, New J. Phys. 16, 063065 (2014)
