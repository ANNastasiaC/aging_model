## Systems-Level Modelling of DNA Damage, Senescence, and Stem Cell Dynamics in Ageing

This repository contains Python code for simulating the dynamics of stem cells, differentiated cells, and senescent cells over the human lifespan. The model explores the effects of stochastic processes, anti-aging interventions, and cancerous cell dynamics.

## Repository Structure

```text
aging-model/
│
├─ tables/
│   ├─ params_table.csv     # Simulation parameters
│   └─ rates_table.docx     # Transition between cell states and process descriptions
├─ baseline_model/
│   ├─ deterministic/       # Deterministic model
│   │   ├─ functions.py
│   │   ├─ params.py
│   │   └─ basic_figures.py
│   └─ stochastic.py        # Stochastic model
├─ interventions_simulation/
│   ├─ cancer_simulation.py
│   ├─ stem_cell_params.py
│   └─ therapy_simulations.py
```

## Installation

This code requires Python 3.11+ and the following packages:

- numpy
- matplotlib
