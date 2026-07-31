# Global Temporal Trends in Carbapenem Resistance Among Bloodstream *Acinetobacter* spp. Isolates Reported to WHO GLASS (2020–2023)

## Overview

This repository contains the complete analytical workflow used for the manuscript:

> **Global Temporal Trends in Carbapenem Resistance Among Bloodstream *Acinetobacter* spp. Isolates Reported to WHO GLASS, 2020–2023**

The project evaluates temporal trends in carbapenem resistance among bloodstream *Acinetobacter* spp. isolates reported to the World Health Organization (WHO) Global Antimicrobial Resistance and Use Surveillance System (GLASS).

The repository provides a fully reproducible workflow beginning with the processed GLASS dataset and reproduces all primary analyses, sensitivity analyses, manuscript tables, and publication-quality figures.

---

## Study Objectives

The study aimed to:

- quantify regional carbapenem resistance among bloodstream *Acinetobacter* spp.
- evaluate temporal resistance trends from 2020–2023
- compare meropenem and imipenem resistance
- assess robustness using a fixed-panel sensitivity analysis restricted to continuously reporting countries

---

## Data Source

Data were obtained from the WHO Global Antimicrobial Resistance and Use Surveillance System (GLASS).

The analyses include bloodstream isolates of *Acinetobacter* spp. tested against:

- Meropenem
- Imipenem

Study years:

- 2020
- 2021
- 2022
- 2023

WHO regions included:

- African Region
- Eastern Mediterranean Region
- European Region
- Region of the Americas
- South-East Asia Region
- Western Pacific Region

---

# Repository Structure

```
acinetobacter-carbapenem-glass-trendsV2/

├── data/
│   └── processed/
│       ├── glass_combined.csv
│       ├── glass_carbapenem_validated.csv
│       ├── carbapenem_pooled_summary.csv
│       ├── logistic_regression_results.csv
│       ├── fixed_panel_regression_results.csv
│       └── fixed_panel_countries.csv
│
├── scripts/
│   ├── 01_load_and_validate_data.R
│   ├── 02_pooled_resistance_wilson_ci.R
│   ├── 03_primary_logistic_regression.R
│   ├── 04_fixed_panel_sensitivity.R
│   ├── 05_figure1_regional_trends.R
│   └── 06_figure2_heatmap.R
│
├── output/
│   ├── figures/
│   └── tables/
│
└── README.md
```

---

# Analytical Workflow

Run the scripts in the following order:

| Script | Purpose |
|---------|----------|
| 01 | Load and validate the processed GLASS dataset |
| 02 | Calculate pooled regional resistance and Wilson confidence intervals (Table 1) |
| 03 | Primary logistic regression analysis (Table 2) |
| 04 | Fixed-panel sensitivity analysis (Supplementary Table S2) |
| 05 | Generate Figure 1 |
| 06 | Generate Figure 2 |

---

# Statistical Methods

The primary analysis uses:

- isolate-weighted pooled resistance estimates
- Wilson 95% confidence intervals
- binomial logistic regression
- odds ratio (OR) for annual temporal change
- Wald confidence intervals
- Benjamini–Hochberg correction for multiple comparisons

Sensitivity analyses were performed using only countries reporting continuously throughout all study years.

---

# Software

Analyses were performed in:

- R (version 4.4 or later)

Primary packages:

- dplyr
- readr
- ggplot2
- tidyr
- writexl
- binom
- viridis

All required packages are installed automatically if absent.

---

# Reproducibility

The repository reproduces:

- Table 1
- Table 2
- Supplementary Table S2
- Figure 1
- Figure 2

using the provided processed dataset.

---

# Results Generated

Running all scripts will generate:

### Tables

- Table 1 – Regional pooled resistance estimates
- Table 2 – Primary logistic regression
- Supplementary Table S2 – Fixed-panel sensitivity analysis

### Figures

- Figure 1 – Regional carbapenem resistance trends
- Figure 2 – Regional heatmap of pooled resistance

---

# Citation

If you use this repository, please cite:

**Waqas M**, *et al.*

*Global Temporal Trends in Carbapenem Resistance Among Bloodstream Acinetobacter spp. Isolates Reported to WHO GLASS, 2020–2023.*

(Manuscript under review.)

---

# License

Code is released under the MIT License.

See the LICENSE file for details.

---

# Contact

**Dr. Muhammad Waqas**

Department of Microbiology

For questions regarding the analyses or repository, please open a GitHub Issue.
