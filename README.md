# Air pollution, socioeconomic position and mental health: genetic correlation analyses

Code for the analyses in:

> Cava C, Sabetian S, Gaetani M. *Socioeconomic stratification explains the genetic overlap between air pollution and mental health.* (submitted)

The analyses use publicly available GWAS summary statistics from the Neale Lab UK Biobank release (round 2, imputed v3, both sexes). No individual-level data are used.

## Overview

| Step | Script | What it does | Main output |
|---|---|---|---|
| 1 | `R/01_download_qc.R` | Screens the Neale Lab manifest for mental health phenotypes; downloads and quality-controls the GWAS of the classified traits, the air pollution measures and four socioeconomic traits | `sumstats_qc/` |
| 2 | `R/02_format_ldsc.R` | Annotates rsIDs and writes LDSC input files (A1 = alternate allele, A2 = reference allele, Z = β/SE) | `sumstats_ldsc/` |
| 3 | `R/03_ldsc_initial.R` | Joint LDSC on all candidate traits (h², Z, intercept, attenuation ratio) and genetic correlations among air pollution measures, used for trait selection | Tables S1–S2 |
| 4 | `R/04_main_analyses.R` | LDSC on the selected traits; genetic correlations; socioeconomic adjustment (partial genetic correlations); subdomain means; LCV | Tables S3–S7 |
| 5 | `python/figures.py` | Figures 2–4 and Supplementary Figure S1 | `figures/` |

`data/trait_classification.csv` records the manual classification of every candidate phenotype: domain, subdomain, inclusion in the main analyses (Z(h²) ≥ 4, one version per phenotype, one measure per pollutant class), eligibility for LCV (Z(h²) ≥ 7), inclusion in the sensitivity set (help-seeking and treatment phenotypes) and reason for exclusion.

## Requirements

- R ≥ 4.3 with `data.table`, `httr`, `MASS` and [`GenomicSEM`](https://github.com/GenomicSEM/GenomicSEM) (analyses run with R 4.5.1 and GenomicSEM 0.0.5)
- Python ≥ 3.9 with `numpy` and `matplotlib` (figures only)

## External files

1. **LD scores**: European LD scores and weights (`eur_w_ld_chr`, HapMap3 SNPs) distributed with LDSC/GenomicSEM.
2. **Variant annotation**: Neale Lab `variants.tsv.bgz` (variant → rsID), placed in the same folder as the LD scores.
3. **LCV**: `MomentFunctions.R` and `RunLCV.R` from the original implementation, https://github.com/lukejoconnor/LCV.

Set the three paths at the top of `config.R`.

## Running

From the repository folder:

```r
source("R/01_download_qc.R")     # several hours, ~150 GWAS files
source("R/02_format_ldsc.R")
source("R/03_ldsc_initial.R")
source("R/04_main_analyses.R")   # LDSC 1-2 h; partial correlations and LCV ~1 h
```

```bash
python python/figures.py <BASE_DIR>/results/main <BASE_DIR>/figures
```

Steps save intermediate results and can be restarted. Results of the partial genetic correlations are based on 2,000 draws from the sampling distribution of the LDSC estimates (`set.seed(2026)`).

## Output tables (`results/main/`)

| File | Content |
|---|---|
| `T1_trait_h2_intercept_ratio.csv` | SNP heritability, LDSC intercept and attenuation ratio per trait |
| `T2_rg_pollution_mental_health.csv` | Genetic correlations, air pollution × mental health (Table S3) |
| `T3_rg_socioeconomic.csv` | Genetic correlations with socioeconomic traits (Table S5) |
| `T4_rg_sensitivity_help_seeking.csv` | Sensitivity set (Table S4) |
| `T5_partial_rg_socioeconomic.csv` | Partial genetic correlations conditional on socioeconomic traits (Table S6) |
| `T6_rg_by_subdomain.csv` | Mean genetic correlation by mental health subdomain |
| `T7_lcv.csv` | Latent causal variable analysis (Table S7) |

## License

MIT (see `LICENSE`).
