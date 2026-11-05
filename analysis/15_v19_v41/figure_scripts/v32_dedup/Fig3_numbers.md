# Fig. 3 - numbers printed on the figure and their sources

| item | value | source |
|---|---|---|
| TCGA-STAD: n/events | 383/156 | figures/data/km_TCGA.csv |
| TCGA-STAD: HR per SD | 0.82 (0.70–0.95) | E8_univariate.csv beta, se (NetImmune) |
| TCGA-STAD: log-rank P | 0.013 | recomputed from figures/data/km_TCGA.csv (median split) |
| GSE84437: n/events | 431/207 | figures/data/km_GSE84437.csv |
| GSE84437: HR per SD | 0.75 (0.65–0.86) | E8_univariate.csv beta, se (NetImmune) |
| GSE84437: log-rank P | 2.9 × 10⁻⁴ | recomputed from figures/data/km_GSE84437.csv (median split) |
| GSE62254 (ACRG): n/events | 283/135 | figures/data/km_GSE62254.csv |
| GSE62254 (ACRG): HR per SD | 0.66 (0.56–0.78) | E8_univariate.csv beta, se (NetImmune) |
| GSE62254 (ACRG): log-rank P | 2.4 × 10⁻⁴ | recomputed from figures/data/km_GSE62254.csv (median split) |
| d NetImmune TCGA-STAD | 0.82 (0.70–0.95) | E8_univariate / E15_per_cohort / E20_per_cohort |
| d NetImmune GSE84437 | 0.75 (0.65–0.86) | E8_univariate / E15_per_cohort / E20_per_cohort |
| d NetImmune GSE62254 | 0.66 (0.56–0.78) | E8_univariate / E15_per_cohort / E20_per_cohort |
| d NetImmune Pooled, benchmark | 0.74 (0.58–0.95), P = 0.036 | E20_pooled.csv (unadjusted, primary) |
| d NetImmune GSE15459 | 0.65 (0.52–0.81) | E8_univariate / E15_per_cohort / E20_per_cohort |
| d NetImmune GSE34942 | 0.92 (0.62–1.36) | E8_univariate / E15_per_cohort / E20_per_cohort |
| d NetImmune GSE26253 (RFS) | 0.73 (0.63–0.86) | E8_univariate / E15_per_cohort / E20_per_cohort |
| d NetImmune Pooled, first three | 0.72 (0.53–1.00), P = 0.049 | E20_pooled.csv (unadjusted, primary) |
| d NetImmune GSE13861 | 0.50 (0.33–0.77) | E8_univariate / E15_per_cohort / E20_per_cohort |
| d NetImmune GSE26899 | 0.58 (0.39–0.85) | E8_univariate / E15_per_cohort / E20_per_cohort |
| d NetImmune GSE26901 | 0.71 (0.55–0.93) | E8_univariate / E15_per_cohort / E20_per_cohort |
| d NetImmune GSE29272 | 0.72 (0.59–0.89) | E8_univariate / E15_per_cohort / E20_per_cohort |
| d NetImmune Pooled, locked 7 | 0.70 (0.62–0.79), P = 2.9 × 10⁻⁴ | E20_pooled.csv (unadjusted, primary) |
| d NetImmune GSE57303 | 0.84 (0.59–1.19) | E8_univariate / E15_per_cohort / E20_per_cohort |
| d NetImmune Pooled, all 10 gastric | 0.72 (0.67–0.78), P = 7.3 × 10⁻⁶ | E20_pooled.csv (unadjusted, primary) |
| d CAF_ECM TCGA-STAD | 1.28 (1.09–1.50) | E8_univariate / E15_per_cohort / E20_per_cohort |
| d CAF_ECM GSE84437 | 1.28 (1.11–1.48) | E8_univariate / E15_per_cohort / E20_per_cohort |
| d CAF_ECM GSE62254 | 1.31 (1.11–1.56) | E8_univariate / E15_per_cohort / E20_per_cohort |
| d CAF_ECM Pooled, benchmark | 1.29 (1.18–1.41), P = 0.0011 | E20_pooled.csv (unadjusted, primary) |
| d CAF_ECM GSE15459 | 1.54 (1.20–1.97) | E8_univariate / E15_per_cohort / E20_per_cohort |
| d CAF_ECM GSE34942 | 0.95 (0.65–1.39) | E8_univariate / E15_per_cohort / E20_per_cohort |
| d CAF_ECM GSE26253 (RFS) | 1.26 (1.05–1.50) | E8_univariate / E15_per_cohort / E20_per_cohort |
| d CAF_ECM Pooled, first three | 1.27 (0.75–2.15), P = 0.19 | E20_pooled.csv (unadjusted, primary) |
| d CAF_ECM GSE13861 | 2.00 (1.26–3.18) | E8_univariate / E15_per_cohort / E20_per_cohort |
| d CAF_ECM GSE26899 | 1.41 (0.97–2.05) | E8_univariate / E15_per_cohort / E20_per_cohort |
| d CAF_ECM GSE26901 | 1.65 (1.22–2.23) | E8_univariate / E15_per_cohort / E20_per_cohort |
| d CAF_ECM GSE29272 | 1.09 (0.89–1.33) | E8_univariate / E15_per_cohort / E20_per_cohort |
| d CAF_ECM Pooled, locked 7 | 1.34 (1.09–1.65), P = 0.014 | E20_pooled.csv (unadjusted, primary) |
| d CAF_ECM GSE57303 | 0.90 (0.64–1.28) | E8_univariate / E15_per_cohort / E20_per_cohort |
| d CAF_ECM Pooled, all 10 gastric | 1.31 (1.17–1.46), P = 3.4 × 10⁻⁴ | E20_pooled.csv (unadjusted, primary) |
| e Benchmark, pooled / 1,000 axes from all genes | rank 1 of 1,001; P ≤ 0.001; observed |log HR| 0.295 | E2_null_joint / E8_biology_null / E20_null |
| e Benchmark, pooled / 1,000 axes from immune and stromal genes | rank 6 of 1,001; P = 0.006; observed |log HR| 0.295 | E2_null_joint / E8_biology_null / E20_null |
| e Locked 7, pooled / 1,000 structure-matched axes | rank 1 of 1,001; P ≤ 0.001; observed |log HR| 0.355 | E2_null_joint / E8_biology_null / E20_null |
| f bench vs CD8A | +0.068 (+0.041 to +0.095), P 6.3 × 10⁻⁷ | E1_paired_dC / E8_h2h_TMEscore / E20_deltaC |
| f bench vs Immune only | +0.046 (+0.023 to +0.069), P 8.9 × 10⁻⁵ | E1_paired_dC / E8_h2h_TMEscore / E20_deltaC |
| f bench vs Barrier | +0.038 (+0.010 to +0.066), P 0.0071 | E1_paired_dC / E8_h2h_TMEscore / E20_deltaC |
| f bench vs CAF_ECM | +0.021 (-0.003 to +0.044), P 0.09 | E1_paired_dC / E8_h2h_TMEscore / E20_deltaC |
| f bench vs ssGSEA scoring | +0.008 (-0.002 to +0.018), P 0.11 | E1_paired_dC / E8_h2h_TMEscore / E20_deltaC |
| f bench vs ESTIMATE | +0.000 (-0.014 to +0.015), P 0.96 | E1_paired_dC / E8_h2h_TMEscore / E20_deltaC |
| f bench vs TMEscore | -0.003 (-0.025 to +0.020), P 0.83 | E1_paired_dC / E8_h2h_TMEscore / E20_deltaC |
| f lock vs CD8A | +0.063 (+0.034 to +0.093), P 2.6 × 10⁻⁵ | E1_paired_dC / E8_h2h_TMEscore / E20_deltaC |
| f lock vs ESTIMATE | +0.015 (-0.002 to +0.032), P 0.082 | E1_paired_dC / E8_h2h_TMEscore / E20_deltaC |
| g median rho at 30% removed | 0.979 | E19_spearman_summary.csv ALL (cohort x draw), drop30 |
| g 5th percentile rho at 30% removed | 0.938 | E19_spearman_summary.csv ALL, drop30 |
| g pooled HR drop10 | median 0.743 (0.733–0.756); direction kept 100% | E19_summary.csv |
| g pooled HR drop20 | median 0.745 (0.730–0.764); direction kept 100% | E19_summary.csv |
| g pooled HR drop30 | median 0.747 (0.727–0.772); direction kept 100% | E19_summary.csv |
| h NetImmune mean held-out C | 0.589 | E6_model_classes.csv NetImmune (fixed), mean of 3 cohorts |
| h Cox on 20 PCs mean C at full size | 0.610 (n = 731) | E6_learning_pcacox.csv |
| h Elastic net mean C at full size | 0.597 (n = 731) | E6_learning_en.csv |
| j TCGA-STAD CAF_ECM | 1.28 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j TCGA-STAD TAM_SUP | 1.12 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j TCGA-STAD APM | 0.97 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j TCGA-STAD IFN | 0.94 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j TCGA-STAD TLS_B | 1.01 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j TCGA-STAD T_CELL | 1.00 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j TCGA-STAD Immune | 0.98 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j TCGA-STAD NetImmune | 0.82 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE84437 CAF_ECM | 1.28 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE84437 TAM_SUP | 0.94 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE84437 APM | 0.88 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE84437 IFN | 0.78 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE84437 TLS_B | 0.90 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE84437 T_CELL | 0.81 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE84437 Immune | 0.81 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE84437 NetImmune | 0.75 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE62254 CAF_ECM | 1.31 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE62254 TAM_SUP | 1.02 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE62254 APM | 0.70 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE62254 IFN | 0.74 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE62254 TLS_B | 0.87 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE62254 T_CELL | 0.89 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE62254 Immune | 0.75 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE62254 NetImmune | 0.66 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE15459 CAF_ECM | 1.54 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE15459 TAM_SUP | 1.58 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE15459 APM | 0.88 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE15459 IFN | 1.08 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE15459 TLS_B | 0.83 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE15459 T_CELL | 0.96 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE15459 Immune | 0.93 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE15459 NetImmune | 0.65 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE34942 CAF_ECM | 0.95 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE34942 TAM_SUP | 0.93 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE34942 APM | 0.95 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE34942 IFN | 0.84 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE34942 TLS_B | 0.96 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE34942 T_CELL | 0.80 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE34942 Immune | 0.88 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE34942 NetImmune | 0.92 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE26253 (RFS) CAF_ECM | 1.26 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE26253 (RFS) TAM_SUP | 1.04 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE26253 (RFS) APM | 0.90 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE26253 (RFS) IFN | 0.91 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE26253 (RFS) TLS_B | 0.86 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE26253 (RFS) T_CELL | 0.82 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE26253 (RFS) Immune | 0.82 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE26253 (RFS) NetImmune | 0.73 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE13861 CAF_ECM | 2.00 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE13861 TAM_SUP | 1.19 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE13861 APM | 0.89 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE13861 IFN | 0.70 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE13861 TLS_B | 0.78 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE13861 T_CELL | 0.79 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE13861 Immune | 0.74 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE13861 NetImmune | 0.50 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE26899 CAF_ECM | 1.41 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE26899 TAM_SUP | 1.34 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE26899 APM | 0.77 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE26899 IFN | 0.68 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE26899 TLS_B | 0.85 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE26899 T_CELL | 0.77 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE26899 Immune | 0.74 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE26899 NetImmune | 0.58 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE26901 CAF_ECM | 1.65 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE26901 TAM_SUP | 1.46 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE26901 APM | 0.96 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE26901 IFN | 1.01 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE26901 TLS_B | 0.93 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE26901 T_CELL | 0.97 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE26901 Immune | 0.96 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE26901 NetImmune | 0.71 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE29272 CAF_ECM | 1.09 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE29272 TAM_SUP | 0.90 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE29272 APM | 0.79 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE29272 IFN | 0.75 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE29272 TLS_B | 0.81 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE29272 T_CELL | 0.73 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE29272 Immune | 0.75 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j GSE29272 NetImmune | 0.72 | E20_estimates_used_for_pooling.csv (unadjusted, primary) |
| j pooled CAF_ECM | 1.31 (1.17–1.46), P 3.4 × 10⁻⁴ | E20_pooled.csv all 10 gastric |
| j pooled TAM_SUP | 1.11 (0.97–1.28), P 0.12 | E20_pooled.csv all 10 gastric |
| j pooled APM | 0.86 (0.79–0.93), P 0.0018 | E20_pooled.csv all 10 gastric |
| j pooled IFN | 0.86 (0.77–0.95), P 0.01 | E20_pooled.csv all 10 gastric |
| j pooled TLS_B | 0.89 (0.83–0.95), P 6.6 × 10⁻⁴ | E20_pooled.csv all 10 gastric |
| j pooled T_CELL | 0.86 (0.80–0.93), P 0.0018 | E20_pooled.csv all 10 gastric |
| j pooled Immune | 0.84 (0.78–0.91), P 5.3 × 10⁻⁴ | E20_pooled.csv all 10 gastric |
| j pooled NetImmune | 0.72 (0.67–0.78), P 7.3 × 10⁻⁶ | E20_pooled.csv all 10 gastric |
| i TCGA-STAD | C 0.662 -> 0.677; LR chi2 8.5, P 0.0035 (n = 314) | revision/results/B4_cindex.csv; E8_lrt.csv |
| i GSE84437 | C 0.675 -> 0.683; LR chi2 9.3, P 0.0023 (n = 431) | revision/results/B4_cindex.csv; E8_lrt.csv |
| i GSE62254 | C 0.729 -> 0.743; LR chi2 19.4, P 1.0 × 10⁻⁵ (n = 281) | revision/results/B4_cindex.csv; E8_lrt.csv |
