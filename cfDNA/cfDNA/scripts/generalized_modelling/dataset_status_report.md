# Dataset Processing and Modeling Status Report

## External Data: Cristiano Cohort
- **Manifest samples:** 459 (Cristiano_manifest.csv)
- **Features computed:** 459/459 samples
- **Feature types (19):** coverage, cposedm, edm, eedm, ends, eoedm, ext_poem_prem, fsd, fsr, iedm, ifs, length, ocf, pfe, poem, poem_prem, prem, wps, wps_compute
  - Each type has up to 4 mapq variants: plain, _mapq15, _mapq30, _mapq45
- **Feature matrices built:** Yes — `matrix/` contains per-feature CSVs (mapq45 variants) plus combined `feature_matrix_with_metadata.csv`
- **SVM configs trained (28):**
  - `svm_linear_cv1_mapq30`
  - `svm_linear_mapq30`
  - `svm_linear_mapq45`
  - `svm_linear_pca0.95`
  - `svm_linear_pca0.95_mapq30`
  - `svm_linear_pca100.0_cv1_mapq30`
  - `svm_linear_pca100.0_gc_cv1_mapq30`
  - `svm_linear_pca150.0_cv1_mapq30`
  - `svm_linear_pca150.0_cv1_mapq30_diag`
  - `svm_linear_pca150.0_gc_cv1_mapq30`
  - `svm_linear_pca150.0_gc_cv1_mapq30_platt`
  - `svm_linear_pca150.0_gc_mapq30`
  - `svm_linear_pca150.0_gc_mapq30_old`
  - `svm_linear_pca150.0_gc_mapq30_sig`
  - `svm_linear_pca150.0_mapq30`
  - `svm_linear_pca200.0_cv1_mapq30`
  - `svm_linear_pca25.0_cv1_mapq30`
  - `svm_linear_pca50.0_cv1_mapq30`
  - `svm_linear_pca50.0_gc_cv1_mapq30`
  - `svm_rbf_cv1`
  - `svm_rbf_cv1_mapq15`
  - `svm_rbf_cv1_mapq30`
  - `svm_rbf_cv1_mapq45`
  - `svm_rbf_cv1_v2`
  - `svm_rbf_pca0.95`
  - `svm_rbf_pca100.0_gc_cv1_mapq30`
  - `svm_rbf_pca150.0_gc_cv1_mapq30`
  - `svm_rbf_pca50.0_gc_cv1_mapq30`
- **CV status:** All 28 configs have a `cv/` directory with per-feature `_cv_probs.csv` files

---

## Internal Data: GenoScan (gs) Dataset
- **Metadata samples:** 38 (internal_metadata_gs.csv)
- **Features computed:** 38/38 samples
- **Feature types (19):** coverage, cposedm, edm, eedm, ends, eoedm, ext_poem_prem, fsd, fsr, iedm, ifs, length, ocf, pfe, poem, poem_prem, prem, wps, wps_compute
  - Each type has 2 mapq variants: plain and _mapq30
- **Feature matrices built:** Yes — `internal_matrix_gs/by_feature/` and `matrix_internal_gs/by_feature/`
- **SVM configs trained (1):**
  - `svm_linear_mapq30`
- **CV status:** `svm_linear_mapq30/cv/` contains `svm_linear_mapq30_edm_cv_probs.csv` (1x10 CV, edm only)

---

## Internal Data: PreveLynch (ly) Dataset
- **Metadata samples:** 689 (internal_metadata_ly.csv)
- **Features computed:** 459/689 samples — **230 samples missing features**
- **Feature types (16):** coverage, cposedm, edm, eedm, ends, eoedm, fsd, fsr, iedm, ifs, length, ocf, pfe, poem, prem, wps
  - Missing compared to gs: `ext_poem_prem`, `poem_prem`, `wps_compute`
  - Each type has 2 mapq variants: plain and _mapq30
- **Feature matrices built:** No — no `internal_matrix_ly/` directory
- **SVM configs trained:** None


The list of final features that will be used and are actually important to this work: 
- length, fsd, fsr, edm, iedm, eedm, eoedm, cposedm, pfe, coverage, ends, ocf, ifs, wps

'poem', 'prem', 'poem_prem', 'ext_poem_prem', and 'wps_compute' are old namings or badly designed features and should be ignored
