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

---

## Thesis Goals and Scope

### Primary contribution
This thesis reimplements the Zhou et al. 2024 cfDNA fragmentomics pipeline entirely from paper descriptions — no source code or repository was available from the authors. The pipeline covers preprocessing (FinaleDB Snakemake), feature extraction, matrix construction, and SVM modelling. On top of this reimplementation, four new end-motif features were designed and evaluated: iEDM, eEDM, eoEDM, and ocEDM (also called cposedm in some scripts).

### Possible thesis goals (to be confirmed)
1. **Reproducibility**: Demonstrate that the Zhou et al. feature pipeline can be faithfully reproduced from paper descriptions alone, achieving comparable AUC scores on the Cristiano cohort.
2. **New feature design**: Evaluate whether the newly designed end-motif features (iedm, eedm, eoedm, ocedm) provide additional discriminative signal beyond EDM.
3. **Ensemble modelling**: Train a stacking ensemble (per-feature SVM base models → meta-learner trained on OOF predictions) and evaluate its performance vs individual feature models.
4. **External validation**: Validate the Cristiano-trained ensemble on the internal GenoScan (gs) and PreveLynch (ly) datasets.
5. **Internal dataset CV**: Perform 10×10 CV and ensemble training on the ly dataset independently; test cross-dataset generalizability (ly → gs).

### Planned modelling workflow
- **Step 1**: 10×10 CV on Cristiano (459 samples) per feature → OOF predictions → meta-matrix → train ensemble → evaluate
- **Step 2**: Apply Cristiano-trained ensemble directly to gs (external validation, no retraining — gs is only 38 samples)
- **Step 3**: 10×10 CV on ly (459 available samples, 230 missing due to FASTQ preprocessing failure) → meta-matrix → train ensemble → test on held-out ly subset (hold out before CV)
- **Step 4**: Apply ly-trained ensemble to gs as cross-dataset test

### Open question
Cancer vs healthy control breakdown for gs and ly datasets not yet confirmed — needed to contextualize validation results and assess class balance.

### Hyperparameter config (fixed for all analyses)
Linear SVM kernel, MAPQ 30 filtering, PCA 150 components, GC correction on. Selected via sequential search on Cristiano with 1×10 CV.
