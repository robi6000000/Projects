---
name: cfDNA project overview
description: Datasets, pipeline stages, directory layout, and current modelling goals for the liquid biopsy cfDNA project
type: project
originSessionId: 3823c645-9b3e-436d-8129-ba9c964f2fd8
---
## Working directory
`/data/projects/liquid_biopsy/Projects/cfDNA/cfDNA/`

## Datasets
- **Cristiano** — 459 samples (EE-prefixed IDs), fragment files in `data/cristiano_features/`, matrices in `data/matrix/by_feature/`, metadata `data/manifest/Cristiano_metadata.csv`
- **PreveLynch (LY)** — Lynch syndrome prevention cohort, features in `data/internal_features_ly/`, matrices in `data/matrix_internal_ly/` (full), `data/matrix_internal_ly_train/` (80%), `data/matrix_internal_ly_test/` (20%), metadata `data/manifest/internal_metadata_ly_filtered*.csv`
- **GenoScan (GS)** — used for external testing only (name from "gene", not "geo")

## Pipeline stages
1. **Feature computation** — `compute_features_mapq.sbatch` (SLURM array), calls `sample_prep_script.sh` → `sample_features_script.py`. `rerun=false` skips already-computed features.
2. **Matrix building** — `matrix_build_array.sbatch` (array) or `matrix_build_single.sbatch` (single feature). Output: `matrix_{feature}_mapq{N}.csv` in `by_feature/`.
3. **SVM CV** — `gen_svm.sbatch` (array 0-13) or `gen_svm_feature_single.sbatch`. Feature index: 0-6 light (8G), 7-12 medium (24G), 13=fsr heavy (64G).
4. **Evaluate** — `python scripts/modelling/evaluate_models.py data/matrix/.../cv`

## Feature index (gen_svm.sbatch)
0=length, 1=edm, 2=iedm, 3=eedm, 4=eoedm, 5=cposedm, 6=fsd, 7=pfe, 8=coverage, 9=ends, 10=ocf, 11=ifs, 12=wps, 13=fsr

## Per-region vs non-per-region features
- **Per-region** (over ~560k open chromatin regions, PCA 150 → SVM): fsr, coverage, ends, pfe, wps, ifs, ocf
- **Non-per-region**: length, fsd, edm, iedm, eedm, eoedm, cposedm

## Best config
`svm_linear_pca150.0_gc_mapq30` — linear SVM, PCA=150, GC correction, MAPQ 30, 10×10 CV

## Naming
- cposedm renamed to **ocedm** (orientation-aware composite end motif) in the thesis. Use `RENAME = {'cposedm': 'ocedm'}` in plots.
- GenoScan = GS (from "gene", not "geo")

## Current status (as of 2026-05-10)
- Cristiano 10×10 CV: **COMPLETE** (all 14 features in `data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30/cv/`)
- PreveLynch CV (cv_ly_med_train, 1717429): med features running ~3d7h at 24G, ~21.8G peak — within normal range, wait until morning
- CRC CV (cv_crc_med, 1864): resubmitted at 32G after memory pressure concerns
- Cristiano per-feature training: light+med done; heavy (fsr) running (1613)
- CRC per-feature training: med (1744) running 6-10h

## Ensemble next steps
Now that Cristiano 10×10 CV is complete, next: build meta-matrix → ensemble CV → ensemble train

## Medium feature memory usage
Peak RSS at 24G: Cristiano ~21.5G (completed OK), LY ~21.8G (running, within range).
pfe at 48G (_safe job) peaked at ~22.1G — 24G is sufficient but tight.
Use 32G for CRC CV med features to avoid pressure.

## Cristiano per-feature training (submitted 2026-05-07)
Config: linear kernel, PCA=150, GC correction, mapq30, cv_repeats=1
Output: `data/matrix/svm_by_feature/svm_linear_pca150.0_gc_cv1_mapq30/models/`
Next: move model files into `svm_linear_pca150.0_gc_mapq30/models/` once all done.

## PreveLynch per-feature training (submitted 2026-05-08)
Config: linear kernel, PCA=150, GC correction, mapq30, cv_repeats=1
Output: `data/matrix_internal_ly_train/svm_by_feature/svm_linear_pca150.0_gc_cv1_mapq30/models/`
Next: same as Cristiano.

## Stage stratification
- Cristiano: stage I/II/III/IV available and used in results.ipynb
- LY: stage column entirely empty — stage AUC not feasible
- CRC: stage AUC cells added to results.ipynb (cells ad109548, c6d614f6), guarded by os.path.exists

## Thesis visualisation strategy (deadline 2026-05-18, 7 days)
- **Cristiano-trained models**: show all available features (14 for CV, as many as complete for external test). Primary comparison to Zhou et al.
- **CRC/LY-trained models**: EDM-feature focus as fallback if med jobs don't finish — restrict to edm, iedm, eedm, eoedm, ocedm (all light features, complete everywhere).

## Visualisation TODO (results.ipynb)
Done: 01-09 hyperparameter/CV plots, 20-23 ensemble plots, 06b CRC stage AUC, feature correlation heatmap.
Missing/needed:
- Score distributions — violin/box plots of predicted probabilities, cancer vs healthy, per feature (Cristiano CV + ensemble)
- Cross-dataset AUC heatmap — feature × condition summary matrix (high priority, can use EDM features as fallback)
- External validation ROC + AUC bars for all training configs (cells exist, just need jobs to finish)

## Ensemble (Cristiano)
- CV AUC: 0.9932 (10×10 CV)
- CV probs: data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30/ensemble/cv/cv_probs.csv
- Model: data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30/ensemble/models/ensemble.pkl
- Meta-matrix: data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30/meta_matrix/meta_matrix.csv
