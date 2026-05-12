# cfDNA pipeline

Working directory: `/data/projects/liquid_biopsy/Projects/cfDNA/cfDNA/`

## Directory structure

```
cfDNA/
├── data/
│   ├── manifest/                           # metadata files
│   │   ├── Cristiano_metadata.csv
│   │   ├── Cristiano_metadata_crc.csv
│   │   ├── internal_metadata_ly_filtered_train.csv
│   │   ├── internal_metadata_ly_filtered_test.csv
│   │   └── internal_metadata_gs.csv
│   ├── processing/                         # reference files built once by static_script.sh
│   │   ├── openchrom_with_id.bed           # OC regions with IDs
│   │   ├── openchrom_with_id_centroid.bed  # OC regions with IDs, with centroid coordinates in place of fragment coordinates
│   │   ├── wps_openchrom_with_id.bed       # OC regions extended ±60bp for WPS 
│   │   └── gc_content_per_region.csv       # GC content % per OC region
│   ├── source/                             # downloaded reference files
│   │   ├── hg38ToHg19.over.chain.gz
│   │   ├── E029/E032/E034-DNase.hotspot.fdr0.01.broad.bed.gz # B, T, Monocyte DNAse peaks
│   │   └── TCGA-ATAC_PanCancer_PeakSet.txt                   # TCGA ATAC-seq multi-cancer peaks
│   ├── ref_genome/hg19/                    # reference genome FASTA
│   │   └── hg19.fa
│   ├── sample_temp/                        # per-sample intermediate BED files, auto-deleted
│   ├── cristiano_features/
│   │   └── {feature}/{feature}_mapq30/
│   ├── internal_features_ly/
│   ├── internal_features_gs/
│   ├── matrix/                             # final matrices, model runs, Cristiano cohort
│   │   ├── by_feature/                     # matrix building output
│   │   │   └── matrix_{feature}_{mapq}.csv
│   │   └── svm_by_feature/                 # specific SVM config runs, named after config parameters of the base models
│   │       └── svm_{kernel}_{pca}_{gc}_{cv}_{mapq}_{run_tag}/  
│   │           ├── cv/                     # per-feature CV probabilities
│   │           ├── models/                 # per-feature trained model .pkl files 
│   │           ├── test/                   # per-feature test predictions
│   │           │   └── {test_dataset}/
│   │           ├── meta_matrix/            # meta-matrix - per-feature CV probabilities merged
│   │           └── ensemble_{run_tag}/     # ensemble model run folder
│   │               ├── cv/                 # ensemble CV probabilities
│   │               ├── models/             # ensemble model .pkl file
│   │               └── test/               # ensemble test predictions
│   ├── matrix_internal_ly_train/
│   ├── matrix_internal_ly_test/
│   ├── matrix_internal_gs/
│   └── matrix_crc/
└── scripts/
    ├── generalized_processing/             # sample processing and matrix building scripts
    │   ├── static_script.sh                # one-time prerequisite file builder
    │   ├── compute_features_mapq.sbatch    # SLURM array job, calls sample_prep_script.sh and sample_features_script.py for each sample
    │   ├── sample_prep_script.sh
    │   ├── sample_features_script.py
    │   ├── matrix_build_array.sbatch       # SLURM array job, calls matrix_build.py for each feature
    │   └── matrix_build.py
    └── generalized_modelling/              # modelling scripts
        ├── gen_svm.sbatch                  # SLURM array job, calls gen_svm_feature.py for CV or training of per-feature models
        ├── gen_svm_feature.py
        ├── gen_svm_test.sbatch             # SLURM array job, calls gen_svm_test.py for testing per-feature models on external cohorts
        ├── gen_svm_test.py
        ├── gen_ensemble_svm.sbatch         # SLURM job, calls gen_ensemble_svm.py for CV, training or testing of the ensemble model
        └── gen_ensemble_svm.py
```

## Step 0 — Snakemake preprocessing (optional)

Only needed if input data is in paired FASTQ format rather than fragment BED files. Located at `scripts/snakemake_finaledb/`. Aligns reads, marks duplicates, and produces gzip-compressed fragment BED files with columns: chrom, start, end, MAPQ, strand. Skip this step if the data is already in that format.

Implementation of the FinaleDB Snakemake Workflow is out of scope for this thesis, but 

## Step 1 — Build static reference files

`scripts/generalized_processing/static_script.sh` — run once per project, no arguments.

```bash
bash scripts/generalized_processing/static_script.sh
```

Requires wget, bedtools, liftOver (downloaded automatically), and python with pandas.

Downloads:
- hg19 reference FASTA (`data/ref_genome/hg19/hg19.fa`)
- UCSC liftover chain hg38→hg19
- Roadmap Epigenomics DNase peaks (cell lines E029, E032, E034)
- TCGA ATAC-seq pan-cancer peak set
Produces in `data/processing/`:
- `openchrom_with_id.bed` — ~560k merged open chromatin regions (DNase + ATAC), 200bp windows centred on peak centroids, with unique integer region IDs
- `openchrom_with_id_centroid.bed` — same with an added centroid coordinate column
- `wps_openchrom_with_id.bed` — regions extended ±60bp on each side for WPS computation
- `gc_content_per_region.csv` — GC fraction per open chromatin region

## Step 2 — Compute per-sample features

### a. SLURM array job for feature vector computation

`scripts/generalized_processing/compute_features_mapq.sbatch` — submits one SLURM task per sample. Each task reads the sample's fragment file path from the metadata CSV by array index and calls `sample_prep_script.sh`.

```bash
sbatch --array=1-459%20 \
    scripts/generalized_processing/compute_features_mapq.sbatch \
    30 \                              # MAPQ_FILTER (int or "none")
    false \                           # RERUN ("true" recomputes existing files)
    data/manifest/Cristiano_metadata.csv \
    data/cristiano_features
```

Arguments (positional):
- `MAPQ_FILTER` — minimum MAPQ score; reads below threshold are excluded (default: none)
- `RERUN` — if false, skips samples whose output files already exist (default: true)
- `METADATA` — sample manifest CSV (default: `data/manifest/internal_metadata.csv`)
- `FEATURES_PATH` — output root directory for feature CSVs (default: `data/internal_features`)

Required metadata columns:
- `sample_id` — unique sample identifier
- `disease` — disease label (e.g. Healthy, Breast cancer)
- `dataset` — cohort name
- `material` — sample material (e.g. pl for plasma)
- `frag_path` — path to gzipped fragment BED file
- `stage` — cancer stage (can be empty)
- `cancer_true` — binary label: 1=cancer, 0=healthy

Fragment BED files are gzip-compressed, BED5 or BED6 format. Chromosomes can be bare (`1`) or prefixed (`chr1`); the script normalises both.

### b. Per-sample preprocessing

`scripts/generalized_processing/sample_prep_script.sh` — called automatically by the array job. Filters fragments to autosomes and creates genomic intersection files against the open chromatin reference. Arguments: `<frag_file> <sample_id> <mapq_filter> <rerun>`.

Intermediate files are written to `data/sample_temp/` and deleted after feature computation:
- `{id}_autosomes_chr_id.bed` — fragments filtered to chr1–chr22, MAPQ ≥ threshold
- `{id}_frag_centroids_openchrom_intersect.bed` — fragment midpoints overlapping open chromatin regions
- `{id}_frag_ends_U_D.bed` — fragment start/end positions with upstream/downstream directionality
- `{id}_frag_ends_openchrom_intersect.bed` — fragment ends overlapping open chromatin regions
- `{id}_frag_wps_intersect.bed` — fragments overlapping WPS-padded regions
- `{id}_frag_ends_with_region.bed` — fragment ends annotated with region IDs
- `{id}_frag_ends_ocf.bed` — fragment ends with OCF position annotation

### c. Feature computation

`scripts/generalized_processing/sample_features_script.py` — computes all 14 features for one sample from the intersection BED files. Arguments: `<sample_id> <centroids_bed> <ends_bed> <ocf_bed> <wps_bed> [mapq_filter] [rerun] [features_path]`.

Output: one CSV per feature in `{features_path}/{feature}_mapq{N}/{sample_id}_{feature}.csv`.

Features (name — shape per sample — description):
- `length` — 22 chroms × 31 bins — fragment length distribution per chromosome (10bp bins, 0–≥300bp)
- `fsd` — 22 chroms × 68 bins — fragment size distribution, 5bp bins (65–405bp)
- `fsr` — ~560k regions × 3 — short/medium/long fragment size ratio per open chromatin region
- `pfe` — ~560k regions × 1 — per-region fragmentation entropy (Shannon entropy of length bins)
- `coverage` — ~560k regions × 1 — fragment count per open chromatin region
- `ends` — ~560k regions × 1 — fragment end count per open chromatin region
- `ocf` — ~560k regions × 1 — open chromatin fragmentation score (end asymmetry relative to centroid)
- `ifs` — ~560k regions × 1 — intra-fragment score (normalised fragment length) per region
- `wps` — ~560k regions × 1 — windowed protection score (120bp window, 120–180bp fragments only)
- `edm` — 22 chroms × 256 — 4-mer end dinucleotide motif frequencies at fragment termini
- `iedm` — 22 chroms × 256 — internal-offset end motifs (positions 4–8 from fragment ends)
- `eedm` — 22 chroms × 256 — external flanking motifs (positions −4 and beyond)
- `eoedm` — 22 chroms × 256 — external-offset flanking motifs (positions −8, +4)
- `cposedm` (→ `ocedm` in thesis) — 22 chroms × 2048 — position-aware composite end motif matrix (8 positions × 256 motifs)

All per-region features are merged against the full open chromatin region list; missing regions are zero-filled. Requires `data/processing/openchrom_with_id.bed` and `data/ref_genome/hg19/hg19.fa`.

## Step 3 — Build feature matrices

### a. SLURM array job for matrix building

`scripts/generalized_processing/matrix_build_array.sbatch` — submits 14 tasks, one per feature.

```bash
sbatch scripts/generalized_processing/matrix_build_array.sbatch \
    30 \                                  # MAPQ_FILTER
    data/cristiano_features \             # FEATURES_FOLDER
    data/matrix/by_feature \             # MATRIX_FOLDER
    data/manifest/Cristiano_metadata.csv  # METADATA_PATH
```

Feature order (array indices 0–13): `length edm iedm eedm eoedm cposedm fsd pfe coverage ends ocf ifs wps fsr`

For a single feature, use `matrix_build_single.sbatch` instead.

### b. Matrix builder

`scripts/generalized_processing/matrix_build.py` — called once per feature by the array job. Aggregates all per-sample feature CSVs into one matrix. Arguments: `<feature> [mapq_filter] [features_folder] [matrix_folder] [metadata_path]`.

Output: `{matrix_folder}/matrix_{feature}_mapq{N}.csv` with leading metadata columns (sample_id, disease, dataset, material, stage, cancer_true) followed by feature columns. Per-region features use `{feature}_region_{i}` naming; matrix features use `{feature}_{row}_{col}`. Only samples present in the metadata CSV are included. Memory: 32G per task; per-region features (~560k columns) are the most intensive.

## Modelling

All modelling scripts are in `scripts/generalized_modelling/`. The best-performing configuration throughout is linear kernel, PCA=150, GC correction, MAPQ 30, 10×10 CV — config name `svm_linear_pca150.0_gc_mapq30`. The config name is auto-derived by `gen_svm_feature.py` from its arguments as `svm_{kernel}[_pca{N}][_gc][_cv{R}][_mapq{M}][_{run_tag}]` (the `_cv{R}` part is omitted when cv_repeats=10). All output lives under `{matrix_folder}/svm_by_feature/{config_name}/`.

### Step 4 — Per-feature cross-validation (10×10)

`gen_svm.sbatch` calls `gen_svm_feature.py`. Runs repeated stratified k-fold CV — 10 folds × 10 repeats = 100 splits — producing out-of-fold predicted probabilities for every sample.

```bash
sbatch --array=0-6  --mem=8G  gen_svm.sbatch cv data/matrix linear true true 150 10 30
sbatch --array=7-12 --mem=32G gen_svm.sbatch cv data/matrix linear true true 150 10 30
sbatch --array=13   --mem=64G gen_svm.sbatch cv data/matrix linear true true 150 10 30
```

Arguments (positional): option (cv/train), matrix_folder, kernel, gc_correction, pca, pca_components, cv_repeats, mapq_filter, run_tag (optional).

GC correction is only applied to per-region features (pfe, coverage, ends, ocf, ifs, wps). Passing `true` for any other feature is silently ignored.

Memory by group:
- light (indices 0–6): length, edm, iedm, eedm, eoedm, cposedm, fsd — 8G
- med (indices 7–12): pfe, coverage, ends, ocf, ifs, wps — 32G
- heavy (index 13): fsr — 64G

Output: `{matrix_folder}/svm_by_feature/{config_name}/cv/{feature}_cv_probs.csv`

### Step 5 — Per-feature training (full dataset)

Same script, `train` as first argument. Fits on the entire matrix with no CV splits.

```bash
sbatch --array=0-6  --mem=8G  gen_svm.sbatch train data/matrix linear true true 150 10 30
sbatch --array=7-12 --mem=32G gen_svm.sbatch train data/matrix linear true true 150 10 30
sbatch --array=13   --mem=64G gen_svm.sbatch train data/matrix linear true true 150 10 30
```

Output: `{matrix_folder}/svm_by_feature/{config_name}/models/{config_name}_{feature}.pkl`

The pickle stores the full preprocessing + model pipeline: feature name, kernel, gc_correction flag, pca flag, pca_components, feature_columns list (used to align test matrices), fitted PCA object, fitted GC scaler, and fitted SVC.

### Step 6 — Build meta-matrix (ensemble prerequisite)

`scripts/generalized_processing/meta_matrix_build.py` — collects all per-feature CV probability files and joins them into a single N_samples × (metadata + F features) matrix. This is the input to the ensemble model.

```bash
python scripts/generalized_processing/meta_matrix_build.py \
    data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30/cv/ \
    data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30/meta_matrix/ \
    data/manifest/Cristiano_metadata.csv
```

The script globs all `*.csv` files in the probs folder and extracts the feature name as `filename.split("_")[-3]`, which expects standard naming (`{feature}_cv_probs.csv`). Non-standard filenames (e.g. old `svm_linear_..._ifs_cv_probs.csv`) will produce wrong column names — rename them before running. Only samples present in both the CV probs and the metadata are included (inner join).

Output: `{output_folder}/meta_matrix.csv`

### Step 7 — Ensemble cross-validation

`gen_ensemble_svm.sbatch` calls `gen_ensemble_svm.py` with option `cv`. Runs 10×10 CV on the meta-matrix with no PCA and no GC correction — the ensemble operates directly on per-feature probabilities.

```bash
sbatch --job-name=ensemble_cv \
    gen_ensemble_svm.sbatch \
    cv \
    data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30/meta_matrix/meta_matrix.csv \
    linear 10 none none none nostd false
```

Arguments (positional): option, meta_matrix_path, kernel, cv_repeats, dataset_tag (none for cv/train), model_path (none = default), test_output_path (none = default), run_tag (sets output dir suffix, e.g. `nostd` → `ensemble_nostd/`), standardize (true/false).

Output directory is derived from the parent of meta_matrix_path: `{config_dir}/ensemble_{run_tag}/`. Output: `{config_dir}/ensemble_{run_tag}/cv/cv_probs.csv`

### Step 8 — Ensemble training

Same script, `train` option. Fits the ensemble SVM on the full meta-matrix.

```bash
sbatch --job-name=ensemble_train \
    gen_ensemble_svm.sbatch \
    train \
    data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30/meta_matrix/meta_matrix.csv \
    linear 10 none none none nostd false
```

Output: `{config_dir}/ensemble_{run_tag}/models/ensemble.pkl`. Needs only 4G — the meta-matrix has F columns (one per feature) and inference is a single linear SVM pass.

### Step 9 — Per-feature external test

`gen_svm_test.sbatch` calls `gen_svm_test.py`. Tests a trained per-feature model on a different cohort's feature matrix. Runs as a SLURM array because loading full feature matrices and replaying fitted PCA/GC transforms is memory-intensive.

```bash
sbatch --array=0-6  --mem=8G  gen_svm_test.sbatch \
    data/matrix_internal_ly_test \
    data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30 \
    prevelynch_test \
    30

sbatch --array=7-12 --mem=32G gen_svm_test.sbatch \
    data/matrix_internal_ly_test \
    data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30 \
    prevelynch_test \
    30
```

Arguments: test_matrix_folder, model_config_dir (contains `models/`), dataset_dir (output subfolder name under `test/`), mapq_filter.

Resolved paths:
- test matrix: `{test_matrix_folder}/by_feature/matrix_{feature}_mapq{N}.csv`
- model pickle: `{model_config_dir}/models/{config_name}_{feature}.pkl`
- output: `{model_config_dir}/test/{dataset_dir}/{feature}_probs.csv`

The pickle's stored `feature_columns` list aligns the test matrix to the training columns. The fitted PCA and GC correction transforms are replayed on the test data before inference.

### Step 10 — Ensemble external test

Two sub-steps: first build a test meta-matrix from the per-feature test predictions, then run ensemble inference on it.

#### a. Build the test meta-matrix:
```bash
python scripts/generalized_processing/meta_matrix_build.py \
    data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30/test/prevelynch_test/ \
    data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30/ensemble_nostd/test/prevelynch_test_meta/ \
    data/manifest/internal_metadata_ly_filtered_test.csv
```

#### b. Run ensemble inference:
```bash
sbatch gen_ensemble_svm.sbatch \
    test \
    .../ensemble_nostd/test/prevelynch_test_meta/meta_matrix.csv \
    linear 10 \
    prevelynch_test \
    .../ensemble_nostd/models/ensemble.pkl \
    none nostd false
```

Output: `{config_dir}/ensemble_{run_tag}/test/{dataset_tag}_preds.csv`. Needs only 4G.

### Summary

```
Step 4:  gen_svm.sbatch cv    (array 0-13, per cohort)   -  cv/{feature}_cv_probs.csv
Step 5:  gen_svm.sbatch train (array 0-13, per cohort)   -  models/{config}_{feature}.pkl
Step 6:  meta_matrix_build.py (cv probs)                 -  meta_matrix/meta_matrix.csv
Step 7:  gen_ensemble_svm.sbatch cv                      -  ensemble_{tag}/cv/cv_probs.csv
Step 8:  gen_ensemble_svm.sbatch train                   -  ensemble_{tag}/models/ensemble.pkl
Step 9:  gen_svm_test.sbatch  (array 0-13, per target)   -  test/{dataset}/{feature}_probs.csv
Step 10a: meta_matrix_build.py (test probs)              -  ensemble_{tag}/test/{dataset}_meta/
Step 10b: gen_ensemble_svm.sbatch test                   -  ensemble_{tag}/test/{dataset}_preds.csv
```

Steps 4–8 are run once per training cohort. Steps 9–10 are run once per training cohort × test cohort combination.

The complete pipeline command sequence from static reference to performing 10x10 cv on the cristiano cohort we used is as follows:
```
<!-- prerequisite files -->
sbatch scripts/generalized_processing/static_script.sbatch

<!-- compute feature vectors per sample -->
sbatch --array=1-459%20 \
    scripts/generalized_processing/compute_features_mapq.sbatch \
    30 false data/manifest/Cristiano_metadata.csv data/cristiano_features

<!-- build feature matrices -->
sbatch --array=0-6 --mem=8G scripts/generalized_processing/matrix_build_array.sbatch 30 data/cristiano_features data/matrix/by_feature data/manifest/Cristiano_metadata.csv
sbatch --array=7-12 --mem=32G scripts/generalized_processing/matrix_build_array.sbatch 30 data/cristiano_features data/matrix/by_feature data/manifest/Cristiano_metadata.csv
sbatch --array=13-13 --mem=64G scripts/generalized_processing/matrix_build_array.sbatch 30 data/cristiano_features data/matrix/by_feature data/manifest/Cristiano_metadata.csv

<!-- perform cross-validation on per-feature models -->
sbatch --array=0-6  --mem=8G  gen_svm.sbatch cv data/matrix linear true true 150 10 30
sbatch --array=7-12 --mem=32G gen_svm.sbatch cv data/matrix linear true true 150 10 30
sbatch --array=13   --mem=64G gen_svm.sbatch cv data/matrix linear true true 150 10 30

<!-- build meta-matrix -->
python scripts/generalized_processing/meta_matrix_build.py data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30/cv/ data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30/meta_matrix/ data/manifest/Cristiano_metadata.csv

<!-- perform cross-validation on the ensemble model -->
sbatch gen_ensemble_svm.sbatch cv data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30/meta_matrix/meta_matrix.csv linear 10 none none none nostd false
```

Additionally, to train the ensemble to use for external testing:
```
<!-- train per-feature models on the full cristiano dataset -->
sbatch --array=0-6  --mem=8G  gen_svm.sbatch train data/matrix linear true true 150 10 30
sbatch --array=7-12 --mem=32G gen_svm.sbatch train data/matrix linear true true 150 10 30
sbatch --array=13   --mem=64G gen_svm.sbatch train data/matrix linear true true 150 10 30

<!-- train the ensemble on the meta-matrix built from CV probs -->
sbatch gen_ensemble_svm.sbatch train data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30/meta_matrix/meta_matrix.csv linear 10 none none none nostd false
```

And to test the trained ensemble on external prevelynch cohort:
```
<!-- per-feature test predictions -->
sbatch --array=0-6  --mem=8G  gen_svm_test.sbatch data/matrix_internal_ly_test data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30 prevelynch_test 30
sbatch --array=7-12 --mem=32G gen_svm_test.sbatch data/matrix_internal_ly_test data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30 prevelynch_test 30
sbatch --array=13   --mem=64G gen_svm_test.sbatch data/matrix_internal_ly_test data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30 prevelynch_test 30

<!-- build test meta-matrix from per-feature test probs -->
python scripts/generalized_processing/meta_matrix_build.py data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30/test/prevelynch_test/ data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30/ensemble_nostd/test/prevelynch_test_meta/ data/manifest/internal_metadata_ly_filtered_test.csv

<!-- ensemble test predictions -->
sbatch gen_ensemble_svm.sbatch test data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30/ensemble_nostd/test/prevelynch_test_meta/meta_matrix.csv linear 10 prevelynch_test data/matrix/svm_by_feature/svm_linear_pca150.0_gc_mapq30/ensemble_nostd/models/ensemble.pkl none nostd false
```