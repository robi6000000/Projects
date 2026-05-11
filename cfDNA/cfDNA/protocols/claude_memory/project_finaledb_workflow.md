---
name: finaledb snakemake workflow
description: Setup, changes made, and MAPQ filtering behaviour for the finaledb snakemake pipeline
type: project
originSessionId: 3823c645-9b3e-436d-8129-ba9c964f2fd8
---
Snakemake workflow used for cfDNA fragment-level feature computation, based on finaledb.

Key changes made during setup:
- MAPQ filtering added as a configurable parameter
- `rerun=false` flag skips samples with already-computed feature files
- Output paths parameterised by mapq threshold (e.g. `_mapq30`)

**Why:** needed to support MAPQ 0/15/30/45 comparisons for hyperparameter selection plots in results.ipynb.
**How to apply:** when adding new MAPQ levels, check the sbatch array script and matrix build scripts for consistent path naming.
