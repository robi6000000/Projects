# Current state of Zhou pipeline for feature extraction and modelling

### files needed:
## data processing
- `static_script.sh` - script for downloading all needed files that are universal to all samples 
(reference genome, annotation, open chromatin regions, ...)

- `sample_prep_script.sh` - script for creating sample-specific files needed for matrix building

- `sample_features_script.py` - script for calculating feature vectores for each sample 
                              
    - called by `sample_prep_script.sh`

- `matrix_build_script.py` - script for building feature matrix for one feature 

- `matrix_build_array.sbatch` - script for building matrices for all samples in parallel on cluster
    - calls `matrix_build_script.py` for each of the 10 features

## modelling
 - `calculate_gc_content.sbatch` and `calculate_gc_content.py` - script for calculating gc content for each region
 - `gc_content_per_region.csv` - output of gc content calculation, used for gc correction in **MatrixProcessor**

 - `matrix_processor.py` - class containing methods for **standardizing**, **gc** correction and **pca**