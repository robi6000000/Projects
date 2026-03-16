# 1. Title
Good morning, my name is Robert Kenderes and today I'll be presenting my bachelor's thesis on computational analysis of cfDNA fragmentation patterns from liquid biopsy.

# 2. Early Cancer Detection
The main problem this thesis addresses is early cancer detection.
- Most cancers are detected in late stages, which significantly reduces survival rates.
- Standard methods — tissue biopsy and imaging — are invasive, expensive, and non-repeatable.
So what's the solution? Liquid biopsy — a simple blood draw.
- Noninvasive, repeatable, scalable — it can be done frequently and cheaply.
- cfDNA carries information about the presence, stage, and tissue of origin of cancer.

# 3. Cell-free DNA and Fragmentation
But what exactly is cfDNA and why does it carry this information?
- DNA released into the bloodstream during cell death — this is cfDNA.
Looking at this diagram:
- Apoptosis/Necrosis — enzymes cut DNA in blood during cell death.
- DNA wrapped around histones forms a nucleosome — this part is protected from cutting.
- DNA in open chromatin regions is exposed and gets cut.
- Fragments reflect the chromatin structure of the source cell and tissue.
- Cancer cells have a distinct chromatin structure — so their fragments look different.

# 4. Thesis Goal and Structure
Now that we have the background — here's what this thesis does.
- Goal: analyze these fragmentation patterns to predict cancer and apply it to internal data.
The thesis has three phases:
- phase 1 - feature extraction and modelling
- following the methodology of Zhou et al 2024
- collect cfdna samples, extract fragmentation patterns and train an ensemble model to predict cancer
- phase 2 - preprocessing pipeline from FinaleDB 
- finaleddb is a database storing cfdna data from which zhou got their samples. They provide a pipeline that turns raw sequencing data in fastq format into a usable bed file format. 
- this is important because our internal data is stored as fastqs
The last stage is applying the preprocessing pipeline to our data and testing the performance of the pre trained model on them
- Raw reads → features → pre-trained model → cancer prediction — applied to our internal colorectal data.

# 5. Zhou et al. 2024 Workflow
Let's zoom into the first phase — the workflow from Zhou et al. 2024.
- Collect cfDNA data — the Cristiano dataset, 459 samples across 8 cancer types.
- Collect relevant open chromatin regions — around 560 thousand regions from the genome.
- Calculate 10 fragmentation features per region per sample — this gives us 10 large matrices.
- Train an SVM model per feature — each model predicts cancer probability from one feature.
- Combine outputs from each model into an ensemble model — the final cancer prediction.

# 6. Two Example Features
Let's look at two of the simpler features to get a feel for what these patterns look like.
- What fraction of fragments are short vs. long? — this is FSR, fragment size ratio, computed per open chromatin region. Cancer cfDNA tends to be shorter, so this ratio shifts.
- Which DNA bases appear at the ends of fragments? — this is EDM, end motif distribution. There are 256 possible 4-base combinations.
- Some 4-mer sequences are more present in cancer — for example the GAAG motif increases in cancer.

# 7. Window Protection Score
Some features get more complex. A good example is WPS — window protection score.
- For each genomic position, how many fragments span it vs. how many end near it?
- High WPS means the position is nucleosome-protected — more fragments span than end there.
- Low WPS means the position is in an open region — many fragment ends nearby.
- Computed across all 560k open chromatin regions — giving a full nucleosome occupancy map per sample.
- With cancer, altered nucleosome positioning produces a distinct WPS profile.

# 8. Zhou et al. 2024 — Conclusions
So what did Zhou find with these 10 features?
- EDM is strongest in cross-validation but unstable across independent datasets.
- Length-only features are weakest overall.
- Features combining length and coverage — like WPS, IFS and OCF — generalize best.
This is why the ensemble matters:
- AUC > 0.90 for 7 out of 8 cancer types.
- Works in early-stage samples — which brings us back to the original problem.
(point to figure) This is the IFP score on the Cristiano dataset — the separation between healthy and cancer is clear across all cancer types.

# 9. Next Steps
That covers phase 1. Phases 2 and 3 are still ahead.
- .fastq to .bed preprocessing pipeline — building and testing the pipeline for raw sequencing data.
- Apply the pre-trained model to internal colorectal cancer data — this is a generalization test.
The internal data comes from an independent source, different lab and platform. Whether the fragmentation signal holds or not — either result tells us something useful.