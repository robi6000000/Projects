# Notes on jiang 2026:

### 'Holistic determination of ends of cfDNA molecules'

- tldr: 
    - reconstructing 3' ends lost during the sequencing process, monitoring prem(pre-end motifs) and poem(post-end motifs).
    - prem has strong correlation with em3 and poem with em5 (but not completely 1) - new information

    - shorter fragments show lower correlations betwwen prem/em3 and poem/3m5:
        - `The correlations between
            PREM and EM3 and between EM5 and POEM decreased to
            0.73 and 0.81, respectively.`
        - `Additionally, the Pearson’s r
            values between EM5 and POEM and between PREM and EM3
            were lower in the relatively short cfDNA population with a size of
            42–70 nt (Pearson’s r, 0.79 and 0.71) (Figure S2A), compared
            with the long cfDNA populations with a size of 70–166 nt
            (Pearson’s r, 0.94 and 0.92) (Figure S2B) and 166–600 nt (Pearson’s
            r, 0.93 and 0.94) (Figure S2C). These observations suggest
            that the **long cfDNA population exhibits stronger correlations between
            EM5 and POEM, as well as between PREM and EM3**,
            consistent with the patterns observed in size ranges associated
            with the 1 st
            and 2 nd
            peaks.`
    - #### Results:
        - EM5 alone: 0.90
        - Adding EM3, PREM, POEM progressively: up to 0.95
        - 4 mer edm is the best: 
            - `In addition, we studied the performance using
1-, 2-, 3-, 4-, and 5-mer end motifs, respectively. As a result, the
4-end motifs gave rise to the best performance (AUC, 0.954)
compared with other k-mer features (AUC range, 0.860–0.951)`

- feature design: 
    - jiang did edm entropy for 5' 4mers in 2020 paper. (MDE)
    - in 2026 : `We further determined an end-motif
ratio (EMR) metric, which was derived by calculating the ratio of
the cumulative frequency of end motifs that exhibited increased
representation in the HCC group to those that exhibited
decreased representation.`
    - `Limitations of the study
...First, the clinical sample size was relatively small, particularly
for the 4-end sequencing dataset, highlighting the need for larger
cohorts to validate these newly identified fragmentomics
markers....Second, sequencing throughput was limited (<1 million
reads per sample)`