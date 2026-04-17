# Modelling results

## Best resul ts so far:
```
    feature       auc  zhou_auc
0  coverage  0.906141    0.9638
1       edm  0.962670    0.9736
2      ends  0.938566    0.9639
3       fsd  0.916962    0.9271
4       fsr  0.772759    0.9441
5       ifs  0.949826    0.9653
6    length  0.911886    0.8741
7       ocf  0.900225    0.9467
8       pfe  0.939158    0.9579
9       wps  0.966144    0.9658
```
## WPS implementation needs to be revisited
- currently only using the precomputed wps values from finaledb

kenderes@gen-manager:/data/projects/liquid_biopsy/reads/original$ ls gsca_*_pl* | wc -l
38
kenderes@gen-manager:/data/projects/liquid_biopsy/reads/original$ ls gspp_*_pl* | wc -l
38