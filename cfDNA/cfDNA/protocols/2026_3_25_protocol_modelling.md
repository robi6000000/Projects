# Modelling results

### rbf kernel, no gc_correction, in-loop standardizing and pca, applied to train and test separately, test data unfortunately using its separate standardization
- This is a methodical improvement over Zhou but yields worse results

```
    feature       auc  zhou_auc
0       ocf  0.811061    0.9467
1    length  0.819936    0.8741
2       ifs  0.928890    0.9653
3       fsd  0.888221    0.9271
4       edm  0.949655    0.9736
5       wps  0.928299    0.9658
6      ends  0.937688    0.9639
7  coverage  0.915779    0.9638
8       pfe  0.863201    0.9579
```

### rbf kernel - standardizing, gc correction and pca in init
- currently trying whole-matrix standardization using preprocessing.scale, gc correction and pca applied in initialization of the model.

