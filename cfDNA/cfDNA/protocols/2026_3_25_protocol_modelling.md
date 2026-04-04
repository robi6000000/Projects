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

## The following model results were only 1x10fold crossvalidated


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

### linear kernel - standardizing, gc correction in init, no pca
```
    feature       auc  zhou_auc
0  coverage  0.929177    0.9638
1      ends  0.932841    0.9639
2       fsd  0.896904    0.9271
3       wps  0.957498    0.9658
4       ifs  0.956201    0.9653
5       pfe  0.957823    0.9579
6    length  0.905435    0.8741
7       edm  0.957384    0.9736
8       ocf  0.952555    0.9467
9       fsr  0.775793    0.9441
```

### fsr is problematic:
- data_leak version (standardizing, gc correction and pca in init):

linear kernel, gc_correction==False, pca==True (n_components=0.95):
```
  feature       auc  zhou_auc
0     fsr  0.493311    0.9441
```
linear kernel, gc_correction==False, pca==True (n_components=0.99):
```
  feature       auc  zhou_auc
0     fsr  0.689645    0.9441
```
linear kernel, gc_correction==True, pca==True (n_components=0.95):
```
  feature       auc  zhou_auc
0     fsr  0.491431    0.9441
```
linear kernel, gc_correction==True, pca==True (n_components=0.99):
```
  feature       auc  zhou_auc
0     fsr  0.690427    0.9441
```
linear kernel, gc_correction==False, pca==False:
```
  feature       auc  zhou_auc
0     fsr  0.776136    0.9441
```
linear kernel, gc_correction==True, pca==False:
```
  feature       auc  zhou_auc
0     fsr  0.775698    0.9441
```

rbf kernel, gc_correction==False, pca==True (n_components=0.95):
```
  feature       auc  zhou_auc
0     fsr  0.746422    0.9441
```

rbf kernel, gc_correction==False, pca==True (n_components=0.99):
```
  feature       auc  zhou_auc
0     fsr  0.714894    0.9441
```
rbf kernel, gc_correction==True, pca==True (n_components=0.95):
```
  feature       auc  zhou_auc
0     fsr  0.745133    0.9441
```
rbf kernel, gc_correction==True, pca==True (n_components=0.99):
```
  feature       auc  zhou_auc
0     fsr  0.716249    0.9441
```
rbf kernel, gc_correction==False, pca==False:
```
  feature       auc  zhou_auc
6     fsr  0.717546    0.9441
```
rbf kernel, gc_correction==True, pca==False:
```
  feature       auc  zhou_auc
0     fsr  0.719464    0.9441
``` 
rbf kernel, gc_correction==False, pca==True (n_components=0.9):
```
  feature       auc  zhou_auc
0     fsr  0.767462    0.9441
```
rbf kernel, gc_correction==True, pca==True (n_components=0.9):
```
  feature       auc  zhou_auc
0     fsr  0.767567    0.9441
```
rbf kernel, gc_correction==False, pca==True (n_components=0.85):
```
  feature       auc  zhou_auc
0     fsr  0.784018    0.9441
```

### Overvall best:
- linear, gc==False, pca==False