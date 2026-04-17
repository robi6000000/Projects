overleaf student - git versions

zettelkasten system

hyperparamter tuning sekcia pre kazdy parameter tuning

mapq, pca, gc, kernel
-zacat s 
mapq - mapq0, mapq15, mapq30, mapq45

first we need to get cv1x10 results before running 10x10 to save time.
the first hyperparameter to choose is test cv1 on mapq, testing baseline model rbf, no gc , no pca
after we find the best performing mapq on rbf kernel, we will test whether linear performs better on that mapq threshold
after we find the best performing kernel we will use that kernel with the found mapq to test multiple pca options
- we will test no pca, 50 pca, 100pca, 150 pca
after we find the best option 
we test it with gc correction and without gc correction





grammarly setup