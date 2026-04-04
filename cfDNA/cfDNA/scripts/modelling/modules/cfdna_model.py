import sys
import os
import pandas as pd
import numpy as np
from sklearn import svm
from sklearn.model_selection import \
    train_test_split, StratifiedShuffleSplit, StratifiedKFold, RepeatedStratifiedKFold, cross_val_score
from sklearn.metrics import roc_auc_score
from sklearn.decomposition import PCA
from modules.matrix_processor import MatrixProcessor


class CFDNAModel:
    metadata_cols = ['seqrun_id','sample_name', 'stage', 'disease', 'tissue']
    gc_correct_features = ['pfe', 'coverage', 'ends', 'ocf', 'ifs','wps']

    def __init__(self, mx: pd.DataFrame, gc_content: pd.DataFrame,  
                 feature: str = None, kernel: str = None, gc_correction: bool = False,
                 pca: bool = False, pca_components: float = 0.95, cv_repeats: int = 10):

        metadata_mx = mx[self.metadata_cols]
        metadata_mx['cancer_true'] = metadata_mx['disease'].apply(lambda x: 0 if x == 'Healthy' else 1)
        self.labels = metadata_mx['cancer_true'].values
        self.feature = feature
        self.gc_correction = gc_correction
        self.gc_content = gc_content
        X = mx.drop(columns=[c for c in self.metadata_cols if c in mx.columns])

        self.processor = MatrixProcessor(X, gc_content)
        self.pca = pca
        self.pca_components = pca_components
        self.matrix = X
        self.sample_ids = X.index.tolist()
        self.features = X.columns.tolist()
        self.kernel = kernel
        self.cv_repeats = cv_repeats
    def _gc_correct(self, X):
        if self.feature == 'fsr':
            fsr_short = X[[c for c in X.columns if c.endswith('65-151')]]
            fsr_medium = X[[c for c in X.columns if c.endswith('151-221')]]
            fsr_long = X[[c for c in X.columns if c.endswith('221-400')]]

            region_ids = range(fsr_short.shape[1])
            gc_aligned = self.gc_content.loc[region_ids].reset_index(drop=True)

            corrected_short = MatrixProcessor(fsr_short, gc_aligned).GC_correction()
            corrected_medium = MatrixProcessor(fsr_medium, gc_aligned).GC_correction()
            corrected_long = MatrixProcessor(fsr_long, gc_aligned).GC_correction()

            X = pd.concat([corrected_short, corrected_medium, corrected_long], axis=1)
        else:
            region_ids = range(X.shape[1])
            gc_aligned = self.gc_content.loc[region_ids].reset_index(drop=True)
            processor2 = MatrixProcessor(X, gc_aligned)
            X = processor2.GC_correction()
        return X
    
    def split(self, test_size=0.2, random_state=1):
        return train_test_split(self.matrix, self.labels, test_size=test_size, random_state=random_state)
    
    def train_svm(self):
        """
        Trains an SVM model using 10x10 repeated stratified cross validation.
        Returns one averaged out-of-fold probability per sample.
        """
        X = self.matrix
        y = self.labels
        rskf = RepeatedStratifiedKFold(n_splits=10, n_repeats=self.cv_repeats, random_state=1)
        prob_sum = np.zeros(len(y), dtype=float)
        prob_count = np.zeros(len(y), dtype=int)
        
        for train_i, test_i in rskf.split(X, y):
            X_train, X_test = X.iloc[train_i], X.iloc[test_i]
            y_train = y[train_i]
            processor_train = MatrixProcessor(X_train, self.gc_content)
            processor_test = MatrixProcessor(X_test, self.gc_content)
            X_train_scaled = processor_train.standardize()
            X_test_scaled = processor_test.standardize(X_train)
            if self.gc_correction:
                X_train_scaled = self._gc_correct(X_train_scaled)
                X_test_scaled = self._gc_correct(X_test_scaled)
            X_train_scaled = X_train_scaled.replace([np.inf, -np.inf], np.nan).fillna(0)
            X_test_scaled = X_test_scaled.replace([np.inf, -np.inf], np.nan).fillna(0)
            if self.pca:
                pca = PCA(n_components=self.pca_components)
                X_train_scaled = pca.fit_transform(X_train_scaled)
                X_test_scaled = pca.transform(X_test_scaled)
            if self.kernel is None:
                model = svm.SVC(probability=True)
            else:
                model = svm.SVC(probability=True, kernel=self.kernel)
            model.fit(X_train_scaled, y_train)
            fold_probs = model.predict_proba(X_test_scaled)[:, 1]
            prob_sum[test_i] += fold_probs
            prob_count[test_i] += 1

        probabilities = np.divide(
            prob_sum,
            prob_count,
            out=np.full(len(y), np.nan, dtype=float),
            where=prob_count > 0,
        )
        
        return {
            'probabilities': probabilities,
            'sample_ids': self.sample_ids,
            'labels': y
        }
