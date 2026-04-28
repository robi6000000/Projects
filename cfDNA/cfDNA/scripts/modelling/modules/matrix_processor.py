import pandas as pd
import numpy as np
from sklearn.preprocessing import scale
from statsmodels.nonparametric.smoothers_lowess import lowess

class MatrixProcessor:
    metadata_cols = ['sample_name', 'stage', 'disease', 'tissue', 'cancer_true']

    def __init__(self, mx: pd.DataFrame, gc_content: pd.DataFrame, X_train: pd.DataFrame = None):
        """
        mx: DataFrame containing the feature matrix with samples as rows and region features as columns
        gc_content: DataFrame containing GC content for each region feature
        X_train: DataFrame containing the training feature matrix for standardization
        """
        
        self.gc_content = gc_content
        self.matrix = mx
        self.sample_ids = self.matrix.index.tolist()
        self.features = self.matrix.columns.tolist()

    def standardize(self, X_train: pd.DataFrame = None):
        ''' 
            X_train is given when self.matrix is the test set 
            - uses train mean and std to standardize test set.

            Otherwise uses preprocessing.scale to standardize the matrix.
        '''
        if X_train is not None:
            self.matrix = self.matrix.replace([np.inf, -np.inf], np.nan).fillna(0)
            X_train = X_train.replace([np.inf, -np.inf], np.nan).fillna(0)
            mean = X_train.mean()
            # fixing division by zero
            std = X_train.std(ddof=0).replace(0, 1)
            self.matrix = (self.matrix - mean) / std
        else:
            self.matrix = self.matrix.replace([np.inf, -np.inf], np.nan).fillna(0)
            self.matrix = self.matrix.apply(scale, axis=0)
        return self.matrix
    
    def GC_correction(self):
        """
        Applies gc correction using lowess regression with a span of 0.75.
        Has to be called on standardized matrix.
        Uses delta-based interpolation to skip redundant evaluations.
        """
        gc_content = self.gc_content['gc_content'].values
        gc_range = gc_content.max() - gc_content.min()
        delta = gc_range / 200.0 if len(gc_content) > 5000 else 0.0
        new_rows = []
        self.lowess_fitted_ = {}
        
        for sample in self.sample_ids:
            coverage = self.matrix.loc[sample].values
            lowess_smoothed = lowess(coverage, gc_content, frac=0.75, delta=delta)
            lowess_fitted = np.interp(gc_content, lowess_smoothed[:, 0], lowess_smoothed[:, 1])
            self.lowess_fitted_[sample] = lowess_fitted
            corrected_coverage = coverage - lowess_fitted
            new_rows.append(corrected_coverage)
        self.matrix = pd.DataFrame(new_rows, index=self.sample_ids, columns=self.features)
        return self.matrix
    
    def transform_gc(self, X_new: pd.DataFrame):
        """Transform new data using stored LOWESS fits"""
        # gc_content = self.gc_content['gc_content'].values
        new_rows = []
        
        for sample in X_new.index:
            coverage = X_new.loc[sample].values
            # use mean lowess fit from training data
            mean_fitted = np.mean(list(self.lowess_fitted_.values()), axis=0)
            corrected_coverage = coverage - mean_fitted
            new_rows.append(corrected_coverage)
        
        return pd.DataFrame(new_rows, index=X_new.index, columns=X_new.columns)
    
    # def PCA(self):
    #     pca = PCA(n_components=0.95)  
    #     pca_result = pca.fit_transform(self.matrix)
    #     pca_df = pd.DataFrame(pca_result, index=self.sample_ids)
    #     return pca_df
    
    def process(self):
        self.standardize()
        self.GC_correction()
        return self.matrix

