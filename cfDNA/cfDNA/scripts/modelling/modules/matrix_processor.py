import pandas as pd
import numpy as np
from sklearn.preprocessing import scale
from statsmodels.nonparametric.smoothers_lowess import lowess

class MatrixProcessor:
    metadata_cols = ['sample_name', 'stage', 'disease', 'tissue', 'cancer_true']

    def __init__(self, mx: pd.DataFrame, gc_content: pd.DataFrame):
        """
        mx: DataFrame containing the feature matrix with samples as rows and region features as columns
        gc_content: DataFrame containing GC content for each region feature
        """
        
        self.gc_content = gc_content
        self.matrix = mx
        self.sample_ids = self.matrix.index.tolist()
        self.features = self.matrix.columns.tolist()

        # self.mean_ = None
        # self.std_ = None
        # self.lowess_fitted_ = None
    
    def standardize(self):
        self.matrix = self.matrix.fillna(0)
        self.matrix = pd.DataFrame(scale(self.matrix), index=self.sample_ids, columns=self.features)
        return self.matrix
    
    def GC_correction(self):
        """
        Applies gc correction using lowess regression with a span of 0.75.
        Has to be called on standardized matrix. 
        """
        gc_content = self.gc_content['gc_content'].values
        new_rows = []
        self.lowess_fitted_ = {}
        
        for sample in self.sample_ids:
            coverage = self.matrix.loc[sample].values
            lowess_smoothed = lowess(coverage, gc_content, frac=0.75)
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

