import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale
from statsmodels.nonparametric.smoothers_lowess import lowess

class MatrixProcessor:
    def __init__(self, mx: pd.DataFrame, gc_file: pd.DataFrame):
        '''
        mx: DataFrame containing the feature matrix with samples as rows and region features as columns
        gc_file: DataFrame containing GC content for each region feature
        '''
        
        self.gc_content = gc_file
        # TODO delete this after removing length from gc generator script
        self.matrix = mx
        self.sample_ids = self.matrix.index.tolist()
        self.features = self.matrix.columns.tolist()
    
    def standardize(self):
        self.matrix = pd.DataFrame(scale(self.matrix), index=self.sample_ids, columns=self.features)
        return self.matrix
    
    def GC_correction(self):
        corrected_matrix = self.matrix.copy()
        
        for feature in self.features:
            gc_values = self.gc_content[feature]
            feature_values = self.matrix[feature]
            
            # Fit lowess model
            lowess_fit = lowess(feature_values, gc_values, frac=0.75)
            fitted_values = lowess_fit[:, 1]
            
            # Correct the feature values by subtracting the fitted values
            corrected_matrix[feature] = feature_values - fitted_values
        self.matrix = corrected_matrix
        return self.matrix
    
    def PCA(self):
        pca = PCA(n_components=0.95)  
        pca_result = pca.fit_transform(self.matrix)
        pca_df = pd.DataFrame(pca_result, index=self.sample_ids)
        return pca_df
    
    def process(self):
        self.standardize()
        self.GC_correction()
        self.PCA()
        return self.matrix

