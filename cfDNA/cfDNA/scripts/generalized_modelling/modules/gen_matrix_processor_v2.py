import pandas as pd
import numpy as np
from sklearn.preprocessing import scale
from statsmodels.nonparametric.smoothers_lowess import lowess

class MatrixProcessor:
    metadata_cols = ['disease', 'dataset', 'material', 'stage', 'cancer_true']

    def __init__(self, mx: pd.DataFrame, gc_content: pd.DataFrame):
        """
        mx: DataFrame containing the feature matrix with samples as rows and region features as columns
        gc_content: DataFrame containing GC content for each region feature
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
    

