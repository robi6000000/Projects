import pandas as pd
import numpy as np
from sklearn.preprocessing import scale
from statsmodels.nonparametric.smoothers_lowess import lowess

class MatrixProcessor:
    def __init__(self, mx_file):
        self.matrix = pd.read_csv(mx_file, index_col=0)
        self.sample_ids = self.matrix.index.tolist()
        self.features = self.matrix.columns.tolist()
    
    def standardize(self):
        self.matrix = pd.DataFrame(scale(self.matrix), index=self.sample_ids, columns=self.features)
        return self.matrix
    
    def GC_correction(self):
        pass
