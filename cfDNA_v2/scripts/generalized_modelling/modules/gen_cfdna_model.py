import gc
import pandas as pd
import numpy as np
from sklearn import svm
from sklearn.model_selection import train_test_split, RepeatedStratifiedKFold
from sklearn.decomposition import PCA
from statsmodels.nonparametric.smoothers_lowess import lowess
import warnings
warnings.filterwarnings("once", category=RuntimeWarning, message="mean of empty slice")


class CFDNAModel:
    '''
    Class for training, testing and perfroming repeated stratified kfold CV with SVM models on cfDNA fragment data.
    '''
    metadata_cols = ['disease', 'dataset', 'material', 'stage', 'cancer_true']
    gc_correct_features = ['pfe', 'coverage', 'ends', 'ocf', 'ifs', 'wps']

    def __init__(self, mx: pd.DataFrame, gc_content: pd.DataFrame = None,
                 feature: str = None, kernel: str = None, gc_correction: bool = False,
                 pca: bool = False, pca_components: float = 0.95, cv_repeats: int = 10,
                 standardize: bool = True):

        metadata_mx = mx[[c for c in self.metadata_cols if c in mx.columns]].copy()
        self.present_metadata_cols = metadata_mx.columns.tolist()
        if 'cancer_true' in metadata_mx.columns:
            self.labels = metadata_mx['cancer_true'].astype(int).values
        else:
            self.labels = metadata_mx['disease'].apply(lambda x: 0 if x in ('Healthy', 'Control', 'ctr') else 1).values
        self.feature = feature
        self.gc_correction = gc_correction
        self.gc_content = gc_content
        X = mx.drop(columns=[c for c in self.metadata_cols if c in mx.columns])

        self.pca = pca
        self.pca_components = float(pca_components)
        if self.pca_components > 1:
            if not self.pca_components.is_integer():
                raise ValueError(f"PCA component count must be an integer when > 1, got {pca_components}")
            self.pca_components = int(self.pca_components)
        self.matrix = X
        self.sample_ids = X.index.tolist()
        self.features = X.columns.tolist()
        self.kernel = kernel
        self.cv_repeats = cv_repeats
        self.standardize = standardize

    def build_estimator(self):
        # create SVM model with specified kernel, with Platt scaling to get probabilities
        if self.kernel is None:
            return svm.SVC(probability=True, random_state=1)
        return svm.SVC(probability=True, kernel=self.kernel, random_state=1)

    def fit_stdizer(self, X_train: np.ndarray):
        # standardize input training matrix and return standardized matrix
        # also returns the training mean and std for standardizing test set
        mean = np.nanmean(X_train, axis=0, dtype=np.float64).astype(np.float32)
        std = np.nanstd(X_train, axis=0, ddof=0, dtype=np.float64).astype(np.float32)
        std[std == 0] = 1.0
        X_train_scaled = np.nan_to_num((X_train - mean) / std, copy=True).astype(np.float32, copy=False)
        return X_train_scaled, {'mean': mean, 'std': std}

    def standardize_vectors(self, X_train: np.ndarray, X_test: np.ndarray):
        mean = np.nanmean(X_train, axis=0, dtype=np.float64).astype(np.float32)
        std = np.nanstd(X_train, axis=0, ddof=0, dtype=np.float64).astype(np.float32)
        std[std == 0] = 1.0
        X_train = np.nan_to_num((X_train - mean) / std, copy=True).astype(np.float32, copy=False)
        X_test = np.nan_to_num((X_test - mean) / std, copy=True).astype(np.float32, copy=False)
        return X_train, X_test

    def gc_correct(self, X):
        '''
        Perform GC correction on input matrix.
        Uses lowess smoothing to fit a curve to the relationship between per-region gc content and coverage, per sample.
        The fitted curve is then subtracted from each sample's coverage to remove gc bias.
        Returns the corrected matrix.
        '''
        region_ids = range(X.shape[1])
        gc_values = self.gc_content.loc[region_ids].reset_index(drop=True)['gc_content'].values
        new_rows = []
        for sample in X.index:
            coverage = X.loc[sample].values
            smoothed = lowess(coverage, gc_values, frac=0.75)
            fitted = np.interp(gc_values, smoothed[:, 0], smoothed[:, 1])
            new_rows.append(coverage - fitted)
        return pd.DataFrame(new_rows, index=X.index, columns=X.columns)

    def split(self, test_size=0.2, random_state=1):
        return train_test_split(self.matrix, self.labels, test_size=test_size, random_state=random_state)

    def cv_svm(self):
        '''
        kx10 repeated stratified CV returning one averaged OOF probability per sample.
        Order within each fold: standardize (train stats) -> GC correct -> PCA -> SVM.
        '''
        X_values = self.matrix.to_numpy(dtype=np.float32, copy=False)
        y = self.labels
        rskf = RepeatedStratifiedKFold(n_splits=10, n_repeats=self.cv_repeats, random_state=1)
        prob_sum = np.zeros(len(y), dtype=float)
        prob_count = np.zeros(len(y), dtype=int)
        total_folds = 10 * self.cv_repeats
        fold_idx = 0

        for train_i, test_i in rskf.split(X_values, y):
            fold_idx += 1
            print(f"Fold {fold_idx}/{total_folds}", flush=True)
            X_train = X_values[train_i].copy()
            X_test = X_values[test_i].copy()
            y_train = y[train_i]

            if self.standardize:
                X_train_scaled, X_test_scaled = self.standardize_vectors(X_train, X_test)
            else:
                X_train_scaled, X_test_scaled = X_train, X_test

            if self.gc_correction:
                X_train_scaled = self.gc_correct(pd.DataFrame(X_train_scaled, columns=self.features)).to_numpy(dtype=np.float32, copy=True)
                X_test_scaled = self.gc_correct(pd.DataFrame(X_test_scaled, columns=self.features)).to_numpy(dtype=np.float32, copy=True)
                X_train_scaled = np.nan_to_num(X_train_scaled, copy=True).astype(np.float32, copy=False)
                X_test_scaled = np.nan_to_num(X_test_scaled, copy=True).astype(np.float32, copy=False)

            if self.pca:
                pca = PCA(n_components=self.pca_components, svd_solver='randomized', random_state=1)
                X_train_scaled = pca.fit_transform(X_train_scaled)
                X_test_scaled = pca.transform(X_test_scaled)

            model = self.build_estimator()
            model.fit(X_train_scaled, y_train)
            prob_sum[test_i] += model.predict_proba(X_test_scaled)[:, 1]
            prob_count[test_i] += 1

            del X_train, X_test, X_train_scaled, X_test_scaled, model
            gc.collect()

        probabilities = np.divide(
            prob_sum, prob_count,
            out=np.full(len(y), np.nan, dtype=float),
            where=prob_count > 0,
        )
        return {'probabilities': probabilities, 'sample_ids': self.sample_ids, 'labels': y}

    def train_model(self):
        '''
        Fit one SVM on the input matrix and return a pickled model.
        Pickle contains the config parameters, standardization/scaler, Pca model and trained SVM.
        '''
        X_values = self.matrix.to_numpy(dtype=np.float32, copy=False)
        y = self.labels

        if self.standardize:
            X_scaled, scaler_state = self.fit_stdizer(X_values.copy())
        else:
            X_scaled = X_values.copy()
            scaler_state = None

        if self.gc_correction:
            X_scaled = self.gc_correct(pd.DataFrame(X_scaled, columns=self.features)).to_numpy(dtype=np.float32, copy=True)

        X_scaled = np.nan_to_num(X_scaled, copy=True).astype(np.float32, copy=False)

        pca_model = None
        if self.pca:
            pca_model = PCA(n_components=self.pca_components)
            X_scaled = pca_model.fit_transform(X_scaled)

        model = self.build_estimator()
        model.fit(X_scaled, y)

        return {
            'model': model,
            'feature': self.feature,
            'kernel': self.kernel,
            'gc_correction': self.gc_correction,
            'pca': self.pca,
            'pca_components': self.pca_components,
            'feature_columns': self.features,
            'sample_ids': self.sample_ids,
            'labels': y,
            'scaler': scaler_state,
            'pca_model': pca_model,
            'metadata_cols': self.present_metadata_cols,
            'classes_': model.classes_.tolist(),
        }

    def predict(self, trained_model: dict) -> pd.DataFrame:
        '''
        Predict probabilities using input trained model dict. 
        If gc_correction was used in training, the same gc_content file must be provided to the testing call.

        Returns a Dataframe with sample_id, predictedprobability, and true label (not predicted label).
        '''
        X = self.matrix.to_numpy(dtype=np.float32, copy=False)
        scaler = trained_model['scaler']
        if scaler is not None:
            X = np.nan_to_num((X - scaler['mean']) / scaler['std'], copy=True).astype(np.float32, copy=False)

        if trained_model.get('gc_correction'):
            if self.gc_content is None:
                raise ValueError("gc_content must be provided to CFDNAModel for predict when the trained model uses gc_correction")
            X_df = pd.DataFrame(X, index=self.sample_ids, columns=self.features)
            X = self.gc_correct(X_df).to_numpy(dtype=np.float32, copy=True)
            X = np.nan_to_num(X, copy=True).astype(np.float32, copy=False)

        if trained_model.get('pca_model') is not None:
            X = trained_model['pca_model'].transform(X)

        probs = trained_model['model'].predict_proba(X)[:, 1]
        return pd.DataFrame({
            'sample_id': self.sample_ids,
            'probability': probs,
            'label': self.labels
        })
