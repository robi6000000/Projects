import gc
import pandas as pd
import numpy as np
from sklearn import svm
from sklearn.model_selection import train_test_split, RepeatedStratifiedKFold
from sklearn.decomposition import PCA
from modules.gen_matrix_processor import MatrixProcessor
import warnings
warnings.filterwarnings("once", category=RuntimeWarning, message="mean of empty slice")


class CFDNAModel:
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

        self.processor = MatrixProcessor(X, gc_content) if gc_content is not None else None
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

    def _build_estimator(self):
        if self.kernel is None:
            return svm.SVC(probability=True, random_state=1)
        return svm.SVC(probability=True, kernel=self.kernel, random_state=1)

    @staticmethod
    def _fit_standardizer(X_train: np.ndarray):
        mean = np.nanmean(X_train, axis=0, dtype=np.float64).astype(np.float32)
        std = np.nanstd(X_train, axis=0, ddof=0, dtype=np.float64).astype(np.float32)
        std[std == 0] = 1.0
        X_train_scaled = np.nan_to_num((X_train - mean) / std, copy=True).astype(np.float32, copy=False)
        return X_train_scaled, {'mean': mean, 'std': std}

    @staticmethod
    def _standardize_arrays(X_train: np.ndarray, X_test: np.ndarray):
        mean = np.nanmean(X_train, axis=0, dtype=np.float64).astype(np.float32)
        std = np.nanstd(X_train, axis=0, ddof=0, dtype=np.float64).astype(np.float32)
        std[std == 0] = 1.0
        X_train = np.nan_to_num((X_train - mean) / std, copy=True).astype(np.float32, copy=False)
        X_test = np.nan_to_num((X_test - mean) / std, copy=True).astype(np.float32, copy=False)
        return X_train, X_test

    def _gc_correct(self, X):
        region_ids = range(X.shape[1])
        gc_aligned = self.gc_content.loc[region_ids].reset_index(drop=True)
        return MatrixProcessor(X, gc_aligned).GC_correction()

    def _fit_gc_correction(self, X_train_scaled: np.ndarray):
        X_train_df = pd.DataFrame(X_train_scaled, columns=self.features)
        region_ids = range(X_train_df.shape[1])
        gc_aligned = self.gc_content.loc[region_ids].reset_index(drop=True)
        processor = MatrixProcessor(X_train_df, gc_aligned)
        corrected_df = processor.GC_correction()
        gc_state = {
            'mean_fitted': np.mean(list(processor.lowess_fitted_.values()), axis=0).astype(np.float32),
        }
        return corrected_df.to_numpy(dtype=np.float32, copy=True), gc_state

    @staticmethod
    def _apply_gc_state(X_scaled: np.ndarray, gc_state: dict) -> np.ndarray:
        return X_scaled - gc_state['mean_fitted']

    def split(self, test_size=0.2, random_state=1):
        return train_test_split(self.matrix, self.labels, test_size=test_size, random_state=random_state)

    def cv_svm(self):
        """
        10x10 repeated stratified CV returning one averaged OOF probability per sample.
        Order within each fold: standardize (train stats) -> GC correct -> PCA -> SVM.
        """
        X_values = self.matrix.to_numpy(dtype=np.float32, copy=False)
        y = self.labels
        rskf = RepeatedStratifiedKFold(n_splits=10, n_repeats=self.cv_repeats, random_state=1)
        prob_sum = np.zeros(len(y), dtype=float)
        prob_count = np.zeros(len(y), dtype=int)

        for train_i, test_i in rskf.split(X_values, y):
            X_train = X_values[train_i].copy()
            X_test = X_values[test_i].copy()
            y_train = y[train_i]

            if self.standardize:
                X_train_scaled, X_test_scaled = self._standardize_arrays(X_train, X_test)
            else:
                X_train_scaled, X_test_scaled = X_train, X_test

            if self.gc_correction:
                X_train_scaled = self._gc_correct(pd.DataFrame(X_train_scaled, columns=self.features)).to_numpy(dtype=np.float32, copy=True)
                X_test_scaled = self._gc_correct(pd.DataFrame(X_test_scaled, columns=self.features)).to_numpy(dtype=np.float32, copy=True)
                X_train_scaled = np.nan_to_num(X_train_scaled, copy=True).astype(np.float32, copy=False)
                X_test_scaled = np.nan_to_num(X_test_scaled, copy=True).astype(np.float32, copy=False)

            if self.pca:
                pca = PCA(n_components=self.pca_components, svd_solver='randomized', random_state=1)
                X_train_scaled = pca.fit_transform(X_train_scaled)
                X_test_scaled = pca.transform(X_test_scaled)

            model = self._build_estimator()
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
        """Fit one SVM on the full dataset and return a serialisable model dict."""
        X_values = self.matrix.to_numpy(dtype=np.float32, copy=False)
        y = self.labels

        if self.standardize:
            X_scaled, scaler_state = self._fit_standardizer(X_values.copy())
        else:
            X_scaled = X_values.copy()
            scaler_state = None

        gc_state = None
        if self.gc_correction:
            X_scaled, gc_state = self._fit_gc_correction(X_scaled)

        X_scaled = np.nan_to_num(X_scaled, copy=True).astype(np.float32, copy=False)

        pca_model = None
        if self.pca:
            pca_model = PCA(n_components=self.pca_components)
            X_scaled = pca_model.fit_transform(X_scaled)

        model = self._build_estimator()
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
            'gc_state': gc_state,
            'pca_model': pca_model,
            'metadata_cols': self.present_metadata_cols,
            'classes_': model.classes_.tolist(),
        }

    def predict(self, trained_model: dict) -> pd.DataFrame:
        """Apply a trained model to this instance's matrix and return probabilities."""
        X = self.matrix.to_numpy(dtype=np.float32, copy=False)
        scaler = trained_model['scaler']
        if scaler is not None:
            X = np.nan_to_num((X - scaler['mean']) / scaler['std'], copy=True).astype(np.float32, copy=False)

        if trained_model.get('gc_correction') and trained_model.get('gc_state'):
            X = self._apply_gc_state(X, trained_model['gc_state'])

        if trained_model.get('pca_model') is not None:
            X = trained_model['pca_model'].transform(X)

        probs = trained_model['model'].predict_proba(X)[:, 1]
        return pd.DataFrame({
            'sample_id': self.sample_ids,
            'probability': probs,
            'label': self.labels
        })
