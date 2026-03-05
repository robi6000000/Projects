import pandas as pd
import numpy as np
from sklearn.preprocessing import scale
from statsmodels.nonparametric.smoothers_lowess import lowess
from pyfaidx import Fasta

class MatrixProcessor:
    def __init__(self, matrix_path, metadata_cols=None):
        """
        Load feature matrix and separate metadata from features.
        
        Args:
            matrix_path: Path to CSV with features and metadata
            metadata_cols: List of metadata column names (default: standard set)
        """
        print(f"Loading matrix from {matrix_path}")
        self.matrix_full = pd.read_csv(matrix_path)
        
        # Default metadata columns
        if metadata_cols is None:
            metadata_cols = ['sample_id', 'seqrun_id', 'sample_name', 'stage', 'disease', 'tissue']
        
        self.metadata_cols = metadata_cols
        self.metadata = self.matrix_full[metadata_cols].copy()
        
        # Feature matrix only (numeric data)
        feature_cols = [col for col in self.matrix_full.columns if col not in metadata_cols]
        self.features = self.matrix_full[feature_cols].copy()
        self.feature_names = feature_cols
        
        print(f"Matrix loaded: {self.features.shape[0]} samples × {self.features.shape[1]} features")
        print(f"Metadata columns: {metadata_cols}")
    
    def standardize(self):
        """
        Z-score standardization using sklearn.preprocessing.scale
        Zhou et al.: "standardized with a Z-score using the preprocessing.scale function"
        """
        print("\nStandardizing features (Z-score)...")
        self.features = pd.DataFrame(
            scale(self.features),
            index=self.features.index,
            columns=self.feature_names
        )
        print("✓ Standardization complete")
        return self
    
    def filter_low_variance(self, threshold=0.01):
        """
        Remove features with very low variance (near-constant across samples)
        Useful before PCA/modeling to reduce dimensionality
        
        Args:
            threshold: Minimum variance threshold (default 0.01)
        """
        print(f"\nFiltering low variance features (threshold={threshold})...")
        variances = self.features.var()
        keep_features = variances[variances > threshold].index
        
        removed = len(self.feature_names) - len(keep_features)
        print(f"Removing {removed} low-variance features ({removed/len(self.feature_names)*100:.1f}%)")
        
        self.features = self.features[keep_features]
        self.feature_names = list(keep_features)
        
        print(f"✓ Remaining features: {len(self.feature_names)}")
        return self
    
    def gc_correction(self, openchrom_bed_path, hg19_fasta_path, span=0.75):
        """
        GC content bias correction using LOWESS regression.
        Zhou et al.: "GC% covariates were regressed from the original fragmentation 
        pattern scores of each open chromatin region using locally weighted smoothing 
        linear regression (Lowess) with a span of 0.75"
        
        NOTE: This is complex and computationally expensive!
        - Requires calculating GC% for 561,414 regions
        - Applies LOWESS regression per feature
        - Takes significant time (~30-60 minutes)
        
        Args:
            openchrom_bed_path: Path to open chromatin regions BED file
            hg19_fasta_path: Path to hg19.fa reference genome
            span: LOWESS smoothing span (default 0.75)
        """
        print("\n⚠ GC correction is computationally expensive and may take 30-60 minutes")
        print("Calculating GC content for 561,414 regions...")
        
        # Load open chromatin regions
        openchrom = pd.read_csv(
            openchrom_bed_path,
            sep="\t",
            header=None,
            names=["chrom", "start", "end", "region_id"]
        )
        
        # Calculate GC content for each region
        genome = Fasta(hg19_fasta_path)
        gc_content = []
        
        for idx, row in openchrom.iterrows():
            if idx % 50000 == 0:
                print(f"  Progress: {idx}/{len(openchrom)} regions...")
            
            try:
                seq = genome[row['chrom']][row['start']:row['end']].seq.upper()
                gc = (seq.count('G') + seq.count('C')) / len(seq) if len(seq) > 0 else 0
                gc_content.append(gc)
            except:
                gc_content.append(0)
        
        openchrom['gc_content'] = gc_content
        
        # Map region_id to GC content
        region_to_gc = dict(zip(openchrom['region_id'], openchrom['gc_content']))
        
        print("\nApplying LOWESS regression to remove GC bias...")
        # For each feature column, regress out GC effect
        # This assumes feature names encode region_id
        
        corrected_features = self.features.copy()
        
        for col_idx, col in enumerate(self.feature_names):
            if col_idx % 10000 == 0:
                print(f"  Progress: {col_idx}/{len(self.feature_names)} features...")
            
            # Extract region_id from feature name (pattern: prefix_region_XXXX)
            # Adjust this parsing based on your actual feature naming
            if 'region_' in col:
                region_id = int(col.split('region_')[1])
                
                if region_id in region_to_gc:
                    gc_val = region_to_gc[region_id]
                    feature_vals = self.features[col].values
                    
                    # Apply LOWESS: fit GC vs feature values, subtract fitted values
                    try:
                        gc_array = np.full(len(feature_vals), gc_val)
                        lowess_fit = lowess(feature_vals, gc_array, frac=span, return_sorted=False)
                        corrected_features[col] = feature_vals - lowess_fit
                    except:
                        pass  # Keep original if LOWESS fails
        
        self.features = corrected_features
        print("✓ GC correction complete")
        return self
    
    def get_processed_matrix(self, include_metadata=True):
        """
        Get the processed feature matrix, optionally with metadata
        
        Args:
            include_metadata: If True, merge metadata back with features
        
        Returns:
            DataFrame with processed features (and metadata if requested)
        """
        if include_metadata:
            return pd.concat([self.metadata, self.features], axis=1)
        else:
            return self.features.copy()
    
    def save(self, output_path, include_metadata=True):
        """
        Save processed matrix to CSV
        
        Args:
            output_path: Where to save the processed matrix
            include_metadata: If True, include metadata columns
        """
        print(f"\nSaving processed matrix to {output_path}")
        matrix_to_save = self.get_processed_matrix(include_metadata=include_metadata)
        matrix_to_save.to_csv(output_path, index=False)
        print(f"✓ Saved: {matrix_to_save.shape[0]} samples × {matrix_to_save.shape[1]} columns")


# Example usage
if __name__ == "__main__":
    # Load and process matrix
    processor = MatrixProcessor("./data/matrix/feature_matrix_with_metadata.csv")
    
    # Apply preprocessing steps
    processor.standardize()
    processor.filter_low_variance(threshold=0.01)
    
    # Optional: GC correction (expensive!)
    # processor.gc_correction(
    #     openchrom_bed_path="./data/processing/openchrom_with_id.bed",
    #     hg19_fasta_path="./data/hg19/hg19.fa"
    # )
    
    # Save processed matrix
    processor.save("./data/matrix/feature_matrix_processed.csv")
    
    # Or get as DataFrame for modeling
    X = processor.get_processed_matrix(include_metadata=False)
    y = processor.metadata['stage']