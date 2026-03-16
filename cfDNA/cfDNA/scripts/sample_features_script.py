import pandas as pd
import os
import sys
import subprocess
import numpy as np

from pyfaidx import Fasta
import pyBigWig
import itertools

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

# take sample id as input
# take './data/processing/openchrom_with_id.bed' (openchrom with ids) as input
# openchrom_path = './data/processing/openchrom_with_id.bed'
# frag_centroids_openchrom_intersect_path = './data/processing/frag_centroids_openchrom_intersect.bed'
# frag_ends_openchrom_intersect_path = './data/processing/frag_ends_openchrom_intersect.bed'
# frag_ends_ocf_path = './data/processing/frag_ends_ocf.bed'
# hg19_fasta_path = './data/hg19/hg19.fa'

class SampleFeatures:
    def __init__(self, sample_id, openchrom_path, frag_centroids_openchrom_intersect_path, 
                 frag_ends_openchrom_intersect_path, frag_ends_ocf_path, hg19_fasta_path, rerun=False):
        
        self.sample_id = sample_id
        self.openchrom_path = openchrom_path
        self.frag_centroids_openchrom_intersect_path = frag_centroids_openchrom_intersect_path
        self.frag_ends_openchrom_intersect_path = frag_ends_openchrom_intersect_path
        self.frag_ends_ocf_path = frag_ends_ocf_path
        self.hg19_fasta_path = hg19_fasta_path

        self.rerun = rerun
        
        self.load()
        print("Number of unique region IDs:", self.df_region_ids.shape[0])
        print("frag_centroids_openchrom_intersect.shape:", self.frag_centroids_openchrom_intersect.shape)
        print("frag_ends_openchrom_intersect.shape:", self.frag_ends_openchrom_intersect.shape)
        print("frag_ends_ocf.shape:", self.frag_ends_ocf.shape)

        
    def load(self):
        self.openchrom_with_id = pd.read_csv(
            self.openchrom_path,
            sep="\t",
            header=None,
            names=["oc_chrom", "oc_start", "oc_end", "region_id"])
        self.df_region_ids = self.openchrom_with_id[['region_id']].drop_duplicates().reset_index(drop=True)

        self.frag_centroids_openchrom_intersect = pd.read_csv(
            self.frag_centroids_openchrom_intersect_path,
            sep="\t",
            header=None,
            names=["f_chrom", "centroid1", "centroid2", "f_start", "f_end", "score", "strand", "oc_chrom", "oc_start", "oc_end", "region_id"])

        self.frag_ends_openchrom_intersect = pd.read_csv(
            self.frag_ends_openchrom_intersect_path,
            sep="\t",
            header=None,
            names=["f_chrom", "end1", "end2", "end_type", "oc_chrom", "oc_start", "oc_end", "region_id"])

        self.frag_ends_ocf = pd.read_csv(
            self.frag_ends_ocf_path,
            sep="\t",
            header=None,
            names=["chrom", "end1", "end2", "end_type", "oc_start", "oc_end", "region_id", "centroid", "rel_pos"]
        )

        self.hg19_fasta = Fasta(self.hg19_fasta_path)
        # create length column ahead because its used multiple features
        self.frag_centroids_openchrom_intersect["length"] = (
            self.frag_centroids_openchrom_intersect["f_end"] - 
            self.frag_centroids_openchrom_intersect["f_start"]
        )

        # check if folders exist for each sample in cristiano_features
        feature_folders = [
            'length', 'pfe', 'fsr', 'fsd', 'coverage', 'ends', 'ocf', 'ifs', 'wps', 'edm'
        ]
        for folder in feature_folders:
            if not os.path.exists(f'./data/cristiano_features/{folder}'):
                os.makedirs(f'./data/cristiano_features/{folder}')

    def calculate_features(self):
        try:
            self.length = self.get_length()
        except Exception as e:
            print(f"Error calculating length: {e}")
            self.length = None
        try:
            self.pfe = self.get_pfe()
        except Exception as e:
            print(f"Error calculating pfe: {e}")
            self.pfe = None
        try:
            self.fsr = self.get_fsr()
        except Exception as e:
            print(f"Error calculating fsr: {e}")
            self.fsr = None
        try:
            self.fsd = self.get_fsd()
        except Exception as e:
            print(f"Error calculating fsd: {e}")
            self.fsd = None
        try:
            self.coverage = self.get_coverage()
        except Exception as e:
            print(f"Error calculating coverage: {e}")
            self.coverage = None
        try:
            self.ends = self.get_ends()
        except Exception as e:
            print(f"Error calculating ends: {e}")
            self.ends = None
        try:
            self.ocf = self.get_ocf()
        except Exception as e:
            print(f"Error calculating ocf: {e}")
            self.ocf = None
        try:
            self.ifs = self.get_ifs()
        except Exception as e:
            print(f"Error calculating ifs: {e}")
            self.ifs = None
        try:
            self.wps = self.get_wps()
        except Exception as e:
            print(f"Error calculating wps: {e}")
            self.wps = None
        try:
            self.edm = self.get_edm()
        except Exception as e:
            print(f"Error calculating edm: {e}")
            self.edm = None

    def make_feature_vector(self):
        self.calculate_features()
        # vector features: df_pfe, df_cov, df_end, df_ocf, df_ifs
        # merge vecotr features on region_id and flatten
        # matrix features: df_length, df_fsr, df_fsd, df_edm
        region_index = self.df_region_ids["region_id"].values

        # (df, value_column, feature_prefix)
        vec_features = []
        if self.pfe is not None:
            vec_features.append((self.pfe, "pfe", "pfe"))
        if self.coverage is not None:
            vec_features.append((self.coverage, "coverage", "cov"))
        if self.ends is not None:
            vec_features.append((self.ends, "end", "end"))
        if self.ocf is not None:
            vec_features.append((self.ocf, "ocf", "ocf"))
        if self.ifs is not None:
            vec_features.append((self.ifs, "IFS", "ifs"))
        if self.wps is not None:
            vec_features.append((self.wps, "wps", "wps"))
        
        mx_features = []
        if self.length is not None:
            mx_features.append((self.length, None, "len"))
        if self.fsr is not None:
            mx_features.append((self.fsr, None, "fsr"))
        if self.fsd is not None:
            mx_features.append((self.fsd, None, "fsd"))
        if self.edm is not None:
            mx_features.append((self.edm, None, "edm"))

        if len(vec_features) == 0 and len(mx_features) == 0:
            print("error: features failed to compute")
            return None


        feature_vectors = []
        feature_names = []

        for df, col, prefix in vec_features:
            vec = df.loc[region_index, col].to_numpy()
            feature_vectors.append(vec)

            feature_names.extend(
                [f"{prefix}_region_{rid}" for rid in df.index]
            )
        # have to double index matrix features
        for mx, col, prefix in mx_features:
            # print(prefix, mx.shape)
            for chrom in mx.index:
                row = mx.loc[chrom]
                vec = row.to_numpy()
                feature_vectors.append(vec)
                # print(f"  {chrom} vec shape: {vec.shape}")
                feature_names.extend(
                    [f"{prefix}_{chrom}_bin_{b}" for b in row.index]
                )

        # final 1D feature vector
        feature_vector = np.concatenate(feature_vectors)
        feature_vector_df = pd.DataFrame([feature_vector], columns=feature_names)
        feature_vector_df.insert(0, "sample_id", self.sample_id)

        print("Feature vector shape:", feature_vector_df.shape)
        return feature_vector_df

    def get_length(self):
        # check if length already exists in feature folder (for given sample)
        filename = f"./data/cristiano_features/length/{self.sample_id}_length.csv"

        if os.path.exists(filename) and self.rerun == False:
            print(f"Length file already exists {self.sample_id}")
            df_length = pd.read_csv(filename, index_col=0)
            return df_length
        

        # create 33 bins 
        bin_edges = list(range(0, 310, 10)) + [np.inf]
        bin_labels = [f"{i}-{i+10}" for i in range(0, 300, 10)] + [">=300"]


        # get length column, cut into bins
        self.frag_centroids_openchrom_intersect["len_bin"] = pd.cut(
            self.frag_centroids_openchrom_intersect["length"], 
            bins=bin_edges, 
            labels=bin_labels, 
            right=False,
            include_lowest=True)
        
        all_bins = self.frag_centroids_openchrom_intersect["len_bin"].cat.categories
        # group by chromosome and get proporotions
        length_matrix = (self.frag_centroids_openchrom_intersect
                    .groupby(["f_chrom", "len_bin"])
                    .size()
                    .unstack(fill_value=0))
        length_matrix = length_matrix.reindex(columns=all_bins, fill_value=0)
    
        chrom_order = [f"chr{i}" for i in range(1, 23)]
        df_length = length_matrix.reindex(chrom_order, fill_value=0)
        print("df_length:", df_length.shape)
        df_length.to_csv(filename)      # save to feature folder
        return df_length
    
    def get_pfe(self):
        filename = f"./data/cristiano_features/pfe/{self.sample_id}_pfe.csv"
        if os.path.exists(filename) and self.rerun == False:
            print(f"PFE file already exists {self.sample_id}")
            df_pfe = pd.read_csv(filename, index_col=0)
            return df_pfe
        
        bins_pfe = ([0, 100] +
                    list(range(100, 260, 10)) +
                    [np.inf])

        bins_pfe = sorted(set(bins_pfe))
        bin_labels_pfe = [f"{bins_pfe[i]}_{bins_pfe[i+1]}" for i in range(len(bins_pfe)-1)]

        # cut into categories
        self.frag_centroids_openchrom_intersect["pfe_bin"] = pd.cut(
                self.frag_centroids_openchrom_intersect["length"],
                bins=bins_pfe,
                labels=bin_labels_pfe,
                right=False,
                include_lowest=True)

        # calculate P_i
        counts = (self.frag_centroids_openchrom_intersect
                .groupby(["region_id", "pfe_bin"])
                .size()
                .unstack(fill_value=0))

        pfe_proportions = counts.div(counts.sum(axis=1), axis=0)
        pfe_proportions.replace(0, np.nan, inplace=True)

        # get actual PFE per refion
        pfe = - (pfe_proportions * np.log2(pfe_proportions)).sum(axis=1)
        df_pfe = pfe.to_frame("pfe")

        # merge pfe result with region ids to fill blanks
        df_pfe = self.df_region_ids.merge(df_pfe, on="region_id", how="left").set_index('region_id')
        df_pfe["pfe"] = df_pfe["pfe"].fillna(0)

        print("df_pfe:", df_pfe.shape)
        df_pfe.to_csv(filename)     
        return df_pfe
    
    def get_fsr(self):
        filename = f"./data/cristiano_features/fsr/{self.sample_id}_fsr.csv"
        if os.path.exists(filename) and self.rerun == False:
            print(f"FSR file already exists {self.sample_id}")
            df_fsr = pd.read_csv(filename, index_col=0)
            return df_fsr
        
        # bin edges
        bins_fsr = [65, 151, 221, 400]
        bin_labels_fsr = [f"{bins_fsr[i]}_{bins_fsr[i+1]}" for i in range(len(bins_fsr)-1)]
        self.frag_centroids_openchrom_intersect["fsr_bin"] = pd.cut(
            self.frag_centroids_openchrom_intersect["length"],
            bins=bins_fsr,
            labels=bin_labels_fsr,
            right=False,
            include_lowest=True
        )
        all_bins = self.frag_centroids_openchrom_intersect["fsr_bin"].cat.categories
        # get proportions
        counts = (self.frag_centroids_openchrom_intersect
                .groupby(["region_id", "fsr_bin"])
                .size()
                .unstack(fill_value=0))
        counts = counts.reindex(columns=all_bins, fill_value=0)
        df_fsr = counts.div(counts.sum(axis=1), axis=0).fillna(0)

        # merge with region ids to fill blanks
        df_fsr = df_fsr.merge(self.df_region_ids, on="region_id", how="right").set_index('region_id')
        print("df_fsr:", df_fsr.shape)
        df_fsr.to_csv(filename)
        return df_fsr
    
    def get_fsd(self):

        filename = f"./data/cristiano_features/fsd/{self.sample_id}_fsd.csv"
        if os.path.exists(filename) and self.rerun == False:
            print(f"FSD file already exists {self.sample_id}")
            df_fsd = pd.read_csv(filename, index_col=0)
            return df_fsd
        # vin 
        bins_fsd = list(range(65, 405, 5)) 
        self.frag_centroids_openchrom_intersect["fsd_bin"] = pd.cut(
            self.frag_centroids_openchrom_intersect["length"],
            bins=bins_fsd,
            right=False,
            include_lowest=True,

        )
        all_bins = self.frag_centroids_openchrom_intersect["fsd_bin"].cat.categories

        # proportions
        counts = (
            self.frag_centroids_openchrom_intersect
            .groupby(["f_chrom", "fsd_bin"])
            .size()
            .unstack(fill_value=0)
        )
        counts = counts.reindex(columns=all_bins, fill_value=0)
        df_fsd = counts.div(counts.sum(axis=1), axis=0).fillna(0)

        # reorder for possible plotting
        chrom_order = [f"chr{i}" for i in range(1, 23)]
        df_fsd = df_fsd.reindex(chrom_order, fill_value=0)

        print("df_fsd:", df_fsd.shape)
        df_fsd.to_csv(filename)
        return df_fsd
    
    def get_coverage(self):

        filename = f"./data/cristiano_features/coverage/{self.sample_id}_coverage.csv"
        if os.path.exists(filename) and self.rerun == False:
            print(f"Coverage file already exists {self.sample_id}")
            df_cov = pd.read_csv(filename, index_col=0)
            return df_cov
        # for coverage we need to merge with the original region id file at the start, since we are just counting occurences
        centroids_intersect = self.df_region_ids.merge(
            self.frag_centroids_openchrom_intersect,
            on="region_id",
            how="left"
        )

        # group by regions and count each region
        centroids_intersect_grouped = centroids_intersect.groupby("region_id")
        df_cov = pd.DataFrame(centroids_intersect_grouped.size())
        df_cov = df_cov.rename(columns={0: "coverage"})

        # fill the empty regions with 0
        df_cov["coverage"] = df_cov["coverage"].fillna(0)
        print("df_cov:", df_cov.shape)
        df_cov.to_csv(filename)
        return df_cov
    
    def get_ends(self):

        filename = f"./data/cristiano_features/ends/{self.sample_id}_ends.csv"
        if os.path.exists(filename) and self.rerun == False:
            print(f"Ends file already exists {self.sample_id}")
            df_end = pd.read_csv(filename, index_col=0)
            return df_end
        ends_intersect = self.df_region_ids.merge(
            self.frag_ends_openchrom_intersect,
            on="region_id",
            how="left"
        )
        counts = (
            ends_intersect
            .groupby("region_id")
            .size()
        )
        df_end = pd.DataFrame(counts).rename(columns={0: "end"})
        df_end["end"] = df_end["end"].fillna(0)
        print("df_end:", df_end.shape)
        df_end.to_csv(filename)
        return df_end
    
    def get_ocf(self):
        filename = f"./data/cristiano_features/ocf/{self.sample_id}_ocf.csv"
        if os.path.exists(filename) and self.rerun == False:
            print(f"OCF file already exists {self.sample_id}")
            df_ocf = pd.read_csv(filename, index_col=0)
            return df_ocf

        # set boundaries of peak regions aroundn the center
        leftmin, leftmax = -70, -50
        rightmin, rightmax = 50, 70
        df = self.frag_ends_ocf
        df["window"] = np.where((df["rel_pos"] >= leftmin) & (df["rel_pos"] < leftmax),
                                "left",
                                np.where((df["rel_pos"] >= rightmin) & (df["rel_pos"] < rightmax),
                                        "right",
                                        None))

        df = df[df["window"].notna()]
        counts = (df
                .groupby(["region_id", "window", "end_type"])
                .size()
                .unstack(fill_value=0))

        df["window"].value_counts()
        counts["left_term"] = (counts["D"] - counts["U"])
        counts["right_term"] = (counts["U"] - counts["D"])

        ocf = (counts
                .groupby("region_id")[["left_term", "right_term"]]
                .sum()
                .sum(axis=1))
        df_ocf = (self.df_region_ids
                .merge(ocf.rename("ocf"), on="region_id", how="left")
                .fillna(0))
        df_ocf.set_index('region_id', inplace=True)
        print("df_ocf:", df_ocf.shape)
        df_ocf.to_csv(filename)
        return df_ocf
    
    def get_ifs(self):
        filename = f"./data/cristiano_features/ifs/{self.sample_id}_ifs.csv"
        if os.path.exists(filename) and self.rerun == False:
            print(f"IFS file already exists {self.sample_id}")
            df_ifs = pd.read_csv(filename, index_col=0)
            return df_ifs
        # count fragments (n), calculate average lengths (l),
        region_counts = self.frag_centroids_openchrom_intersect.groupby("region_id").size()
        region_avg_length = self.frag_centroids_openchrom_intersect.groupby("region_id")["length"].mean()

        # average length per chromosome (L):
        chrom_avg_length = self.frag_centroids_openchrom_intersect.groupby("f_chrom")["length"].mean()

        # add chrom avg to main df
        self.frag_centroids_openchrom_intersect["chrom_avg"] = self.frag_centroids_openchrom_intersect["f_chrom"].map(chrom_avg_length)
        df_ifs = pd.DataFrame({
            "region_id": region_counts.index,
            "n": region_counts.values,
            "l": region_avg_length.values
        })

        # map chromosome for L
        region_chrom = self.frag_centroids_openchrom_intersect.groupby("region_id")["f_chrom"].first()
        df_ifs["chrom"] = df_ifs["region_id"].map(region_chrom)
        df_ifs["L"] = df_ifs["chrom"].map(chrom_avg_length)

        # compute IFS and fill the missing regions
        df_ifs["IFS"] = df_ifs["n"] * (1 + df_ifs["l"] / df_ifs["L"])
        df_ifs = self.df_region_ids.merge(df_ifs[["region_id", "IFS"]], on="region_id", how="left").fillna(0)
        df_ifs.set_index('region_id', inplace=True)
        print('df_ifs:', df_ifs.shape)
        df_ifs.to_csv(filename)
        return df_ifs
    
    def get_wps(self):
        """
        (WPS), which is the number of DNA fragments completely spanning a 120 bp window centered at a given genomic coordinate
        minus the number of fragments with an endpoint within that
        same window (Figure 2A). As expected, the WPS correlates
        with the locations of nucleosomes within strongly positioned arrays, 
        """

        """
        We define the windowed protection
        score (WPS) of a window of size k as the number of molecules spanning
        the window minus those with an endpoint within the window. We assign
        the determined WPS to the center of the window. For 35–80 bp fragments
        (short fraction, S-WPS), k = 16; for 120–180 bp fragments (long fraction,
        L-WPS), k = 120.
        """
        
        """
        For WPS,[44] according to the genomic coordinate position of each
        cfDNA fragment, a window of 120 bp was slid at 1 bp intervals, and the
        likelihood of each base pair being covered at the whole genome level, fully
        covered (+1), and partially covered (−1), was counted. The mean value of
        all loci within each open chromatin region was calculated.
        """
        filename = f"./data/cristiano_features/wps/{self.sample_id}_wps.csv"
        if os.path.exists(filename) and self.rerun == False:
            print(f"WPS file already exists {self.sample_id}")
            df_wps = pd.read_csv(filename, index_col=0)
            return df_wps
        # load bigwig file for sample (this is created in download script so no need to check and create)
        wps_sample = f"./data/source/cristiano_features/cristiano_wps_download/{self.sample_id}.hg19.wps.mapq30.bw"
        bigwig = pyBigWig.open(wps_sample)
        wps_mean_values = []
        for idx, row in self.openchrom_with_id.iterrows():
            chrom = row['oc_chrom']
            start = row['oc_start']
            end = row['oc_end']
            region_id = row['region_id']
            try:
                wps_mean = bigwig.stats(chrom, start, end, type="mean")[0]
                wps_mean_values.append((region_id, wps_mean))
            except Exception as e:
                # handle none values - no coverage
                wps_mean_values.append((region_id, np.nan))
        bigwig.close()
        df_wps = pd.DataFrame(wps_mean_values, columns=["region_id", "wps"])
        df_wps = self.df_region_ids.merge(df_wps, on="region_id", how="left").set_index('region_id')
        df_wps["wps"] = df_wps["wps"].fillna(0)
        print("df_wps:", df_wps.shape)
        df_wps.to_csv(filename)
        return df_wps

    def _get_edm_motif(self, df, k=4):
        # helper function to get 4 base motif at 5' end positions
        # df input is a groupby dataframe
        chrom = df.name  
        
        if chrom not in self.hg19_fasta:
            return pd.Series([None] * len(df), index=df.index)
        
        seq = self.hg19_fasta[chrom]
        motifs = []
        
        for pos in df['pos'].values:
            if pos >= 0 and pos + k <= len(seq):
                motif = str(seq[pos:pos+k].seq).upper()
                motifs.append(motif)
            else:
                motifs.append(None)
        
        return pd.Series(motifs, index=df.index)


    def get_edm(self):
        filename = f"./data/cristiano_features/edm/{self.sample_id}_edm.csv"
        if os.path.exists(filename) and self.rerun == False:
            print(f"EDM file already exists {self.sample_id}")
            df_edm = pd.read_csv(filename, index_col=0)
            return df_edm
        
        # get end positions
        ends_5 = self.frag_centroids_openchrom_intersect[["f_chrom", "f_start"]].copy()
        ends_5 = ends_5.rename(columns={"f_start": "pos"})
                
        # get motifs at each position
        ends_5['motif'] = (
            ends_5
            .groupby('f_chrom', group_keys=False)
            .apply(self._get_edm_motif, k=4)
        )
        ends_5 = ends_5[ends_5['motif'].notna()]
        ends_5 = ends_5[ends_5['motif'].str.match('^[ACGT]{4}$')]
                
        # count and convert to proportions
        motif_counts = (
            ends_5
            .groupby(['f_chrom', 'motif'])
            .size()
            .unstack(fill_value=0)
        )
        # have tp reindex to include all motifs
        all_motifs = [''.join(p) for p in itertools.product('ACGT', repeat=4)]
        motif_counts = motif_counts.reindex(columns=all_motifs, fill_value=0)

        df_edm = motif_counts.div(motif_counts.sum(axis=1), axis=0)
        
        # Reorder to chr1-chr22
        chrom_order = [f"chr{i}" for i in range(1, 23)]
        df_edm = df_edm.reindex(chrom_order, fill_value=0)
        
        print(f"df_edm: {df_edm.shape}")
        df_edm.to_csv(filename)
        return df_edm
        


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python sample_features_script.py sample_id frag_centroids_openchrom_intersect_path frag_ends_openchrom_intersect_path frag_ends_ocf_path hg19_fasta_path")
        sys.exit(1)
    openchrom_path = './data/processing/openchrom_with_id.bed'
    sample_id = sys.argv[1]
    frag_centroids_openchrom_intersect_path = sys.argv[2]
    frag_ends_openchrom_intersect_path = sys.argv[3]
    frag_ends_ocf_path = sys.argv[4]
    hg19_fasta_path = './data/hg19/hg19.fa'

    sample_features = SampleFeatures(
        sample_id,
        openchrom_path,
        frag_centroids_openchrom_intersect_path,
        frag_ends_openchrom_intersect_path,
        frag_ends_ocf_path,
        hg19_fasta_path,
        rerun=False
    )
    feature_vector_df = sample_features.make_feature_vector()
    data_temp_path = f'./data/rows_sample_temp/{sample_id}_features.csv'
    feature_vector_df.to_csv(data_temp_path, index=False)
    