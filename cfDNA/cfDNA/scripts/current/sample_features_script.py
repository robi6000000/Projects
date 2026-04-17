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
                 frag_ends_openchrom_intersect_path, frag_ends_ocf_path, frag_wps_intersect_path, hg19_fasta_path,
                 rerun=False, rerun_features: list = None, mapq_filter: int = None, mapq_filter_features: list = None):
        
        self.sample_id = sample_id
        self.openchrom_path = openchrom_path
        self.frag_centroids_openchrom_intersect_path = frag_centroids_openchrom_intersect_path
        self.frag_ends_openchrom_intersect_path = frag_ends_openchrom_intersect_path
        self.frag_ends_ocf_path = frag_ends_ocf_path
        self.frag_wps_intersect_path = frag_wps_intersect_path
        self.hg19_fasta_path = hg19_fasta_path

        self.rerun = rerun
        self.rerun_features = rerun_features
        self.mapq_filter = mapq_filter
        self.mapq_filter_features = mapq_filter_features
        
        self.load()
        print("Number of unique region IDs:", self.df_region_ids.shape[0])
        print("frag_centroids_openchrom_intersect.shape:", self.frag_centroids_openchrom_intersect.shape)
        print("frag_ends_openchrom_intersect.shape:", self.frag_ends_openchrom_intersect.shape)
        print("frag_ends_ocf.shape:", self.frag_ends_ocf.shape)

    def _use_saved_file(self, filename, featurename):
        if os.path.exists(filename):
            if not self.rerun:
                return True
            elif self.rerun_features is None:
                return False
            elif featurename.split('_')[0] in self.rerun_features:
                print(f'recomputing {featurename}')
                return False
            else:
                return True
        else:
            return False

    def _use_mapq_filter(self, feature_name):
        """
        Helper to check whether to use mapq filter for a given feature based on user input.
        """
        if self.mapq_filter is None:
            return False
        elif self.mapq_filter_features is None:
            return True
        else:
            return feature_name in self.mapq_filter_features

    def _filtered(self, df, feature_name):
        """
        Filter df by set mapq score
        """
        if self._use_mapq_filter(feature_name):
            return df[df['score'] >= self.mapq_filter].copy()
        return df

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
            names=["f_chrom", "centroid1", "centroid2", "f_start", "f_end", "score", "strand", "frag_id", "oc_chrom", "oc_start", "oc_end", "region_id"])

        self.frag_ends_openchrom_intersect = pd.read_csv(
            self.frag_ends_openchrom_intersect_path,
            sep="\t",
            header=None,
            names=["f_chrom", "end1", "end2", "end_type", "score", "strand", "frag_id", "oc_chrom", "oc_start", "oc_end", "region_id"])

        self.frag_ends_ocf = pd.read_csv(
            self.frag_ends_ocf_path,
            sep="\t",
            header=None,
            names=["chrom", "end1", "end2", "end_type", "score", "strand", "frag_id", 
                "oc_start", "oc_end", "region_id", "centroid", "rel_pos"]
        )

        self.frag_wps_openchrom_intersect = pd.read_csv(
            self.frag_wps_intersect_path,
            sep="\t",
            header=None,
            names=["f_chrom", "f_start", "f_end", "score", "strand", "frag_id","oc_chrom", "oc_start", "oc_end", "region_id"]
        )

        print("Loading hg19 reference fasta")
        self.hg19_fasta = Fasta(self.hg19_fasta_path)
        # create length column ahead because its used multiple features
        self.frag_centroids_openchrom_intersect["length"] = (
            self.frag_centroids_openchrom_intersect["f_end"] - 
            self.frag_centroids_openchrom_intersect["f_start"]
        )
        self.frag_wps_openchrom_intersect["length"] = (
            self.frag_wps_openchrom_intersect["f_end"] - 
            self.frag_wps_openchrom_intersect["f_start"]
        )

        # check if folders exist for each sample in cristiano_features
        base_features = ['length', 'pfe', 'fsr', 'fsd', 'coverage', 'ends', 'ocf', 'ifs', 'wps', 'edm', 'poem', 'prem', 'iedm', 'eedm']
        folders_to_create = list(base_features)
        if self.mapq_filter is not None:
            mapq_features = self.mapq_filter_features if self.mapq_filter_features is not None else base_features
            folders_to_create += [f'{f}_mapq{self.mapq_filter}' for f in mapq_features]
        for folder in folders_to_create:
            os.makedirs(f'./data/cristiano_features/{folder}', exist_ok=True)

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
        # try:
        #     self.wps = self.get_wps()
        # except Exception as e:
        #     print(f"Error calculating wps: {e}")
        #     self.wps = None
        try:
            self.edm = self.get_edm()
        except Exception as e:
            print(f"Error calculating edm: {e}")
            self.edm = None
        try:
            self.iedm = self.get_iedm()
        except Exception as e:
            print(f"Error calculating iedm: {e}")
            self.iedm = None
        try:
            self.eedm = self.get_eedm()
        except Exception as e:
            print(f"Error calculating eedm: {e}")
            self.eedm = None
        try:
            self.poem = self.get_poem()
        except Exception as e:
            print(f"Error calculating poem: {e}")
            self.poem = None
        try:
            self.prem = self.get_prem()
        except Exception as e:
            print(f"Error calculating prem: {e}")
            self.prem = None
        try:
            self.wps = self.get_wps()
        except Exception as e:
            print(f"Error calculating wps: {e}")
            self.wps = None


    def get_length(self):
        # check if length already exists in feature folder (for given sample)
        featurename = "length"
        if self._use_mapq_filter("length"):
            featurename = f"length_mapq{self.mapq_filter}"
        filepath = f"./data/cristiano_features/{featurename}/{self.sample_id}_{featurename}.csv"
        if self._use_saved_file(filepath, featurename):
            print(f"file already exists {self.sample_id}")
            df = pd.read_csv(filepath, index_col=0)
            return df
        print(f"calculating {featurename}")

        frags = self._filtered(self.frag_centroids_openchrom_intersect, "length")
        # create 33 bins 
        bin_edges = list(range(0, 310, 10)) + [np.inf]
        bin_labels = [f"{i}-{i+10}" for i in range(0, 300, 10)] + [">=300"]

        # get length column, cut into bins
        frags["len_bin"] = pd.cut(
            frags["length"], 
            bins=bin_edges, 
            labels=bin_labels, 
            right=False,
            include_lowest=True)
        
        all_bins = frags["len_bin"].cat.categories
        # group by chromosome and get proporotions
        length_matrix = (frags
                    .groupby(["f_chrom", "len_bin"])
                    .size()
                    .unstack(fill_value=0))
        length_matrix = length_matrix.reindex(columns=all_bins, fill_value=0)
    
        chrom_order = [f"chr{i}" for i in range(1, 23)]
        df_length = length_matrix.reindex(chrom_order, fill_value=0)
        print("df_length:", df_length.shape)
        df_length.to_csv(filepath)      # save to feature folder
        return df_length
    
    def get_pfe(self):
        featurename = "pfe"
        if self._use_mapq_filter("pfe"):
            featurename = f"pfe_mapq{self.mapq_filter}"
        filename = f"./data/cristiano_features/{featurename}/{self.sample_id}_{featurename}.csv"
        if self._use_saved_file(filename, featurename):
            print(f"file already exists {self.sample_id}")
            df = pd.read_csv(filename, index_col=0)
            return df
        print(f"calculating {featurename}")

        frags = self._filtered(self.frag_centroids_openchrom_intersect, "pfe")
        bins_pfe = ([0, 100] +
                    list(range(100, 260, 10)) +
                    [np.inf])

        bins_pfe = sorted(set(bins_pfe))
        bin_labels_pfe = [f"{bins_pfe[i]}-{bins_pfe[i+1]}" for i in range(len(bins_pfe)-1)]

        # cut into categories
        frags["pfe_bin"] = pd.cut(
            frags["length"],
                bins=bins_pfe,
                labels=bin_labels_pfe,
                right=False,
                include_lowest=True)

        # calculate P_i
        counts = (frags
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
        featurename = 'fsr'
        if self._use_mapq_filter("fsr"):
            featurename = f"fsr_mapq{self.mapq_filter}"
        filename = f"./data/cristiano_features/{featurename}/{self.sample_id}_{featurename}.csv"
        if self._use_saved_file(filename, featurename):
            print(f"file already exists {self.sample_id}")
            df = pd.read_csv(filename, index_col=0)
            return df

        print(f"calculating {featurename}")
        frags = self._filtered(self.frag_centroids_openchrom_intersect, "fsr")
        # bin edges
        bins_fsr = [64, 150, 220, 400]
        bin_labels_fsr = ["65-150", "151-220", "221-400"]
        frags["fsr_bin"] = pd.cut(
            frags["length"],
            bins=bins_fsr,
            labels=bin_labels_fsr,
            right=True,
            include_lowest=True
        )
        all_bins = frags["fsr_bin"].cat.categories
        # get proportions
        counts = (frags     
                .groupby(["region_id", "fsr_bin"])
                .size()
                .unstack(fill_value=0))
        counts = counts.reindex(columns=all_bins, fill_value=0)
        df_fsr = counts.div(counts.sum(axis=1), axis=0).fillna(0)

        # merge with region ids to fill the regions with no fragments
        df_fsr = df_fsr.merge(self.df_region_ids, on="region_id", how="right").set_index('region_id')
        df_fsr = df_fsr.fillna(0)
        print("df_fsr:", df_fsr.shape)
        df_fsr.to_csv(filename)
        return df_fsr
    
    

    def get_fsd(self):

        featurename = 'fsd'
        if self._use_mapq_filter("fsd"):
            featurename = f"fsd_mapq{self.mapq_filter}"
        filename = f"./data/cristiano_features/{featurename}/{self.sample_id}_{featurename}.csv"
        if self._use_saved_file(filename, featurename):
            print(f"file already exists {self.sample_id}")
            df = pd.read_csv(filename, index_col=0)
            return df

        print(f"calculating {featurename}")
        frags = self._filtered(self.frag_centroids_openchrom_intersect, "fsd")
        # bin 
        bins_fsd = list(range(65, 405, 5)) 
        bin_labels_fsd = [f"{bins_fsd[i]}-{bins_fsd[i+1]}" for i in range(len(bins_fsd)-1)]

        frags["fsd_bin"] = pd.cut(
            frags["length"],
            labels=bin_labels_fsd,
            bins=bins_fsd,
            right=False,
            include_lowest=True,
        )
        all_bins = frags["fsd_bin"].cat.categories

        # proportions
        counts = (
            frags
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

        featurename = 'coverage'
        if self._use_mapq_filter("coverage"):
            featurename = f"coverage_mapq{self.mapq_filter}"
        filename = f"./data/cristiano_features/{featurename}/{self.sample_id}_{featurename}.csv"
        if self._use_saved_file(filename, featurename):
            print(f"file already exists {self.sample_id}")
            df = pd.read_csv(filename, index_col=0)
            return df
        
        print(f"calculating {featurename}")

        frags = self._filtered(self.frag_centroids_openchrom_intersect, "coverage")
        # for coverage we need to merge with the original region id file at the start, since we are just counting occurences
        centroids_intersect = self.df_region_ids.merge(
            frags,
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

        featurename = 'ends'
        if self._use_mapq_filter("ends"):
            featurename = f"ends_mapq{self.mapq_filter}"
        filename = f"./data/cristiano_features/{featurename}/{self.sample_id}_{featurename}.csv"
        if self._use_saved_file(filename, featurename):
            print(f"file already exists {self.sample_id}")
            df = pd.read_csv(filename, index_col=0)
            return df

        print(f"calculating {featurename}")
        ends = self._filtered(self.frag_ends_openchrom_intersect, "ends")
        ends_intersect = self.df_region_ids.merge(
            ends,
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
        featurename = 'ocf'
        if self._use_mapq_filter("ocf"):
            featurename = f"ocf_mapq{self.mapq_filter}"
        filename = f"./data/cristiano_features/{featurename}/{self.sample_id}_{featurename}.csv"
        if self._use_saved_file(filename, featurename):
            print(f"file already exists {self.sample_id}")
            df = pd.read_csv(filename, index_col=0)
            return df
        print(f"calculating {featurename}")

        # set boundaries around peaks
        leftmin, leftmax = -70, -50
        rightmin, rightmax = 50, 70
        df = self._filtered(self.frag_ends_ocf, "ocf")

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

        # calculate D-U and U-D for both windows, use D-U for left window and U-D for right window
        counts["left_term"] = (counts.get("D", 0) - counts.get("U", 0))
        counts["right_term"] = (counts.get("U", 0) - counts.get("D", 0))

        counts['left_window'] = np.where(
            counts.index.get_level_values('window') == "left",
            counts["left_term"], 0)
        counts['right_window'] = np.where(
            counts.index.get_level_values('window') == "right",
            counts["right_term"], 0)

        ocf_left = counts.groupby("region_id")["left_window"].sum()
        ocf_right = counts.groupby("region_id")["right_window"].sum()

        ocf = ocf_left + ocf_right

        print(f"OCF describe:")
        print(ocf.describe())

        print(ocf.head(20))

        df_ocf = (self.df_region_ids
                .merge(ocf.rename("ocf"), on="region_id", how="left")
                .fillna(0))
        df_ocf.set_index('region_id', inplace=True)
        print("df_ocf:", df_ocf.shape)
        df_ocf.to_csv(filename)
        return df_ocf
    
    def get_ifs(self):
        featurename = 'ifs'
        if self._use_mapq_filter("ifs"):
            featurename = f"ifs_mapq{self.mapq_filter}"
        filename = f"./data/cristiano_features/{featurename}/{self.sample_id}_{featurename}.csv"
        if self._use_saved_file(filename, featurename):
            print(f"file already exists {self.sample_id}")
            df = pd.read_csv(filename, index_col=0)
            return df

        print(f"calculating {featurename}")
        frags = self._filtered(self.frag_centroids_openchrom_intersect, "ifs")
        # count fragments (n), calculate average lengths (l),
        region_counts = frags.groupby("region_id").size()
        region_avg_length = frags.groupby("region_id")["length"].mean()

        # average length per chromosome (L):
        chrom_avg_length = frags.groupby("f_chrom")["length"].mean()

        # add chrom avg to main dataframe
        frags["chrom_avg"] = frags["f_chrom"].map(chrom_avg_length)
        df_ifs = pd.DataFrame({
            "region_id": region_counts.index,
            "n": region_counts.values,
            "l": region_avg_length.values
        })

        # map chromosome for L
        region_chrom = frags.groupby("region_id")["f_chrom"].first()
        df_ifs["chrom"] = df_ifs["region_id"].map(region_chrom)
        df_ifs["L"] = df_ifs["chrom"].map(chrom_avg_length)

        # compute IFS and fill the missing regions
        df_ifs["IFS"] = df_ifs["n"] * (1 + df_ifs["l"] / df_ifs["L"])
        df_ifs = self.df_region_ids.merge(df_ifs[["region_id", "IFS"]], on="region_id", how="left").fillna(0)
        df_ifs.set_index('region_id', inplace=True)
        print('df_ifs:', df_ifs.shape)
        df_ifs.to_csv(filename)
        return df_ifs
    
    def get_wps_old(self):
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
        featurename = "wps"
        if self._use_mapq_filter("wps"):
            featurename = f"wps_mapq{self.mapq_filter}"
        filename = f"./data/cristiano_features/{featurename}/{self.sample_id}_{featurename}.csv"
        if self._use_saved_file(filename, featurename):
            print(f"file already exists {self.sample_id}")
            df = pd.read_csv(filename, index_col=0)
            return df
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

    def get_wps(self):
        """
        compute wps from fragment ends instead of bigwig file
        uses frag_wps_openchrom_intersect.
        """
        featurename = "wps"
        if self._use_mapq_filter("wps"):
            featurename = f"wps_mapq{self.mapq_filter}"
        filename = f"./data/cristiano_features/{featurename}/{self.sample_id}_{featurename}.csv"
        if self._use_saved_file(filename, featurename):
            print(f"file already exists {self.sample_id}")
            df = pd.read_csv(filename, index_col=0)
            return df
        
        print(f"calculating {featurename}")
        half_window = 60
        wps_means = []

        wps_frags = self._filtered(self.frag_wps_openchrom_intersect, "wps")
        filtered_frags = wps_frags[
            (wps_frags['length'] >= 120) &
            (wps_frags['length'] <= 180)
        ]
        frags_in_region = {}
        for region_id, group in filtered_frags.groupby('region_id'):
            starts = group['f_start'].values  
            ends = group['f_end'].values     
            frags_in_region[region_id] = (starts, ends)

        for idx, row in self.openchrom_with_id.iterrows():
            chrom = row['oc_chrom']
            start = row['oc_start']
            end = row['oc_end']
            region_id = row['region_id']

            if region_id not in frags_in_region:
                wps_means.append((region_id, 0.0))
                continue
            
            f_starts = frags_in_region[region_id][0]
            f_ends = frags_in_region[region_id][1]
            wps_scores_region = []

            for pos in range(start, end):
                window_start = pos - half_window
                window_end = pos + half_window

                starts_in_window = (f_starts >= window_start) & (f_starts < window_end)
                ends_in_window = (f_ends > window_start) & (f_ends <= window_end)
                fragments_ending = np.sum(starts_in_window | ends_in_window)
                fragments_spanning = np.sum((f_starts < window_start) & (f_ends > window_end))

                wps_scores_region.append(fragments_spanning - fragments_ending)
            
            mean_wps = np.mean(wps_scores_region) if wps_scores_region else 0.0
            wps_means.append((region_id, mean_wps))

        df_wps = pd.DataFrame(wps_means, columns=["region_id", "wps"])
        df_wps = self.df_region_ids.merge(df_wps, on="region_id", how="left").set_index('region_id')
        df_wps["wps"] = df_wps["wps"].fillna(0)
        
        print("df_wps:", df_wps.shape)
        df_wps.to_csv(filename)
        return df_wps
        
    def _reverse_complement(self, seq):
        """
        Takes non-reversed end sequence of the 3' downstream end as input.
        Returns reverse complement.
        """
        map = {"A": "T", "T": "A", 
               "C": "G", "G": "C"}
        res = ""
        for char in seq:
            res += map.get(char, 'N') 
        return res[::-1]
        

    def _get_motif(self, df, offset=0, k=4, rev_complement=False):
        """
        Extract a 4-mer end motif at window given the startin gcoordinate [pos+offset, pos+offset+k) 
        df: groupby dataframe
        rev_complement: reverse complement the extracted sequence

        Offset/rc reference for k=4:
          edm D:      offset=0,   rev_complement=False,  [pos,   pos+4)   
          edm U:      offset=-4,  rev_complement=True,   [pos-4, pos)    
          poem:       offset=+4,  rev_complement=False,  [pos+4, pos+8)   
          prem:       offset=-8,  rev_complement=True,   [pos-8, pos-4)   
          ext prem:   offset=-4,  rev_complement=False,   [pos-4, pos)     
          ext poem:   offset=0,   rev_complement=True,  [pos,   pos+4)   
        """
        chrom = df.name
        if chrom not in self.hg19_fasta:
            return pd.Series([None] * len(df), index=df.index)

        seq = self.hg19_fasta[chrom]
        motifs = []

        for pos in df['pos'].values:
            start = pos + offset
            end = start + k
            try:
                motif_seq = str(seq[start:end].seq).upper()
                if rev_complement:
                    motif_seq = self._reverse_complement(motif_seq)
                motifs.append(motif_seq)
            except Exception as e:
                print(f'Error: {chrom}, {start}, {end}, {e}')
                motifs.append(None)

        return pd.Series(motifs, index=df.index, name="motif")


    def get_edm(self):
        featurename = "edm"
        if self._use_mapq_filter("edm"):
            featurename = f"edm_mapq{self.mapq_filter}"
        filepath = f"./data/cristiano_features/{featurename}/{self.sample_id}_{featurename}.csv"
        if self._use_saved_file(filepath, featurename):
            print(f"file already exists {self.sample_id}")
            df = pd.read_csv(filepath, index_col=0)
            return df

        print(f"calculating {featurename}")
        frags = self._filtered(self.frag_centroids_openchrom_intersect, "edm")
        if len(frags) == 0:
            print(f"Warning - No fragments after MAPQ filtering for {featurename}")
            all_motifs = [''.join(p) for p in itertools.product('ACGT', repeat=4)]
            chrom_order = [f"chr{i}" for i in range(1, 23)]
            df_edm = pd.DataFrame(0, index=chrom_order, columns=all_motifs)
            df_edm.to_csv(filepath)
            return df_edm
        # get end positions
        ends_start = frags[["f_chrom", "f_start"]].copy()
        ends_start = ends_start.rename(columns={"f_start": "pos"})
        ends_end = frags[["f_chrom", "f_end"]].copy()
        ends_end = ends_end.rename(columns={"f_end": "pos"})

        # get motifs at each position
        edm_start = (
            ends_start
            .groupby("f_chrom")
            .apply(lambda df: self._get_motif(df, offset=0, rev_complement=False))
            .rename("motif")
            .reset_index()
            .drop(columns=["level_1"])
        )
        edm_end = (
            ends_end
            .groupby("f_chrom")
            .apply(lambda df: self._get_motif(df, offset=-4, rev_complement=True))
            .rename("motif")
            .reset_index()
            .drop(columns=["level_1"])
        )

        df_motif = pd.concat([edm_start, edm_end], ignore_index=True)
        motif_counts = (
            df_motif
            .groupby(["f_chrom", "motif"])
            .size()
            .unstack(fill_value=0)
        )
        print(f"motif_counts: {motif_counts.shape}")
        # have to reindex to include all motifs
        all_motifs = [''.join(p) for p in itertools.product('ACGT', repeat=4)]
        motif_counts = motif_counts.reindex(columns=all_motifs, fill_value=0)

        df_edm = motif_counts.div(motif_counts.sum(axis=1), axis=0)
        
        # Reorder to chr1-chr22
        chrom_order = [f"chr{i}" for i in range(1, 23)]
        df_edm = df_edm.reindex(chrom_order, fill_value=0)
        
        print(f"df_edm: {df_edm.shape}")
        df_edm.to_csv(filepath)
        return df_edm
    
    def get_poem(self):
        """
        get post-end 4mer motif, 4 bases after frag_start
        """
        featurename = 'poem'
        if self._use_mapq_filter("poem"):
            featurename = f"poem_mapq{self.mapq_filter}"
        filename = f"./data/cristiano_features/{featurename}/{self.sample_id}_{featurename}.csv"
        if self._use_saved_file(filename, featurename):
            print(f"file already exists {self.sample_id}")
            df = pd.read_csv(filename, index_col=0)
            return df

        print(f"calculating {featurename}")
        frags = self._filtered(self.frag_centroids_openchrom_intersect, "poem")
        if len(frags) == 0:
            print(f"Warning - No fragments after MAPQ filtering for {featurename}")
            all_motifs = [''.join(p) for p in itertools.product('ACGT', repeat=4)]
            chrom_order = [f"chr{i}" for i in range(1, 23)]
            df_poem = pd.DataFrame(0, index=chrom_order, columns=all_motifs)
            df_poem.to_csv(filename)
            return df_poem
        # get end positions
        ends_start = frags[["f_chrom", "f_start"]].copy()
        ends_start = ends_start.rename(columns={"f_start": "pos"})

        # get motifs 4 bp after the normal downstream edm
        poem_start = (
            ends_start
            .groupby("f_chrom")
            .apply(lambda df: self._get_motif(df, offset=4, rev_complement=False))
            .rename("motif")
            .reset_index()
            .drop(columns=["level_1"])
        )

        # count and convert to proportions
        motif_counts = (
            poem_start
            .groupby(["f_chrom", "motif"])
            .size()
            .unstack(fill_value=0)
        )
        print(f"motif_counts: {motif_counts.shape}")
        # have to reindex to include all motifs
        all_motifs = [''.join(p) for p in itertools.product('ACGT', repeat=4)]
        motif_counts = motif_counts.reindex(columns=all_motifs, fill_value=0)

        df_poem = motif_counts.div(motif_counts.sum(axis=1), axis=0)
        
        # Reorder to chr1-chr22
        chrom_order = [f"chr{i}" for i in range(1, 23)]
        df_poem = df_poem.reindex(chrom_order, fill_value=0)
        
        print(f"df_poem: {df_poem.shape}")
        df_poem.to_csv(filename)
        return df_poem

    def get_prem(self):
        """
        get pre-end 4mer motif, 4 bases before frag_end, reverse complemented
        """
        featurename = 'prem'
        if self._use_mapq_filter("prem"):
            featurename = f"prem_mapq{self.mapq_filter}"
        filename = f"./data/cristiano_features/{featurename}/{self.sample_id}_{featurename}.csv"
        if self._use_saved_file(filename, featurename):
            print(f"file already exists {self.sample_id}")
            df = pd.read_csv(filename, index_col=0)
            return df

        print(f"calculating {featurename}")
        frags = self._filtered(self.frag_centroids_openchrom_intersect, "prem")
        if len(frags) == 0:
            print(f"Warning - No fragments after MAPQ filtering for {featurename}")
            all_motifs = [''.join(p) for p in itertools.product('ACGT', repeat=4)]
            chrom_order = [f"chr{i}" for i in range(1, 23)]
            df_prem = pd.DataFrame(0, index=chrom_order, columns=all_motifs)
            df_prem.to_csv(filename)
            return df_prem
        # get end positions
        ends_end = frags[["f_chrom", "f_end"]].copy()
        ends_end = ends_end.rename(columns={"f_end": "pos"})

        # get motifs 4 bp before the upstream edm, reverse complemented
        prem_end = (
            ends_end
            .groupby("f_chrom")
            .apply(lambda df: self._get_motif(df, offset=-8, rev_complement=True))
            .rename("motif")
            .reset_index()
            .drop(columns=["level_1"])
        )

        # count and convert to proportions
        motif_counts = (
            prem_end
            .groupby(["f_chrom", "motif"])
            .size()
            .unstack(fill_value=0)
        )
        print(f"motif_counts: {motif_counts.shape}")

        # have to reindex to include all motifs
        all_motifs = [''.join(p) for p in itertools.product('ACGT', repeat=4)]
        motif_counts = motif_counts.reindex(columns=all_motifs, fill_value=0)

        df_prem = motif_counts.div(motif_counts.sum(axis=1), axis=0)

        # Reorder to chr1-chr22
        chrom_order = [f"chr{i}" for i in range(1, 23)]
        df_prem = df_prem.reindex(chrom_order, fill_value=0)

        print(f"df_prem: {df_prem.shape}")
        df_prem.to_csv(filename)
        return df_prem
    
    def get_iedm(self):
        """
        get combined poem and prem motif, 4 bases after frag_start and 4 bases before frag_end (reverse complemented)
        """
        featurename = 'iedm'
        if self._use_mapq_filter("iedm"):
            featurename = f"iedm_mapq{self.mapq_filter}"
        filename = f"./data/cristiano_features/{featurename}/{self.sample_id}_{featurename}.csv"
        if self._use_saved_file(filename, featurename):
            print(f"file already exists {self.sample_id}")
            df = pd.read_csv(filename, index_col=0)
            return df

        print(f"calculating {featurename}")
        frags = self._filtered(self.frag_centroids_openchrom_intersect, "iedm")
        if len(frags) == 0:
            print(f"Warning - No fragments after MAPQ filtering for {featurename}")
            all_motifs = [''.join(p) for p in itertools.product('ACGT', repeat=4)]
            chrom_order = [f"chr{i}" for i in range(1, 23)]
            df_poem = pd.DataFrame(0, index=chrom_order, columns=all_motifs)
            df_poem.to_csv(filename)
            return df_poem
        # get end positions
        ends_start = frags[["f_chrom", "f_start"]].copy()
        ends_start = ends_start.rename(columns={"f_start": "pos"})
        ends_end = frags[["f_chrom", "f_end"]].copy()
        ends_end = ends_end.rename(columns={"f_end": "pos"})

        # keep same compact style as get_edm: groupby-apply motifs then combine
        poem = ends_start.groupby("f_chrom").apply(
            lambda df: self._get_motif(df, offset=4, rev_complement=False)
        )
        prem = ends_end.groupby("f_chrom").apply(
            lambda df: self._get_motif(df, offset=-8, rev_complement=True)
        )

        # count and convert to proportions
        motif_counts = (
            pd.concat([poem, prem])
            .groupby(level=0)
            .value_counts()
            .unstack(fill_value=0)
        )
        print(f"motif_counts: {motif_counts.shape}")
        # have to reindex to include all motifs
        all_motifs = [''.join(p) for p in itertools.product('ACGT', repeat=4)]
        motif_counts = motif_counts.reindex(columns=all_motifs, fill_value=0)

        df_iedm = motif_counts.div(motif_counts.sum(axis=1), axis=0)
        
        # Reorder to chr1-chr22
        chrom_order = [f"chr{i}" for i in range(1, 23)]
        df_iedm = df_iedm.reindex(chrom_order, fill_value=0)
        
        print(f"df_iedm: {df_iedm.shape}")
        df_iedm.to_csv(filename)
        return df_iedm
    
    
    def get_eedm(self):
        """
        get combined external pre/post-end motif:
        - external pre-end: 4 bases before fragment start [f_start-4, f_start), reverse complemented
        - external post-end: 4 bases after fragment end [f_end, f_end+4), not reverse complemented
        """
        featurename = 'eedm'
        if self._use_mapq_filter("eedm"):
            featurename = f"eedm_mapq{self.mapq_filter}"
        filename = f"./data/cristiano_features/{featurename}/{self.sample_id}_{featurename}.csv"
        if self._use_saved_file(filename, featurename):
            print(f"file already exists {self.sample_id}")
            df = pd.read_csv(filename, index_col=0)
            return df

        print(f"calculating {featurename}")

        frags = self._filtered(self.frag_centroids_openchrom_intersect, "eedm")
        if len(frags) == 0:
            print(f"Warning - No fragments after MAPQ filtering for {featurename}")
            all_motifs = [''.join(p) for p in itertools.product('ACGT', repeat=4)]
            chrom_order = [f"chr{i}" for i in range(1, 23)]
            df_eedm = pd.DataFrame(0, index=chrom_order, columns=all_motifs)
            df_eedm.to_csv(filename)
            return df_eedm
        # get end positions
        ends_start = frags[["f_chrom", "f_start"]].copy()
        ends_start = ends_start.rename(columns={"f_start": "pos"})
        ends_end = frags[["f_chrom", "f_end"]].copy()
        ends_end = ends_end.rename(columns={"f_end": "pos"})

        ext_pre = ends_start.groupby("f_chrom").apply(
            lambda df: self._get_motif(df, offset=-4, rev_complement=False)
        )
        ext_post = ends_end.groupby("f_chrom").apply(
            lambda df: self._get_motif(df, offset=0, rev_complement=True)
        )

        # count and convert to proportions
        motif_counts = (
            pd.concat([ext_pre, ext_post])
            .groupby(level=0)
            .value_counts()
            .unstack(fill_value=0)
        )
        print(f"motif_counts: {motif_counts.shape}")
        # have to reindex to include all motifs
        all_motifs = [''.join(p) for p in itertools.product('ACGT', repeat=4)]
        motif_counts = motif_counts.reindex(columns=all_motifs, fill_value=0)

        df_eedm = motif_counts.div(motif_counts.sum(axis=1), axis=0)
        
        # Reorder to chr1-chr22
        chrom_order = [f"chr{i}" for i in range(1, 23)]
        df_eedm = df_eedm.reindex(chrom_order, fill_value=0)
        
        print(f"df_eedm: {df_eedm.shape}")
        df_eedm.to_csv(filename)
        return df_eedm
    

    
        


if __name__ == "__main__":
    if len(sys.argv) not in [6, 7, 8]:
        print("use: python sample_features_script.py sample_id frag_centroids_openchrom_intersect_path frag_ends_openchrom_intersect_path frag_ends_ocf_path frag_wps_intersect_path [mapq_filter] [rerun]")
        sys.exit(1)
    
    openchrom_path = './data/processing/openchrom_with_id.bed'
    sample_id = sys.argv[1]
    frag_centroids_openchrom_intersect_path = sys.argv[2]
    frag_ends_openchrom_intersect_path = sys.argv[3]
    frag_ends_ocf_path = sys.argv[4]
    frag_wps_intersect_path = sys.argv[5]
    hg19_fasta_path = './data/ref_genome/hg19/hg19.fa'
    
    mapq_filter = None
    if len(sys.argv) >= 7:
        mapq_arg = sys.argv[6].lower()
        if mapq_arg not in ['none', '']:
            mapq_filter = int(sys.argv[6])

    rerun = True
    if len(sys.argv) == 8:
        rerun_arg = sys.argv[7].strip().lower()
        if rerun_arg in ['true', '1', 'yes', 'y']:
            rerun = True
        elif rerun_arg in ['false', '0', 'no', 'n']:
            rerun = False
        else:
            print(f"Invalid rerun value: {sys.argv[7]}. Use true/false.")
            sys.exit(1)
    
    sample_features = SampleFeatures(
        sample_id,
        openchrom_path,
        frag_centroids_openchrom_intersect_path,
        frag_ends_openchrom_intersect_path,
        frag_ends_ocf_path,
        frag_wps_intersect_path,
        hg19_fasta_path,
        rerun=rerun,
        rerun_features=None,  
        mapq_filter=mapq_filter,
        mapq_filter_features=None  # Filter all features
    )
    sample_features.calculate_features()