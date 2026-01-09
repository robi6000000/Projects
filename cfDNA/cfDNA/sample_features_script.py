import pandas as pd
import sys
import subprocess
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pyfaidx import Fasta
import sys

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

# take sample id as input
# take './data/processing/openchrom_with_id.bed' (openchrom with ids) as input
# openchrom_path = './data/processing/openchrom_with_id.bed'
# frag_centroids_openchrom_intersect_path = './data/processing/frag_centroids_openchrom_intersect.bed'
# frag_ends_openchrom_intersect_path = './data/processing/frag_ends_openchrom_intersect.bed'
# frag_ends_ocf_path = './data/processing/frag_ends_ocf.bed'


def main():
    if len(sys.argv) != 5:
        print("Usage: python calc_features.py <SAMPLE_ID> <FRAG_CENTROIDS_OPENCHROM_INTERSECT_PATH> <FRAG_ENDS_OPENCHROM_INTERSECT_PATH> <FRAG_ENDS_OCF_PATH>")
        sys.exit(1)

    sample_id = sys.argv[1]
    frag_centroids_openchrom_intersect_path = sys.argv[2]
    frag_ends_openchrom_intersect_path = sys.argv[3]
    frag_ends_ocf_path = sys.argv[4]
    openchrom_path = './data/processing/openchrom_with_id.bed'
    print(f"Using sample: {sample_id}")


    openchrom_with_id = pd.read_csv(
        openchrom_path,
        sep="\t",
        header=None,
        names=["oc_chrom", "oc_start", "oc_end", "region_id"])
    df_region_ids = openchrom_with_id[['region_id']].drop_duplicates().reset_index(drop=True)

    frag_centroids_openchrom_intersect = pd.read_csv(
        frag_centroids_openchrom_intersect_path,
        sep="\t",
        header=None,
        names=["f_chrom", "centroid1", "centroid2", "f_start", "f_end", "score", "strand", "oc_chrom", "oc_start", "oc_end", "region_id"])

    frag_ends_openchrom_intersect = pd.read_csv(
        frag_ends_openchrom_intersect_path,
        sep="\t",
        header=None,
        names=["f_chrom", "end1", "end2", "end_type", "oc_chrom", "oc_start", "oc_end", "region_id"])

    print("Number of unique region IDs:", df_region_ids.shape[0])
    print("frag_centroids_openchrom_intersect.shape:", frag_centroids_openchrom_intersect.shape)
    print("frag_ends_openchrom_intersect.shape:", frag_ends_openchrom_intersect.shape)


    ### LENGTH
    bin_edges = list(range(0, 310, 10)) + [np.inf] 
    bin_labels = [f"{i}-{i+10}" for i in range(0, 300, 10)] + [">=300"]

    # fragment lengths
    frag_centroids_openchrom_intersect["length"] = frag_centroids_openchrom_intersect["f_end"] - frag_centroids_openchrom_intersect["f_start"]
    frag_centroids_openchrom_intersect["len_bin"] = pd.cut(
        frag_centroids_openchrom_intersect["length"], 
        bins=bin_edges, 
        labels=bin_labels, 
        right=False,
        include_lowest=True)

    # group by chromosome and bin and get counts
    length_matrix = (frag_centroids_openchrom_intersect
                .groupby(["f_chrom", "len_bin"], observed=True)
                .size()
                .unstack(fill_value=0))

    # convert to proportions per chromosome
    chrom_order = [f"chr{i}" for i in range(1, 23)]
    df_length = length_matrix
    df_length = df_length.reindex(chrom_order)
    print("df_length:", df_length.shape)

        
    ### **PFE**
    bins_pfe = ([0, 100] +
                list(range(100, 260, 10)) +
                [np.inf])

    bins_pfe = sorted(set(bins_pfe))
    bin_labels_pfe = [f"{bins_pfe[i]}-{bins_pfe[i+1]}" for i in range(len(bins_pfe)-1)]

    # create length column for later binning
    frag_centroids_openchrom_intersect["length"] = (
        frag_centroids_openchrom_intersect["f_end"] -
        frag_centroids_openchrom_intersect["f_start"])

    # cut into categories
    frag_centroids_openchrom_intersect["pfe_bin"] = pd.cut(
            frag_centroids_openchrom_intersect["length"],
            bins=bins_pfe,
            labels=bin_labels_pfe,
            right=False,
            include_lowest=True)

    # calculate P_i
    counts = (frag_centroids_openchrom_intersect
            .groupby(["region_id", "pfe_bin"], observed=True)
            .size()
            .unstack(fill_value=0))

    pfe_proportions = counts.div(counts.sum(axis=1), axis=0)
    pfe_proportions.replace(0, np.nan, inplace=True)

    # get actual PFE per refion
    pfe = - (pfe_proportions * np.log2(pfe_proportions)).sum(axis=1)
    df_pfe = pfe.to_frame("pfe")

    # merge pfe result with previously stored region ids to fill blanks
    df_pfe = df_region_ids.merge(df_pfe, on="region_id", how="left").set_index('region_id')
    df_pfe["pfe"] = df_pfe["pfe"].fillna(0)
    print("df_pfe:", df_pfe.shape)



    ### **FSR**
    bins_fsr = [65, 151, 221, 400]
    bin_labels_fsr = [f"{bins_fsr[i]}-{bins_fsr[i+1]}" for i in range(len(bins_fsr)-1)]

    frag_centroids_openchrom_intersect["fsr_bin"] = pd.cut(
        frag_centroids_openchrom_intersect["length"],
        bins=bins_fsr,
        labels=bin_labels_fsr,
        right=False,
        include_lowest=True
    )
    frag_centroids_openchrom_intersect_grouped = frag_centroids_openchrom_intersect.groupby("region_id")

    counts = (frag_centroids_openchrom_intersect
            .groupby(["region_id", "fsr_bin"], observed=True)
            .size()
            .unstack(fill_value=0))
    df_fsr = counts.div(counts.sum(axis=1), axis=0)
    df_fsr = df_fsr.merge(df_region_ids, on="region_id", how="right").set_index('region_id')
    print("df_fsr:", df_fsr.shape)



    ### **FSD**
    bins_fsd = list(range(65, 405, 5)) 
    frag_centroids_openchrom_intersect["fsd_bin"] = pd.cut(
        frag_centroids_openchrom_intersect["length"],
        bins=bins_fsd,
        right=False,
        include_lowest=True
    )
    counts = (
        frag_centroids_openchrom_intersect
        .groupby(["f_chrom", "fsd_bin"], observed=True)
        .size()
        .unstack(fill_value=0)
    )

    df_fsd = counts.div(counts.sum(axis=1), axis=0)
    df_fsd = df_fsd.reindex(chrom_order)
    print("df_fsd:", df_fsd.shape)



    ### **coverage**
    # for coverage we need to merge with the original region id file at the start, since we are just counting occurences
    centroids_intersect = df_region_ids.merge(
        frag_centroids_openchrom_intersect,
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
    print("region with max coverage:", df_cov['coverage'].idxmax())



    ### END
    ends_intersect = df_region_ids.merge(
        frag_ends_openchrom_intersect,
        on="region_id",
        how="left"
    )
    counts = (
        ends_intersect
        .groupby("region_id")
        .size()
    )
    df_end = pd.DataFrame(counts).rename(columns={0: "end"})
    print("df_end:", df_end.shape)
    print("region with max ends:", df_end['end'].idxmax())



    ###  OCF
    frag_ends_ocf = pd.read_csv(
        frag_ends_ocf_path,
        sep="\t",
        header=None,
        names=["chrom", "end1", "end2", "end_type", "oc_start", "oc_end", "region_id", "centroid", "rel_pos"]
    )
    leftmin, leftmax = -70, -50
    rightmin, rightmax = 50, 70
    df = frag_ends_ocf
    df["window"] = np.where((df["rel_pos"] >= leftmin) & (df["rel_pos"] < leftmax),
                            "left",
                            np.where((df["rel_pos"] >= rightmin) & (df["rel_pos"] < rightmax),
                                    "right",
                                    None))

    df = df[df["window"].notna()]
    counts = (df
            .groupby(["region_id", "window", "end_type"], observed=True)
            .size()
            .unstack(fill_value=0))

    df["window"].value_counts()
    counts["left_term"] = (counts["D"] - counts["U"])
    counts["right_term"] = (counts["U"] - counts["D"])

    ocf = (counts
            .groupby("region_id")[["left_term", "right_term"]]
            .sum()
            .sum(axis=1))
    df_ocf = (df_region_ids
            .merge(ocf.rename("ocf"), on="region_id", how="left")
            .fillna(0))
    df_ocf.set_index('region_id', inplace=True)
    print("df_ocf:", df_ocf.shape)



    ### IFS
    # count fragments (n), calculate average lengths (l),
    region_counts = frag_centroids_openchrom_intersect.groupby("region_id").size()
    region_avg_length = frag_centroids_openchrom_intersect.groupby("region_id")["length"].mean()

    # average length per chromosome (L):
    chrom_avg_length = frag_centroids_openchrom_intersect.groupby("f_chrom")["length"].mean()

    # add chrom avg to main df
    frag_centroids_openchrom_intersect["chrom_avg"] = frag_centroids_openchrom_intersect["f_chrom"].map(chrom_avg_length)
    df_ifs = pd.DataFrame({
        "region_id": region_counts.index,
        "n": region_counts.values,
        "l": region_avg_length.values
    })

    # map chromosome for L
    region_chrom = frag_centroids_openchrom_intersect.groupby("region_id")["f_chrom"].first()
    df_ifs["chrom"] = df_ifs["region_id"].map(region_chrom)
    df_ifs["L"] = df_ifs["chrom"].map(chrom_avg_length)

    # compute IFS and fill the missing regions
    df_ifs["IFS"] = df_ifs["n"] * (1 + df_ifs["l"] / df_ifs["L"])
    df_ifs = df_region_ids.merge(df_ifs[["region_id", "IFS"]], on="region_id", how="left").fillna(0)
    df_ifs.set_index('region_id', inplace=True)
    print('df_ifs:', df_ifs.shape)
    # df_ifs.plot(kind='line', title='Distribution of IFS values')



    ### WPS
    #TODO


    ### EDM
    hg19 = Fasta('./data/source_data/hg19.fa')
    # 5'ends
    ends_5 = frag_centroids_openchrom_intersect[['f_chrom','f_start','region_id']].copy()
    ends_5['pos'] = ends_5['f_start']
    ends_5['end_type'] = '5'

    # 3' ends
    ends_3 = frag_centroids_openchrom_intersect[['f_chrom','f_end','region_id']].copy()
    ends_3['pos'] = ends_3['f_end'] - 4
    ends_3['end_type'] = '3'
    ends = pd.concat([ends_5, ends_3], ignore_index=True)

    def get_motif(df, genome, k=4):
        chrom = df.name
        seq = genome[chrom]
        motifs = [
            str(seq[pos:pos+k].seq).upper()
            if pos >= 0 else None
            for pos in df['pos'].values
        ]
        return pd.Series(motifs, index=df.index)

    # filter out non-acgt letters (N)
    ends['motif'] = (
        ends
        .groupby('f_chrom', group_keys=False)
        .apply(get_motif, genome=hg19, k=4, include_groups=False)
    )
    ends.loc[~ends['motif'].str.match('^[ACGT]{4}$'), 'motif'].head()
    ends = ends[ends['motif'].str.match('^[ACGT]{4}$')]

    # now we calculate proportions of each motif at each end per chromosome - 256x22
    motif_counts = (
        ends
        .groupby(['f_chrom', 'motif'], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    df_motif = motif_counts.div(motif_counts.sum(axis=1), axis=0)
    df_motif = df_motif.reindex(chrom_order)
    print("df_motif.shape:", df_motif.shape)
    print("motif with max proportion:", df_motif.max().idxmax())


if __name__ == "__main__":
    main()
