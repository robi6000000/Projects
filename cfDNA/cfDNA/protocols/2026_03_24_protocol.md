# Feature matrix description

## length:
Dimensions - samples (459) x (metadata_cols (5) + length_bins (31) x autosomes (22))
- feature column names: length_chr1_bin_0-10

## pfe:
Dimensions - samples (459) x (metadata_cols (5) + openchrom_regions (561 414) x pfe_value (1))
- feature column names: pfe_regionid

## fsd:
Dimensions - samples (459) x (metadata_cols (5) + length_bins (67) x autosomes (22))
- feature column names: fsd_chr1_bin_65-70

# fsr:
Dimensions - samples (459) x (metadata_cols (5) + openchrom_regions (561 414) x fsr_bin (3))
- feature column names: fsr_regionid_bin_60-120

# coverage:
Dimensions - samples (459) x (metadata_cols (5) + openchrom_regions (561 414) x coverage (1))
- feature column names: coverage_regionid

# ends:
Dimensions - samples (459) x (metadata_cols (5) + openchrom_regions (561 414) x end_count (1))
- feature column names: ends_regionid

# ifs:
Dimensions - samples (459) x (metadata_cols (5) + openchrom_regions (561 414) x ifs_value (1))
- feature column names: ifs_regionid

# ocf:
Dimensions - samples (459) x (metadata_cols (5) + openchrom_regions (561 414) x ocf_value (1))
- feature column names: ocf_regionid

# wps:
Dimensions - samples (459) x (metadata_cols (5) + openchrom_regions (561 414) x wps_value (1))
- feature column names: wps_regionid

# edm:
Dimensions - samples (459) x (metadata_cols (5) + endmotif (256) x autosomes (22))
- feature column names: edm_chr1_bin_AAAA
