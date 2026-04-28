gc_per_region = []
gc_df = pd.DataFrame(gc_per_region)
os.makedirs(os.path.dirname(args.output), exist_ok=True)
gc_df.to_csv(args.output, index=False)

import pandas as pd
from pyfaidx import Fasta
import numpy as np
import sys
import os

def gc_content(sequence):
    """Calculate GC content percentage of a DNA sequence."""
    sequence = sequence.upper()
    gc_count = sequence.count("G") + sequence.count("C")
    total_count = len(sequence)
    if total_count == 0:
        return 0.0
    return (gc_count / total_count) * 100

def load_regions(path):
    return pd.read_csv(
        path,
        sep='\t',
        header=None,
        names=['chrom', 'start', 'end', 'region_id']
    )

if __name__ == "__main__":
    # Usage: python calculate_gc_content_general.py <regions_bed> <output_csv> <genome_fasta>
    # Defaults previously used:
    #   regions_bed: ./data/processing/openchrom_with_id.bed
    #   output_csv: ./data/processing/gc_content_per_region.csv
    #   genome_fasta: ./data/hg19/hg19.fa

    if len(sys.argv) < 4:
        print("Usage: python calculate_gc_content_general.py <regions_bed> <output_csv> <genome_fasta>")
        print("Example: python calculate_gc_content_general.py ./data/processing/openchrom_with_id.bed ./data/processing/gc_content_per_region.csv ./data/hg19/hg19.fa")
        sys.exit(1)

    regions_path = sys.argv[1]
    output_path = sys.argv[2]
    genome_path = sys.argv[3]

    print(f"Calculating GC content per region\nRegions file: {regions_path}\nOutput file: {output_path}\nGenome fasta: {genome_path}")

    regions = load_regions(regions_path)
    genome = Fasta(genome_path)
    print(f"Regions: {len(regions)}")

    gc_per_region = []
    for i, row in regions.iterrows():
        if i % 50000 == 0:
            print(f"processing {i}")
        chrom = row['chrom']
        start = row['start']
        end = row['end']
        region_id = row['region_id']
        try:
            seq = str(genome[chrom][start:end].seq).upper()
            gc_percent = gc_content(seq)
            length = len(seq)
            gc_per_region.append({
                'region_id': region_id,
                'gc_content': gc_percent,
                'length': length
            })
        except Exception as e:
            print(f"Error on {region_id}: {e}")
            gc_per_region.append({
                'region_id': region_id,
                'gc_content': np.nan,
                'length': 0
            })

    gc_df = pd.DataFrame(gc_per_region)
    nan_count = gc_df['gc_content'].isna().sum()
    if nan_count > 0:
        print(f"{nan_count} regions with missing GC content")

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    gc_df.to_csv(output_path, index=False)
    print(f"GC content file finished - {output_path}")
