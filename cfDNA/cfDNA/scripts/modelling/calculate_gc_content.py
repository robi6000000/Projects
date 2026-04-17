import pandas as pd
from pyfaidx import Fasta
import numpy as np

print("Calculating GC content per region")

# Load regions and reference
regions = pd.read_csv(
    './data/processing/openchrom_with_id.bed',
    sep='\t',
    header=None,
    names=['chrom', 'start', 'end', 'region_id']
)
genome = Fasta('./data/hg19/hg19.fa')
print(f"Regions: {len(regions)}")


def gc_content(sequence):
    """Calculate GC content percentage of a DNA sequence."""
    sequence = sequence.upper()
    gc_count = sequence.count("G") + sequence.count("C")
    total_count = len(sequence)
    if total_count == 0:
        return 0.0
    return (gc_count / total_count) * 100

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
        
        # Calculate GC percentage
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
# check nan values
nan_count = gc_df['gc_content'].isna().sum()
if nan_count > 0:
    print(f"{nan_count} regions with missing GC content")

file_path = './data/processing/gc_content_per_region.csv'
gc_df.to_csv(file_path, index=False)

print(f"GC content file finished - {file_path}")

