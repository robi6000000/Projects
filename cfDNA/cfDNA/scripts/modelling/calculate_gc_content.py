import pandas as pd
from pyfaidx import Fasta
import numpy as np

print("Calculating GC content per region")

# Load regions
print("Loading open chromatin regions...")
regions = pd.read_csv(
    './data/processing/openchrom_with_id.bed',
    sep='\t',
    header=None,
    names=['chrom', 'start', 'end', 'region_id']
)
print(f"Loaded {len(regions)} regions")

# Load reference genome
print("Loading reference genome...")
genome = Fasta('./data/hg19/hg19.fa')
print("Genome loaded")

# Calculate GC% for each region
print("Calculating GC content...")
gc_content = []

for idx, row in regions.iterrows():
    if idx % 50000 == 0:
        print(f"Processing region {idx+1}/{len(regions)}...")
    
    chrom = row['chrom']
    start = row['start']
    end = row['end']
    region_id = row['region_id']
    
    try:
        # Get sequence
        seq = str(genome[chrom][start:end].seq).upper()
        
        # Calculate GC percentage
        gc_count = seq.count('G') + seq.count('C')
        total = len(seq)
        gc_percent = gc_count / total if total > 0 else 0
        
        gc_content.append({
            'region_id': region_id,
            'gc_content': gc_percent,
            'length': total
        })
    
    except Exception as e:
        # Handle edge cases (invalid chromosomes, etc.)
        gc_content.append({
            'region_id': region_id,
            'gc_content': np.nan,
            'length': 0
        })

# Save results
print("Saving results...")
gc_df = pd.DataFrame(gc_content)

# Check for missing values
nan_count = gc_df['gc_content'].isna().sum()
if nan_count > 0:
    print(f"Warning: {nan_count} regions with missing GC content")

output_path = './data/processing/gc_content_per_region.csv'
gc_df.to_csv(output_path, index=False)

print(f"Saved GC content for {len(gc_df)} regions to {output_path}")

# Summary statistics
print("\nSummary statistics:")
print(f"Mean GC: {gc_df['gc_content'].mean():.3f}")
print(f"Median GC: {gc_df['gc_content'].median():.3f}")
print(f"Std GC: {gc_df['gc_content'].std():.3f}")
print(f"Range: {gc_df['gc_content'].min():.3f} - {gc_df['gc_content'].max():.3f}")