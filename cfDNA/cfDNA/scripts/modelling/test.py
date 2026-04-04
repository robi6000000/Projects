import pandas as pd
sample = pd.read_csv('./data/sample_temp/EE87788_frag_centroids_openchrom_intersect.bed', 
                     sep='\t', header=None)
sample['length'] = sample[4] - sample[3]

print("Fragment counts by length:")
print(sample['length'].describe())
print(f"\nFragments exactly 400bp: {(sample['length'] == 400).sum()}")
print(f"Fragments in [65,150]: {((sample['length'] >= 65) & (sample['length'] <= 150)).sum()}")
print(f"Fragments in [151,220]: {((sample['length'] >= 151) & (sample['length'] <= 220)).sum()}")
print(f"Fragments in [221,399]: {((sample['length'] >= 221) & (sample['length'] < 400)).sum()}")
print(f"Fragments exactly 400: {(sample['length'] == 400).sum()}")