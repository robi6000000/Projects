import pandas as pd

gc_df = pd.read_csv('./data/processing/gc_content_per_region.csv')

print(gc_df['gc_content'].describe())



# Check extremes
print("\nVery low GC regions (< 0.3):")
print(len(gc_df[gc_df['gc_content'] < 0.3]))

print("Very high GC regions (> 0.7):")
print(len(gc_df[gc_df['gc_content'] > 0.7]))