import pandas as pd

gc_df = pd.read_csv('./data/processing/gc_content_per_region.csv')

# Look at distribution
print(gc_df['gc_content'].describe())

# Plot histogram
import matplotlib.pyplot as plt
plt.hist(gc_df['gc_content'], bins=50)
plt.xlabel('GC content')
plt.ylabel('Number of regions')
plt.title('GC content distribution across open chromatin regions')
plt.savefig('gc_distribution.png')
plt.show()

# Check extremes
print("\nVery low GC regions (< 0.3):")
print(len(gc_df[gc_df['gc_content'] < 0.3]))

print("Very high GC regions (> 0.7):")
print(len(gc_df[gc_df['gc_content'] > 0.7]))