import pandas as pd
import os

print("="*60)
print("HDF5 FEATURE MATRIX VALIDATION")
print("="*60)

matrix_dir = "./data/matrix/by_feature"

# Expected dimensions (samples, features + 6 metadata cols)
expected = {
    'length': (459, 688),      # 682 + 6
    'pfe': (459, 561420),      # 561414 + 6
    'fsr': (459, 1684248),     # 1684242 + 6
    'fsd': (459, 1480),        # 1474 + 6
    'coverage': (459, 561420),
    'ends': (459, 561420),
    'ocf': (459, 561420),
    'ifs': (459, 561420),
    'wps': (459, 561420),
    'edm': (459, 5638)         # 5632 + 6
}

print("\nValidating all matrices...")

all_ok = True
for feature in ['length', 'pfe', 'fsr', 'fsd', 'coverage', 'ends', 'ocf', 'ifs', 'wps', 'edm']:
    h5_path = f"{matrix_dir}/matrix_{feature}.h5"
    
    if not os.path.exists(h5_path):
        print(f"✗ {feature}: FILE NOT FOUND")
        all_ok = False
        continue
    
    # Load and check shape (FAST with HDF5!)
    df = pd.read_hdf(h5_path, key='data')
    size_mb = os.path.getsize(h5_path) / (1024**2)
    
    exp_shape = expected[feature]
    status = "✓" if df.shape == exp_shape else "⚠"
    
    print(f"{status} {feature:8s}: {df.shape} ({size_mb:6.1f} MB)")
    
    if df.shape != exp_shape:
        print(f"   Expected: {exp_shape}")
        all_ok = False

print("\n" + "="*60)
if all_ok:
    print("✓ ALL MATRICES VALIDATED!")
    print("Ready for Zhou et al. ensemble modeling!")
else:
    print("⚠ Some matrices need attention")
print("="*60)

# Quick load test
print("\nTesting load speed...")
import time

start = time.time()
df_length = pd.read_hdf(f"{matrix_dir}/matrix_length.h5", key='data')
elapsed = time.time() - start

print(f"✓ Loaded LENGTH matrix ({df_length.shape}) in {elapsed:.2f} seconds")
print("="*60)