import pandas as pd
import glob
import os

def fix_fsd_headers():
    """Fix FSD headers from '[65, 70)' to '65-70'"""
    print("Fixing FSD headers...")
    files = glob.glob("./data/cristiano_features/fsd/*.csv")
    
    for i, file_path in enumerate(files):
        if i % 50 == 0:
            print(f"  Processing {i+1}/{len(files)}...")
        
        df = pd.read_csv(file_path, index_col=0)
        
        # Rename columns: remove brackets, parentheses, replace ', ' with '-'
        new_columns = {}
        for col in df.columns:
            clean = str(col).replace('[', '').replace(']', '').replace('(', '').replace(')', '').replace(', ', '-').replace(' ', '')
            new_columns[col] = clean
        
        df.rename(columns=new_columns, inplace=True)
        df.to_csv(file_path)
    
    print(f"  Fixed {len(files)} FSD files")

def fix_fsr_headers():
    """Fix FSR headers from '65_151' to '65-151'"""
    print("Fixing FSR headers...")
    files = glob.glob("./data/cristiano_features/fsr/*.csv")
    
    for i, file_path in enumerate(files):
        if i % 50 == 0:
            print(f"  Processing {i+1}/{len(files)}...")
        
        df = pd.read_csv(file_path, index_col=0)
        
        # Check if needs fixing (has underscores)
        if any('_' in str(col) and '-' not in str(col) for col in df.columns):
            new_columns = {col: str(col).replace('_', '-') for col in df.columns}
            df.rename(columns=new_columns, inplace=True)
            df.to_csv(file_path)
    
    print(f"  Fixed {len(files)} FSR files")

if __name__ == "__main__":
    print("Starting header renaming...")
    
    fix_fsd_headers()
    fix_fsr_headers()
    
    print("\nDone! All headers now use dash format (e.g., 65-70)")
    
    # Verify with a sample
    print("\nVerifying FSD sample:")
    sample_fsd = pd.read_csv(glob.glob("./data/cristiano_features/fsd/*.csv")[0], index_col=0)
    print(f"  Columns: {list(sample_fsd.columns)[:5]}")
    
    print("\nVerifying FSR sample:")
    sample_fsr = pd.read_csv(glob.glob("./data/cristiano_features/fsr/*.csv")[0], index_col=0)
    print(f"  Columns: {list(sample_fsr.columns)}")