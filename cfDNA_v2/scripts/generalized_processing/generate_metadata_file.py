import os
import pandas as pd
import sys

'''
Script to generate metadata for PreveLynch and GenoScan datasets, based on their filenames
'''
def generate_metadata(ds_name):
    rows = []
    for fname in os.listdir("frag/"):
        if not fname.endswith(".GRCh37.frag.bed.gz") or not fname.startswith(ds_name):
            continue
        sample_id = fname.replace(".GRCh37.frag.bed.gz", "")
        metadata = sample_id.split("_")
        dataset = metadata[0]  
        cancer_true = 1 if dataset[-2:] in ("cc", "ca") else 0
        disease = metadata[1]
        material = metadata[-1]
        if material != "pl":
            continue 
        rows.append({
            "sample_id": sample_id,
            "disease": disease,
            "dataset": dataset,
            "material": material,
            "frag_path": f"frag/{fname}",
            "stage": None,
            "cancer_true": cancer_true
        })

    df = pd.DataFrame(rows)
    out_path = f"data/metadata/internal_metadata_{ds_name}.csv"
    df.to_csv(out_path, index=False)
    return df

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usecase: python generate_metadata_file.py [dataset_prefix]")
        print("  e.g. gs, ly, gsca, lycc")
        sys.exit(1)
    ds_name = sys.argv[1]
    df = generate_metadata(ds_name)
    print(f"saved to data/metadata/internal_metadata_{ds_name}.csv")