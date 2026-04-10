import os
import pandas as pd

rows = []
for fname in os.listdir("frag/"):
    if not fname.endswith(".GRCh37.frag.bed.gz"):
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
df.to_csv("data/manifest/internal_metadata.csv", index=False)