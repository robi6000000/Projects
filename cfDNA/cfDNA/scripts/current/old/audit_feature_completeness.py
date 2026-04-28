#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def load_manifest_samples(manifest_path: Path):
    sample_to_index = {}
    with manifest_path.open(newline="") as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader, start=1):
            seqrun_id = (row.get("seqrun_id") or "").strip()
            if not seqrun_id:
                continue
            sample_id = f"EE{seqrun_id}"
            sample_to_index[sample_id] = i
    return sample_to_index


def samples_in_feature_folder(folder: Path, feature_name: str):
    suffix = f"_{feature_name}.csv"
    found = set()
    for fp in folder.glob("*.csv"):
        name = fp.name
        if name.endswith(suffix):
            found.add(name[: -len(suffix)])
    return found


def main():
    parser = argparse.ArgumentParser(
        description="Audit feature folder completeness and prepare missing mapq30 backfill lists."
    )
    parser.add_argument(
        "--features-dir",
        default="data/cristiano_features",
        help="Base folder that contains feature subfolders.",
    )
    parser.add_argument(
        "--manifest",
        default="data/manifest/Cristiano_manifest.csv",
        help="Manifest CSV with seqrun_id column.",
    )
    parser.add_argument(
        "--reference-feature",
        default="coverage",
        help="Feature folder used as the canonical processed sample set.",
    )
    parser.add_argument(
        "--focus-suffix",
        default="_mapq30",
        help="Only union missing samples from features ending with this suffix.",
    )
    parser.add_argument(
        "--report",
        default="log/cristiano_feature_missing_vs_reference.tsv",
        help="Output TSV report path.",
    )
    parser.add_argument(
        "--missing-ids-out",
        default="log/missing_mapq30_ids.txt",
        help="Output newline-delimited sample IDs missing in any focused feature.",
    )
    parser.add_argument(
        "--array-indices-out",
        default="log/missing_mapq30_array_indices.txt",
        help="Output newline-delimited SLURM array indices for missing sample IDs.",
    )
    args = parser.parse_args()

    features_dir = Path(args.features_dir)
    manifest_path = Path(args.manifest)
    report_path = Path(args.report)
    missing_ids_path = Path(args.missing_ids_out)
    array_indices_path = Path(args.array_indices_out)

    if not features_dir.exists():
        raise FileNotFoundError(f"Features folder not found: {features_dir}")
    if not manifest_path.exists():
        raise FileNotFoundError(f"Manifest not found: {manifest_path}")

    sample_to_index = load_manifest_samples(manifest_path)
    manifest_samples = set(sample_to_index.keys())

    reference_folder = features_dir / args.reference_feature
    if not reference_folder.exists():
        raise FileNotFoundError(f"Reference feature folder not found: {reference_folder}")

    reference_samples = samples_in_feature_folder(reference_folder, args.reference_feature)
    if not reference_samples:
        raise RuntimeError(
            f"Reference feature folder {reference_folder} has no files matching '*_{args.reference_feature}.csv'"
        )

    report_rows = []
    focused_missing_union = set()

    for folder in sorted([p for p in features_dir.iterdir() if p.is_dir()]):
        feature = folder.name
        found = samples_in_feature_folder(folder, feature)
        missing_vs_ref = sorted(reference_samples - found)
        extra_vs_ref = sorted(found - reference_samples)

        report_rows.append(
            {
                "feature": feature,
                "found": len(found),
                "missing_vs_reference": len(missing_vs_ref),
                "extra_vs_reference": len(extra_vs_ref),
                "missing_ids": ",".join(missing_vs_ref),
                "extra_ids": ",".join(extra_vs_ref),
            }
        )

        if feature.endswith(args.focus_suffix) and missing_vs_ref:
            focused_missing_union.update(missing_vs_ref)

    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            delimiter="\t",
            fieldnames=[
                "feature",
                "found",
                "missing_vs_reference",
                "extra_vs_reference",
                "missing_ids",
                "extra_ids",
            ],
        )
        writer.writeheader()
        writer.writerows(report_rows)

    missing_ids = sorted(focused_missing_union)
    missing_ids_path.parent.mkdir(parents=True, exist_ok=True)
    missing_ids_path.write_text("\n".join(missing_ids) + ("\n" if missing_ids else ""))

    array_indices = [sample_to_index[sid] for sid in missing_ids if sid in sample_to_index]
    array_indices_path.parent.mkdir(parents=True, exist_ok=True)
    array_indices_path.write_text(
        "\n".join(str(i) for i in array_indices) + ("\n" if array_indices else "")
    )

    print(f"manifest samples: {len(manifest_samples)}")
    print(f"reference set ({args.reference_feature}) samples: {len(reference_samples)}")
    print(f"focused missing union ({args.focus_suffix}): {len(missing_ids)}")
    print(f"report: {report_path}")
    print(f"missing ids: {missing_ids_path}")
    print(f"array indices: {array_indices_path}")


if __name__ == "__main__":
    main()
