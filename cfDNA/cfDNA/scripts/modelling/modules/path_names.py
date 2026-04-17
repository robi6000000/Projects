from pathlib import Path


features = (
    "eedm",
    "iedm",
    "length",
    "edm",
    "fsd",
    "fsr",
    "pfe",
    "coverage",
    "ends",
    "ocf",
    "ifs",
    "wps"
)


def normalize_mapq_filter(mapq_filter):
    if mapq_filter is None:
        return None
    if str(mapq_filter).lower() == "none":
        return None
    return str(mapq_filter)


def build_run_prefix(
    kernel,
    pca=False,
    pca_components=0.95,
    gc_correction=False,
    cv_repeats=10,
    mapq_filter=None,
):
    parts = [f"svm_{kernel}"]
    if pca:
        parts.append(f"pca{pca_components}")
    if gc_correction:
        parts.append("gc")
    if cv_repeats != 10:
        parts.append(f"cv{cv_repeats}")

    mapq_filter = normalize_mapq_filter(mapq_filter)
    if mapq_filter is not None:
        parts.append(f"mapq{mapq_filter}")

    return "_".join(parts)


def build_probs_filename(feature, **kwargs):
    prefix = build_run_prefix(**kwargs)
    return f"{prefix}_{feature}_probs.csv"


def build_probs_path(base_dir, feature, **kwargs):
    prefix = build_run_prefix(**kwargs)
    return Path(base_dir) / prefix / build_probs_filename(feature, **kwargs)


def split_probs_filename(file_name, known_features=None):
    feature = known_features or features
    if not file_name.endswith("_probs.csv"):
        raise ValueError(f"Not a probability output file: {file_name}")

    stem = file_name[: -len("_probs.csv")]
    for feature in sorted(known_features, key=len, reverse=True):
        suffix = f"_{feature}"
        if stem.endswith(suffix):
            return stem[: -len(suffix)], feature

    raise ValueError(f"Could not parse feature from file name: {file_name}")


def parse_run_prefix(prefix):
    if not prefix.startswith("svm_"):
        raise ValueError(f"Invalid run prefix: {prefix}")

    parts = prefix.split("_")
    if len(parts) < 2:
        raise ValueError(f"Invalid run prefix: {prefix}")

    parsed = {
        "kernel": parts[1],
        "pca": False,
        "pca_components": None,
        "gc_correction": False,
        "cv_repeats": 10,
        "mapq_filter": None,
    }

    for token in parts[2:]:
        if token == "gc":
            parsed["gc_correction"] = True
        elif token.startswith("pca"):
            parsed["pca"] = True
            parsed["pca_components"] = float(token[3:])
        elif token.startswith("cv"):
            parsed["cv_repeats"] = int(token[2:])
        elif token.startswith("mapq"):
            parsed["mapq_filter"] = token[5:]
        else:
            raise ValueError(f"Unknown run prefix token: {token}")

    return parsed


def parse_probs_filename(file_name, known_features=None):
    prefix, feature = split_probs_filename(file_name, known_features=known_features)
    parsed = parse_run_prefix(prefix)
    parsed["prefix"] = prefix
    parsed["feature"] = feature
    return parsed