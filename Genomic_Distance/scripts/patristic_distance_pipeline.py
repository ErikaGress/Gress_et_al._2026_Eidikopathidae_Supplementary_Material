# Gress et al. 2026. Supplementary Material. Python script to recompute patristic-distance matrices and family summaries from the tree.

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo


def normalize(value):
    return str(value).strip().replace("-", "_")


def bootstrap_ci(values, n=500, rng=None):
    values = np.asarray(values, dtype=float)
    values = values[~np.isnan(values)]

    if values.size < 2:
        return np.nan, np.nan

    rng = np.random.default_rng() if rng is None else rng
    means = np.empty(n, dtype=float)
    for idx in range(n):
        sample = rng.choice(values, size=values.size, replace=True)
        means[idx] = np.mean(sample)

    return tuple(np.percentile(means, [2.5, 97.5]))


def remove_outliers(values):
    values = np.asarray(values, dtype=float)
    values = values[~np.isnan(values)]

    if values.size < 4:
        return values

    q1, q3 = np.percentile(values, [25, 75])
    iqr = q3 - q1
    filtered = values[(q1 - 1.5 * iqr <= values) & (values <= q3 + 1.5 * iqr)]
    return filtered if filtered.size else values


def read_csv(path, **kwargs):
    with open(path, newline="") as handle:
        skiprows = 0
        first_line = handle.readline().strip()
        if first_line.startswith("Gress et al. 2026. Supplementary Material."):
            skiprows = 1
    return pd.read_csv(path, comment="#", skip_blank_lines=True, skiprows=skiprows, **kwargs)


def read_tree(path):
    raw = Path(path).read_text()
    lines = [line for line in raw.splitlines() if not line.lstrip().startswith("[")]
    return Phylo.read_from_string("\n".join(lines), "newick")


def load_unresolved_species(path):
    unresolved_df = read_csv(path)
    column = "species" if "species" in unresolved_df.columns else unresolved_df.columns[0]
    series = unresolved_df[column].dropna().astype(str).str.strip()
    return {normalize(value) for value in series if value}


def load_taxa_table(path, excluded_families):
    taxa_df = read_csv(path)
    required = {"family", "species"}
    missing = required.difference(taxa_df.columns)
    if missing:
        raise ValueError(f"Taxa file is missing required columns: {sorted(missing)}")

    taxa_df = taxa_df.loc[:, ["family", "species"]].copy()
    taxa_df["family"] = taxa_df["family"].astype(str).str.strip()
    taxa_df["species"] = taxa_df["species"].astype(str).map(normalize)
    taxa_df = taxa_df[
        taxa_df["family"].ne("")
        & taxa_df["species"].ne("")
        & ~taxa_df["family"].str.lower().eq("nan")
        & ~taxa_df["species"].str.lower().eq("nan")
        & ~taxa_df["family"].isin(excluded_families)
    ]

    duplicate_species = taxa_df[taxa_df.duplicated("species", keep=False)]
    if not duplicate_species.empty:
        conflicting = duplicate_species.groupby("species")["family"].nunique()
        conflicting = conflicting[conflicting > 1]
        if not conflicting.empty:
            raise ValueError(
                "Some species map to multiple families in the taxa file: "
                + ", ".join(conflicting.index.tolist())
            )
        taxa_df = taxa_df.drop_duplicates(subset=["species"], keep="first")

    return taxa_df


def validate_tree_terminals(tree):
    terminals = [terminal.name for terminal in tree.get_terminals()]
    if not terminals:
        raise ValueError("The tree has no terminal taxa.")
    if any(name is None or str(name).strip() == "" for name in terminals):
        raise ValueError("All terminal taxa in the tree must be named.")

    normalized = pd.Series(terminals).map(normalize)
    duplicated = normalized[normalized.duplicated()].unique().tolist()
    if duplicated:
        raise ValueError("Tree contains duplicated terminal labels after normalization: " + ", ".join(duplicated))

    return terminals


def main():
    parser = argparse.ArgumentParser(description="Compute patristic distances and family-level summaries.")
    parser.add_argument("-tree", required=True, help="Input tree in Newick format.")
    parser.add_argument("-out", required=True, help="Output prefix.")
    parser.add_argument("-taxa", required=True, help="CSV containing family and species columns.")
    parser.add_argument("-unres", required=True, help="CSV listing unresolved species to exclude from family analyses.")
    parser.add_argument(
        "--exclude-family",
        action="append",
        default=["Ameripathidae"],
        help="Family to exclude before analyses. Defaults to Ameripathidae; repeat for multiple families.",
    )
    args = parser.parse_args()

    tree = read_tree(args.tree)
    terminals = validate_tree_terminals(tree)
    unresolved = load_unresolved_species(args.unres)
    taxa_df = load_taxa_table(args.taxa, set(args.exclude_family))

    allowed_species = set(taxa_df["species"])
    terminals = [terminal for terminal in terminals if normalize(terminal) in allowed_species]
    if len(terminals) < 2:
        raise ValueError("Fewer than two taxa remain after applying taxa and family filters.")

    matrix = pd.DataFrame(index=terminals, columns=terminals, dtype=float)
    for idx, taxon1 in enumerate(terminals):
        matrix.loc[taxon1, taxon1] = 0.0
        for taxon2 in terminals[idx + 1 :]:
            try:
                distance = tree.distance(taxon1, taxon2)
            except Exception as exc:
                warnings.warn(f"Could not compute patristic distance for {taxon1} vs {taxon2}: {exc}")
                distance = np.nan
            matrix.loc[taxon1, taxon2] = distance
            matrix.loc[taxon2, taxon1] = distance

    all_values = matrix.to_numpy(dtype=float).flatten()
    all_values = all_values[~np.isnan(all_values)]

    var_mat = pd.DataFrame(np.nan, index=terminals, columns=terminals, dtype=float)
    ci_low = pd.DataFrame(np.nan, index=terminals, columns=terminals, dtype=float)
    ci_high = pd.DataFrame(np.nan, index=terminals, columns=terminals, dtype=float)

    if all_values.size:
        global_var = float(np.var(all_values))
        ci_l, ci_h = bootstrap_ci(all_values)
        var_mat.loc[:, :] = global_var
        ci_low.loc[:, :] = ci_l
        ci_high.loc[:, :] = ci_h

    matrix.to_csv(args.out + "_cophenetic_matrix.csv")
    var_mat.to_csv(args.out + "_var.csv")
    ci_low.to_csv(args.out + "_ci_low.csv")
    ci_high.to_csv(args.out + "_ci_high.csv")

    taxa_df_family = taxa_df[~taxa_df["species"].isin(unresolved)].copy()
    species_to_family = dict(zip(taxa_df_family["species"], taxa_df_family["family"]))
    valid_species = [tip for tip in matrix.index if normalize(tip) in species_to_family]

    if len(valid_species) < 2:
        raise ValueError("Fewer than two taxa remain for family analysis after unresolved-species filtering.")

    matrix_filt = matrix.loc[valid_species, valid_species]
    families = sorted({species_to_family[normalize(species)] for species in valid_species})
    family_to_taxa = {
        family: [tip for tip in valid_species if species_to_family[normalize(tip)] == family]
        for family in families
    }

    fam_matrix = pd.DataFrame(np.nan, index=families, columns=families, dtype=float)
    summary = []

    for idx, family1 in enumerate(families):
        fam_matrix.loc[family1, family1] = 0.0

        for family2 in families[idx + 1 :]:
            values = []
            for taxon1 in family_to_taxa[family1]:
                for taxon2 in family_to_taxa[family2]:
                    distance = matrix_filt.loc[taxon1, taxon2]
                    if not np.isnan(distance):
                        values.append(float(distance))

            if not values:
                continue

            clean_values = remove_outliers(values)
            mean_value = float(np.mean(clean_values))
            ci_l, ci_h = bootstrap_ci(clean_values)

            fam_matrix.loc[family1, family2] = mean_value
            fam_matrix.loc[family2, family1] = mean_value
            summary.append(
                {
                    "family1": family1,
                    "family2": family2,
                    "mean": mean_value,
                    "min": float(np.min(values)),
                    "max": float(np.max(values)),
                    "CI_low": ci_l,
                    "CI_high": ci_h,
                    "n": int(len(values)),
                }
            )

    fam_matrix.to_csv(args.out + "_family_matrix.csv")
    pd.DataFrame(summary).to_csv(args.out + "_family_summary.csv", index=False)


if __name__ == "__main__":
    main()
