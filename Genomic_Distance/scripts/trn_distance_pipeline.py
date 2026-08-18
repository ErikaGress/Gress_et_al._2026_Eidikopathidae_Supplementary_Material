# Gress et al. 2026. Supplementary Material. Python script to recompute Tamura-Nei genomic-distance matrices and family summaries.

import argparse
import csv
import itertools
import random
from collections import defaultdict
from math import log
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO


def normalize_id(value):
    return str(value).strip().replace("-", "_")


def uncommented_lines(path):
    with open(path, newline="") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("Gress et al. 2026. Supplementary Material."):
                continue
            if not line.lstrip().startswith("#"):
                yield line


def read_taxa(path, excluded_families):
    species_list = []
    species_to_family = {}

    reader = csv.DictReader(uncommented_lines(path))
    for row in reader:
        family = row["family"].strip()
        species = row["species"].strip()
        if family in excluded_families:
            continue
        species_list.append(species)
        species_to_family[species] = family

    return species_list, species_to_family


def read_unresolved(path):
    if path is None:
        return set()

    reader = csv.DictReader(uncommented_lines(path))
    return {normalize_id(row["species"]) for row in reader if row.get("species")}


def trn_distance(seq1, seq2):
    p1 = p2 = q = count = 0
    freq = {"A": 0, "C": 0, "G": 0, "T": 0}

    for base1, base2 in zip(seq1.upper(), seq2.upper()):
        if base1 not in "ACGT" or base2 not in "ACGT":
            continue

        count += 1
        freq[base1] += 1
        freq[base2] += 1

        if base1 == base2:
            continue
        if base1 in "AG" and base2 in "AG":
            p1 += 1
        elif base1 in "CT" and base2 in "CT":
            p2 += 1
        else:
            q += 1

    if count == 0:
        return np.nan, 0

    p1 /= count
    p2 /= count
    q /= count

    total = count * 2
    pi_a = freq["A"] / total + 1e-10
    pi_c = freq["C"] / total + 1e-10
    pi_g = freq["G"] / total + 1e-10
    pi_t = freq["T"] / total + 1e-10

    base_sum = pi_a + pi_c + pi_g + pi_t
    pi_a /= base_sum
    pi_c /= base_sum
    pi_g /= base_sum
    pi_t /= base_sum

    purine = pi_a + pi_g
    pyrimidine = pi_c + pi_t
    if purine < 1e-9 or pyrimidine < 1e-9:
        return np.nan, count

    a_term = 2 * pi_a * pi_g / purine
    b_term = 2 * pi_c * pi_t / pyrimidine
    c_term = 2 * purine * pyrimidine

    try:
        w1 = 1 - p1 / (2 * a_term) - q / (2 * c_term) if a_term > 1e-9 else 1 - q / (2 * c_term)
        w2 = 1 - p2 / (2 * b_term) - q / (2 * c_term) if b_term > 1e-9 else 1 - q / (2 * c_term)
        w3 = 1 - q / (2 * c_term)

        if w1 <= 0 or w2 <= 0 or w3 <= 0:
            return np.nan, count

        distance = -a_term * log(w1) - b_term * log(w2) - c_term * log(w3)
        return distance, count
    except (FloatingPointError, ValueError, ZeroDivisionError):
        return np.nan, count


def bootstrap_ci(values, n=500):
    if len(values) < 2:
        return np.nan, np.nan

    means = []
    for _ in range(n):
        sample = [random.choice(values) for _ in range(len(values))]
        means.append(np.mean(sample))

    return np.percentile(means, [2.5, 97.5])


def remove_outliers(values):
    if len(values) < 4:
        return values

    q1, q3 = np.percentile(values, [25, 75])
    iqr = q3 - q1
    filtered = [value for value in values if (q1 - 1.5 * iqr) <= value <= (q3 + 1.5 * iqr)]
    return filtered if filtered else values


def main():
    parser = argparse.ArgumentParser(description="Compute TrN genomic distances from aligned UCE loci.")
    parser.add_argument("-dir", required=True, help="Folder containing Nexus or FASTA alignments.")
    parser.add_argument("-taxa", required=True, help="CSV with family and species columns.")
    parser.add_argument("-unres", default=None, help="CSV listing unresolved species to exclude from family analyses.")
    parser.add_argument("-out", required=True, help="Output prefix.")
    parser.add_argument("-fmt", default="nexus", choices=["fasta", "nexus"], help="Alignment format.")
    parser.add_argument(
        "--exclude-family",
        action="append",
        default=["Ameripathidae"],
        help="Family to exclude before analyses. Defaults to Ameripathidae; repeat for multiple families.",
    )
    args = parser.parse_args()

    alignment_folder = Path(args.dir)
    species_list, species_to_family = read_taxa(args.taxa, set(args.exclude_family))
    unresolved = read_unresolved(args.unres)
    normalized_species = {normalize_id(species) for species in species_list}

    pairs_dist = defaultdict(list)
    pairs_weights = defaultdict(list)
    files = sorted(path for path in alignment_folder.iterdir() if path.suffix.lower() == f".{args.fmt}")

    for path in files:
        records = list(SeqIO.parse(str(path), args.fmt))
        taxa_seq = {}

        for record in records:
            record_id = normalize_id(record.id)
            if record_id not in normalized_species:
                continue

            for species in species_list:
                if normalize_id(species) == record_id:
                    taxa_seq[species] = str(record.seq)
                    break

        for taxon1, taxon2 in itertools.combinations(taxa_seq.keys(), 2):
            distance, sites = trn_distance(taxa_seq[taxon1], taxa_seq[taxon2])
            if not np.isnan(distance) and sites > 0:
                key = tuple(sorted([taxon1, taxon2]))
                pairs_dist[key].append(distance)
                pairs_weights[key].append(sites)

    species_sorted = sorted(species_list)
    mean_mat = pd.DataFrame(index=species_sorted, columns=species_sorted)
    var_mat = pd.DataFrame(index=species_sorted, columns=species_sorted)
    ci_low = pd.DataFrame(index=species_sorted, columns=species_sorted)
    ci_high = pd.DataFrame(index=species_sorted, columns=species_sorted)
    np.fill_diagonal(mean_mat.values, 0)

    for idx, taxon1 in enumerate(species_sorted):
        for taxon2 in species_sorted[idx + 1 :]:
            key = tuple(sorted([taxon1, taxon2]))
            values = pairs_dist[key]
            weights = pairs_weights[key]

            if values:
                mean = np.average(values, weights=weights)
                var = np.var(values)
                low, high = bootstrap_ci(values)
            else:
                mean = var = low = high = np.nan

            mean_mat.loc[taxon1, taxon2] = mean
            mean_mat.loc[taxon2, taxon1] = mean
            var_mat.loc[taxon1, taxon2] = var
            var_mat.loc[taxon2, taxon1] = var
            ci_low.loc[taxon1, taxon2] = low
            ci_low.loc[taxon2, taxon1] = low
            ci_high.loc[taxon1, taxon2] = high
            ci_high.loc[taxon2, taxon1] = high

    mean_mat.to_csv(args.out + "_mean.csv")
    var_mat.to_csv(args.out + "_var.csv")
    ci_low.to_csv(args.out + "_ci_low.csv")
    ci_high.to_csv(args.out + "_ci_high.csv")

    species_to_family_norm = {normalize_id(species): family for species, family in species_to_family.items()}
    valid_species = []
    for species in mean_mat.index:
        species_norm = normalize_id(species)
        if species_norm in species_to_family_norm and species_norm not in unresolved:
            if mean_mat.loc[species].notna().sum() > 1:
                valid_species.append(species)

    mean_mat = mean_mat.loc[valid_species, valid_species]
    families = sorted({species_to_family_norm[normalize_id(species)] for species in valid_species})
    family_to_taxa = {
        family: [species for species in valid_species if species_to_family_norm[normalize_id(species)] == family]
        for family in families
    }

    fam_matrix = pd.DataFrame(index=families, columns=families)
    summary = []

    for family1 in families:
        for family2 in families:
            if family1 == family2:
                fam_matrix.loc[family1, family2] = 0
                continue

            values = []
            for taxon1 in family_to_taxa[family1]:
                for taxon2 in family_to_taxa[family2]:
                    distance = mean_mat.loc[taxon1, taxon2]
                    if not pd.isna(distance):
                        values.append(float(distance))

            if not values:
                fam_matrix.loc[family1, family2] = np.nan
                continue

            clean_values = remove_outliers(values)
            ci_l, ci_h = bootstrap_ci(clean_values)
            fam_matrix.loc[family1, family2] = np.mean(clean_values)
            summary.append(
                {
                    "family1": family1,
                    "family2": family2,
                    "mean": np.mean(clean_values),
                    "min": np.min(values),
                    "max": np.max(values),
                    "CI_low": ci_l,
                    "CI_high": ci_h,
                    "n": len(values),
                }
            )

    fam_matrix.to_csv(args.out + "_family_mean.csv")
    pd.DataFrame(summary).to_csv(args.out + "_family_summary.csv", index=False)


if __name__ == "__main__":
    main()
