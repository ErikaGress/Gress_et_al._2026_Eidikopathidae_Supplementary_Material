# Gress et al. 2026: Eidikopathidae fam. nov. Supplementary Material

This repository contains supplementary files, analysis scripts, and result tables associated with the genomic-distance analyses for *Eidikopathidae* fam. nov.

## Contents

```text
Supplementary_Genomic_Distance/
+-- Supplementary_Table_1.xlsx
+-- Supplementary_Table_2.xlsx
+-- Genomic_Distance/
|   +-- lists/
|   |   +-- taxaID.csv
|   |   +-- unresolved_species.csv
|   +-- results/
|   |   +-- Genomic_distance_Heatmaps.docx
|   |   +-- Patristic/
|   |   +-- TrN/
|   +-- scripts/
|       +-- patristic_distance_pipeline.py
|       +-- trn_distance_pipeline.py
+-- Trees/
    +-- Gress_et_al_2026_ASTRAL.tre
    +-- Gress_et_al_2026_iqtree.treefile
```

## File Description

### Supplementary Tables

- `Supplementary_Table_1.xlsx`  
  Specimen collection details for *Eidikopathidae* fam. nov.

- `Supplementary_Table_2.xlsx`  
  Results of UCE analyses from newly collected specimens.

### Genomic Distance

Located in `Genomic_Distance/`.

This folder contains the lookup tables, Python scripts, and output files used for the genomic-distance analyses.

### Lists

- `Genomic_Distance/lists/taxaID.csv`  
  Taxon and family lookup table used for the genomic-distance analyses.

- `Genomic_Distance/lists/unresolved_species.csv`  
  List of unresolved specimens excluded from family- and genus-level analyses.

### Scripts

- `Genomic_Distance/scripts/trn_distance_pipeline.py`  
  Python script used to calculate genomic distances from aligned UCE loci under the Tamura-Nei (TrN) model and summarize distances among species, genera, and families.

- `Genomic_Distance/scripts/patristic_distance_pipeline.py`  
  Python script used to calculate patristic distances from a phylogenetic tree and summarize distances among species, genera, and families.

### Results

Located in `Genomic_Distance/results/`.

- `Genomic_distance_Heatmaps.docx`  
  Supplementary heatmap figures summarizing genomic-distance results.

#### TrN Distance Results

Located in `Genomic_Distance/results/TrN/`.

- `TrN_mean.csv`: mean pairwise TrN genomic-distance matrix.
- `TrN_var.csv`: variance estimates for pairwise TrN distances.
- `TrN_ci_low.csv` and `TrN_ci_high.csv`: lower and upper confidence intervals for TrN distances.
- `TrN_family_mean.csv`: family-level mean TrN distance matrix.
- `TrN_family_summary.csv`: summary statistics for family-level TrN distances.
- `family_min_max_distances.csv`: minimum and maximum family-level TrN distances.
- `TrN_between_genus_distances.csv`: TrN distances between genera.
- `TrN_within_genus_distances.csv`: TrN distances within genera.
- `Eidikopathidae_interspecific.csv`: interspecific TrN distances within *Eidikopathidae*.
- `Eidikopathidae_intraspecific.csv`: intraspecific TrN distances within *Eidikopathidae*.

#### Patristic Distance Results

Located in `Genomic_Distance/results/Patristic/`.

- `Patristic_mean.csv`: mean pairwise patristic-distance matrix.
- `Patristic_var.csv`: variance estimates for pairwise patristic distances.
- `Patristic_ci_low.csv` and `Patristic_ci_high.csv`: lower and upper confidence intervals for patristic distances.
- `Patristic_family_mean.csv`: family-level mean patristic-distance matrix.
- `Patristic_family_summary.csv`: summary statistics for family-level patristic distances.
- `Patristic_between_genus_distances.csv`: patristic distances between genera.
- `Patristic_within_genus_distances.csv`: patristic distances within genera.
- `Eidikopathidae_interspecific.csv`: interspecific patristic distances within *Eidikopathidae*.
- `Eidikopathidae_intraspecific.csv`: intraspecific patristic distances within *Eidikopathidae*.

### Trees

Located in `Trees/`.

- `Gress_et_al_2026_ASTRAL.tre`  
  ASTRAL species tree used in the supplementary analyses.

- `Gress_et_al_2026_iqtree.treefile`  
  IQ-TREE phylogenetic tree file used in the supplementary analyses.

## Citation

Gress, E., Li, Y.-X., Opresko, D. M., & Qiu, J.-W. 2026. Integrating genomic distance analyses in the description of a new family, genus, and species of sponge-associated antipatharians (black corals). Molecular Phylogenetics and Evolution. https://doi.org/10.1016/j.ympev.2026.108718
