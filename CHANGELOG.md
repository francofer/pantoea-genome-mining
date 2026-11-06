# Changelog

## v1.1 — 2026-11

Fixes in the functional classification of `Pantoea_minning.py`. The categories
assigned by v1.0 are affected; `results_genes_of_interest.csv` must be
regenerated.

### Fixed

- **Evaluation order of `CATEGORY_MAP`.** "Secretion system" was the first
  entry and its pattern included `type (I|II|III|IV|VI)`. Products such as
  *SymE family type I addiction module toxin* and *type II toxin-antitoxin
  system VapB family antitoxin* were therefore assigned to the secretion
  category. Toxin-antitoxin is now evaluated first, and the secretion pattern
  requires the word *secretion*.
- **Gene-symbol patterns without word boundaries.** All patterns are matched
  with `re.IGNORECASE`, so character classes written in upper case also match
  lower case. `r"Rel[BE]"` matched *peptide chain **rele**ase factor 2* and
  *lipoprotein-**rele**asing ABC transporter*; `r"Fli[A-Z]"` matched
  *ArnE/ArnF **fli**ppase*. Word boundaries (`\b`) were added.
- **`murein` in the virulence pattern.** *MpaA murein peptide amidase* was
  classified as a virulence factor. It now falls under cell wall / murein.
- **Duplicated `r"toxin"`** in `INTERESTING_PATTERNS`.

### Added

- Categories that previously had to be curated by hand: `LPS modification`,
  `Cell wall/murein`, `Transport`, `Metabolism`, `Translation`, `Defense`,
  `Mobile genetic element`.
- Patterns for fimbrial genes (`fimbri`, `Cfa[A-E]`, `Cbl[A-D]`), which left
  23 fimbrial proteins in `Other` in v1.0.
- `group` column with the eight-category regrouping used in the figures and in
  Table S2.

## v1.0

Version used for the first submission of the manuscript.

## v1.2 — 2026-11

Version used for the results reported in the manuscript.

### Fixed

- **Over-matching in `INTERESTING_PATTERNS`.** v1.1 widened `r"Tra[YH]"` to
  `r"Tra[A-Z]"`, which under `re.IGNORECASE` also matches *transporter*,
  *transferase*, *transcription* and *translation*; the table grew from 469 to
  1,718 entries. Word boundaries were added to all gene-symbol patterns.
- **Genes lost in v1.1.** The autotransporter assembly module (TamA, TamB), the
  type III restriction–modification endonuclease and a protein annotated simply
  as "toxin" were no longer captured. Patterns were added for all three.

### Added

- Patterns for prepilin peptidases, fimbrial usher proteins and T6SS effectors.
  An Ssp4-family T6SS effector, not detected by v1.0, is now recovered and
  classified under secretion systems.

### Changed

- Mobile genetic elements (transposases, integrases, recombinases) are **not**
  included. v1.1 added them, which brought 32 elements of plasmid p3 into the
  table and shifted the focus away from virulence and adaptation. The complete
  set of coding sequences remains available in `results_all_CDS.csv`.

### Scripts added

- `effectors_table_fig.py`, `hivir_alt_scan.py` and `make_hivir_fig.py`.
