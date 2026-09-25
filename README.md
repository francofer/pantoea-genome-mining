# Pantoea ananatis 31npA-Arg1-genome mining

Python script to parse a PGAP-annotated GenBank file and extract
genes of interest (virulence, secretion systems, toxin–antitoxin
modules, resistance, motility, conjugation, biofilm, stress response)
from *Pantoea ananatis* strain 31npA-Arg1.

## Requirements

- Python ≥ 3.8
- biopython

Install:

    pip install biopython

## Usage

    python analisis.py -i genome.gb -o results

Arguments:
- `-i` / `--input`: input GenBank file (PGAP annotation)
- `-o` / `--output`: output prefix (default: `results`)
- `-t` / `--threads`: number of threads (accepted for compatibility)

## Outputs

- `results.txt` — human-readable report with per-replicon statistics
- `results_all_CDS.csv` — all coding sequences
- `results_genes_of_interest.csv` — filtered genes of interest

## Data availability

Genome assembly and annotation: NCBI BioProject **PRJNA1446712**.


## License

MIT
