#!/usr/bin/env python3
"""
Pantoea_minning.py  (v1.2)
Parse a PGAP-annotated GenBank file for Pantoea ananatis.
Generates a text report and CSV tables.

v1.2 — see CHANGELOG.md for the full history.

v1.1 — fixes in the functional classification:
  * CATEGORY_MAP reordered: toxin-antitoxin is evaluated BEFORE secretion, so
    "type I addiction module" and "type II toxin-antitoxin system" are no
    longer captured by the secretion pattern.
  * word boundaries added to gene-symbol patterns: with re.IGNORECASE,
    r"Rel[BE]" matched "release" and r"Fli[A-Z]" matched "flippase".
  * "murein" moved from the virulence pattern to a cell-wall category.
  * duplicated r"toxin" removed from INTERESTING_PATTERNS.
  * categories added that previously had to be curated by hand:
    LPS modification, Cell wall/murein, Transport, Metabolism, Translation,
    Defense, Mobile genetic element.
  * new "group" column with the 8-category regrouping used in the figures.

Usage: python analisis.py -i genome.gb -t 4 -o results
"""

import argparse
import csv
import re
import sys
from collections import OrderedDict
from Bio import SeqIO

# Patterns for genes of interest (case-insensitive)
INTERESTING_PATTERNS = [
    # Virulence and secretion
    r"virulence", r"toxin", r"antitoxin", r"secretion", r"effector",
    r"\bTss[A-Z]\b", r"VgrG", r"Hcp", r"\bTag[HF]\b", r"\bDot[ADHMIK]\b", r"IcmK",
    r"\bTra[A-Z]\b", r"\bTrb[A-Z]\b", r"\bGsp[A-O]\b", r"type (I|II|III|IV|VI) secretion",
    # Resistance
    r"resistance", r"efflux", r"multidrug", r"antibiotic", r"tellurite",
    r"arsenate", r"bleomycin", r"\bTer[BC]\b", r"TehB", r"MarB", r"\bEmr[AB]\b", r"HlyD",
    # Toxin-antitoxin
    r"antitoxin", r"addiction module", r"VENN motif", r"\btoxin\b", r"\bHip[AB]\b", r"\bMaz[EF]\b", r"\bRel[BE]F?\b", r"\bParD?E\b",
    r"\bVap[BC]\b", r"PrlF", r"YhaV", r"SymE", r"RatA", r"\bGho[ST]\b", r"CbtA",
    r"YeeU", r"HicB", r"HigA", r"VagC", r"BsmA", r"YmgB", r"LrgA", r"TisB", r"VENN",
    # Flagellum and motility
    r"flagell", r"\bFlg[A-N]\b", r"\bFli[A-Z]\b", r"\bFlh[A-D]\b", r"\bMot[AB]\b", r"chemotaxis", r"methyl-accepting",
    # Biofilm and stress
    r"biofilm", r"acid resistance", r"peroxide", r"oxidative",
    # Other virulence factors
    r"\bSrf[BC]\b", r"MsgA", r"TspB", r"BrkB", r"murein", r"hemolysin",
    r"adhesin", r"invasion", r"capsule", r"autotransporter", r"\bTam[AB]\b",
    # Defence and mobile elements
    r"restriction-modification", r"restriction endonuclease", r"CRISPR",
    r"prepilin", r"\bT6SS\b",
    # Conjugation
    r"conjug", r"pil(us|in)", r"\bPil[A-Z]\b", r"fimbri", r"\bCfa[A-E]\b", r"\bCbl[A-D]\b",
]

CATEGORY_MAP = [
    (r"biofilm|\b(BsmA|YmgB|AriR|CidA|LrgA)\b", "Biofilm/stress"),
    (r"antitoxin|addiction module|VENN motif|pre-toxin|^toxin$|\btoxin\b.*\bsystem\b|"
     r"\b(Hip[AB]|Maz[EF]|Rel[BE]F?|ParDE|ParE|Vap[BC]|PrlF|YhaV|SymE|RatA|"
     r"Gho[ST]|CbtA|YeeU|HicB|HigA|VagC|TisB)\b", "Toxin-antitoxin"),
    (r"conjugat|prepilin|\bpil(us|in|T|W|N|Q)\b|fimbri|adhesin|relaxase|mobilization protein|usher|"
     r"\b(Tra[A-Z]|Trb[A-Z]|Cfa[A-E]|Cbl[A-D])\b", "Conjugation/pili"),
    (r"secretion|\bT6SS\b|autotransporter|\bTam[AB]\b|\b(Tss[A-M]|VgrG|Hcp|Tag[A-Z]|Dot[A-Z]|Icm[A-Z]|Gsp[A-O]|HlyB|HlyD|PulD|Ssp[0-9])\b",
     "Secretion system"),
    (r"flagell|chemotaxis|methyl-accepting|\b(Flg[A-N]|Fli[A-Z]|Flh[A-D]|Mot[AB]|Che[ABRWYZ])\b",
     "Motility/flagellum"),
    (r"lipopolysaccharide|\bLPS\b|O-antigen|flippase|\b(Arn[A-F]|Lpx[A-Z]|Waa[A-Z])\b",
     "LPS modification"),
    (r"murein|peptidoglycan|\bMpaA\b|cell wall|D-alanyl|\bamidase\b", "Cell wall/murein"),
    (r"restriction|CRISPR|\bCas[0-9]|abortive infection|anti-phage", "Defense"),
    (r"resistance|efflux|multidrug|antibiotic|tellurite|arsenate|bleomycin|beta-lactamase|"
     r"\b(Ter[A-F]|TehB|Mar[ABR]|Emr[ABKY]|Acr[AB]|TolC)\b", "Resistance"),
    (r"peroxide|oxidative|catalase|superoxide|universal stress|cold shock|heat shock|"
     r"acid resistance|chaperone", "Biofilm/stress"),
    (r"transposase|integrase|recombinase|insertion sequence|\bIS[0-9]", "Mobile genetic element"),
    (r"hemolysin|invasion|capsule|virulence|\b(Srf[BC]|MsgA|TspB|BrkB|GfcC|SirB)\b",
     "Virulence factor"),
    (r"transporter|permease|symporter|antiporter|\bABC\b", "Transport"),
    (r"dehydratase|dehydrogenase|\bsynthase\b|\bsynthetase\b|isomerase|\btransferase\b|"
     r"\bkinase\b|\blyase\b|\bhydrolase\b|\breductase\b|\bmutase\b|carbonic anhydrase|"
     r"aminopeptidase|dehalogenase|cyclo-ligase", "Metabolism"),
    (r"release factor|ribosom|\btRNA\b|translation", "Translation"),
]



# Reagrupacion de las categorias finas en 8 grupos para figuras y tablas
GROUP = {
    "Secretion system": "Secretion systems",
    "Motility/flagellum": "Motility and chemotaxis",
    "Toxin-antitoxin": "Toxin-antitoxin systems",
    "Conjugation/pili": "Adhesion and conjugation",
    "Mobile genetic element": "Adhesion and conjugation",
    "Resistance": "Resistance and defence",
    "Defense": "Resistance and defence",
    "Cell wall/murein": "Cell envelope and LPS",
    "LPS modification": "Cell envelope and LPS",
    "Biofilm/stress": "Biofilm and stress",
    "Virulence factor": "Other virulence-associated",
    "Metabolism": "Other virulence-associated",
    "Transport": "Other virulence-associated",
    "Translation": "Other virulence-associated",
    "Other": "Other virulence-associated",
}


def is_hypothetical(product):
    """Return True if the product description matches hypothetical keywords."""
    if not product:
        return False
    return bool(re.search(r"hypothetical|uncharacterized", product, re.I))


def categorize(product):
    """Assign a functional category based on the product description."""
    if not product:
        return "Other"
    for pattern, cat in CATEGORY_MAP:
        if re.search(pattern, product, re.I):
            return cat
    return "Other"


def is_interesting(product, note=""):
    """Return True if the product or note matches any pattern of interest."""
    text = (product + " " + note).lower()
    for pat in INTERESTING_PATTERNS:
        if re.search(pat, text, re.I):
            return True
    return False


def gc_content(seq):
    """Return GC content (%) of a sequence, rounded to 2 decimals."""
    seq = seq.upper()
    gc = (seq.count('G') + seq.count('C')) / len(seq) * 100
    return round(gc, 2)


def main():
    parser = argparse.ArgumentParser(
        description="Parse a PGAP-annotated GenBank file for Pantoea ananatis"
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Input GenBank file")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of threads (not used, kept for compatibility)")
    parser.add_argument("-o", "--output", default="results",
                        help="Output prefix")
    args = parser.parse_args()

    records = list(SeqIO.parse(args.input, "genbank"))
    if not records:
        sys.exit("No records found in the GenBank file.")

    # Sort by size: the largest record is assumed to be the chromosome
    records_sorted = sorted(records, key=lambda r: len(r), reverse=True)
    replicon_info = {}
    for i, rec in enumerate(records_sorted):
        if i == 0:
            rep_type = "chromosome"
            rep_name = "chromosome"
        else:
            rep_type = "plasmid"
            rep_name = f"p{i}"
        replicon_info[rec.id] = {"type": rep_type, "name": rep_name, "record": rec}

    stats = OrderedDict()
    all_cds = []
    interesting_genes = []

    for rec in records:
        info = replicon_info[rec.id]
        rep_name = info["name"]
        rep_type = info["type"]
        length = len(rec)
        gc = gc_content(rec.seq)

        counts = {
            "CDS": 0,
            "pseudogene": 0,
            "rRNA": 0,
            "tRNA": 0,
            "hypothetical": 0,
            "other": 0,
        }
        for feature in rec.features:
            ft = feature.type
            if ft == "CDS":
                product = feature.qualifiers.get("product", [""])[0]
                locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                gene = feature.qualifiers.get("gene", [""])[0]
                note = feature.qualifiers.get("note", [""])[0]
                is_pseudo = bool(feature.qualifiers.get("pseudogene"))
                hypothetical = is_hypothetical(product)
                start = int(feature.location.start) + 1
                end = int(feature.location.end)
                strand = "+" if feature.location.strand == 1 else "-"
                length_bp = len(feature.location)
                if is_pseudo:
                    counts["pseudogene"] += 1
                else:
                    counts["CDS"] += 1
                    if hypothetical:
                        counts["hypothetical"] += 1
                cds_entry = {
                    "replicon": rep_name,
                    "type": rep_type,
                    "contig": rec.id,
                    "locus_tag": locus_tag,
                    "gene": gene,
                    "product": product,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "length_bp": length_bp,
                    "hypothetical": "yes" if hypothetical else "no",
                    "pseudogene": "yes" if is_pseudo else "no",
                    "note": note,
                    "category": categorize(product) if not hypothetical else "Hypothetical",
                }
                cds_entry["group"] = GROUP.get(cds_entry["category"], "Other virulence-associated")
                all_cds.append(cds_entry)
                if not is_pseudo and not hypothetical and is_interesting(product, note):
                    interesting_genes.append(cds_entry)
            elif ft == "pseudogene":
                counts["pseudogene"] += 1
            elif ft == "rRNA":
                counts["rRNA"] += 1
            elif ft == "tRNA":
                counts["tRNA"] += 1
            else:
                counts["other"] += 1

        stats[rep_name] = {
            "type": rep_type,
            "contig": rec.id,
            "length": length,
            "gc": gc,
            "features": counts,
        }

    # Text report
    txt_file = f"{args.output}.txt"
    with open(txt_file, "w") as f:
        f.write("Genome Annotation Report for Pantoea ananatis\n")
        f.write("=" * 70 + "\n\n")
        total_cds = sum(s["features"]["CDS"] for s in stats.values())
        total_pseudo = sum(s["features"]["pseudogene"] for s in stats.values())
        total_rrna = sum(s["features"]["rRNA"] for s in stats.values())
        total_trna = sum(s["features"]["tRNA"] for s in stats.values())
        total_hypo = sum(s["features"]["hypothetical"] for s in stats.values())
        total_length = sum(s["length"] for s in stats.values())
        f.write("GENERAL STATISTICS\n")
        f.write("-" * 70 + "\n")
        f.write(f"Total genome size (bp): {total_length:,}\n")
        f.write(f"Number of replicons: {len(stats)}\n")
        f.write(f"Total CDSs: {total_cds:,}\n")
        f.write(f"  - Hypothetical proteins: {total_hypo:,} ({total_hypo/total_cds*100:.1f}%)\n")
        f.write(f"Total pseudogenes: {total_pseudo:,}\n")
        f.write(f"Total rRNA genes: {total_rrna:,}\n")
        f.write(f"Total tRNA genes: {total_trna:,}\n\n")

        f.write("PER-REPLICON BREAKDOWN\n")
        f.write("-" * 70 + "\n")
        for rep_name, s in stats.items():
            feats = s["features"]
            hypo_pct = (feats["hypothetical"] / feats["CDS"] * 100) if feats["CDS"] > 0 else 0
            f.write(f"\nReplicon: {rep_name} ({s['type']})\n")
            f.write(f"  Contig: {s['contig']}\n")
            f.write(f"  Length (bp): {s['length']:,}\n")
            f.write(f"  GC content (%): {s['gc']}\n")
            f.write(f"  CDSs: {feats['CDS']:,}\n")
            f.write(f"    - Hypothetical: {feats['hypothetical']:,} ({hypo_pct:.1f}%)\n")
            f.write(f"  Pseudogenes: {feats['pseudogene']:,}\n")
            f.write(f"  rRNA genes: {feats['rRNA']:,}\n")
            f.write(f"  tRNA genes: {feats['tRNA']:,}\n")

        f.write("\n" + "=" * 70 + "\n")
        f.write("SELECTED FEATURES OF INTEREST\n")
        f.write("-" * 70 + "\n")
        if interesting_genes:
            for g in interesting_genes:
                f.write(f"\nReplicon: {g['replicon']} ({g['type']})\n")
                f.write(f"  Locus tag: {g['locus_tag']}\n")
                f.write(f"  Gene: {g['gene']}\n")
                f.write(f"  Product: {g['product']}\n")
                f.write(f"  Category: {g['category']}\n")
                f.write(f"  Location: {g['start']}..{g['end']} ({g['strand']})\n")
                if g['note']:
                    f.write(f"  Note: {g['note']}\n")
        else:
            f.write("No genes of interest were found.\n")

        f.write("\n" + "=" * 70 + "\n")
        f.write("End of report\n")

    # CSV: all CDSs
    cds_csv = f"{args.output}_all_CDS.csv"
    if all_cds:
        with open(cds_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=all_cds[0].keys())
            writer.writeheader()
            writer.writerows(all_cds)

    # CSV: genes of interest
    interest_csv = f"{args.output}_genes_of_interest.csv"
    if interesting_genes:
        with open(interest_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=interesting_genes[0].keys())
            writer.writeheader()
            writer.writerows(interesting_genes)

    print(f"Report written to {txt_file}")
    print(f"All CDS table written to {cds_csv}")
    print(f"Genes-of-interest table written to {interest_csv}")


if __name__ == "__main__":
    main()
