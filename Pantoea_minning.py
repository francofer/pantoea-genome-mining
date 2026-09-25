#!/usr/bin/env python3
"""
analisis.py
Analiza un archivo GenBank de PGAP para Pantoea ananatis.
Genera un reporte de texto y tablas CSV.
Uso: python analisis.py -i archivo.gb -t 4 -o results
"""

import argparse
import csv
import re
import sys
from collections import OrderedDict
from Bio import SeqIO

# Patrones para genes de interés (case-insensitive)
INTERESTING_PATTERNS = [
    # Virulencia y secreción
    r"virulence", r"toxin", r"antitoxin", r"secretion", r"effector",
    r"Tss[A-Z]", r"VgrG", r"Hcp", r"Tag[HF]", r"Dot[ADHMIK]", r"IcmK",
    r"Tra[YH]", r"GspE", r"type I", r"type II", r"type III", r"type IV", r"type VI",
    # Resistencia
    r"resistance", r"efflux", r"multidrug", r"antibiotic", r"tellurite",
    r"arsenate", r"bleomycin", r"Ter[BC]", r"TehB", r"MarB", r"Emr[AB]", r"HlyD",
    # Toxina-antitoxina
    r"toxin", r"antitoxin", r"Hip[AB]", r"Maz[EF]", r"Rel[BE]", r"ParE",
    r"Vap[BC]", r"PrlF", r"YhaV", r"SymE", r"RatA", r"Gho[ST]", r"CbtA",
    r"YeeU", r"HicB", r"HigA", r"VagC", r"BsmA", r"YmgB", r"LrgA", r"TisB", r"VENN",
    # Flagelo y motilidad
    r"flagell", r"Flg[A-N]", r"Fli[A-Z]", r"Flh[A-D]", r"Mot[AB]", r"chemotaxis",
    # Biofilm y estrés
    r"biofilm", r"acid resistance", r"peroxide", r"oxidative",
    # Otros factores de virulencia
    r"Srf[BC]", r"MsgA", r"TspB", r"BrkB", r"murein", r"hemolysin",
    r"adhesin", r"invasion", r"capsule",
    # Conjugación
    r"conjug", r"pilus", r"PilN", r"DotD", r"TraH",
]

CATEGORY_MAP = [
    (r"secretion|Tss|VgrG|Hcp|Tag|Dot|Icm|GspE|type (I|II|III|IV|VI)", "Secretion system"),
    (r"toxin|antitoxin|Hip|Maz|Rel|ParE|Vap|PrlF|YhaV|SymE|RatA|Gho|CbtA|YeeU|HicB|HigA|VagC|BsmA|YmgB|LrgA|TisB|VENN", "Toxin-antitoxin"),
    (r"resistance|efflux|multidrug|antibiotic|tellurite|arsenate|bleomycin|Ter|TehB|MarB|Emr|HlyD", "Resistance"),
    (r"flagell|Flg|Fli|Flh|Mot|chemotaxis", "Motility/flagellum"),
    (r"biofilm|acid resistance|peroxide|oxidative", "Biofilm/stress"),
    (r"conjug|pilus|PilN|DotD|TraH", "Conjugation/pili"),
    (r"virulence|Srf|MsgA|TspB|BrkB|murein|hemolysin|adhesin|invasion|capsule", "Virulence factor"),
]

def is_hypothetical(product):
    if not product:
        return False
    return bool(re.search(r"hypothetical|uncharacterized", product, re.I))

def categorize(product):
    if not product:
        return "Other"
    for pattern, cat in CATEGORY_MAP:
        if re.search(pattern, product, re.I):
            return cat
    return "Other"

def is_interesting(product, note=""):
    text = (product + " " + note).lower()
    for pat in INTERESTING_PATTERNS:
        if re.search(pat, text, re.I):
            return True
    return False

def gc_content(seq):
    seq = seq.upper()
    gc = (seq.count('G') + seq.count('C')) / len(seq) * 100
    return round(gc, 2)

def main():
    parser = argparse.ArgumentParser(description="Analiza GenBank de PGAP para Pantoea ananatis")
    parser.add_argument("-i", "--input", required=True, help="Archivo GenBank de entrada")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Número de hilos (no usado, por compatibilidad)")
    parser.add_argument("-o", "--output", default="results", help="Prefijo de salida")
    args = parser.parse_args()

    records = list(SeqIO.parse(args.input, "genbank"))
    if not records:
        sys.exit("No se encontraron registros en el archivo GenBank.")

    # Ordenar por tamaño: el más grande es el cromosoma
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
            "other": 0
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
                    "category": categorize(product) if not hypothetical else "Hypothetical"
                }
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
            "features": counts
        }

    # Reporte de texto
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
            f.write("No se encontraron genes de interés.\n")

        f.write("\n" + "=" * 70 + "\n")
        f.write("End of report\n")

    # CSV: todos los CDS
    cds_csv = f"{args.output}_all_CDS.csv"
    if all_cds:
        with open(cds_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=all_cds[0].keys())
            writer.writeheader()
            writer.writerows(all_cds)

    # CSV: genes de interés
    interest_csv = f"{args.output}_genes_of_interest.csv"
    if interesting_genes:
        with open(interest_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=interesting_genes[0].keys())
            writer.writeheader()
            writer.writerows(interesting_genes)

    print(f"Reporte escrito en {txt_file}")
    print(f"Tabla de todos los CDS en {cds_csv}")
    print(f"Tabla de genes de interés en {interest_csv}")

if __name__ == "__main__":
    main()