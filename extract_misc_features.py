#!/usr/bin/env python3

import argparse
import os
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation

# --------------------------------------------------
# GENE ALIAS TABLE (BUILT IN)
# --------------------------------------------------
gene_aliases = {
    "ATP6": ["ATP6","ATP SYNTHASE F0 SUBUNIT 6","APT6","ATP SYNTHASE A0 SUBUNIT 6",
             "ATP SYNTHASE SUBUNIT 6","ATP SYNTHASE FO SUBUNIT 6","ATPASE6","ATPASE SUBUNIT 6"],
    "ATP8": ["ATP8","ATP SYNTHASE F0 SUBUNIT 8","APT8","ATP SYNTHASE A0 SUBUNIT 8",
             "ATP SYNTHASE SUBUNIT 8","ATP SYNTHASE FO SUBUNIT 8","ATPASE8","ATPASE SUBUNIT 8"],
    "COX1": ["COX1","CYTOCHROME C OXIDASE SUBUNIT 1","CYTOCHROME OXIDASE SUBUNIT I",
             "CYTOCHROME C OXIDASE SUBUNIT I","COXI","CO1","COI","CYTOCHROME COXIDASE SUBUNIT I",
             "CYTOCHROME OXIDASE SUBUNIT 1","CYTOCHROME OXYDASE SUBUNIT 1","CYTOCHROME OXIDASE I",
             "MT-CO1","CYTOCHOME OXIDASE SUBUNIT I","COX-1","C01","COX I",
             "CYTOCHROME OXIDASE C SUBUNIT I","CYTOCHROME OXIDASE 1",
             "CYTCHROME OXIDASE SUBUNIT I","CYTOMCHROME C OXIDASE SUBUNIT 1",
             "CO I","COX 1", "COl"],
    "COX2": ["COX2","CYTOCHROME C OXIDASE SUBUNIT 2","CYTOCHROME OXIDASE SUBUNIT II",
             "CYTOCHROME C OXIDASE SUBUNIT II","COXII","CO2","COII",
             "CYTOCHROME COXIDASE SUBUNIT II","CYTOCHROME OXIDASE SUBUNIT 2",
             "CYTOCHROME OXYDASE SUBUNIT 2","CYTOCHROME C OXIDASE II",
             "CYTOCHROME OXIDASE (CO) II","CO II","COX II","COX II"],
    "COX3": ["COX3","CYTOCHROME C OXIDASE SUBUNIT 3","CYTOCHROME OXIDASE SUBUNIT III",
             "CYTOCHROME C OXIDASE SUBUNIT III","COXIII","CO3","COIII",
             "CYTOCHROME COXIDASE SUBUNIT III","CYTOCHROME OXIDASE SUBUNIT 3",
             "CYTOCHROME OXYDASE SUBUNIT 3","CYTOCHROME OXIDASE III","CO III"],
    "CYTB": ["CYTB","CYTOCHROME B","CYB","COB","COB / CYTB","CYTOCHROME B (COB) GENE",
             "CYT B","COB/CYTB","CB"],
    "ND1":  ["ND1","NAD1","NSD1","NADH1","NADH DEHYDROGENASE SUBUNIT I",
             "NADH DEHYDROGENASE SUBUNIT 1","NADH DESHYDROGENASE SUBUNIT 1",
             "NAD1-0","NDI","NADHI"],
    "ND2":  ["ND2","NAD2","NSD2","NADH2","NADH DEHYDROGENASE SUBUNIT II",
             "NADH DEHYDROGENASE SUBUNIT 2","NADH DESHYDROGENASE SUBUNIT 2","NAD2-0"],
    "ND3":  ["ND3","NAD3","NSD3","NADH3","NADH DEHYDROGENASE SUBUNIT III",
             "NADH DEHYDROGENASE SUBUNIT 3","NADH DESHYDROGENASE SUBUNIT 3","NAD3-0"],
    "ND4":  ["ND4","NAD4","NSD4","NADH4","NADH DEHYDROGENASE SUBUNIT IV",
             "NADH DEHYDROGENASE SUBUNIT 4","NADH DESHYDROGENASE SUBUNIT 4","ND4-0"],
    "ND4L": ["ND4L","NAD4L","NSD4L","NADH4L","NADH DEHYDROGENASE SUBUNIT IVL",
             "NADH DEHYDROGENASE SUBUNIT 4L","NADH DESHYDROGENASE SUBUNIT 4L","NAD4L-0"],
    "ND5":  ["ND5","NAD5","NSD5","NADH5","NADH DEHYDROGENASE SUBUNIT V",
             "NADH DEHYDROGENASE SUBUNIT 5","NADH DESHYDROGENASE SUBUNIT 5","NAD5-0",
             "NADH DEHYDROGENASE SUBUNIT 5 (ND5)", "MTND5"],
    "ND6":  ["ND6","NAD6","NSD6","NADH6","NADH DEHYDROGENASE SUBUNIT VI",
             "NADH DEHYDROGENASE SUBUNIT 6","NADH DESHYDROGENASE SUBUNIT 6","NAD6-0"],
    "LSU": ["28S","28S RRNA","28S RIBOSOMAL RNA","28S-RRNA","28S RNA","LSU",
            "LARGE SUBUNIT RRNA","28S LARGE SUBUNIT RIBOSOMAL RNA",
            "LARGE RIBOSOMAL SUBUNIT RNA","NUCLEAR LARGE SUBUNIT RRNA",
            "NUCLEAR LSU RRNA","LARGE SUBUNIT RIBOSOMAL RNA",
            "RRNL","RNL","28SRRNA","28S-RIBOSOMAL RNA","28S RIBOSOME RNA",
            "28S LARGE SUBUNIT","28S rDNA","28S RIBOSOMAL DNA",
            "28 SRNA","28S-RNA","28S r-rna","28S ribo RNA"],
    "SSU": ["18S","18S RRNA","18S RIBOSOMAL RNA","18S-RRNA","18S RNA","SSU",
            "SMALL SUBUNIT RRNA","18S SMALL SUBUNIT RIBOSOMAL RNA",
            "SMALL RIBOSOMAL SUBUNIT RNA","NUCLEAR SMALL SUBUNIT RRNA",
            "NUCLEAR SSU RRNA","SMALL SUBUNIT RIBOSOMAL RNA",
            "RRNS","RNS","18SRRNA","18S-RIBOSOMAL RNA","18S RIBOSOME RNA",
            "18S rDNA","18S RIBOSOMAL DNA","18 SRNA","18S-RNA"],
    "12S": ["12S","12S RRNA","12S RIBOSOMAL RNA","12S-RRNA","12S RNA",
            "12S SMALL SUBUNIT RIBOSOMAL RNA","SSU","RRNS","12SRRN","12S-RNA"],
    "16S": ["16S","16S RRNA","16S RIBOSOMAL RNA","16S-RRNA","16S RNA",
            "16S LARGE SUBUNIT RIBOSOMAL RNA","LSU","RRNL","16SRRN","16S-RNA"],
    "tRNA": ["TRNA", "TRANSFER_RNA"]
}

# build reverse lookup
alias_lookup = {}
for canon, arr in gene_aliases.items():
    for a in arr:
        alias_lookup[a.upper()] = canon


def normalise_name(name):
    if not name:
        return None
    key = name.strip().upper()
    return alias_lookup.get(key)


def extract_seq(seqrecord, feature):
    loc = feature.location
    if isinstance(loc, FeatureLocation):
        return loc.extract(seqrecord).seq
    parts = [p.extract(seqrecord).seq for p in loc.parts]
    return sum(parts, seqrecord.seq[:0])


def main():
    ap = argparse.ArgumentParser(description="Extract only alias-matched genes from CDS-free/rRNA-free/tRNA-free records.")
    ap.add_argument("-g","--genbank", nargs="+", required=True)
    ap.add_argument("-o","--output", required=True)
    ap.add_argument("-n","--organism", action="store_true")
    ap.add_argument("-f","--filter")
    args = ap.parse_args()

    outdir = os.path.join(args.output, "MISC_LABELS")
    os.makedirs(outdir, exist_ok=True)

    filterlist = set()
    if args.filter:
        with open(args.filter) as f:
            filterlist = {l.strip() for l in f if l.strip()}

    for gb in args.genbank:
        for rec in SeqIO.parse(gb, "genbank"):
            if args.filter and rec.name not in filterlist:
                continue

            # rule 1: skip any record containing CDS / rRNA / tRNA
            types = {f.type.upper() for f in rec.features}
            if any(t in types for t in ["CDS","RRNA","TRNA"]):
                continue

            org = rec.annotations.get("organism","UNKNOWN").replace(" ","_")
            prefix = org if args.organism else rec.name

            seen = {}  # canonical_name → sequence string

            for f in rec.features:
                # Attempt to read gene-like qualifiers
                candidates = []
                for q in ["gene","product","note","gene_synonym"]:
                    v = f.qualifiers.get(q, [])
                    for item in v:
                        candidates.append(item)

                canon = None
                for c in candidates:
                    canon = normalise_name(c)
                    if canon:
                        break
                if not canon:
                    continue  # skip everything not in alias table

                seq = str(extract_seq(rec, f))
                if not seq:
                    continue

                # enforce "only first unique extraction per canonical gene"
                if canon in seen:
                    if seen[canon] == seq:
                        continue  # identical, skip
                else:
                    seen[canon] = seq

                    with open(os.path.join(outdir, f"{canon}.fasta"), "a") as w:
                        w.write(f">{prefix}\n{seq}\n")


if __name__ == "__main__":
    main()
