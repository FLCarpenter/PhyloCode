#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Extract CDS, rRNA, and tRNA sequences from GenBank files,
grouping them by standardized gene names.
Outputs are split into subfolders: CDS/, rRNA/, tRNA/
"""

import argparse
import textwrap as _textwrap
import os
import sys
import datetime
from collections import defaultdict
from Bio import SeqIO

# -----------------------------------------------------------------------------
# Logging
# -----------------------------------------------------------------------------
def init_logger(output_dir):
    logpath = os.path.join(output_dir, "run.log")
    fh = open(logpath, "w")
    def log(msg):
        timestamp = datetime.datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")
        line = f"{timestamp} {msg}\n"
        fh.write(line)
        fh.flush()
        sys.stderr.write(line)
    return log

# -----------------------------------------------------------------------------
# Custom help formatter
# -----------------------------------------------------------------------------
class MultilineFormatter(argparse.HelpFormatter):
    def _fill_text(self, text, width, indent):
        text = self._whitespace_matcher.sub(' ', text).strip()
        paragraphs = text.split('|n ')
        multiline_text = ''
        for paragraph in paragraphs:
            formatted_paragraph = _textwrap.fill(
                paragraph,
                width,
                initial_indent=indent,
                subsequent_indent=indent
            ) + '\n\n'
            multiline_text += formatted_paragraph
        return multiline_text

# -----------------------------------------------------------------------------
# Gene alias definitions
# -----------------------------------------------------------------------------
gene_aliases = {
    "ATP6": ["ATP6","ATP SYNTHASE F0 SUBUNIT 6","APT6","ATP SYNTHASE A0 SUBUNIT 6",
             "ATP SYNTHASE SUBUNIT 6","ATP SYNTHASE FO SUBUNIT 6","ATPASE6","ATPASE SUBUNIT 6"],
    "ATP8": ["ATP8","ATP SYNTHASE F0 SUBUNIT 8","APT8","ATP SYNTHASE A0 SUBUNIT 8",
             "ATP SYNTHASE SUBUNIT 8","ATP SYNTHASE FO SUBUNIT 8","ATPASE8","ATPASE SUBUNIT 8"],
    "COX1": ["COX1","CYTOCHROME C OXIDASE SUBUNIT 1","CYTOCHROME OXIDASE SUBUNIT I",
             "CYTOCHROME C OXIDASE SUBUNIT I","COXI","CO1","COI","CYTOCHROME COXIDASE SUBUNIT I",
             "CYTOCHROME OXIDASE SUBUNIT 1","CYTOCHROME OXYDASE SUBUNIT 1","CYTOCHROME OXIDASE I",
             "MT-CO1","CYTOCHOME OXIDASE SUBUNIT I","COX-1","C01","COX I","CYTOCHROME OXIDASE C SUBUNIT I",
             "CYTOCHROME OXIDASE 1","CYTCHROME OXIDASE SUBUNIT I","CYTOMCHROME C OXIDASE SUBUNIT 1","CO I","COX 1", "COl"],
    "COX2": ["COX2","CYTOCHROME C OXIDASE SUBUNIT 2","CYTOCHROME OXIDASE SUBUNIT II",
             "CYTOCHROME C OXIDASE SUBUNIT II","COXII","CO2","COII","CYTOCHROME COXIDASE SUBUNIT II",
             "CYTOCHROME OXIDASE SUBUNIT 2","CYTOCHROME OXYDASE SUBUNIT 2","CYTOCHROME C OXIDASE II",
             "CYTOCHROME OXIDASE (CO) II","CO II","COX II","COX II"],
    "COX3": ["COX3","CYTOCHROME C OXIDASE SUBUNIT 3","CYTOCHROME OXIDASE SUBUNIT III",
             "CYTOCHROME C OXIDASE SUBUNIT III","COXIII","CO3","COIII","CYTOCHROME COXIDASE SUBUNIT III",
             "CYTOCHROME OXIDASE SUBUNIT 3","CYTOCHROME OXYDASE SUBUNIT 3","CYTOCHROME OXIDASE III", "CO III"],
    "CYTB": ["CYTB","CYTOCHROME B","CYB","COB","COB / CYTB","CYTOCHROME B (COB) GENE","CYT B","COB/CYTB","CB"],
    "ND1":  ["ND1","NAD1","NSD1","NADH1","NADH DEHYDROGENASE SUBUNIT I",
             "NADH DEHYDROGENASE SUBUNIT 1","NADH DESHYDROGENASE SUBUNIT 1","NAD1-0","NDI","NADHI"],
    "ND2":  ["ND2","NAD2","NSD2","NADH2","NADH DEHYDROGENASE SUBUNIT II",
             "NADH DEHYDROGENASE SUBUNIT 2","NADH DESHYDROGENASE SUBUNIT 2","NAD2-0"],
    "ND3":  ["ND3","NAD3","NSD3","NADH3","NADH DEHYDROGENASE SUBUNIT III",
             "NADH DEHYDROGENASE SUBUNIT 3","NADH DESHYDROGENASE SUBUNIT 3","NAD3-0"],
    "ND4":  ["ND4","NAD4","NSD4","NADH4","NADH DEHYDROGENASE SUBUNIT IV",
             "NADH DEHYDROGENASE SUBUNIT 4","NADH DESHYDROGENASE SUBUNIT 4","NAD4-0"],
    "ND4L": ["ND4L","NAD4L","NSD4L","NADH4L","NADH DEHYDROGENASE SUBUNIT IVL",
             "NADH DEHYDROGENASE SUBUNIT 4L","NADH DESHYDROGENASE SUBUNIT 4L","NAD4L-0"],
    "ND5":  ["ND5","NAD5","NSD5","NADH5","NADH DEHYDROGENASE SUBUNIT V",
             "NADH DEHYDROGENASE SUBUNIT 5","NADH DESHYDROGENASE SUBUNIT 5","NAD5-0","NADH DEHYDROGENASE SUBUNIT 5 (ND5)", "MTND5"],
    "ND6":  ["ND6","NAD6","NSD6","NADH6","NADH DEHYDROGENASE SUBUNIT VI",
             "NADH DEHYDROGENASE SUBUNIT 6","NADH DESHYDROGENASE SUBUNIT 6","NAD6-0"],
    "LSU": ["28S","28S RRNA","28S RIBOSOMAL RNA","28S-RRNA","28S RNA","LSU",
            "LARGE SUBUNIT RRNA", "28S LARGE SUBUNIT RIBOSOMAL RNA",
            "LARGE RIBOSOMAL SUBUNIT RNA","NUCLEAR LARGE SUBUNIT RRNA","NUCLEAR LSU RRNA",
            "LARGE SUBUNIT RIBOSOMAL RNA","RRNL","RNL","28SRRNA","28S-RIBOSOMAL RNA",
            "28S RIBOSOME RNA", "28S LARGE SUBUNIT","28S rDNA","28S RIBOSOMAL DNA",
            "28 SRNA","28S-RNA","28S r-rna","28S ribo RNA","28S RIBOSOML RNA"],
    "SSU": ["18S","18S RRNA","18S RIBOSOMAL RNA","18S-RRNA","18S RNA","SSU",
            "SMALL SUBUNIT RRNA", "18S SMALL SUBUNIT RIBOSOMAL RNA","SMALL RIBOSOMAL SUBUNIT RNA",
            "NUCLEAR SMALL SUBUNIT RRNA","NUCLEAR SSU RRNA", "SMALL SUBUNIT RIBOSOMAL RNA",
            "RRNS","RNS","18SRRNA","18S-RIBOSOMAL RNA","18S RIBOSOME RNA",
            "18S rDNA","18S RIBOSOMAL DNA","18 SRNA","18S-RNA","18S r-rna","18S ribo RNA",
            "18S RIBOSOML RNA", "18S small subunit"],
    "12S": ["12S","12S RRNA","12S RIBOSOMAL RNA","12S-RRNA","12S RNA",
            "12S SMALL SUBUNIT RIBOSOMAL RNA","SSU","RRNS","12SRRN","12S-RNA",
            "12S RIBOSOMAL DNA", "12S small ribosomal RNA", "12SrRNA", "12srRNA"],
    "16S": ["16S","16S RRNA","16S RIBOSOMAL RNA","16S-RRNA","16S RNA",
            "16S LARGE SUBUNIT RIBOSOMAL RNA","LSU","RRNL","16SRRN","16S-RNA",
            "16S RIBOSOMAL DNA", "16S rRNA gene", "16S ribosomal RNA, large subunit","16SrRNA"],
    "tRNA": ["TRNA", "TRANSFER_RNA"]
}

# -----------------------------------------------------------------------------
# Load gene variants
# -----------------------------------------------------------------------------
def loadnamevariants(report=False):
    output = {}      # variant → canonical name
    fullparse = {}   # canonical name → {type, variants}
    alltypes = set() # set of all annotation types used

    rRNA_KEYS = {"12S", "16S", "LSU", "SSU"}
    tRNA_KEYS = {"TRNA"}

    for name, variants in gene_aliases.items():
        key = name.upper()
        if key in rRNA_KEYS:
            annotype = "rRNA"
        elif key in tRNA_KEYS:
            annotype = "tRNA"
        else:
            annotype = "CDS"

        alltypes.add(annotype)
        fullvariants = []
        for v in variants:
            v = v.upper()
            for suffix in ['', ' GENE', f' {annotype.upper()}']:
                variant = v + suffix
                fullvariants.append(variant)
                output[variant] = name
        fullparse[name] = {'type': annotype, 'variants': fullvariants}

        if report:
            formatted = _textwrap.fill(', '.join(fullvariants), width=80,
                                       initial_indent='\t', subsequent_indent='\t')
            print(f"Standard name = {name}, type = {annotype}:\n{formatted}")

    return output, alltypes, fullparse

# -----------------------------------------------------------------------------
# Extract feature name
# -----------------------------------------------------------------------------
def get_feat_name(feat):
    featname = "unknown"
    nametags = ['gene', 'product', 'label', 'standard_name']
    for t in nametags:
        if t in feat.qualifiers:
            featname = feat.qualifiers[t][0].strip().upper()
            break
    return featname

# -----------------------------------------------------------------------------
# CLI arguments
# -----------------------------------------------------------------------------
def getcliargs(arglist=None, knowngenes=None, knowntypes=None):
    parser = argparse.ArgumentParser(
        description="Extract CDS, rRNA, tRNA from GenBank files",
        formatter_class=MultilineFormatter)

    parser.add_argument("-g", "--genbank", type=str, metavar='PATH', required=True, nargs='+')
    parser.add_argument("-o", "--output", type=str, metavar='PATH', required=True)
    parser.add_argument("-m", "--mingenes", type=int, default=0)
    parser.add_argument("-q", "--reqgenes", type=str, nargs='*', choices=knowngenes)
    parser.add_argument("-f", "--filter", type=str)
    parser.add_argument("-n", "--organism", action='store_true')
    parser.add_argument("-w", "--writeunknowns", action='store_true')
    parser.add_argument("-k", "--keepframe", action='store_true')
    parser.add_argument("-s", "--showgenes", action='store_true')
    parser.add_argument("-v", "--version", action='store_true')

    args = parser.parse_args(arglist) if arglist else parser.parse_args()
    return args

# -----------------------------------------------------------------------------
# Priority dictionary
# -----------------------------------------------------------------------------
PRIORITY = {"CDS": 3, "rRNA": 3, "tRNA": 3, "gene": 2, "misc_feature": 1, "unknown": 0}

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    nameconvert, annotypes, namevariants = loadnamevariants()
    args = getcliargs(None, nameconvert.keys(), annotypes)

    if args.showgenes:
        loadnamevariants(report=True)
        sys.exit(0)

    # Create output folder and subfolders
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    subfolders = {t: os.path.join(args.output, t) for t in ["CDS","rRNA","tRNA"]}
    for f in subfolders.values():
        os.makedirs(f, exist_ok=True)

    log = init_logger(args.output)
    log("=== Extraction run started ===")
    log(f"Input files: {args.genbank}")

    unrecgenes = defaultdict(list)
    warnings = []
    nrejected = 0
    outfh = {t:{} for t in ["CDS","rRNA","tRNA"]}

    filterlist = []
    if args.filter:
        with open(args.filter) as fh:
            filterlist = [line.strip() for line in fh]

    for gbpath in args.genbank:
        total_records = sum(1 for _ in SeqIO.parse(gbpath, "genbank"))
        processed = 0
        log(f"Processing file: {gbpath} ({total_records} records)")

        for seqrecord in SeqIO.parse(gbpath, "genbank"):
            processed += 1
            if processed % 1000 == 0 or processed == total_records:
                pct = (processed / total_records) * 100
                log(f"Progress: {processed}/{total_records} records ({pct:.1f}%)")

            if args.filter and seqrecord.name not in filterlist:
                nrejected += 1
                continue

            seqname = seqrecord.name
            outname = seqrecord.organism.replace(' ', '_') if args.organism and hasattr(seqrecord, "organism") else seqname

            foundgenes = {}       # gene -> sequence
            foundgenes_type = {}  # gene -> canonical type

            for feat in seqrecord.features:
                feat_type = feat.type
                name = get_feat_name(feat)

                if name not in nameconvert:
                    unrecgenes[name].append(seqname)
                    if not args.writeunknowns:
                        continue
                    stdname = name
                    anotype = "unknown"
                else:
                    stdname = nameconvert[name]
                    anotype = namevariants[stdname]['type']

                # Extract feature sequence
                featseq = feat.extract(seqrecord.seq)

                # Apply codon_start frame only for CDS
                frame_applied = False
                if anotype == "CDS" and args.keepframe and 'codon_start' in feat.qualifiers:
                    codon_start = int(feat.qualifiers['codon_start'][0])
                    featseq = featseq[(codon_start - 1):]
                    frame_applied = True

                # ---- Priority-based storage (applied to all types) ----
                existing_type = foundgenes_type.get(stdname)
                existing_seq = foundgenes.get(stdname)

                # Adjust priority if frame was applied to CDS
                current_priority = PRIORITY.get(anotype, 0)
                if frame_applied:
                    current_priority += 1  # ensure frame-adjusted CDS always replaces unadjusted

                existing_priority = PRIORITY.get(existing_type, 0) if existing_type else -1

                # Store sequence if higher priority or first time
                if existing_type is None or current_priority > existing_priority:
                    foundgenes[stdname] = str(featseq)
                    foundgenes_type[stdname] = anotype
                elif current_priority == existing_priority:
                    # Warn only if sequences differ
                    if str(featseq) != existing_seq:
                        msg = f"{seqname} has multiple distinct {anotype} annotations of {stdname}"
                        warnings.append(msg)
                        log(f"  Warning: {msg}")



            # Write output
            if len(foundgenes) >= args.mingenes:
                if not args.reqgenes or all(g in foundgenes for g in args.reqgenes):
                    for gene, seq in foundgenes.items():
                        typ = foundgenes_type[gene]
                        if typ not in ["CDS","rRNA","tRNA"]:
                            typ = "CDS"  # fallback
                        if gene not in outfh[typ]:
                            outfh[typ][gene] = open(os.path.join(subfolders[typ], f"{gene}.fasta"), 'w')
                        outfh[typ][gene].write(f">{outname}\n{seq}\n")
                else:
                    warnings.append(f"{seqname} missing required genes, skipped")
            else:
                warnings.append(f"{seqname} has too few annotated genes, skipped")

    # Close all files
    for subdict in outfh.values():
        for fh in subdict.values():
            fh.close()

    print("\n" + "="*40)
    print("Warnings Summary")
    print("="*40 + "\n")

    if unrecgenes:
        print("Unrecognized genes present:")
        for gene, entries in unrecgenes.items():
            print(f"  {gene} - {', '.join(sorted(set(entries)))}")

    if nrejected > 0:
        print(f"{nrejected} entries rejected by filter list")

    if warnings:
        print("\nOther warnings during processing:")
        for w in warnings:
            print(f"  {w}")
    else:
        print("No warnings.")

    log("=== Extraction run finished ===")
