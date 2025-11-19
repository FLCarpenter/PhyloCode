#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Extracts CDS, rRNA, and tRNA sequences from GenBank-format files,
grouping them by standardized gene names."""

# Imports
import argparse
import textwrap as _textwrap
import urllib
import os
import sys
from collections import defaultdict
from Bio import SeqIO

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
# Global variables
# -----------------------------------------------------------------------------
version = '1.1-integrated'

# -----------------------------------------------------------------------------
# Gene variant definitions (local fallback)
# -----------------------------------------------------------------------------
gene_aliases = {
    "ATP6": ["ATP6","ATP SYNTHASE F0 SUBUNIT 6","APT6","ATP SYNTHASE A0 SUBUNIT 6",
             "ATP SYNTHASE SUBUNIT 6","ATP SYNTHASE FO SUBUNIT 6","ATPASE6","ATPASE SUBUNIT 6"],
    "ATP8": ["ATP8","ATP SYNTHASE F0 SUBUNIT 8","APT8","ATP SYNTHASE A0 SUBUNIT 8",
             "ATP SYNTHASE SUBUNIT 8","ATP SYNTHASE FO SUBUNIT 8","ATPASE8","ATPASE SUBUNIT 8"],
    "COX1": ["COX1","CYTOCHROME C OXIDASE SUBUNIT 1","CYTOCHROME OXIDASE SUBUNIT I",
             "CYTOCHROME C OXIDASE SUBUNIT I","COXI","CO1","COI","CYTOCHROME COXIDASE SUBUNIT I",
             "CYTOCHROME OXIDASE SUBUNIT 1","CYTOCHROME OXYDASE SUBUNIT 1"],
    "COX2": ["COX2","CYTOCHROME C OXIDASE SUBUNIT 2","CYTOCHROME OXIDASE SUBUNIT II",
             "CYTOCHROME C OXIDASE SUBUNIT II","COXII","CO2","COII","CYTOCHROME COXIDASE SUBUNIT II",
             "CYTOCHROME OXIDASE SUBUNIT 2","CYTOCHROME OXYDASE SUBUNIT 2"],
    "COX3": ["COX3","CYTOCHROME C OXIDASE SUBUNIT 3","CYTOCHROME OXIDASE SUBUNIT III",
             "CYTOCHROME C OXIDASE SUBUNIT III","COXIII","CO3","COIII","CYTOCHROME COXIDASE SUBUNIT III",
             "CYTOCHROME OXIDASE SUBUNIT 3","CYTOCHROME OXYDASE SUBUNIT 3"],
    "CYTB": ["CYTB","CYTOCHROME B","CYB","COB","COB / CYTB"],
    "ND1":  ["ND1","NAD1","NSD1","NADH1","NADH DEHYDROGENASE SUBUNIT I",
             "NADH DEHYDROGENASE SUBUNIT 1","NADH DESHYDROGENASE SUBUNIT 1","NAD1-0"],
    "ND2":  ["ND2","NAD2","NSD2","NADH2","NADH DEHYDROGENASE SUBUNIT II",
             "NADH DEHYDROGENASE SUBUNIT 2","NADH DESHYDROGENASE SUBUNIT 2","NAD2-0"],
    "ND3":  ["ND3","NAD3","NSD3","NADH3","NADH DEHYDROGENASE SUBUNIT III",
             "NADH DEHYDROGENASE SUBUNIT 3","NADH DESHYDROGENASE SUBUNIT 3","NAD3-0"],
    "ND4":  ["ND4","NAD4","NSD4","NADH4","NADH DEHYDROGENASE SUBUNIT IV",
             "NADH DEHYDROGENASE SUBUNIT 4","NADH DESHYDROGENASE SUBUNIT 4","NAD4-0"],
    "ND4L": ["ND4L","NAD4L","NSD4L","NADH4L","NADH DEHYDROGENASE SUBUNIT IVL",
             "NADH DEHYDROGENASE SUBUNIT 4L","NADH DESHYDROGENASE SUBUNIT 4L","NAD4L-0"],
    "ND5":  ["ND5","NAD5","NSD5","NADH5","NADH DEHYDROGENASE SUBUNIT V",
             "NADH DEHYDROGENASE SUBUNIT 5","NADH DESHYDROGENASE SUBUNIT 5","NAD5-0"],
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
# Function: loadnamevariants
# -----------------------------------------------------------------------------
def loadnamevariants(report=False):
    """Load or define gene name variants, distinguishing CDS vs rRNA vs tRNA."""
    output = {}
    fullparse = {}
    alltypes = set()

    url = "https://raw.githubusercontent.com/tjcreedy/constants/master/gene_name_variants.txt"
    success = False

    try:
        for line in urllib.request.urlopen(url):
            line = line.decode('utf-8').strip()
            description, variants = line.split(":")
            name, annotype, fullname = description.split(";")
            variants = variants.split(',')
            variants.extend([name, fullname.upper()])
            fullvariants = []
            for v in variants:
                for g in ['', ' ']:
                    v = v.replace(g, '')
                    for s in ['', ' GENE', ' ' + annotype.upper()]:
                        var = v + s
                        fullvariants.append(var)
                        output[v + s] = name
            alltypes.add(annotype)
            fullparse[name] = {'type': annotype, 'variants': fullvariants}
            if report:
                fullvariants = [v.replace(' ', '\u00A0') if len(v) < 12 else v for v in fullvariants]
                fullvariants = _textwrap.fill(', '.join(fullvariants), width=80,
                                              initial_indent='\t', subsequent_indent='\t')
                print(f"Standard name = {name}, type = {annotype}, full name = {fullname}:\n{fullvariants}")
        success = True
    except Exception as e:
        sys.stderr.write(f"Warning: could not load gene_name_variants.txt from GitHub ({e}). Using local aliases.\n")

    if not success:
        # fallback from the embedded dictionary
        for name, variants in gene_aliases.items():
            if "rRNA" in name.upper():
                annotype = "rRNA"
            elif "TRNA" in name.upper():
                annotype = "tRNA"
            else:
                annotype = "CDS"
            alltypes.add(annotype)

            fullvariants = []
            for v in variants:
                v = v.upper()
                for s in ['', ' ' + annotype.upper(), ' GENE']:
                    fullvariants.append(v + s)
                    output[v + s] = name
            fullparse[name] = {'type': annotype, 'variants': fullvariants}
            if report:
                fullvariants = _textwrap.fill(', '.join(fullvariants), width=80,
                                              initial_indent='\t', subsequent_indent='\t')
                print(f"Standard name = {name}, type = {annotype}:\n{fullvariants}")

    if not report:
        return output, alltypes, fullparse

# -----------------------------------------------------------------------------
# Function: get_feat_name
# -----------------------------------------------------------------------------
def get_feat_name(feat):
    featname = "unknown"
    nametags = ['gene', 'product', 'label', 'standard_name']
    for t in nametags:
        if t in feat.qualifiers:
            featname = feat.qualifiers[t][0].upper()
            break
    return featname

# -----------------------------------------------------------------------------
# Function: getcliargs
# -----------------------------------------------------------------------------
def getcliargs(arglist=None, knowngenes=None, knowntypes=None):
    parser = argparse.ArgumentParser(
        description="""
    This script finds and extracts sequences corresponding to CDS, rRNA and/or tRNA annotations 
    from a set of sequences from one or more genbank-format flat files. Specify which regions to 
    extract using the -r/--regiontypes argument (default is CDS)
    |n
    The region name is identified based on the /gene, /label or /product tag in the genbank-format 
    flat file. Two methods are used to ensure that naming variants are correctly identified as the 
    same gene:
        1. the gene name is converted to uppercase (so that "atp8" is the same as "ATP8")
        2. the script removes semicolons, underscores, hyphens or spaces, and any characters 
           following these, from gene names ("ATP8-0" -> "ATP8")
    |n
    Optionally, you can print the known gene name variants with -s/--showgenes.
    """,
        formatter_class=MultilineFormatter)

    parser.add_argument("-g", "--genbank", type=str, metavar='PATH', required=True, nargs='+')
    parser.add_argument("-o", "--output", type=str, metavar='PATH', required=True)
    parser.add_argument("-m", "--mingenes", type=int, metavar='N', default=0)
    parser.add_argument("-q", "--reqgenes", type=str, metavar='GENE', nargs='*',
                        choices=knowngenes)
    parser.add_argument("-r", "--genetypes", type=str, metavar='TYPE', nargs='+',
                        choices=knowntypes, default=list(knowntypes))
    parser.add_argument("-f", "--filter", type=str, metavar='PATH')
    parser.add_argument("-n", "--organism", action='store_true')
    parser.add_argument("-w", "--writeunknowns", action='store_true')
    parser.add_argument("-k", "--keepframe", action='store_true')
    parser.add_argument("-p", "--presence", type=str, metavar='PATH')
    parser.add_argument("-s", "--showgenes", action='store_true')
    parser.add_argument("-v", "--version", action='store_true')

    args = parser.parse_args(arglist) if arglist else parser.parse_args()

    if args.version:
        print(version)
        sys.exit(0)
    if args.mingenes < 0:
        parser.error(f"{args.mingenes} is not a valid minimum gene count")

    return args

# -----------------------------------------------------------------------------
# Main script
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    nameconvert, annotypes, namevariants = loadnamevariants()
    args = getcliargs(None, nameconvert.keys(), annotypes)

    if args.showgenes:
        loadnamevariants(report=True)
        sys.exit(0)

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    unrecgenes = defaultdict(list)

    if args.presence:
        args.presence = open(args.presence, 'w')

    filterlist = []
    if args.filter:
        with open(args.filter) as fh:
            filterlist = [line.strip() for line in fh]

    nrejected = 0
    outfh = {}

    for gbpath in args.genbank:
        for seqrecord in SeqIO.parse(gbpath, "genbank"):
            if args.filter and seqrecord.name not in filterlist:
                nrejected += 1
                continue

            seqname = seqrecord.name
            outname = seqrecord.organism.replace(' ', '_') if args.organism and hasattr(seqrecord, "organism") else seqname

            foundgenes = defaultdict(list)

            for feat in seqrecord.features:
                if feat.type not in args.genetypes:
                    continue

                name = get_feat_name(feat)
                if name in nameconvert:
                    stdname = nameconvert[name]
                else:
                    unrecgenes[name].append(seqname)
                    if args.writeunknowns:
                        stdname = name
                    else:
                        continue

                featsequence = feat.extract(seqrecord.seq)
                if args.keepframe and 'codon_start' in feat.qualifiers:
                    featsequence = featsequence[(int(feat.qualifiers['codon_start'][0]) - 1):]
                foundgenes[stdname].append(featsequence)

            if len(foundgenes) >= args.mingenes:
                if not args.reqgenes or all(g in foundgenes for g in args.reqgenes):
                    for gene, seqs in foundgenes.items():
                        if len(seqs) > 1:
                            sys.stderr.write(
                                f"Warning: {seqname} has multiple annotations of {gene}\n"
                            )
                        if gene not in outfh:
                            outfh[gene] = open(os.path.join(args.output, f"{gene}.fasta"), 'w')
                        for seq in seqs:
                            outfh[gene].write(f">{outname}\n{seq}\n")
                    if args.presence:
                        args.presence.write(f"{outname},{','.join(foundgenes.keys())}\n")
                else:
                    sys.stderr.write(f"Warning: {seqname} missing required genes, skipped\n")
            else:
                sys.stderr.write(f"Warning: {seqname} has too few annotated genes, skipped\n")

    for fh in outfh.values():
        fh.close()
    if args.presence:
        args.presence.close()

    if unrecgenes:
        sys.stderr.write("Warning: unrecognised genes present in these entries:\n")
        for gene, entries in unrecgenes.items():
            sys.stderr.write(f"\t{gene} - {', '.join(sorted(set(entries)))}\n")

    if nrejected > 0:
        sys.stderr.write(f"Warning: {nrejected} entries rejected by filter list\n")
