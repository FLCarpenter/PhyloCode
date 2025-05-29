# -*- coding: utf-8 -*-
import argparse
import re
import os
import pandas as pd
from Bio import SeqIO

def read_partition_file(partition_file):
    partition_info = {}
    with open(partition_file, 'r') as pf:
        for line in pf:
            match = re.match(r'(\S+)\s+=\s+(\d+)-(\d+)', line.strip())
            if match:
                gene_name = match.group(1)
                start_pos = int(match.group(2)) - 1  # 0-based indexing
                end_pos = int(match.group(3))       # end inclusive for Python slicing
                partition_info[gene_name] = (start_pos, end_pos)
    return partition_info

def read_fasta(supermatrix_file):
    taxa = []
    sequences = []
    with open(supermatrix_file, 'r') as sm:
        for record in SeqIO.parse(sm, "fasta"):
            taxa.append(record.id)
            sequences.append(str(record.seq))
    return taxa, sequences

def gene_presence_matrix(partition_info, taxa, sequences):
    matrix = []

    for taxon, seq in zip(taxa, sequences):
        row = {'Taxon': taxon}
        for gene_name, (start, end) in partition_info.items():
            gene_seq = seq[start:end]
            gene_present = any(base.upper() not in ['-', 'N', '?'] for base in gene_seq)
            base_gene_name = os.path.splitext(os.path.basename(gene_name))[0]
            row[base_gene_name] = gene_present
        matrix.append(row)

    return pd.DataFrame(matrix)

def main():
    parser = argparse.ArgumentParser(description="Create gene presence/absence matrix from a supermatrix and partition file.")
    parser.add_argument('partition_file', help="Partition file with gene regions.")
    parser.add_argument('supermatrix_file', help="FASTA file with the supermatrix alignment.")
    parser.add_argument('--output_csv', help="Optional: Output CSV filename.", default=None)

    args = parser.parse_args()

    partition_info = read_partition_file(args.partition_file)
    taxa, sequences = read_fasta(args.supermatrix_file)
    df = gene_presence_matrix(partition_info, taxa, sequences)

    if args.output_csv:
        df.to_csv(args.output_csv, index=False)
    else:
        print(df)

if __name__ == "__main__":
    main()

