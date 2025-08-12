#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sys import argv

# IUPAC nucleotide codes (DNA) + gap/missing symbols
valid_nucleotides = set("ACGTURYKMSWBDHVN?-")

def clean_sequence(seq):
    """Replace any character not in valid_nucleotides with 'N'."""
    return ''.join(base if base.upper() in valid_nucleotides else 'N' for base in seq)

def print_usage():
    print("Usage: ./replace_non_iupac.py input.fasta output.fasta")

if __name__ == "__main__":
    if len(argv) != 3:
        print_usage()
        exit(1)

    infile = argv[1]
    outfile = argv[2]

    with open(infile, 'r') as fin, open(outfile, 'w') as fout:
        for line in fin:
            line = line.strip('\n')
            if line.startswith('>'):
                fout.write(line + '\n')
            else:
                fout.write(clean_sequence(line) + '\n')

    print(f"Cleaned FASTA written to {outfile}")
