# -*- coding: utf-8 -*-
import argparse
import re
from Bio import SeqIO
import os

def read_partition_file(partition_file):
    """
    Reads the partition file and returns a dictionary where keys are gene names
    and values are tuples of start and end column indices (1-based).
    """
    partition_info = {}
    with open(partition_file, 'r') as pf:
        for line in pf:
            # Match lines that look like:
            # 2_nt_renamed_for_supermatrix/12S.fasta = 1-1699
            match = re.match(r'(\S+)\s+=\s+(\d+)-(\d+)', line.strip())
            if match:
                gene_name = match.group(1)
                start_pos = int(match.group(2)) - 1  # Convert to 0-based indexing
                end_pos = int(match.group(3)) - 1    # Convert to 0-based indexing
                partition_info[gene_name] = (start_pos, end_pos)
    return partition_info

def read_fasta(supermatrix_file):
    """
    Reads the FASTA supermatrix file and returns the taxa names and the sequence matrix.
    """
    taxa = []
    sequences = []
    with open(supermatrix_file, 'r') as sm:
        # Use BioPython to parse the FASTA file
        for record in SeqIO.parse(sm, "fasta"):
            taxa.append(record.id)  # Taxon name is the record id
            sequences.append(str(record.seq))  # Sequence as a string

    return taxa, sequences

def split_alignment(partition_info, taxa, sequences):
    """
    Splits the alignment based on the partition file and writes each gene alignment
    into separate FASTA files. Ensures the output directory exists.
    """
    # Create the output directory for all the gene files
    output_dir = "split_alignment"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for gene_name, (start_pos, end_pos) in partition_info.items():
        # Extract the columns corresponding to the current gene
        gene_alignment = []
        for seq in sequences:
            gene_alignment.append(seq[start_pos:end_pos+1])  # Slice the gene region

        # Extract the base gene name without path and extension
        base_gene_name = os.path.splitext(os.path.basename(gene_name))[0]  # Get the base name without extension

        # Set the output filename to be stored in the split_alignment directory
        output_filename = os.path.join(output_dir, f'{base_gene_name}_split.fasta')

        # Write the gene-specific alignment to a FASTA file
        with open(output_filename, 'w') as out_file:
            for taxon, seq in zip(taxa, gene_alignment):
                out_file.write(f">{taxon}\n")
                out_file.write(f"{seq}\n")

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Split a FASTA supermatrix based on a partition file.")
    parser.add_argument('partition_file', help="The partition file that defines gene regions.")
    parser.add_argument('supermatrix_file', help="The FASTA file containing the supermatrix alignment.")

    # Parse the command line arguments
    args = parser.parse_args()

    # Step 1: Read the partition file
    partition_info = read_partition_file(args.partition_file)

    # Step 2: Read the FASTA supermatrix file
    taxa, sequences = read_fasta(args.supermatrix_file)

    # Step 3: Split the alignment and output gene-specific alignments
    split_alignment(partition_info, taxa, sequences)

if __name__ == "__main__":
    main()

