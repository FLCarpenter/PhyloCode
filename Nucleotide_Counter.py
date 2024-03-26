# -*- coding: utf-8 -*-


def count_nucleotides(sequence):
    counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    for base in sequence:
        if base in counts:
            counts[base] += 1
    return counts

def read_fasta(file_path):
    sequences = {}
    current_sequence = None
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                current_sequence = line[1:]
                sequences[current_sequence] = ''
            else:
                sequences[current_sequence] += line
    return sequences

def main():
    fasta_file = "sequences.fasta"  # Change this to your FASTA file path
    sequences = read_fasta(fasta_file)
    
    for seq_name, seq in sequences.items():
        print(f"Counts for sequence '{seq_name}':")
        counts = count_nucleotides(seq)
        for base, count in counts.items():
            print(f"{base}: {count}")
        print()

if __name__ == "__main__":
    main()
