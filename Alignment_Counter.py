import sys
import argparse

def count_nucleotides(sequence):
    """Counts occurrences of A, T, G, C in a sequence."""
    counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    for base in sequence:
        if base in counts:
            counts[base] += 1
    return counts

def calculate_gap_percentage(sequence):
    """Calculates the percentage of gaps ('-') in a sequence."""
    total_length = len(sequence)
    gap_count = sequence.count('-')
    return (gap_count / total_length) * 100 if total_length > 0 else 0

def read_fasta(file_path):
    """Reads sequences from a FASTA file."""
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
    """Main function to process command-line arguments and execute tasks."""
    parser = argparse.ArgumentParser(description="Analyze FASTA sequences.")
    parser.add_argument("fasta_file", help="Path to the FASTA file.")
    parser.add_argument("-g", "--gaps", action="store_true", help="Calculate gap percentage instead of nucleotide counts.")
    parser.add_argument("-a", "--alignment", action="store_true", help="Compute counts or gap percentage for the entire alignment instead of per sequence.")

    args = parser.parse_args()

    sequences = read_fasta(args.fasta_file)

    if args.alignment:
        # Compute totals for the whole alignment
        total_counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        total_length = 0
        total_gaps = 0

        for seq in sequences.values():
            total_length += len(seq)
            total_gaps += seq.count('-')
            counts = count_nucleotides(seq)
            for base in total_counts:
                total_counts[base] += counts[base]

        print("Results for entire alignment:")
        if args.gaps:
            gap_percentage = (total_gaps / total_length) * 100 if total_length > 0 else 0
            print(f"Gap Percentage: {gap_percentage:.2f}%")
        else:
            for base, count in total_counts.items():
                print(f"{base}: {count}")

    else:
        # Compute per-sequence results
        for seq_name, seq in sequences.items():
            print(f"Results for sequence '{seq_name}':")

            if args.gaps:
                gap_percentage = calculate_gap_percentage(seq)
                print(f"Gap Percentage: {gap_percentage:.2f}%")
            else:
                counts = count_nucleotides(seq)
                for base, count in counts.items():
                    print(f"{base}: {count}")

            print()

if __name__ == "__main__":
    main()

