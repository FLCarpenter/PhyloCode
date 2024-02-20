# -*- coding: utf-8 -*-


from sys import argv

valid_nucleotides = set("ACGTURYKMSWBDHVN?-")
amino_acids = set('ABCDEFGHIKLMNPQRSTUVWYZX?-')

def Is_Valid_Sequence(sequence):
#    Returns True  if the sequence is composed of Nucleotides only
    sequence = set(sequence.upper())
    if all([base in valid_nucleotides for base in sequence]):
        return True
    elif all([base in amino_acids for base in sequence]):
        return False
#        print ("The input file contains amino acid sequences.\n")
    else:
        return False
#       print  ("The input file contains invalid sequences")

def RY_Recoding_N(Sequence):
    recoded_seq = ""
    for base in Sequence:
        if base in 'AGagRr':
            recoded_seq += 'R'
        elif base in 'TCtcYy':
            recoded_seq += 'Y'
        else:
            recoded_seq += '-'
    return recoded_seq

def RY_Recoding_Codon_NT(Sequence, Position):
    Pos = int(Position)
    Index = 1 + (3 - Pos)
    recoded_seq = ''
    for base in Sequence:
        if (Index + 3) % 3 == 0:
            if base in 'AGagRr':
                recoded_seq += 'R'
            elif base in 'TCtcYy':
                recoded_seq += 'Y'
            else:
                recoded_seq += '-'
        else:
            recoded_seq += base
        Index += 1
    return recoded_seq

def RY_Recoding_Codon_NT3R(Sequence, Position):
    Pos = int(Position)
    Index = 1 + (2 - Pos)
    recoded_seq = ''
    for base in Sequence:
        if (Index + 2) % 2 == 0:
            if base in 'AGagRr':
                recoded_seq += 'R'
            elif base in 'TCtcYy':
                recoded_seq += 'Y'
            else:
                recoded_seq += '-'
        else:
            recoded_seq += base
        Index += 1
    return recoded_seq

def print_usage():
    print("Usage: python RY_Recoder.py input.fasta datatype position")
    print("  - input.fasta: Input FASTA file containing nucleotide sequences.")
    print("  - datatype: Specify if the input data is nucleotide all 3 codons (nt), or nucleotide with 3rd position removed (nt3r).")
    print("  - position: Specify the codon position to recode (1, 2, 3) or 'N' to recode the whole sequence.")
    
if __name__ == "__main__":
    if len(argv) != 4:
        print("Error: Incorrect number of arguments.")
        print_usage()
    else:
        Script = argv[0]
        Codon = argv[3]
        Targets = [argv[1]]
        Option = argv[2]

        if Codon not in ['1', '2', '3', 'N']:
            print("Error: Invalid position. Use 'N', '1', '2', or '3'.")
            print_usage()
        else:
            for File in Targets:
                with open(File, 'r') as F:
                    if Codon == 'N':
                        OutFName = f'{File.split(".")[0]}_RY.{File.split(".")[1]}'
                    else:
                        OutFName = f'{File.split(".")[0]}_Position{Codon}_RY.{File.split(".")[1]}'
                    Out = open(OutFName, 'w')
                    for Line in F:
                        Line = Line.strip('\n')
                        if Line.startswith('>'):
                            Out.write(Line + '\n')
                        else:
                            Line = Line.upper()
                            if Is_Valid_Sequence(Line):
                                if Option == 'nt':
                                    if Codon == 'N':
                                        Out.write(RY_Recoding_N(Line) + '\n')
                                    else:
                                        Out.write(RY_Recoding_Codon_NT(Line, Codon) + '\n')
                                elif Option == 'nt3r': 
                                    if Codon == 'N':
                                        Out.write(RY_Recoding_N(Line) + '\n')
                                    else:
                                        Out.write(RY_Recoding_Codon_NT3R(Line, Codon) + '\n')
                            else:
                                print("Warning: Record skipped with non-nucleotide characters.")
