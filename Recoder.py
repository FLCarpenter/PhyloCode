# -*- coding: utf-8 -*-
"""
Created on Tues Apr 14 14:59:18 2023

@author: Fiona
"""

import argparse

valid_nucleotides = set("ACGTURYKMSWBDHVN?-")
amino_acids = set('ABCDEFGHIKLMNPQRSTUVWYZX?-')


def Is_Valid_Sequence(sequence):
    sequence = set(sequence.upper())
    if all([base in valid_nucleotides for base in sequence]):
        return True
    elif all([base in amino_acids for base in sequence]):
        return False
    else:
        return False


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


def Binary_Recoding_N(Sequence):
    recoded_seq = ""
    for base in Sequence:
        if base in 'AGag':
            recoded_seq += '0'
        elif base in 'TCtc':
            recoded_seq += '1'
        else:
            recoded_seq += '-'
    return recoded_seq


def Binary_Recoding_Codon_NT(Sequence, Position):
    Pos = int(Position)
    Index = 1 + (3 - Pos)
    recoded_seq = ''
    for base in Sequence:
        if (Index + 3) % 3 == 0:
            if base in 'AGag':
                recoded_seq += '0'
            elif base in 'TCtc':
                recoded_seq += '1'
            else:
                recoded_seq += '-'
        else:
            recoded_seq += base
        Index += 1
    return recoded_seq


def Binary_Recoding_Codon_NT3R(Sequence, Position):
    Pos = int(Position)
    Index = 1 + (2 - Pos)
    recoded_seq = ''
    for base in Sequence:
        if (Index + 2) % 2 == 0:
            if base in 'AGag':
                recoded_seq += '0'
            elif base in 'TCtc':
                recoded_seq += '1'
            else:
                recoded_seq += '-'
        else:
            recoded_seq += base
        Index += 1
    return recoded_seq


# ---------------- CLI ---------------- #

def parse_args():
    parser = argparse.ArgumentParser(
        description="Recodes nucleotide sequences into RY or Binary formats."
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input FASTA file"
    )

    parser.add_argument(
        "-d", "--datatype",
        required=True,
        choices=["nt", "nt3r"],
        help="Input datatype: nt or nt3r"
    )

    parser.add_argument(
        "-p", "--position",
        required=True,
        choices=["1", "2", "3", "N"],
        help="Codon position: 1, 2, 3, or N"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        choices=["RY", "Binary"],
        help="Output format: RY or Binary"
    )

    return parser.parse_args()


# ---------------- MAIN ---------------- #

if __name__ == "__main__":
    args = parse_args()

    File = args.input
    Option = args.datatype
    Codon = args.position
    OutData = args.output

    with open(File, 'r') as F:
        if OutData == 'RY':
            if Codon == 'N':
                OutFName = f'{File.split(".")[0]}_RY.{File.split(".")[1]}'
            else:
                OutFName = f'{File.split(".")[0]}_Position{Codon}_RY.{File.split(".")[1]}'
        elif OutData == 'Binary':
            if Codon == 'N':
                OutFName = f'{File.split(".")[0]}_Binary.{File.split(".")[1]}'
            else:
                OutFName = f'{File.split(".")[0]}_Position{Codon}_Binary.{File.split(".")[1]}'

        Out = open(OutFName, 'w')

        for Line in F:
            Line = Line.strip('\n')

            if Line.startswith('>'):
                Out.write(Line + '\n')
            else:
                Line = Line.upper()

                if Is_Valid_Sequence(Line):
                    if OutData == 'RY':
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

                    elif OutData == 'Binary':
                        if Option == 'nt':
                            if Codon == 'N':
                                Out.write(Binary_Recoding_N(Line) + '\n')
                            else:
                                Out.write(Binary_Recoding_Codon_NT(Line, Codon) + '\n')
                        elif Option == 'nt3r':
                            if Codon == 'N':
                                Out.write(Binary_Recoding_N(Line) + '\n')
                            else:
                                Out.write(Binary_Recoding_Codon_NT3R(Line, Codon) + '\n')

                else:
                    print("Warning: Record skipped with non-nucleotide characters.")

        Out.close()
