# -*- coding: utf-8 -*-

import sys
import re

def calculate_first_element(gene_name, start_number):
    if gene_name.endswith('_1') and start_number != 1:
        return "{}".format(((start_number - 1) // 3) * 2 + 1)
    elif gene_name.endswith('_2') and start_number != 2:
        return "{}".format(((start_number - 2) // 3) * 2 + 2)
    else:
        # Handle other cases if needed
        return "{}".format(start_number)

def calculate_second_part(end_number):
    return (end_number // 3) * 2

#Gene+Codon

def process_line_gene_codon(line):
    match = re.match(r'\s*charset (\w+) = (\d+)-(\d+)\\(\d+);', line)
    if match:
        gene_name = match.group(1)
        start_number = int(match.group(2))
        end_number = int(match.group(3))
        base_number = int(match.group(4))

        modified_start = calculate_first_element(gene_name, start_number)
        modified_second = calculate_second_part(end_number)

        # Use 2 as the base number in the modified line
        modified_line = "charset {} = {}-{}\\2;".format(gene_name, modified_start, modified_second)
        return modified_line
    else:
        return line

#Strand+Codon

def process_line_strand_codon(line):
    match = re.match(r'\s*charset (\w+) = (.+?);', line)
    if match:
        gene_name = match.group(1)
        positions = match.group(2).split()

        modified_positions = []
        for position in positions:
            start_number, end_number, base_number = map(int, re.split(r'[-\\]', position))
            modified_start = calculate_first_element(gene_name, start_number)
            modified_end = calculate_second_part(end_number)

            # Use 2 as the base number in the modified positions
            modified_positions.append("{}-{}\\2".format(modified_start, modified_end))

        modified_line = "charset {} = {};".format(gene_name, ' '.join(modified_positions))
        return modified_line
    else:
        return line

#Codon Only

def process_line_codon(line):
    match = re.match(r'\s*charset (\w+) = (.+?);', line)
    if match:
        gene_name = match.group(1)
        positions = match.group(2).split()

        modified_positions = []
        for position in positions:
            start_number, end_number, base_number = map(int, re.split(r'[-\\]', position))
            modified_start = calculate_first_element(gene_name, start_number)
            modified_end = calculate_second_part(end_number)

            # Use 2 as the base number in the modified positions
            modified_positions.append("{}-{}\\2".format(modified_start, modified_end))

        modified_line = "charset {} = {};".format(gene_name, ' '.join(modified_positions))
        return modified_line
    else:
        return line



if len(sys.argv) != 5 or sys.argv[1] != '-f':
    print("Usage: python Repartitioner.py -f <gc, dc, or c> input.txt output.txt")
    sys.exit(1)


format_option = sys.argv[2]
input_filename = sys.argv[3]
output_filename = sys.argv[4]

try:
    with open(input_filename, 'r') as file:
        input_data = file.read()
except FileNotFoundError:
    print(f"Error: Input file '{input_filename}' not found.")
    sys.exit(1)

# Process each line based on the specified format
if format_option == 'gc':
    process_line = process_line_gene_codon
elif format_option == 'dc':
    process_line = process_line_strand_codon
elif format_option == 'c':
    process_line = process_line_codon
else:
    print("Error: Invalid format option. Use gc, dc, or c.")
    sys.exit(1)

# Process each line
modified_lines = [process_line(line) for line in input_data.split('\n')]

# Join the modified lines back into a string
modified_input = '\n'.join(modified_lines)

# Write the modified content to the output file
with open(output_filename, 'w') as file:
    file.write(modified_input)

print(f"Output written to {output_filename}")
