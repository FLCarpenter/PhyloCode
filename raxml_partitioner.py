import re
import sys

def parse_nexus(file_path):
    charsets = {}
    charpartition = {}

    with open(file_path, 'r') as file:
        nexus_content = file.read()

    # Extract charsets
    charset_matches = re.finditer(r'charset (\w+) = (.+?);', nexus_content, re.DOTALL)
    for match in charset_matches:
        charset_name = match.group(1)
        ranges = match.group(2).split()
        charsets[charset_name] = ranges

    # Extract charpartition with whole model name
    charpartition_match = re.search(r'charpartition mymodels =(.+?);', nexus_content, re.DOTALL)
    if charpartition_match:
        charpartition_text = charpartition_match.group(1)
        # Adjusted pattern to capture optional semicolon
        partition_matches = re.finditer(r'(\w+\+\w+\+\w+): (\w+),?', charpartition_text)
        for match in partition_matches:
            model_name = match.group(1)
            charset_name = match.group(2)
            if model_name not in charpartition:
                charpartition[model_name] = []
            charpartition[model_name].append(charset_name)

    return charsets, charpartition

def format_output(charsets, charpartition):
    output_lines = []

    for model_name, charset_names in charpartition.items():
        for charset_name in charset_names:
            formatted_ranges = ', '.join(charsets.get(charset_name, []))
            output_lines.append(f"{model_name}, {charset_name} = {formatted_ranges}")

    return output_lines

def write_to_file(output_lines, output_file):
    with open(output_file, 'w') as file:
        file.write('\n'.join(output_lines))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python raxml_partitioner.py input.nex output.txt")
        sys.exit(1)

    input_file_path = sys.argv[1]
    output_file_path = sys.argv[2]

    charsets, charpartition = parse_nexus(input_file_path)
    output_lines = format_output(charsets, charpartition)
    write_to_file(output_lines, output_file_path)

    print(f"Output written to {output_file_path}")
