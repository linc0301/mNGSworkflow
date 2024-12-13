#!/usr/bin/env python3
import argparse

def process_quality(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # dealing with the 4th line
            if line[0] == '+':
                quality_line = infile.readline().strip()
                new_quality = ''.join(['I' if q != 'I' else q for q in quality_line])
                outfile.write(line)
                outfile.write(new_quality + '\n')
            else:
                outfile.write(line)

def main():
    parser = argparse.ArgumentParser(description="Process quality values in a FASTQ file.")
    parser.add_argument('-i', '--input', required=True, help="Input FASTQ file.")
    parser.add_argument('-o', '--output', required=True, help="Output FASTQ file.")
    
    args = parser.parse_args()
    process_quality(args.input, args.output)

if __name__ == "__main__":
    main()