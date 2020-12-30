#!/usr/bin/env python3

# GFF Feature finder
#
# Author: Tyler Aprati
#
# This program finds genes between a coordinate range within a gff file
# Example usage:
# featurefinder.py -i input.gff -o chromosome -s start_index -e end_index
#

import sys, argparse
from argparse import RawTextHelpFormatter

print()

parser = argparse.ArgumentParser(description=
                                 'This program is written for execution in the Python (version 3) language.\n'+
                                 'This program filters a gff file by criteria.\n',
                                 formatter_class=RawTextHelpFormatter)

requiredParam = parser.add_argument_group('required parameters')
requiredParam.add_argument('-i', type=str, metavar='infile', required=True, help='Name of input file')
requiredParam.add_argument('-o', type=str, metavar='outfile', required=True, help='Name for output file')
requiredParam.add_argument('-c', type=str, metavar='chromosome', required=True, help='Chromosome of interest')
requiredParam.add_argument('-t', type=str, metavar='feature_type', required=True, help='Feature Type')
requiredParam.add_argument('-s', type=int, metavar='start_index', required=True, help='Start index for search range')
requiredParam.add_argument('-e', type=int, metavar='end_index', required=True, help='End index for search range')

args = parser.parse_args()

# Load in other inputs
filename = args.i
output = args.o
chrom = args.c
type = args.t
s_index = args.s
e_index = args.e

# Open input and output files
file = open(filename,'r')
output = open(output,'w')

# Print update
print("Searching through",filename)

# Find features
count = 0
while True:
    l = file.readline()
    if l == '':
        break
    line = l.split("\t")
    if line[0] == chrom:
        if line[2] == type:
            if (int(line[3]) > s_index) and (int(line[3]) < e_index):
                output.write(l)
                count += 1

print(count,"Features Found")

output.close()

