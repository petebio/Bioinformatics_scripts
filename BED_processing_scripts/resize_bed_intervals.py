#!/usr/bin/env python
import argparse
import math

###########################################################################################################################################

parser = argparse.ArgumentParser(description = 'Resize BED intervals to a given size')
parser.add_argument('infile', type = str, help = 'Input BED file')
parser.add_argument('outfile', type = str, help = 'Output BED file')
parser.add_argument('-l', '--length', dest = 'length', type = int, default = 400, help = 'Length of output features. Default = 400')

args = parser.parse_args()

###########################################################################################################################################

extend = int(math.ceil(args.length / 2))

input_bed = open(args.infile, 'r')
output_bed = open(args.outfile, 'w')

for line in input_bed:
	fields = line.strip().split('\t')
	mid_pos = int(math.ceil((float(fields[1]) + float(fields[2])) / 2))

	start = mid_pos - extend - 1
	end = mid_pos + extend

	if start < 0:
		continue

	fields[1] = str(start)
	fields[2] = str(end)

	output_bed.write('\t'.join(fields) + '\n')

input_bed.close()
output_bed.close()

###########################################################################################################################################
