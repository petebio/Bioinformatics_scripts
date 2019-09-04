#!/usr/bin/env python
import argparse
import sys

#################################################################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('outfile', type = str, help = 'Output file')
parser.add_argument('-a', type = str, default = 'Sample_A', help = 'Sample ID for sample A, Default = Sample_A')
parser.add_argument('-b', type = str, default = 'Sample_B', help = 'Sample ID for sample B, Default = Sample_B')

args = parser.parse_args()

#################################################################################################################################

header = 'Chr\tStart\tEnd\tClosest_Gene\tDist_to_TSS\t' + args.a + '\t' + args.b

out = open(args.outfile, 'w')
out.write(header + '\n')

for num,line in enumerate(sys.stdin):
	if num == 0:
		continue

	fields = line.strip().split('\t')
	
	chrom = fields[1]
	start = str(int(fields[2]) - 1)
	end = fields[3]

	gene_id = fields[15]

	if gene_id == '':
		gene_id = 'No_Gene'

	dist_to_tss = fields[9].lstrip('-')

	if dist_to_tss == '':
		dist_to_tss = 'NA'

	tag_count_a = fields[19]
	if tag_count_a == '':
		tag_count_a = '0'

	if len(fields) == 21:
		tag_count_b = fields[20]
	else:
		tag_count_b = '0'

	output = '\t'.join([chrom, start, end, gene_id, dist_to_tss, tag_count_a, tag_count_b])
	out.write(output + '\n')

out.close()

#################################################################################################################################
