#!/usr/bin/env python
import argparse
import sys

#################################################################################################################################

parser = argparse.ArgumentParser(description = 'Separate peaks into Promoter/Distal elements based on Homer output.')
parser.add_argument('-d', '--dist', type = int, default = 1500, help = 'Promoter distance cut-off. Default = 1500bp')
parser.add_argument('-p', '--prefix', dest = 'prefix', type = str, required = False, help = 'Optional prefix for output files')

args = parser.parse_args()

#################################################################################################################################

promoter_out = 'Promoter.bed'
distal_out = 'Distal.bed'

if args.prefix:
	promoter_out = args.prefix + '_' + promoter_out
	distal_out = args.prefix + '_' + distal_out

promoter = open(promoter_out, 'w')
distal = open(distal_out, 'w')

for num,line in enumerate(sys.stdin):
	if num == 0:
		continue

	fields = line.strip().split('\t')

	chrom = fields[1]
	start = str(int(fields[2]) - 2)
	end = fields[3]
	gene_id = fields[15]

	if gene_id == '':
		gene_id = 'No_Gene'

	dist_to_tss = abs(int(fields[9]))

	bed_line = '\t'.join([chrom, start, end, str(dist_to_tss), gene_id])

	if dist_to_tss < args.dist:
		promoter.write(bed_line + '\n')
	else:
		distal.write(bed_line + '\n')

promoter.close()
distal.close()

#################################################################################################################################
