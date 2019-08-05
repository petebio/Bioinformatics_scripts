#!/usr/bin/env python
import pybedtools
import argparse
import math
import sys

#####################################################################################################################################################

parser = argparse.ArgumentParser(description = 'Creake a union of two BED summit files')
parser.add_argument('bed_file_a', type = str, help = 'BED file A')
parser.add_argument('bed_file_b', type = str, help = 'BED file B')
parser.add_argument('outfile', type = str, help = 'Output BED file')
parser.add_argument('-d', '--dist', dest = 'dist', type = int, default = 400, help = 'Maximum distance between summits to merge. Default = 400bp')

args = parser.parse_args()

#####################################################################################################################################################

def prepare_BED_file(bed_file):
	input_bed = open(bed_file, 'r')
	output_bed = list()

	extend = int(math.ceil(float(args.dist) / 2))

	for line in input_bed:
		fields = line.strip().split('\t')

		chrom = fields[0]
		mid_pos = int(math.ceil((float(fields[1]) + float(fields[2])) / 2))

		start = mid_pos - extend
		end = mid_pos + extend

		if start < 0:
			continue

		bed_line = '\t'.join([chrom, str(start), str(end), str(mid_pos)])
		output_bed.append(bed_line)

	output_bed = pybedtools.BedTool(output_bed)
	return output_bed

def extract_unique_sites(bed_file_x, bed_file_y):
	unique_sites = bed_file_x.intersect(bed_file_y, v = True)
	unique_site_list = list()

	for line in unique_sites:
		chrom = line[0]
		summit_pos = line[3]
		start = int(summit_pos) - 1

		bed_line = '\t'.join([chrom, str(start), summit_pos])
		unique_site_list.append(bed_line)

	return unique_site_list

#####################################################################################################################################################

bed_sites_a = prepare_BED_file(args.bed_file_a)
bed_sites_b = prepare_BED_file(args.bed_file_b)

sys.stderr.write('\n------------------------------\n')
sys.stderr.write('Read %d sites from ' % (bed_sites_a.count()) + args.bed_file_a + '\n')
sys.stderr.write('Read %d sites from ' % (bed_sites_b.count()) + args.bed_file_b + '\n')
sys.stderr.write('------------------------------\n')

common_bed_sites = bed_sites_a.intersect(bed_sites_b, wo = True)

sys.stderr.write('Merging %d sites\n' % (common_bed_sites.count()))
sys.stderr.write('------------------------------\n')

merged_sites = list()
for line in common_bed_sites:
	chrom = line[0]
	summit_pos_a = line[3]
	summit_pos_b = line[7]

	new_summit_pos = int(math.ceil((float(summit_pos_a) + float(summit_pos_b)) / 2))
	start = new_summit_pos - 1

	bed_line = '\t'.join([chrom, str(start), str(new_summit_pos)])
	merged_sites.append(bed_line)

unique_sites_a = extract_unique_sites(bed_sites_a, bed_sites_b)
unique_sites_b = extract_unique_sites(bed_sites_b, bed_sites_a)

sys.stderr.write('Found %d unique sites in ' % (len(unique_sites_a)) + args.bed_file_a + '\n')
sys.stderr.write('Found %d unique sites in ' % (len(unique_sites_b)) + args.bed_file_b + '\n')
sys.stderr.write('------------------------------\n')

combined_bed = merged_sites + unique_sites_a + unique_sites_b
combined_bed = pybedtools.BedTool(combined_bed).sort()

sys.stderr.write('Writing %d sites to ' % (combined_bed.count()) + args.outfile + '\n')
sys.stderr.write('------------------------------\n\n')

combined_bed.saveas(args.outfile)

#####################################################################################################################################################
