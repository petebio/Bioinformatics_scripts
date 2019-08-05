#!/usr/bin/env python
import pybedtools
import argparse
import numpy
import os

###############################################################################################################################################################

parser = argparse.ArgumentParser(description = 'Align a motif PWM matrix to a BED file using Homer, also removing duplicated (e.g. palindromic) motifs')
parser.add_argument('bed_file', type = str, help = 'BED file to align motifs to')
parser.add_argument('motif_file', type = str, help = 'Motif file to provide to Homer')
parser.add_argument('genome', type = str, help = 'Genome version to use with Homer')
parser.add_argument('outfile', type = str, help = 'Output file')

args = parser.parse_args()

###############################################################################################################################################################

tmp_dir = 'tmp_' + str(numpy.random.randint(1e6))
os.mkdir(tmp_dir)

homer_command = 'annotatePeaks.pl ' + args.bed_file + ' ' + args.genome + ' -noann '
homer_command = homer_command + ' -m ' + args.motif_file + ' -mbed ' + tmp_dir + '/Motif_raw.bed > ' + tmp_dir + '/Homer_stdout.txt'
os.system(homer_command)

bed_input = pybedtools.BedTool(tmp_dir + '/Motif_raw.bed').sort()
bed_output = open(args.outfile, 'w')

prev_line = list()

for line in bed_input:
	is_duplicate = False

	chrom = line[0]
	start = line[1]
	end = line[2]
	strand = line[5]

	mid_pos = int(numpy.ceil((float(start) + float(end)) / 2))

	if len(prev_line) != 0:
		prev_chrom, prev_mid_pos, prev_strand = prev_line

		if chrom == prev_chrom and strand != prev_strand:
			dist = abs(mid_pos - prev_mid_pos)

			if dist < 2:
				is_duplicate = True


	if not is_duplicate:
		bed_line = '\t'.join([chrom, start, end])
		bed_output.write(bed_line + '\n')
	
	prev_line = [chrom, mid_pos, strand]

bed_output.close()

os.system('rm -r ' + tmp_dir)

###############################################################################################################################################################
