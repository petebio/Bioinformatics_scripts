#!/usr/bin/env python
from clint.textui import progress
import argparse
import pyDNase
import numpy
import math
import sys

#################################################################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('bed_file', type = str, help = 'BED file containing regions to plot')
parser.add_argument('bam_file', type = str, help = 'BAM file containing reads to plot')
parser.add_argument('outfile', type = str, help = 'Output file (.tsv)')
parser.add_argument('-w', '--window', dest = 'w', type = int, default = 200, help = 'Window size to plot. Default = 200bp')

args = parser.parse_args()

#################################################################################################################################

# Read BAM file
reads = pyDNase.BAMHandler(args.bam_file)

# Calculate the distance to extend footprints by (window size / 2)
extend = int(math.ceil(args.w / 2))

# Get regions from BED file
regions = pyDNase.GenomicIntervalSet(args.bed_file)

# Keep track of number of forward and reverse reads
fwd_cut_tracking = dict()
rev_cut_tracking = dict()

sys.stderr.write('Counting cuts in regions...\n')

for site in progress.bar(regions):
	# Get chromosome, strand, start and end positions for site
	chrom = site.chromosome
	start = site.startbp
	end = site.endbp
	strand = site.strand

	# Calculate the center position for the site
	# and extend by +/- (window size / 2)
	center = int(math.ceil((start + end) / 2))
	window_start = center - extend
	window_end = center + extend

	for n,i in enumerate(range(window_start, window_end + 1), 1):
		# Initialize list of cut values in fwd and rev tracking dictionaries
		if n not in fwd_cut_tracking:
			fwd_cut_tracking[n] = list()
			rev_cut_tracking[n] = list()

		# Create a genomic interval compatible with pyDNase and extract number of cuts
		genome_pos = ','.join([chrom, str(i), str(i + 1), strand])		
	
		cuts = reads[genome_pos]
		fwd_cut = cuts['+'][0]
		rev_cut = cuts['-'][0]
	
		fwd_cut_tracking[n].append(fwd_cut)
		rev_cut_tracking[n].append(rev_cut)

#################################################################################################################################

sys.stderr.write('Caculating average profiles...\n')

out = open(args.outfile, 'w')
out.write('bin\tfwd_cuts\trev_cuts\n')

for n in sorted(fwd_cut_tracking):
	fwd_cut_average = numpy.average(fwd_cut_tracking[n])
	rev_cut_average = numpy.average(rev_cut_tracking[n])
	out.write(str(n) + '\t' + str(fwd_cut_average) + '\t' + str(rev_cut_average) + '\n')

out.close()

sys.stderr.write('Done!\n')

#################################################################################################################################
