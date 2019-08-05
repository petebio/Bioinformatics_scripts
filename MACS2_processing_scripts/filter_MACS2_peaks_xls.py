#!/usr/bin/env python
import pybedtools
import argparse
import sys
import re

#####################################################################################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('xls_file', type = str, help = 'MACS2 peaks xls file to filter')
parser.add_argument('height', type = int, help = 'Peak height cutoff')
parser.add_argument('outfile', type = str, help = 'Output file')
parser.add_argument('-b', '--blacklist', dest = 'b', required = False, help = 'An optional BED file to filter peaks against (e.g. blacklist)')

args = parser.parse_args()

#####################################################################################################################################################

xls = open(args.xls_file, 'r')
bed = list()

total_count = 0

for line in xls:
	line = line.strip()

	if line.startswith('#') or line == '' or re.search('start', line):
		continue

	fields = line.split('\t')

	chrom = fields[0]
	summit = fields[4]
	height = int(fields[5].rstrip('.00'))

	if height >= args.height:
		start = str(int(summit) - 1)

		bed_line = '\t'.join([chrom, start, summit])
		bed.append(bed_line)

	total_count = total_count + 1

xls.close()

bed = pybedtools.BedTool(bed).sort()

removed_count = total_count - bed.count()

sys.stderr.write('Read %d peaks\n' % total_count)
sys.stderr.write('Removed %d peaks with height < %d\n' % (removed_count, args.height))

#####################################################################################################################################################

if args.b:
	before_count = bed.count()

	blacklist = pybedtools.BedTool(args.b)
	bed = bed.intersect(blacklist, v = True)

	after_count = bed.count()
	bl_removed_count = before_count - after_count

	sys.stderr.write('Removed %d peaks with blacklist\n' % bl_removed_count)

#####################################################################################################################################################

sys.stderr.write('Writing %d peaks to ' % bed.count() + args.outfile + '\n')
bed.saveas(args.outfile)

#####################################################################################################################################################
