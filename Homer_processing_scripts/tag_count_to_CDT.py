#!/usr/bin/env python
import argparse
import numpy
import math
import os

#################################################################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('bed_file', type = str, help = 'Input BED file')
parser.add_argument('bedgraph_file', type = str, help = 'BedGraph file')
parser.add_argument('genome', type = str, help = 'Genome version to use with homer')
parser.add_argument('outfile', type = str, help = 'Output file')
parser.add_argument('-s', '--size', dest = 'size', type = int, default = 200, help = 'Window size. Default = 200bp')

args = parser.parse_args()

#################################################################################################################################

tmp_dir = './tmp_' + str(numpy.random.randint(1e6))
os.mkdir(tmp_dir)

bed = open(args.bed_file, 'r')
ext_bed = open(tmp_dir + '/Resized.bed', 'w')

extend = int(math.ceil(args.size / 2))

for line in bed:
	fields = line.strip().split('\t')

	chrom = fields[0]
	pos = int(math.ceil((float(fields[1]) + float(fields[2])) / 2))

	start = pos - extend
	end = pos + extend

	bed_out = '\t'.join([chrom, str(start), str(end)])
	ext_bed.write(bed_out + '\n')

bed.close()
ext_bed.close()

#################################################################################################################################

homer_command = 'annotatePeaks.pl ' + tmp_dir + '/Resized.bed ' + args.genome + ' -bedGraph ' + args.bedgraph_file + ' '
homer_command = homer_command + '-size given -noann > ' + tmp_dir + '/Homer_output.tsv'

os.system(homer_command)

#################################################################################################################################

homer = open(tmp_dir + '/Homer_output.tsv', 'r')
tag_table = dict()

for num,line in enumerate(homer):
	if num == 0:
		continue

	fields = line.strip().split('\t')

	chrom = fields[1]
	start = str(int(fields[2]) - 1)
	end = fields[3]

	pos = chrom + ':' + start + '-' + end

	tag_count = '0'
	if len(fields) == 20:
		tag_count = fields[19]

	tag_table[pos] = tag_count

homer.close()

#################################################################################################################################

bed = open(tmp_dir + '/Resized.bed', 'r')
out = open(args.outfile, 'w')

out.write('GID\tCOORD\tNAME\tGWEIGHT\tFC\n')

for num,line in enumerate(bed, 1):
	fields = line.strip().split('\t')
	
	chrom = fields[0]
	start = fields[1]
	end = fields[2]

	pos = chrom + ':' + start + '-' + end

	tag_count = '0'
	if pos in tag_table:
		tag_count = tag_table[pos]

	cdt = '\t'.join([str(num), pos, pos, '1', tag_count])
	out.write(cdt + '\n')

bed.close()
out.close()

os.system('rm -r ' + tmp_dir)

#################################################################################################################################
