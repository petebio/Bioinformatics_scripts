#!/usr/bin/env python
import sys
import re

#################################################################################################################################

read_id_tracking = dict()

for line in sys.stdin:
	is_duplicate = False

	if not re.search('^@', line):
		fields = line.strip().split('\t')

		read_id = fields[0]
		chrom_a = fields[2]

		if read_id not in read_id_tracking:
			read_id_tracking[read_id] = 0

		if not chrom_a == '*':
			read_id_tracking[read_id] = 1

			if re.search('XS:i:\d+', line):
				is_duplicate = True
			else:
				read_id_tracking[read_id] = 2

	if not is_duplicate:
		sys.stdout.write(line)

#################################################################################################################################

total_read_count = len(read_id_tracking)

aligned_read_count = 0
unique_read_count = 0

for value in read_id_tracking.values():
	if value != 0:
		aligned_read_count = aligned_read_count + 1

	if value == 2:
		unique_read_count = unique_read_count + 1

perc_aligned = float(aligned_read_count) / float(total_read_count) * 100
perc_unique_aligned = float(unique_read_count) / float(aligned_read_count) * 100

sys.stderr.write('Read %d reads from SAM file\n' % (total_read_count))
sys.stderr.write('Of these, %d (%.2f%%) were aligned\n' %(aligned_read_count, perc_aligned))
sys.stderr.write('%d (%.2f%%) of alignments were unique\n' % (unique_read_count, perc_unique_aligned))

#################################################################################################################################		
