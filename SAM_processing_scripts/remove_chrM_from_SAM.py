#!/usr/bin/env python
import sys
import re

#################################################################################################################################

for line in sys.stdin:
	keep_line = True

	if re.search('^@SQ', line):
		chrom = re.search('SN:(\S+)', line).group(1)

		if chrom == 'chrM':
			keep_line = False

	elif not re.search('^@', line):
		fields = line.strip().split('\t')
		
		chrom_a = fields[2]
		chrom_b = fields[6]

		if chrom_a == 'chrM' or chrom_b == 'chrM':
			keep_line = False

	if keep_line:
		sys.stdout.write(line)

#################################################################################################################################
