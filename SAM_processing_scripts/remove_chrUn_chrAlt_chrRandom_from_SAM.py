#!/usr/bin/env python
import sys
import re

#################################################################################################################################

for line in sys.stdin:
	keep_line = True

	if re.search('^@SQ', line):
		chrom = re.search('SN:(\S+)', line).group(1)
		
		if re.search('random', chrom) or re.search('alt', chrom) or re.search('chrUn', chrom):
			keep_line = False

	elif not re.search('^@', line):
		fields = line.strip().split('\t')
		
		chrom_a = fields[2]
		chrom_b = fields[6]

		if re.search('random', chrom_a) or re.search('random', chrom_b):
			keep_line = False
		elif re.search('alt', chrom_a) or re.search('alt', chrom_b):
			keep_line = False
		elif re.search('chrUn', chrom_a) or re.search('chrUn', chrom_b):
			keep_line = False

	if keep_line:
		sys.stdout.write(line)

#################################################################################################################################
