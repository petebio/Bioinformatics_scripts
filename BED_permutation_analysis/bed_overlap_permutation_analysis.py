#!/usr/bin/env python
import pybedtools
import argparse
import numpy
import sys

###########################################################################################################################################################################################################################

parser = argparse.ArgumentParser(description = 'Carries out a permutation test to determine the significance of an overlap between a set of BED files')
parser.add_argument('fg_bed', type = str, help = 'BED file of foreground sites')
parser.add_argument('bg_bed', type = str, help = 'BED file of background sites to sample from')
parser.add_argument('comp_bed', type = str, help = 'BED file of sites to compare the foreground and randomly sampled sites to')
parser.add_argument('out', type = str, help = 'Output file. A list of overlap counts from N random samples')
parser.add_argument('-n', type = int, default = 1000, help = 'Number of random samples. Default = 1000')

args = parser.parse_args()

###########################################################################################################################################################################################################################

# Read foreground and comparison sites using pybedtools
fg_bed = pybedtools.BedTool(args.fg_bed).sort()
comp_bed = pybedtools.BedTool(args.comp_bed).sort()

# Count the number of sites that intersect between foreground and comparison sites
fg_intersect_count = fg_bed.intersect(comp_bed, wa = True, u = True).count()

sys.stderr.write('Read %d foreground sites from %s\n' % (fg_bed.count(), args.fg_bed))
sys.stderr.write('Comparing these to %d sites from %s\n' % (comp_bed.count(), args.comp_bed))
sys.stderr.write('Found %d common sites\n' % fg_intersect_count)

###########################################################################################################################################################################################################################

# Read background sites into a list
bg_bed_file = open(args.bg_bed, 'r')
bg_bed = bg_bed_file.read().splitlines()
bg_bed_file.close()

sys.stderr.write('Read %d background sites from %s\n' % (len(bg_bed), args.bg_bed))

# Keep track of intersect counts between random samples and the comparison sites
sample_site_counts = list()

# Record the number of times the intersect of the random samples is greater than that of the actual foreground sites
sample_gt_count = 0

n_fg_sites = fg_bed.count()

sys.stderr.write('Beginning random sampling\n')

for i in range(0,args.n):
	sys.stderr.write('Iteration %d of %d\r' % (i + 1, args.n))

	# Take n_fg_sites sites from background without replacement
	sample = list(numpy.random.choice(bg_bed, n_fg_sites, replace = False))
	sample_bed = pybedtools.BedTool(sample).sort()

	# Count the number of sites that intersect between random sample and comparison sites
	sample_intersect_count = sample_bed.intersect(comp_bed, wa = True, u = True).count()
	sample_site_counts.append(sample_intersect_count)

	# Check if random count is greater than foregound count
	if sample_intersect_count > fg_intersect_count:
		sample_gt_count = sample_gt_count + 1

sys.stderr.write('\nDone!!\n')

# Calculate the p-value - this is the proportion of times the random sample count is greater than the actual foreground count
# Add a count of 1 to the numerator and denominator to avoid miscalculation of the p-value
# See Phipson & Smyth (2010). Permutation P-values Should Never Be Zero: Calculating Exact P-values When Permutations Are Randomly Drawn (PMID: 21044043) for more details
p = (sample_gt_count + 1) / (args.n + 1)
sys.stderr.write('p = %f\n' % p)

# Calculate the mean and standard deviation of the random sample count distribution - use these to calculate the Z-Score
mu = numpy.mean(sample_site_counts)
sigma = numpy.std(sample_site_counts)

z = (fg_intersect_count - mu) / sigma

sys.stderr.write('z = %f\n' % z)

###########################################################################################################################################################################################################################

sys.stderr.write('Writing results to %s\n' % args.out)

out = open(args.out, 'w')

# The first line of the output file contains all summary statistics from the test
out.write('# fg = %s; sample = %s; n_fg_sites = %d; mu = %f; sigma = %f; z = %f; p = %f\n' % (args.fg_bed, args.comp_bed, fg_intersect_count, mu, sigma, z, p))

# Write all counts from the random samples - useful if you want to plot the null distribution later
for count in sample_site_counts:
	out.write('%d\n' % count)

out.close()

###########################################################################################################################################################################################################################