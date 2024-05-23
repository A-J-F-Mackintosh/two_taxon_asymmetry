#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Usage: estimate_Am.py -v <STR> -a <STR> -b <STR> -s <INT> [-h -j <INT> -r <STR>]

  [Options]
    -v, --vcf <STR>                             VCF file
    -a, --a_samples <STR>                       A comma delimited list of samples from population A, e.g. sample_1,sample_2,sample_3
    -b, --b_samples <STR>                       A comma delimited list of samples from population B
    -j, --jackknife_blocksize <INT>             Number of SNPs in each block (aim for >20 blocks) [default: 100_000]
    -r, --rollover <STR>                        Allow jackknife blocks to roll over between sequences [default: False]
    -s, --span <INT>                            Look for het' mutations within this many bases of a hetAB
    -h, --help                                  Show this message

"""

import sys
from docopt import docopt
import itertools
import allel
import numpy as np
import math
import copy


def S_from_het_prime_list(het_prime_list):
	a = sum([entry[0] for entry in het_prime_list])
	b = sum([entry[1] for entry in het_prime_list])
	if a + b == 0:
		return 0
	else:
		S = (a - b) / (a + b)
		return S

def block_jackknife(het_prime_counts):
	overall_estimate = S_from_het_prime_list(het_prime_counts)
	n = len(het_prime_counts)
	if n < 2:
		return [overall_estimate, overall_estimate, overall_estimate]
	else:
		pseudo_estimates = []
		for i in range(0, n):
			retained_blocks = []
			for j in range(0, n):
				if j != i:
					retained_blocks.append(het_prime_counts[j])
			retained_blocks_estimate = S_from_het_prime_list(retained_blocks)
			pseudo_estimate = (n * overall_estimate) - ((n - 1) * retained_blocks_estimate)
			pseudo_estimates.append(pseudo_estimate)
		pseudovalue_variance = sum([(e - overall_estimate)**2 for e in pseudo_estimates]) / (n - 1)
		lower_CI = overall_estimate - (1.96 * math.sqrt(pseudovalue_variance / n))
		higher_CI = overall_estimate + (1.96 * math.sqrt(pseudovalue_variance / n))
		return [lower_CI, overall_estimate, higher_CI]


def collect_prime(het_pos, hetAB_pos, fixed_pos, span):

	FGVtrim = False

	het_prime = 0

	hetAB_pos_fl = np.copy(hetAB_pos)
	fixed_pos_fl = np.copy(fixed_pos)

	hetAB_to_be_clipped = np.array([], dtype=np.int8)
	fixed_to_be_clipped = np.array([], dtype=np.int8)

	for het in het_pos:

		#print("considering het: ", het)

		# remove hetAB and fixed from arrays if they are trailing the current het position
		if len(hetAB_to_be_clipped) > 0:
			hetAB_pos_fl = np.delete(hetAB_pos_fl, hetAB_to_be_clipped)
			hetAB_to_be_clipped = np.array([], dtype=np.int8)
		if len(fixed_to_be_clipped) > 0:
			fixed_pos_fl = np.delete(fixed_pos_fl, fixed_to_be_clipped)
			fixed_to_be_clipped = np.array([], dtype=np.int8)

		# now look at the next hetAB
		for hetAB in range(0, len(hetAB_pos_fl)):

			primed = False

			#print("considering hetAB: ", hetAB_pos_fl[hetAB])

			# mark trailing hetAB for removal
			if het - hetAB_pos_fl[hetAB] > span:
				hetAB_to_be_clipped = np.append(hetAB_to_be_clipped, np.array([hetAB]))

				#print("marked hetAB {} for removal".format(hetAB_pos_fl[hetAB]))

			# if nearby the het, check whether there is also a fixed nearby
			elif abs(het - hetAB_pos_fl[hetAB]) <= span:

				if FGVtrim:

					het_hetAB_dist = abs(het - hetAB_pos_fl[hetAB])

					#print("hetAB {} may yield a het_prime".format(hetAB_pos_fl[hetAB]))

					for fixed in range(0, len(fixed_pos_fl)):

						#print("considering fixed: ", fixed_pos_fl[fixed])

						# mark trailing fixed for removal
						if het - fixed_pos_fl[fixed] > 0 and hetAB_pos_fl[hetAB] - fixed_pos_fl[fixed] > 0:

							if het - fixed_pos_fl[fixed] > span:
								fixed_to_be_clipped = np.append(fixed_to_be_clipped, np.array([fixed]))

								#print("marked fixed {} for removal".format(fixed_pos_fl[fixed]))

						# if there is no fixed between (reached a fixed that is beyond), then add a het_prime
						elif fixed_pos_fl[fixed] > het and fixed_pos_fl[fixed] > hetAB_pos_fl[hetAB]:
							het_prime += 1
							#print("found a het_prime")
							primed = True
							break # stop looping over fixed

						# this means the fixed must be inbetween
						else:

							#print("fixed {} is inbetween".format(fixed_pos_fl[fixed]))

							if fixed_pos_fl[fixed] - het > 0: # this means no more hetABs will help find a prime
								primed = True
								#print("and is also beyond the het")
							break

				else:
					het_prime += 1
					primed = True

				# break from hetAB loop if you have already found a het_prime, or will never find one
				if primed:
					break

			# if the next hetAB is beyond, then stop looping over hetABs
			elif hetAB_pos_fl[hetAB] - het > span:
				#print("hetAB {} is beyond".format(hetAB))
				break

	return het_prime


def get_SNP_arrays(hetA_het_array, hetB_het_array, hetA_homref_array, 
	hetB_homref_array, hetA_homalt_array, hetB_homalt_array):

	hetA_hom_array = np.logical_or(hetA_homref_array, hetA_homalt_array)
	hetB_hom_array = np.logical_or(hetB_homref_array, hetB_homalt_array)

	hetAB_sites = np.logical_and(hetA_het_array, hetB_het_array)
	hetA_sites = np.logical_and(hetA_het_array, hetB_hom_array)
	hetB_sites = np.logical_and(hetA_hom_array, hetB_het_array)
	fixed_sites = np.logical_or(np.logical_and(hetA_homref_array, hetB_homalt_array), 
		np.logical_and(hetA_homalt_array, hetB_homref_array))

	hetAB_pos = positions[hetAB_sites > 0]
	hetA_pos = positions[hetA_sites > 0]
	hetB_pos = positions[hetB_sites > 0]
	fixed_pos = positions[fixed_sites > 0]

	return hetAB_pos, hetA_pos, hetB_pos, fixed_pos


if __name__ == '__main__':
	__version__ = '0.1'
	args = docopt(__doc__)

	# deal with args
	vcf_f = args["--vcf"]
	block_size = int(args["--jackknife_blocksize"])
	rollover = args["--rollover"]
	if rollover == "True":
		rollover = True
	else:
		rollover = False
	span = int(args["--span"])
	a_samples = args["--a_samples"].split(",")
	b_samples = args["--b_samples"].split(",")

	# collect sequences from vcf header
	sequences = []
	header = allel.read_vcf_headers(vcf_f)
	for x in header[0]:
		if x.startswith("##contig="):
			contig = x.split("=")[2]
			contig = contig.split(",")[0]
			sequences.append(contig)

	# collect the indices of a_samples and b_samples
	a_sample_keys = []
	b_sample_keys = []
	samples = allel.read_vcf_headers(vcf_f)[-1]
	for i in range(0, len(samples)):
		if samples[i] in a_samples:
			a_sample_keys.append(i)
		elif samples[i] in b_samples:
			b_sample_keys.append(i)

	# collect het_prime here
	het_prime_counts = []

	# block counter, just for printing
	block = 0

	# keep track for rollingover
	SNPs_so_far = 0

	# loop through the sequences
	for sequence in sequences:

		#print(sequence)

		gt_key, pos_key = "calldata/GT", "variants/POS"
		fields, samples, header, chunks = allel.iter_vcf_chunks(vcf_f, fields=[gt_key, pos_key], 
			region=sequence, tabix="tabix", chunk_length=int(round(block_size/10, 0))) # small contigs will be ignored

		for chunk in chunks:

			# if new block, then reset het_prime counters
			if SNPs_so_far == 0:
				block_het_prime = [0, 0] # block specific values that all samples contribute to

			genotypes = allel.GenotypeArray(chunk[0][gt_key])
			positions = chunk[0][pos_key]

			if rollover:
				SNPs_so_far += len(positions)

			for combo in itertools.product(a_sample_keys, b_sample_keys): # analyse all pairwise combos of samples
				index_A, index_B = combo

				# currently assuming that SNPs are biallelic

				hetA_het_array = genotypes[:, index_A].is_het()
				hetB_het_array = genotypes[:, index_B].is_het()

				hetA_homref_array = genotypes[:, index_A].is_hom_ref()
				hetB_homref_array = genotypes[:, index_B].is_hom_ref()

				hetA_homalt_array = genotypes[:, index_A].is_hom_alt()
				hetB_homalt_array = genotypes[:, index_B].is_hom_alt()

				hetAB_pos, hetA_pos, hetB_pos, fixed_pos = get_SNP_arrays(hetA_het_array, hetB_het_array, 
					hetA_homref_array, hetB_homref_array, hetA_homalt_array, hetB_homalt_array)

				hetA_prime = collect_prime(hetA_pos, hetAB_pos, fixed_pos, span)
				hetB_prime = collect_prime(hetB_pos, hetAB_pos, fixed_pos, span)

				block_het_prime[0] += hetA_prime
				block_het_prime[1] += hetB_prime

			if rollover and SNPs_so_far >= block_size or not rollover:

				#print(SNPs_so_far)
				SNPs_so_far = 0

				if block_het_prime[0] + block_het_prime[1] == 0:
						print("[!] WARNING: there were zero het' sites in this block. Consider increasing the jackknife block size")

				else:

					het_prime_counts.append(block_het_prime)

					CIs = block_jackknife(het_prime_counts)

					print("[=] A_m for block_{}: {}".format(block, 
						round(S_from_het_prime_list([block_het_prime]), 6)))
					print("[=] Current estimate of A_m: {} ({}, {})".format(round(CIs[1], 6), round(CIs[0], 6), round(CIs[2], 6)))

					block += 1

	print("[=] Final estimate of A_m: {} ({}, {})".format(round(CIs[1], 6), round(CIs[0], 6), round(CIs[2], 6)))

