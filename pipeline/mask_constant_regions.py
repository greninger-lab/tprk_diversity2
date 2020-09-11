# Masks constant regions in full length sequences from unwrap.fasta, which can be generated from the following command: 
# awk 'BEGIN {RS=">";FS="\n";OFS=""} NR>1 {print ">"$1; $1=""; print}' Isolates_aa_filt_fullORFs.aln.fasta > unwrap.fasta

# Constant regions are found with fuzzy match based on length, and replaced with X's.

import subprocess
import argparse
import sys
import os 
import regex
from itertools import chain
import numpy as np
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pandas as pd

constant_regions = ["MIDPSATSRYGSPRLVSNGFRHRRKVVYQRVGHRRFSLIFFFVVVLGRSPRLWAQVSFTPDIEGYAELAW",
"GFKTTTDFKIVFPIVAKKDFKYRGEGNVYAEINVKALKLSLESNGGAKFDTKGSAKTIEATLHCYGAYLTIGKNPDFKSTFAVLWEPWTANGDYKSKGDKPVYEPGFEGAGGKLGYKQTDIAGTGLTFDIAFKFASNTD",
"GDILFGWERTREDGVQEYIKVELTGNS","KLWGLCALAA","GADALLTSGYRWFSAGGYFAS","KLETKGSDPDTSFLEGLDLGVDVRTYM","YFPVYGKVWGSYRHDMGEYGWVKVYANL","ECGVVVSPLEKVEIRLSWEQGKLQENSNVVIEKNVTERWQFVGACRLIW"]

fasta_file = "unwrap.fasta"
output_file = open("masked_constant_regions.fasta","w+")

output_names = []
output_read_seqs = []

for line_num, line in enumerate(open(fasta_file)):
	## Read name
	if (line_num % 2 == 0):
		read_name = line.rstrip()
		read_count = int(read_name.split("_")[2])

		#output_names.append(read_name)
		#output_file.write(read_name+"\n")
	
	## Read sequence
	else:
		read_seq = line.rstrip()
		masked_read_seq = read_seq
		num_masks = 0

		## Loop through constant regions
		for constant_seq in constant_regions:
			## First find exact match and replace with Ns 
			exact_match = regex.search(constant_seq,read_seq)
			if exact_match:
				masked_read_seq = masked_read_seq.replace(exact_match[0],'X'*len(exact_match[0]))
				num_masks +=1
			
			## Now we look for fuzzy matches
			else:
				## Fuzzy matching for shorter constant regions will be <3 substitutions
				if(len(constant_seq) < 48):
					fuzzy_match = regex.search(r"(?b)("+constant_seq + "){s<=3}", read_seq)
				else:
					fuzzy_match = regex.search(r"(?b)("+constant_seq + "){s<=5}", read_seq)
				
				## Replace fuzzy match with Ns
				if fuzzy_match:
					masked_read_seq = masked_read_seq.replace(fuzzy_match[0],'X'*len(fuzzy_match[0]))
					num_masks+=1

		## First sequence, just add to list
		if output_names == []:
			output_names.append(read_name)
			output_read_seqs.append(masked_read_seq)

		else:
			## Otherwise look for existing read seq
			existing_match = False
			for index_num, existing_read_seq in enumerate(output_read_seqs):
				if masked_read_seq == existing_read_seq:
					prev_name = output_names[index_num]
					prev_count = int(prev_name.split("_")[2])

					sample_name = read_name.split("_")[0]
					prev_sample_name = prev_name.split("_")[0]

					# only collapse if they are within the same sample
					if sample_name == prev_sample_name:
						new_count = prev_count + read_count
						new_name = prev_name.split("_")[0]+"_"+prev_name.split("_")[1]+"_"+ str(new_count)

						output_names[index_num] = new_name
						print(output_names[index_num])
						print("collapsing ",prev_name," and ",read_name," into ",new_name)

						existing_match = True
					else:
						print("Identical matches between different samples found. Not collapsing.")
		
			if not existing_match:
				output_names.append(read_name)
				output_read_seqs.append(masked_read_seq)


for index, name in enumerate(output_names):
	output_file.write(name + "\n")
	output_file.write(output_read_seqs[index] + "\n")