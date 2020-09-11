# Converts masked_constant_regions.fasta into a summary .csv similar to Table_allAAfilt_fullORFs.tsv,
# but with the masked and colapsed regions.

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

fasta_file = "masked_constant_regions.fasta"

table = open("full_length.csv","w+")

# Translates a string of nucleotides into amino acids.
def translate_nucs(read_seq):
	coding_dna = Seq(read_seq, generic_dna)
	translation = str(coding_dna.translate())
	return translation

read_total = 0

for line_num, line in enumerate(open(fasta_file)):
	if (line_num % 2 == 0):
		read_name = line.rstrip()
		read_count = int(read_name.split("_")[2])
		read_name = read_name.replace(">","")
		sample_name = read_name.split("_")[0]
		sample_name = sample_name.replace(">","")
	else:
		read_seq = line.rstrip()
		#read_seq_rev = str((Seq(read_seq, generic_dna).reverse_complement()))
		#table.write(read_name + "," + translate_nucs(read_seq_rev) + "," + 
		#	str(read_count)+"," + "\n")
		table.write(read_name+","+read_seq+","+ sample_name+"," + str(read_count) + "\n")
		read_total += read_count

print(read_total)