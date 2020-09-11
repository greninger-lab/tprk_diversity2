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

sequence = "TGTCGGGGCTAAGGTGAGTATGA"

fasta_file = "all-variable_unwrap.fasta"

total_v3_count = 0
total_count = 0

relative_counts = []
v3_sum = 759886

for line_num, line in enumerate(open(fasta_file)):
	## Read name
	relative_count = 0
	if (line_num % 2 == 0):
		read_name = line.rstrip()
		read_count = int(read_name.split("_")[2])
	
	## Read sequence
	else:
		read_seq = line.rstrip()
		
		exact_match = regex.search(sequence,read_seq)
		if(exact_match):
			total_v3_count+= read_count
			relative_count = read_count / v3_sum
			relative_counts.append(relative_count)

		if("V3" in read_name):
			total_count += read_count

num_samples = 0
relative_count_sum = 0

for relative_count in relative_counts:
	num_samples +=1
	relative_count_sum += relative_count

avg = relative_count_sum / num_samples

