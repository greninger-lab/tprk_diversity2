# Syph_r, but for full length. Makes everything into "V1" (but is actually full length).

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

V1_list = list()

# Takes in the input file and goes through line by line, defining the read name and read sequence.
# Puts the read sequence through region_seq to match each read to the variable region/
def find_region(sample, input_format, is_pacbio, current_dir, strain_name):
	read_name = ""
	read_seq = ""
	region = ""
	sample_name = sample.split("." + input_format)[0]
	num_input_reads = 0

	# Fastq files are in blocks of 4: with first line being name, second line being sequence,
	# whereas fasta files are in blocks of 2.
	if input_format == "fasta":
		line_block = 2
	else:
		line_block = 4
		
	for line_num, line in enumerate(open(current_dir + "/" + sample)):
		# Finds read name
		if (line_num % line_block == 0):
			read_name = line
			# Illumina reads are read as one read = one count. PacBio reads are assumed to be RAD, and counts are
			# found after the underscore in the read name.
			if is_pacbio==False:
				# read_name = read_name.split('-')[0]
				read_count = 1
			# Includes the count in PacBio names after FAD/RAD, may be [1] or [2] depending on file
			else:
				read_count = int(float(read_name.split('_')[2].rstrip()))
			num_input_reads += read_count
		# Finds read sequence
		elif (line_num % line_block == 1):
			read_seq = str((line)).rstrip()
			# Matches each read seq to the region
			region_seq(read_seq, read_name, sample_name, is_pacbio, current_dir, num_input_reads)


# String matches a read to a variable region (V1-V7).
def region_seq(read_seq, read_name, sample_name, is_pacbio, current_dir, num_input_reads):
	# Different beginning portions for the variable regions.
	# Illumina reads are read as one read = one count. PacBio reads are assumed to be RAD, and counts are
	# found after the underscore in the read name.
	if is_pacbio==False:
	# read_name = read_name.split('-')[0]
		read_count = 1
	# Includes the count in PacBio names after FAD/RAD, may be [1] or [2] depending on file
	else:
		read_count = int(float(read_name.split('_')[2].rstrip()))

	
	v_list = V1_list
	
	region_seq = read_seq
	# If the sequence is already in a global list of sequences for this region,
	# which is in [read sequence, count] format, increase the count by 1
	if (region_seq in chain(*v_list)):
		for sublist in v_list:
			if (sublist[0] == region_seq):
				sublist[1] = sublist[1] + read_count
	# Otherwise, add this sequence to the global list of sequences for this region
	else:
		v_list.append([region_seq, read_count])

# Finds the total number of reads that mapped to a region.
def total_count(region_list):
	total = 0
	for read_seq, count in region_list:
		total = total + count
	return total

# Translates a string of nucleotides into amino acids.
def translate_nucs(read_seq):
	coding_dna = Seq(read_seq, generic_dna)
	translation = str(coding_dna.translate())
	return translation

# Makes the final data table in csv format, with columns Region, Read, RelativeFreq, Count.
def make_table(strain_name, current_dir):
	Vlist_of_aas = [V1_list]
	variable_regions = ["V1"]

	table = open(current_dir + "/" + strain_name + "_final_data.csv", "w+")
	table.write("Region,Read,RelativeFreq,Count" + "\n")
	for index, v_list in enumerate(Vlist_of_aas):
		total = total_count(v_list)
		for read_seq, count in v_list:
			table.write(variable_regions[index] + "," + read_seq + "," + 
				str(((count / total) * 100)) + "," + str(count) + "\n")

	# Filters out the lines with greater than 1% to a separate _final_data_fitered.csv.
	#subprocess.call("awk -F\"[,|\\(]\" \'($3+0)>=1{print}\' " + strain_name + "_final_data.csv > " + strain_name + "_final_data_filtered.csv", shell=True)

if __name__ == '__main__': 
	parser = argparse.ArgumentParser(description='tprK project. Example usage: syph.py -i fasta -p pacbio')
	parser.add_argument('-i', '--input_format', required=True,
		help='Specify either fasta or fastq. Will take all the files in this folder with the specified extension. ')
	parser.add_argument('-pacbio', action='store_true', help='Write this flag to specify this file is a pacbio file. '
		'Do not use with -illumina') 
	parser.add_argument('-illumina', action='store_true', help='Write this flag to specify this file is an illumina file. ' 
		'Do not use with -pacbio')
	parser.add_argument('-s', '--strain_name', required=True,
		help='Specify what strain name.')
	parser.add_argument('-d', '--directory', help='Pass directory (used when passed in from within R.')
	
	# Checks for argument sanity.
	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(1)
	
	current_dir = args.directory

	input_format = args.input_format.lower()
	is_pacbio = args.pacbio
	is_illumina = args.illumina
	file_name = args.strain_name

	# Right now the flag for input format is only to distinguish where the read name and read is.
	# Generally, fasta files are in a chunk of 2, fastq in a chunk of 4.
	if input_format in ["fasta", "fastq"]:
		print("Input format is " + input_format + ".")
	else:
		print("ERROR: Please specify either fasta or fastq for -i.")
		sys.exit()

	# Specifying whether the file is PacBio and Illumina changes how the program counts for each read.
	# Currently we assume PacBio sequences will be RAD-ified and in a specific format (i.e. >seq1_274) where
	# the reads have been clustered and the read count will be after one underscore.
	if is_pacbio and is_illumina:
		print("ERROR: -pacbio and -illumina are not compatible. Please select only one. ")
		sys.exit()
	if (is_pacbio == False) and (is_illumina==False):
		print("ERROR: Please select either -pacbio or -illumina.")
		sys.exit()

	# The all assignments file is a list of all the reads and what regions they mapped to, along with both
	# the nucleotide sequence and the amino acid sequence.
	all_assignments = open("all_assignments.csv", "w+")

	input_extension = input_format
	if input_format == "fasta" and is_pacbio:
		input_extension = "RAD.nolines.fix.fasta"

	# Extracts strain name from file name.
	strain_name = file_name.split("." + input_format)[0]
	# Get rid of stupidly long name in PacBio.
	if "RAD" in strain_name:
		strain_name = strain_name.split(".noprimers.filtered.RAD.nolines.fix.fasta")[0]
	
	# Checks if file exists already, and skips.
	file_tobemade = strain_name + "_final_data.csv"
	if os.path.isfile(file_tobemade):
		print(file_tobemade," already exists. Skipping making frequency tables...")
	else:
		print(file_tobemade,"")

		# Matches each read to a region and starts building a list.
		find_region(file_name, input_format, is_pacbio, current_dir, strain_name)
		
		# Builds the frequency final_table.csv for each file.
		make_table(strain_name, current_dir)
		V1_list = list()
		V2_list = list()	
		V3_list = list()	
		V4_list = list()	
		V5_list = list()	
		V6_list = list()	
		V7_list = list()
		V1_dna = list()
		V2_dna = list()
		V3_dna = list()
		V4_dna = list()
		V5_dna = list()
		V6_dna = list()
		V7_dna = list()