# Converts Illumina_validated_Illumina_allreads.csv to allreads_Illumina.csv, switching from 
# DNA to amino acid level and retallying if necessary.

import subprocess
import argparse
import sys
import os 
import regex
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pandas as pd
from itertools import chain

# Translates a string of nucleotides into amino acids.
def translate_nucs(read_seq):
	coding_dna = Seq(read_seq, generic_dna)
	translation = str(coding_dna.translate())
	return translation

if __name__ == '__main__': 
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-i', '--illumina_check', help='Illumina file to turn into allreads.')

	# Checks for argument sanity.
	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(1)

	illumina_data = args.illumina_check
	all_reads = open("allreads.csv","w+")

	row_list = []
	translated_read_seqs = []

	for line_num, line in enumerate(open(illumina_data)):
		if (line_num == 0):
			header = line.rstrip()
			print("header", header)
			all_reads.write(header+"\n")
		else:
			line = line.rstrip()
			line_parts = line.split(",")
			read_seq = line_parts[1]
			read_seq_aa = translate_nucs(read_seq)

			translated_line = line.rstrip().replace(read_seq,read_seq_aa)
			translated_line = translated_line.rstrip()

			if (read_seq_aa in translated_read_seqs):
				index_num = translated_read_seqs.index(read_seq_aa)
				#print(read_seq_aa)
				#print(index_num)
				#print(row_list[index_num])
				#print(translated_read_seqs[index_num] + "\n")
				previous_line = row_list[index_num]
				line_parts2 = previous_line.split(",")
				new_line = ""
				print(previous_line)
				print(translated_line)

				for part_num, part in enumerate(line_parts2):
					part = part.rstrip()
					if (part_num == 1):
						new_line = line_parts2[0] + "," + part
					elif (part_num > 1):
						num1 = 0
						num2 = 0
						if (part != "NA"):
							num1 = float(part)
						if (line_parts[part_num] != "NA"):
							num2 = float(line_parts[part_num])

						new_num = num1 + num2
						if(new_num == 0):
							new_num = "NA"
						new_line = new_line + "," + str(new_num)
				row_list[index_num] = new_line
				print(row_list[index_num]+"\n")

			else:
				translated_read_seqs.append(read_seq_aa)
				row_list.append(translated_line)

	for line in row_list:
		all_reads.write(line)
		all_reads.write("\n")