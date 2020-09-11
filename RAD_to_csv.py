import subprocess
import argparse
import sys
import os 
import regex
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pandas as pd
from itertools import chain

if __name__ == '__main__': 
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-i', '--illumina_check', help='Illumina^2 file to validate against.')
	parser.add_argument('-s', '--sample', help='RAD fasta to turn into csv.')

	
	# Checks for argument sanity.
	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(1)

	rad_fasta = args.sample
	rad_fasta_filtered = rad_fasta.split(".noprimers.filtered.RAD.nolines.fix")[0] + "_illumina_validated.fasta"
	rad_fasta_thrown_away = rad_fasta.split(".noprimers.filtered.RAD.nolines.fix")[0] + "_illumina_validated_thrown_away.fasta"

	#pacbio_fastq = args.sample
	#pacbio_fastq_filtered = pacbio_fastq.split(".Q20.noprimers.filtered.fastq")[0] + "_illumina_validated.fasta"
	#pacbio_fastq_thrown_away = pacbio_fastq.split(".Q20.noprimers.filtered.fastq")[0] + "_illumina_validated_thrown_away.fasta"

	#rad_fasta = pacbio_fastq
	#rad_fasta_filtered = pacbio_fastq_filtered
	#rad_fasta_thrown_away = pacbio_fastq_thrown_away

	rad_fasta_filtered = open(rad_fasta_filtered, "w+")
	rad_fasta_thrown_away = open(rad_fasta_thrown_away, "w+")
	rad_fasta_stats = open(args.sample.split(".Q20.noprimers.filtered.fastq")[0]+"_stats.txt", "w+")

	illumina_data = pd.read_csv(args.illumina_check)

	illumina_reads = list(illumina_data.Read)
	illumina_regions = list(illumina_data.Region)

	num_kept = 0
	count_kept = 0

	num_filtered = 0
	count_filtered = 0

	total_count = 0

	for line_num, line in enumerate(open(rad_fasta)):
		if (line_num % 2 == 0):
			read_name = line.rstrip()
			#read_name = read_name.replace("@",">")
			read_count = int(read_name.split("_")[1])
			#read_count = 1
		elif (line_num % 2 == 1):
			num_matches = 0

			v1_matched = False 
			v2_matched = False 
			v3_matched = False 
			v4_matched = False 
			v5_matched = False 
			v6_matched = False 
			v7_matched = False 
			
			read_seq = line.rstrip()
			read_seq = str((Seq(read_seq, generic_dna).reverse_complement()))
			for index_num, illumina_verified_read in enumerate(illumina_reads):
				exact_match = regex.search(illumina_verified_read,read_seq)
				if exact_match:
					num_matches += 1
					region_match = illumina_regions[index_num]
					if region_match == "V1":
						v1_matched = True
					elif region_match == "V2":
						v2_matched = True
					elif region_match == "V3":
						v3_matched = True
					elif region_match == "V4":
						v4_matched = True
					elif region_match == "V5":
						v5_matched = True
					elif region_match == "V6":
						v6_matched = True
					elif region_match == "V7":
						v7_matched = True

			if v1_matched and v2_matched and v3_matched and v4_matched and v5_matched and v6_matched and v7_matched:
				rad_fasta_filtered.write(read_name + "\n")
				rad_fasta_filtered.write(read_seq + "\n")
				num_kept += 1
				count_kept += read_count
			else:
				num_filtered += 1
				count_filtered += read_count
				rad_fasta_thrown_away.write(read_name + "\n")
				rad_fasta_thrown_away.write(read_seq + "\n")

			# if(num_matches != 7):
			# 	rad_fasta_filtered.write(read_name + "\n")
			# 	rad_fasta_filtered.write(read_seq + "\n")
			# 	num_kept += 1
			# 	count_kept += read_count

			# if (num_matches < 7):
			# 	num_filtered += 1
			# 	count_filtered += read_count
			# 	rad_fasta_thrown_away.write(read_name + "\n")
			# 	rad_fasta_thrown_away.write(read_seq + "\n")

			# if (num_matches > 7):
				print("whoa > 7")

			total_count += read_count

	rad_fasta_stats.write("Threw out " + str(num_filtered) + " sequences totaling " + str(count_filtered) + " reads out of " + str(total_count) + " for not matching 7 Illumina-twice verified reads." + "\n")
	rad_fasta_stats.write("Kept " + str(num_kept) + " sequences totaling " + str(count_kept) + " reads out of " + str(total_count) + " for matching 7 Illumina-twice verified reads.")