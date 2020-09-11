import subprocess
import argparse
import sys
import os 
import regex
from itertools import chain

sample_name = "80-20-cloneamp1"
fasta_file = sample_name+".fasta"

AS10_major_read_list = ["GIASETGGAGALKH","WEGKSNTGVPAGVTPSKYGLG","TLSSGYAPVPANDILWDVGAKVSM","TDVGRKKDGANGDI","KASNVFEGVFLARNIAMREHDCAAYI","PVYYFAARARPGAAADDI","YGGTNKKNDAAAAAPAPATKWKAEYCGYY"]
MI04_major_read_list = ["GIAYENGGAIKH","WEGKPNGNVPAGVTHSKYGLG","TLSSGYAQAAGAAAAAAVNNAILWDVGAKVSM","TDVGHKKENAANGDI","KASNVFKDVFLTTPMLQHDCAAYI","PVYYFAARARAGAAVPAI","YGGTNKKAAALTKWKAEYCGYY"]

output_file = open(sample_name+"_chimeric_matches.txt","w+")

AS10_sum = 0
MI04_sum = 0

for line_num, line in enumerate(open(fasta_file)):
	if (line_num % 2 == 0):
		read_name = line.rstrip()
		read_count = int(read_name.split("_")[2])
	else:
		read_seq = line.rstrip()

		AS10_matches = 0
		MI04_matches = 0

		for region_seq in AS10_major_read_list:
			exact_match = regex.search(region_seq,read_seq)
			if exact_match:
				match_seq = exact_match[0].rstrip()
				AS10_matches +=1
			# Allow for 1 substitution
			else:
				fuzzy_match = regex.search(r"(?b)("+region_seq + "){s<=1}", read_seq)
				if fuzzy_match:
					#print("Fuzzy match found.")
					#print(fuzzy_match)
					#print("Actual sequence was actually "+region_seq)
					AS10_matches +=1


		for region_seq in MI04_major_read_list:
			exact_match = regex.search(region_seq,read_seq)
			if exact_match:
				match_seq = exact_match[0].rstrip()
				MI04_matches +=1
			# Allow for 1 substitution
			else:
				fuzzy_match = regex.search(r"(?b)("+region_seq + "){s<=1}", read_seq)
				if fuzzy_match:
					#print("Fuzzy match found.")
					#print(fuzzy_match)
					#print("Actual sequence was actually "+region_seq)
					MI04_matches +=1


		match_line = read_name+" matched AS10 "+str(AS10_matches)+" and matched MI04 "+str(MI04_matches)+"\n"
		output_file.write(match_line)
		# print(match_line)

		if((AS10_matches + MI04_matches) < 7):
			output_file.write("YO NOT UP TO 7\n")
			#print("not up to 7")
		#output_file.write(read_seq+"\n")
