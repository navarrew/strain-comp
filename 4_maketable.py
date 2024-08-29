#!/usr/bin/env python
'''
This script takes info from the strainlist.txt file (made by process_ncbi.py) and 
the cluster_header_info.tab file (made by mmseqcluser.py) to generate both the 
cluster_table.tab file and a stripped down version of the cluster_table. tab file
called 'cluster_hit_counts.tab' that will be used later for hierarchical clustering
of strains and proteins.
'''

'''
┌----------------------------------------┐
| import modules						 |
└----------------------------------------┘
'''

import re
from datetime import datetime
import argparse
import os

'''
┌----------------------------------------┐
| define functions						 |
└----------------------------------------┘
'''

def get_locus_tag(loc):
	"""
	Returns the locus tag at the front of the input locus ID.
	There are two general formats of tags:
	(1) Alphanumeric characters preceding an underscore and a string of numbers.
	(2) Alphabet characters preceding a string of numbers.
	"""

	if loc is None:
		return
	if '_' in loc:
		locus_tag = str(loc.split('_')[0])	# For locus tags formatted as alphanumeric characters before an underscore
		return locus_tag
	else:
		locus = re.match(r'\D+', loc)  # For locus tags formatted as letters before a string of numbers
		locus_tag = locus.group(0)
		return locus_tag


def convert_fasta_header_format(fasta_header):
	"""
	Takes an NCBI formatted FASTA header and converts it into a format appropriate for our tab files.

	Example input:

	NZ_LR861808.1_cds_WP_000502119.1_1973 [gene=tnpA] [locus_tag=JMT79_RS10105] \
	[protein=IS200/IS605-like element IS200F family transposase] [protein_id=WP_000502119.1] \
	[location=complement(2043708..2044166)] [pctGC=45.9] [gbkey=CDS]

	Example output:

	NZ_LR861808.1_cds_WP_000502119.1_1973|JMT79_RS10105|WP_000502119.1
	"""
	items_in_header = fasta_header.split(' [')
	accession = items_in_header[0][1:]	# Take first item in list and remove the first > character
	fasta_locus_tag = 'none'
	fasta_protein_id = 'none'
	for fasta_item in items_in_header:
		if 'locus_tag=' in fasta_item:
			fasta_locus_tag = fasta_item[10:-1]
		if 'protein_id=' in fasta_item:
			fasta_protein_id = fasta_item[11:-1]
	output_string = fasta_locus_tag + '|' + fasta_protein_id + '|' + accession
	return output_string


def get_matches_from_cluster(str_loc_pre, members_list):
	"""
	Takes a strain locus prefix and a list of members in a cluster from a single line of "cluster_header_info.tab".
	Outputs a string and list of protein matches for the cluster in the given strain.
	"""
	
	match_list = []
	for member in members_list:
#		mem_loc_pre = get_locus_tag(member.split('[locus_tag=')[1].split('] ')[0])
#		if mem_loc_pre == str_loc_pre:
		str_loc_prefix = str_loc_pre + "_"
		if str_loc_prefix in member:
			member_form = convert_fasta_header_format(member)
			# Append the hit to a growing list of hits for this particular strain
			match_list.append(member_form)
			# Remove the match from the list to shorten future searches through it, saving time
			members_list.remove(member)
	# If there are one or more matches in the list,
	# Add the number of matches in the list and join all the matches together
	if len(match_list) >= 1:
		matches = '[' + str(len(match_list)) + ']|' + str(', '.join(match_list))
		return matches, match_list
	# If there are no matches in the list, use "*" as a placeholder
	else:
		matches = '*'
		return matches, match_list


def get_cluster_annotations(cluster_members_list):
	"""
	Takes list of members in a cluster from a single line of "cluster_header_info.tab".
	Outputs various annotation data.
	"""
	GC_counts = []
	annotation_set = set()
	trunc_members_list = []
	slippage_members_list = []
	flag = ""

	for member in cluster_members_list:
		protein_length = 0
		#put GC contents of each member into a list
		member_GC_content = float(member.split('[pctGC=')[1].split('] ')[0])
		GC_counts.append(member_GC_content)

		
		#Get length of proteins and find oddballs
		memberx = member.replace('location=complement(','SPLIT_ME')
		memberx = memberx.replace('location=','SPLIT_ME')
		memberx = memberx.replace(')','')
		member_nucleotide_length = memberx.split('SPLIT_ME')[1].split('] ')[0]
		
		#Put NCBI annotations into a set
		member_annotation = annotation_set.add(member.split('[protein=')[1].split('] ')[0])

		
		if (">" in member_nucleotide_length or "<" in member_nucleotide_length):
			trunc_members_list.append(member)
		elif "join" in member_nucleotide_length:
			slippage_members_list.append(member)
		else: 
			gene_ends = member_nucleotide_length.split('..')
			gene_length = int(gene_ends[1]) - int(gene_ends[0]) + 1
			protein_length = str(gene_length // 3)
			
	if (len(trunc_members_list) != 0 or len(slippage_members_list) != 0):
		if len(trunc_members_list) == len(cluster_members_list):
			flag = 'truncated: ' + str(len(trunc_members_list)) +' of '+ str(len(cluster_members_list))
		if len(slippage_members_list) == len(cluster_members_list):
			flag = 'frameshifts: ' + str(len(slippage_members_list)) +' of '+ str(len(cluster_members_list))
		else:
			flag = str(len(trunc_members_list)) +' trunc; '+ str(len(slippage_members_list)) +  ' frmshft out of ' + str(len(cluster_members_list))

	#put it together for export

	#get the 
	GCavg = "???"
	if len(GC_counts) != 0:
		GCavg = sum(GC_counts) / len(GC_counts)
		GCavg = str('{:.4}'.format(GCavg))	

	#calculate the spread of GC contents by taking the top and bottom GC values from a sorted list
	GC_counts.sort()
	if len(GC_counts) >= 2:
		GC_diff = str('{:.4}'.format(abs(GC_counts[0] - GC_counts[-1])))
	else: GC_diff = "0.00"

	cluster_annotations = ', '.join(annotation_set)

	return GCavg, GC_diff, cluster_annotations, protein_length, flag



'''
┌----------------------------------------┐
| main program							 |
└----------------------------------------┘
'''


if __name__ == '__main__':

	print("\nPipeline for genomic comparisons",
		  "Navarre lab updated August, 2024. UPDATED",
		  "\nPlease have 'cluster_hit_count_table.tab' and 'strainlist.txt' in the working directory.",
		  sep='\n')

	instructions = "tablemaker.py\n" \
					"This will reorder and rename your protein clusters by abundance.\n" \
					

'''
┌------------------------------------------------┐
| set up directories and check for key files	 |
└------------------------------------------------┘
'''


# Set up tab and report directories, if not already present
	home_directory_list = []
	for entry in os.scandir():
		if entry.is_dir():
			home_directory_list.append(entry.name)
	if "tables" not in home_directory_list:
		os.mkdir("tables")

	data_directory_list = []
	for entry in os.scandir('data'):
		if entry.is_dir():
			data_directory_list.append(entry.name)
	if "report" not in data_directory_list:
		os.mkdir("data/report")

#this allows a help message if the user put -h in the command line			
	parser = argparse.ArgumentParser(add_help=True, description=instructions)
	
	# Write the date and time of execution to the report file
	report = open('data/report/report.txt', 'a')
	report.write("\n\n***TABLEMAKER.PY***\nDATE AND TIME: " + str(datetime.now()))


'''
┌---------------------------------------------┐
| parse through files and join to make tables |
└---------------------------------------------┘
'''

	print("\nMaking tables of clusters and their hits in each strain...")
	timestart = datetime.now()
	print("Started at: " + timestart.strftime("%Y-%m-%d %H:%M:%S"))

	# Take the file "cluster_header_info.tab" and parse it into  different tab files:
	# (1) "cluster_table.tab" shows the full information for all protein members belonging to a cluster in a strain
	# (2) "cluster_hit_count_table.tab" only shows the number of hits for a cluster in each strain

	clus_tab = open('tables/cluster_table.tab', 'w')
	hitc_tab = open('tables/cluster_hit_count_table.tab', 'w')

	# Read in the strain list - this will provide the headers for the right half of the table
	with open('strainlist.txt') as f:
		str_list = [line.rstrip() for line in f]
		# Count the total number of strains analyzed
		tot_strains = len(str_list)
		# Populate a new strain list with the locus prefixes
		str_loc_pre_list = []
		for str_info in str_list:
			# First section of each line in the strain list file is the locus id
			str_loc_pre_list.append(str_info.split(' | ')[0])

	print("Clustering " + str(len(str_list)) + " strains.")

	# Write out the header line for each of the table files
	clus_tab.write('CLUSTER\tNCBI gene names\tNCBI annotations\tGCpct\tGC spread\tprotein length\tflags\ttotal count\tstrain count\t' + '\t'.join(str_list) + '\n')
	hitc_tab.write('CLUSTER\ttotal count\tstrain count\t' + '\t'.join(str_list) + '\n')

	# Open "cluster_header_info.tab" file which has each cluster followed by its members in every line
	# The front part gets basic information about the cluster to put in the first few columns of our output files
	# Then iterate through each of the hits in each cluster to place them in the appropriate column for their strain

	with open('data/mmseq_output/cluster_header_info.tab', 'r') as f:
		for line in [line.rstrip() for line in f.readlines()]:

			# Generate the front part for a line in each output file, which has basic information about the cluster
			terms = line.split('\t')
			front_hitc_tab = terms[0] + ' - ' + terms[2] + '\t' + str(terms[3]) #this is the lead info for the hit-count tab file
			#makes a new field for hit count table CLUSTER_X-definition tab 'total counts'

			# Generate lists of the names of hit loci and the counts of hit loci
			members_list = terms[4:]
			cluster_GCpct, GC_spread, cluster_annotations, prot_length, flags = get_cluster_annotations(members_list)
			hit_loci_list = []
			hit_counts_list = []
			for str_loc_pre in str_loc_pre_list:
				matches, match_list = get_matches_from_cluster(str_loc_pre, members_list)
				hit_loci_list.append(matches)
				hit_counts_list.append(str(len(match_list)))
			#the first two items in the line include the cluster name and the NCBI gene ID.
			#grab them here and start a new list called "front cluster info list"
			front_cluster_info_list = terms[0:2]
			#another item in the line is how many total hits there are.  Grab it into a new variable 'total_cluster_hits'
			total_cluster_hits = terms[3]
			#now we make the first few columns as it will appear in the cluster_table.tab file. 
			add_on_list = [cluster_annotations, cluster_GCpct, GC_spread, str(prot_length), flags, total_cluster_hits]
			add_on = '\t'.join(add_on_list)
			front_cluster_info_list.append(add_on)
			front_clus_tab = '\t'.join(front_cluster_info_list)
			
			# Calculate the number of strains that had a hit
			num_str_w_hits = len(hit_loci_list) - hit_loci_list.count('*')

			# Get the names of hit loci and their counts as strings
			hit_loci = '\t'.join(hit_loci_list)
			hit_counts = '\t'.join(hit_counts_list)

			# Write line that joins the preliminary basic info, number of strains with hits, and hit loci or hit counts
			clus_tab.write(front_clus_tab + '\t' + str(num_str_w_hits) + '\t' + hit_loci + '\n')
			hitc_tab.write(front_hitc_tab + '\t' + str(num_str_w_hits) + '\t' + hit_counts + '\n')

	timestart = datetime.now()
	print("Finished at: " + timestart.strftime("%Y-%m-%d %H:%M:%S"))


	# Close the output files
	clus_tab.close()
	hitc_tab.close()
