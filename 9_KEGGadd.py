#!/usr/bin/env python

import os
import argparse

#	  PROGRAM DESCRIPTION
#
#	  Takes a detailed-tsv KEGG file, grabs the significant hits marked with '*' at the first character of the line,
#	  then grabs the cluster ID and adds KO numbers and KEGG definitions in two adjacent columns of the tab file.
#
#	  INPUTS
#
#	  "species_detail.tab": Output from kofamscan, a tab-separated file with KEGG annotation of each protein cluster.
#	  Column 1: Gene name, in our case cluster number.
#	  Column 2: KEGG Orthology (KO) ID.
#	  Column 3: Threshold score for significance of the cluster match to a KO-annotated HMM.
#	  Column 4: Actual score of the cluster match to a KO-annotated HMM.
#	  Column 5: E-value, the probability that a random protein would match to the KO-annotated HMM
#				with a score that is as good or better than the cluster.
#	  Column 6: KO definition.
#
#	  "cluster_table.tab": Output from the main mmseqs2 clustering pipeline, with the following format.
#	  Column 1: Cluster number, in the form 'CLUSTER_N'.
#	  Column 2: Gene names for proteins in the cluster, if any.
#	  Column 3: Gene functions for proteins in the cluster, if any.
#	  Column 4: Total count of all proteins in the cluster.
#	  Column 5: Number of strains with at least one protein in the cluster.
#	  Subsequent columns: Represent each strain, cells contain number and names of proteins in a given cluster.


# Set up directories, if not already present
directory_list = []
for entry in os.scandir():
	if entry.is_dir():
		directory_list.append(entry.name)
if 'tables' not in directory_list:
	print("\nCannot find the tables directory!  Are you sure you're in the right place?")
	exit()
if 'report' not in directory_list:
	os.mkdir('report')
	
data_directory_list = []
for entry in os.scandir('data'):
	if entry.is_dir():
		data_directory_list.append(entry.name)
if 'kegg' not in data_directory_list:
	os.mkdir('data/kegg')
	
outtab1_file = open('data/kegg/cluster_to_kegg.tab', 'w')
outtab2_file = open('tables/cluster_table_KEGG.tab', 'w')

#set the default file names and locations
kegg_input_file = 'species_detail.tab'
input_cluster_table = 'tab/cluster_table.tab'


with open(kegg_input_file, 'r') as f:

	kegg_lines = f.readlines()
	kegg_clus_dict = {}
	KO_info = ""
	KO_list = []
	KO_def_list = []
	prev_clus_name = 'START'
	clus_name = ""
	set_of_significant_cluster_names = set()
	
# Skip the header lines
	for line in kegg_lines:
# If the line is not empty and the first character is a '*', indicating a significant hit
		if line != '' and line [0] != '#' and line[0] == '*':

			info_list = line.split('\t')
			#break the line into pieces
			clus_name = ('_'.join(info_list[1].split('_')[:2]))
			set_of_significant_cluster_names.add(clus_name)
			#rejoin the two terms to get CLUSTER_X and call it 'clus_name'
			# Check if new cluster has been reached.  If so then write out the info for the previous lines.
			if clus_name != prev_clus_name:
				# Add cluster information to dictionary, write to first output file
				KO_info = ';'.join(KO_list) + '\t*' + '; '.join(KO_def_list)
				if prev_clus_name != "START":
					kegg_clus_dict[prev_clus_name] = KO_info
					outtab1_file.write(prev_clus_name + '\t' + KO_info + '\n')
				# Reset for new line of the file
				KO_list = []
				KO_def_list = []
				prev_clus_name = clus_name

			KO_list.append(info_list[2]) #add KO number to a list of KO numbers
			KO_def_list.append(info_list[6].strip('"\n'))

	outtab1_file.write("###INSIGNIFICANT BEST GUESS HITS###\n")
	for line in kegg_lines:
			if line != '' and line[0] != '#' and line[0] != "*":
				info_list_insignificant = line.split('\t')
				insignificant_clus_name = ('_'.join(info_list_insignificant[1].split('_')[:2]))
				if insignificant_clus_name not in set_of_significant_cluster_names:
					if insignificant_clus_name != prev_clus_name:
						insignificant_cluster_info = str(info_list_insignificant[2]+ '\t' + info_list_insignificant[6].strip('"\n') + " (evalue=" + info_list_insignificant[5]+(")"))
						if (float(info_list_insignificant[5]) < 1e-7):
							kegg_clus_dict[insignificant_clus_name] = insignificant_cluster_info
							outtab1_file.write(insignificant_clus_name + "\t" + insignificant_cluster_info + '\n')
					prev_clus_name = insignificant_clus_name


with open(input_cluster_table, 'r') as f:

	intab_lines = f.readlines()

	# Writing out the new header
	header_list = intab_lines[0].split('\t')
	header_list.insert(3, 'KO_number\tKEGG_definition')
	outtab2_file.write('\t'.join(header_list))

	for line in intab_lines[1:]:
		cluster_list = line.split('\t')
		cluster_id = cluster_list[0]
		if cluster_id in kegg_clus_dict:
			cluster_list.insert(3, kegg_clus_dict[cluster_id])
		else:
			cluster_list.insert(3, '*\tnone')
		outtab2_file.write('\t'.join(cluster_list))

outtab1_file.close()
outtab2_file.close()
