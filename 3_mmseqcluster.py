#!/usr/bin/env python

"""
mmseqcluster.py
This script takes input from an 'faa' folder filled with amino acid FASTA files
for different strains.	It concatenates these to make a single large file and 
performs clustering on them using the mmseq2 program.
"""

'''
┌----------------------------------------┐
| import modules						 |
└----------------------------------------┘
'''

import os
from datetime import datetime
from Bio import SeqIO
import argparse

'''
┌----------------------------------------┐
| define functions						 |
└----------------------------------------┘
'''


def line_format(seq):
	"""
	Returns a FASTA sequence with 80 characters per line.
	"""

	newline_list = []
	while seq:
		newline = seq[0:80]
		newline_list.append(newline)
		seq = seq[80:]
	output = '\n'.join(newline_list)
	return output


def pseudoremover(infile):
	"""
	Takes concatenated .faa protein sequences as input,
	removes pseudogenes to output "data/faa/pseudofree.faa".
	"""

	# Create a list of sequence records, where each protein is larger than 30 AA
	seq_rec_list = []
	for seq_rec in SeqIO.parse(infile, 'fasta'):
		if len(seq_rec.seq) > 30:
			seq_rec_list.append(seq_rec)

	# Make a new file to write out to.
	with open('data/faa/pseudofree.faa', 'w') as outfile:
		for seq_rec in seq_rec_list:
			if 'PSEUDOGENE' not in seq_rec.description:
				if '[pseudo=true]' not in seq_rec.description:
					outfile.write('>' + seq_rec.description + '\n')
					outfile.write(line_format(str(seq_rec.seq) + '\n'))


def get_cluster_list(infile):
	"""
	Given a FASTA-formatted set of protein clusters from mmseqs2,
	produces a nested list of clusters and their members,
	where the outermost list is sorted by cluster size.
	"""

	with open(infile) as f:
		lines = [line.rstrip() for line in f]

	i = 1  # Start an index, skip the very first line in the file.
	member_list, cluster_list = [], []
	while i < len(lines):
		firstline = lines[i]
		secondline = lines[i + 1]
		if secondline[0] == '>':  # If there are two > in a row it means a new cluster has started
			i += 1	# Advance the "frame" by one so you get back on track
			cluster_list.append(member_list)  # Add previously assembled cluster to a growing list of clusters
			member_list = []  # Reset the member_list to get ready to add new sequences to
		else:  # Otherwise the first line should be the definition line and the second line should be the sequence
			member_list.append(firstline + '\n' + line_format(secondline))	# Append the sequence to the cluster
			i += 2	# Skip ahead two lines
	cluster_list.append(member_list)  # The last cluster in the queue now needs to be added to the cluster list

	# Sort the cluster_list by size, those with the greatest number of members appear first in the list
	cluster_list = sorted(cluster_list, key=len, reverse=True)

	return cluster_list


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


def get_gene_prot_names(header):
	"""
	Takes the header from each member within a cluster and extracts the gene and protein names.
	"""

	member_fragments = header.split(' [')
	gene_set = set()
	protein_set = set()

	for fragment in member_fragments:
		if 'gene=' in fragment:
			gene_set.add(str(fragment[5:-1].split('_')[0]))
		if 'protein=' in fragment:
			protein_set.add(str(fragment[8:-1]))

	if len(gene_set) > 1:
		gene_set.discard('none')
	if len(gene_set) == 0:
		gene_set.add('none')

	gene_names = ', '.join(gene_set)
	protein_names = ', '.join(protein_set)

	return gene_names, protein_names
	
	
'''
┌----------------------------------------┐
| main program							 |
└----------------------------------------┘
'''


if __name__ == '__main__':

	print("\nPipeline for genomic comparisons",
		  "Navarre lab updated summer 2024",
		  sep='\n')
		  
 # Set up directories, if not already present
	directory_list = []
	for entry in os.scandir('data'):
		if entry.is_dir():
			directory_list.append(entry.name)
	for dir in ['mmseqdb', 'mmseq_output', 'report']:
		if dir not in directory_list:
			os.mkdir('data/'+dir)
		if 'faa' not in directory_list:
			print ("\nCannot find the 'faa' directory (FASTA protein files) in this directory.\n")
			exit()

	#get variable values and command line help:
	instructions = "mmseqcluster.py\n" \
					"This will use mmseq2 based clustering to group proteins together." \
					"You use the -p flag to set percent identity for clusters (default 80)." \
					"You use the -l flag to set length of mutual overlap (default = 90)"
					
	parser = argparse.ArgumentParser(add_help=True, description=instructions)
	parser.add_argument('-p', '--pctid', action='store', dest='clustering_pct', help = "Pct identity for clustering (30-99, default = 80).")
	parser.add_argument('-l', '--covlength', action='store', dest='coverage_length', help = "Pct mutual coverage of protein length (30-99, default = 90).")
	parser.add_argument('-n', '--name', action='store', dest='clustername', help = "Name prefix for clusters (DEFAULT = CLUSTER.")
	
	args = parser.parse_args()

	if args.clustering_pct:
		clustering_pct = int(args.clustering_pct)
	else: clustering_pct = 80

	if args.coverage_length:
		coverage_length = int(args.coverage_length)
	else: coverage_length = 90

	if args.clustername:
		cluster_prefix = args.clustername + "_"
	else: cluster_prefix = "CLUSTER_"

	# Write the date and time of execution, and the % identity cutoff to the report file
	report = open('data/report/report.txt', 'a')
	report.write("\n\n***MMSEQCLUSTER.PY***\nDATE AND TIME: "  + str(datetime.now()) +
				 "\n\nMMSEQS2 PERCENT IDENTITY CUTOFF: " + str(clustering_pct) + "%" +
				 "\nCLUSTER PREFIX: " + str(cluster_prefix[:-1]))

	'''
	┌-------------------------------------------------------------------------------┐
	| 2. Clustering proteins by mmseq2.                                             |
	└-------------------------------------------------------------------------------┘
	'''

	# Concatenate all protein sequences in faa files, generating "temp.faa"
	
	os.system('cat data/faa/*.faa > data/faa/temp.faa')

	# Remove pseudogenes from concatenated protein sequences, generating 'pseudofree.faa'
	pseudoremover('data/faa/temp.faa')

	# Running mmseqs2 to obtain protein clusters
	print("\nClustering with mmseqs2...\n")

	cmd1 = 'mmseqs createdb data/faa/pseudofree.faa data/mmseqdb/DB'
	cmd2 = 'mmseqs cluster data/mmseqdb/DB data/mmseqdb/clusteredDB data/mmseqdb/tmp --min-seq-id ' + str(float(clustering_pct) / 100) + ' -c ' + str(float(coverage_length) / 100) + ' --cov-mode 0'
	cmd3 = 'mmseqs createseqfiledb data/mmseqdb/DB data/mmseqdb/clusteredDB data/mmseqdb/clu_seq'
	cmd4 = 'mmseqs result2flat data/mmseqdb/DB data/mmseqdb/DB data/mmseqdb/clu_seq data/mmseq_output/clustered_sequences.fasta'

	for cmd in cmd1, cmd2, cmd3, cmd4:
		os.system(cmd)

	# Write the commands used for mmseqs2 into the report file, then close it
	report.write("\nMMSEQS2 FULL COMMANDS:\n" + cmd1 + "\n" + cmd2 + "\n" + cmd3 + "\n" + cmd4)


	# Output is a file called "clustered_sequences.fasta"

	# Remove the temporary faa files
	os.system('rm data/faa/temp.faa && rm data/faa/pseudofree.faa')
	
	'''
	┌-------------------------------------------------------------------------------┐
	| 3. Assigning each protein cluster a name in order of their number of members. |
	└-------------------------------------------------------------------------------┘
	'''

	print("\nAssigned each cluster a name in order of abundance.")
	print("Each protein cluster will have the prefix: " + cluster_prefix[:-1])

	# Get a nested list of clusters and their members, ordered by cluster size
	cluster_list = get_cluster_list('data/mmseq_output/clustered_sequences.fasta')

	# Iterate through the list of clusters and rename them, produce three output files
	summ_file = open('data/mmseq_output/cluster_summary_with_sequences.faa', 'w')
	repr_file = open('data/mmseq_output/cluster_representative_sequences.faa', 'w')
	head_file = open('data/mmseq_output/cluster_header_info.tab', 'w')
	cluster_metadata_file = open('data/mmseq_output/cluster_metadata.tab', 'w')
	cluster_metadata_file.write('LOCUS_ID\tCLUSTER_ID\tPROTEIN_ID\tACCESSION\tNCBI_ANNOTATION\n')

	# Set index for cluster number
	j = 1

	summ_file.write('#Number of clusters: ' + str(len(cluster_list)) + '\n')
	summ_file.write('#Cluster prefix: ' + cluster_prefix + '\n')

	for cluster in cluster_list:

		# Write the cluster number, number of members, and headers and protein sequences of all members in each cluster
		cluster_index = str(j).zfill(6)
		summ_file.write('#' + cluster_prefix + cluster_index +
						' [members= ' + str(len(cluster)) + ']\n' +
						(str('\n'.join(cluster)) + '\n'))

		# Write the header and protein sequence of the most representative member of each cluster
		repr_file.write(str('>' + cluster_prefix + cluster_index + '_' + cluster[0][1:] + '\n'))

		# Iterate through each member within a cluster and extract the gene and protein names from the header
		# Write the cluster number, gene and protein names, number of members, and headers of all members
		header_list = []
		for member in cluster:
			header = member.split('\n')[0]
			header_list.append(header)
			gene_names, protein_names = get_gene_prot_names(header)
		head_file.write(cluster_prefix + cluster_index + '\t' +
						gene_names + '\t' +
						protein_names + '\t' +
						str(len(cluster)) + '\t' +
						str('\t'.join(header_list)) + '\n')

		for header2 in header_list:
			accession_in_header2 = "lcl|" + header2[0:header2.find(" ")] #find gives the location of the first space and you trim the string from zero to the first space.
			locus_tag_in_header2_front = header2[(header2.find("locus_tag=")+10):]
			locus_tag_in_header2 = locus_tag_in_header2_front[:locus_tag_in_header2_front.find("] ")]
			WP_tag_in_header2_front = header2[(header2.find("protein_id=")+11):]
			WP_tag_in_header2 = WP_tag_in_header2_front[:WP_tag_in_header2_front.find("] ")]
			info_tag_in_header2_front = header2[(header2.find("protein=")+8):]
			info_tag_in_header2 = info_tag_in_header2_front[:info_tag_in_header2_front.find("] ")]
			cluster_metadata_file.write(locus_tag_in_header2 + '\t' + cluster_prefix + cluster_index + '\t' + WP_tag_in_header2 + '\t' + accession_in_header2 + '\t' + info_tag_in_header2 + "\n")
			
		# Increase the number of clusters processed by 1
		j += 1

	# Close the three output files
	summ_file.close()
	repr_file.close()
	head_file.close()
	cluster_metadata_file.close()

	# Write the total number of clusters to the report file
	report.write("\n\nTotal number of protein clusters: " + str(j))
	report.close()

	print("\nFinished!\nTotal number of protein clusters: " + str(j)+"\n")

