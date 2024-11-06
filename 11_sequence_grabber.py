#!/usr/bin/env python

"""
sequence_grabber.py
This script takes input from a file or command line of clusterIDs or WPs or 
locusIds or accessions or even NCBI annotations and spits out the particular genes
that meet the criteria both in fna and faa format.

These will appear in a 
"""

'''
┌----------------------------------------┐
| import modules						 |
└----------------------------------------┘
'''
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
import os

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

def create_directories(base_dir='genegrab'):
	# Check if the base directory exists
	if not os.path.exists(base_dir):
		# Create the base directory if it doesn't exist
		os.makedirs(base_dir)
		print(f"Created directory: {base_dir}")
	else:
		print(f"Directory '{base_dir}' already exists.")
'''
┌----------------------------------------┐
| main program start					 |
└----------------------------------------┘
'''	
	
if __name__ == '__main__':

	metadata_filename = "data/mmseq_output/cluster_metadata.tab"
	search_filename = 'searchterms.txt'
	output_foldername = 'genegrab'
	instructions = "Put genes you want to grab in a file (default name 'searchterms.txt')\n" \
				"or use the -s flag and separate your terms with a comma and enclose\n" \
				"in quotations if there's spaces in one of your search terms.\n"\
				"e.g. -s 'class C sortase,CLUSTER_000001'"

#parse the command line for user input
	parser = argparse.ArgumentParser(add_help=True, description=instructions)
	parser.add_argument('-i', '--in', action='store', dest='search_term_file', help = "Name and path of file with search terms (default = 'searchterms.txt').")
	parser.add_argument('-s', '--search', action='store', dest='search_term', help = "A specific search term fed by user")
	parser.add_argument('-o', '--out', action='store', dest='output', help = "Output folder name (default = 'genegrab').")
	parser.add_argument('-m', '--meta', action='store', dest='metadata_file', help = "Name and path of input metatdata file (default = data/mmseq_output/cluster_metadata.tab).")

	args = parser.parse_args()

	if args.metadata_file:
		metadata_filename = args.metadata_file

	if args.search_term:
		user_defined_search_term = args.search_term
		print(user_defined_search_term)
		
	if args.search_term_file:
		search_filename = args.search_term_file

	if args.output:
		output_foldername = args.output

	create_directories(output_foldername)

#now create a list of the terms we will search for in the cluster_metadata.tab file
	search_terms = []
	if args.search_term:
		search_terms = user_defined_search_term.split(",")

	elif search_filename:
		with open(search_filename , 'r') as f:
			for line in f:
				search_terms.append(line.rstrip()) 

	input_table = pd.read_csv('data/mmseq_output/cluster_metadata.tab', sep='\t')	
	output_table = input_table[input_table.isin(search_terms).any(axis=1)]

# Set the initial metadata output filename
	index = 0
	metadata_output_filename = os.path.join(output_foldername, f"{output_foldername}_metadata.tab")
# Increment the filename until an available one is found
	while os.path.exists(metadata_output_filename):
	    index += 1
	    metadata_output_filename = os.path.join(output_foldername, f"{output_foldername}_metadata({index}).tab")

	output_table.to_csv(metadata_output_filename, sep='\t', index=False)

	accession_list = list(output_table['ACCESSION'])

# Set the initial protein output filename
	index = 0
	protein_output_filename = os.path.join(output_foldername, f"{output_foldername}_protein.faa")
# Increment the filename until an available one is found
	while os.path.exists(protein_output_filename):
	    index += 1
	    protein_output_filename = os.path.join(output_foldername, f"{output_foldername}_protein({index}).faa")

# Set the initial nucleotide output filename
	index = 0
	nucleotide_output_filename = os.path.join(output_foldername, f"{output_foldername}_nucleotide.fna")
# Increment the filename until an available one is found
	while os.path.exists(nucleotide_output_filename):
	    index += 1
	    nucleotide_output_filename = os.path.join(output_foldername, f"{output_foldername}_nucleotide({index}).fna")

	protein_output_file = open(protein_output_filename, 'w')
	nucleotide_output_file = open(nucleotide_output_filename, 'w')

	input_seq_iterator = SeqIO.parse("data/fna/all.fna", "fasta")
	for record in input_seq_iterator:
		if record.id in accession_list:
			nucleotide_output_file.write(record.format("fasta"))
			protein_output_file.write(">" + record.description[4:] + "\n")
			protein_output_sequence = line_format(str(record.seq.translate()[:-1]))
			protein_output_file.write(protein_output_sequence + "\n")
