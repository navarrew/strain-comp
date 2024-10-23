#!/usr/bin/env python

"""
grabber.py
This script takes input from a list of clusterIDs or WPs or locusIds or accessions and gives back what you want.
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

'''
┌----------------------------------------┐
| main program start					 |
└----------------------------------------┘
'''	
	
if __name__ == '__main__':
	instructions = "TO BE FILLED IN LATER"

	parser = argparse.ArgumentParser(add_help=True, description=instructions)
	parser.add_argument('-m', '--met', action='store', dest='metadata_file', help = "Name and path of metatdata file (default = data/mmseq_output/cluster_metadata.tab).")
	parser.add_argument('-i', '--in', action='store', dest='search_terms', help = "Name and path of file with search terms (default = 'searchterms.txt').")
# 	parser.add_argument('-o', '--out', action='store', dest='output', help = "Output file names.")
	
	args = parser.parse_args()

	if args.metadata_file:
		metadata_filename = args.metadata_file
	else: metadata_filename = "data/mmseq_output/cluster_metadata.tab"

	if args.search_terms:
		search_filename = args.search_terms
	else: search_filename = 'searchterms.txt'

# 	if args.output:
# 		output_filename = args.output
# 	else: output_filename = "XXXX"

	search_terms = []
	with open(search_filename, 'r') as f:
		for line in f:
			search_terms.append(line.rstrip()) 

	input_table = pd.read_csv('data/mmseq_output/cluster_metadata.tab', sep='\t')	
	output_table = input_table[input_table.isin(search_terms).any(axis=1)]
	output_table.to_csv('searchoutput.tab', sep='\t')
	accession_list = list(output_table['ACCESSION'])

	protein_output_file = open('search_proteins.faa', 'w')
	nucleotide_output_file = open('search_nucleotide.fna', 'w')

	input_seq_iterator = SeqIO.parse("data/fna/all.fna", "fasta")
	for record in input_seq_iterator:
		if record.id in accession_list:
			nucleotide_output_file.write(record.format("fasta"))
			protein_output_file.write(">" + record.description[4:] + "\n")
			protein_output_sequence = line_format(str(record.seq.translate()[:-1]))
			protein_output_file.write(protein_output_sequence + "\n")
