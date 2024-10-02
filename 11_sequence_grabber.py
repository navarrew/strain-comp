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

# 	parser = argparse.ArgumentParser(add_help=True, description=instructions)
# 	parser.add_argument('-d', '--data', action='store', dest='metadata_file', help = "Pct identity for clustering (30-99, default = 80).")
# 	parser.add_argument('-q', '--query', action='store', dest='search_terms', help = "Pct mutual coverage of protein length (30-99, default = 90).")
# 	parser.add_argument('-o', '--out', action='store', dest='output', help = "Name prefix for clusters (DEFAULT = CLUSTER.")
# 	
# 	args = parser.parse_args()
# 
# 	if args.metadata_file:
# 		clustering_pct = int(args.clustering_pct)
# 	else: clustering_pct = 80
# 
# 	if args.coverage_length:
# 		coverage_length = int(args.coverage_length)
# 	else: coverage_length = 90
# 
# 	if args.clustername:
# 		cluster_prefix = args.clustername + "_"
# 	else: cluster_prefix = "CLUSTER_"

	search_terms = []
	with open('searchterms.txt', 'r') as f:
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
