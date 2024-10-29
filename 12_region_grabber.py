#!/usr/bin/env python

"""
region_grabber.py
This script takes input from a list of clusterIDs or WPs or Locus Ids or accessions 
and gives back a range of nucleotides including the query gene and a specified 
number of bases upstream and downstream of the gene as well.  The output is in FASTA
nucleotide format and genbank flat-file format.	 The genes are oriented in their 5' to
3' direction by default.
"""

'''
┌----------------------------------------┐
| import modules						 |
└----------------------------------------┘
'''
import argparse
import os
import numpy as np
import pandas as pd
from Bio import SeqIO, Entrez

'''
┌------------------┐
| define functions |
└------------------┘
'''

def create_directories(base_dir='grab', sub_dirs=['full_gbff', 'gbff', 'fna'], suffix=''):
	# Check if the base directory exists
	directory_list = []
	if not os.path.exists(base_dir):
		# Create the base directory if it doesn't exist
		os.makedirs(base_dir)
		print(f"Created directory: {base_dir}")
	else:
		print(f"Directory '{base_dir}' already exists.")

	# Create the full_gbff directory without incrementing
	full_gbff_dir = os.path.join(base_dir, 'full_gbff')
	directory_list.append(full_gbff_dir)
	if not os.path.exists(full_gbff_dir):
		os.makedirs(full_gbff_dir)
		print(f"Created directory: {full_gbff_dir}")
	else:
		print(f"Directory '{full_gbff_dir}' already exists.")

	# Loop through gbff and fna directories and create them with increment

	for sub_dir in sub_dirs[1:]:  # Skip full_gbff
		new_dir_name = f"{sub_dir}{suffix}"
		index = 1

		# Find an available name with an index (i.e., gbff_suffix (1), fna_suffix (1), etc.)
		while os.path.exists(os.path.join(base_dir, new_dir_name)):
			new_dir_name = f"{sub_dir}{suffix}({index})"
			index += 1

		# Create the subdirectory
		os.makedirs(os.path.join(base_dir, new_dir_name))
		directory_list.append(f"{base_dir}/{new_dir_name}")
		print(f"Created directory: {new_dir_name}")

	return(directory_list)

'''
┌----------------------------------------┐
| main program start					 |
└----------------------------------------┘
''' 

if __name__ == '__main__':

# set up constants
	metadata_filename = "data/mmseq_output/cluster_metadata.tab"
	search_filename = 'searchterms.txt'
	output_foldername = 'regiongrab'
	add_to_gene_fiveprime = 0
	add_to_gene_threeprime = 0
	Entrez.email = "william.navarre@utoronto.ca"  # Always tell NCBI who you are
	instructions = "Put genes you want to grab in a file (default name 'searchterms.txt')\n" \
			"Use the -r flag to specify amount to add to 5' and 3' ends separated by a colon."

#parse the command line for user input
	parser = argparse.ArgumentParser(add_help=True, description=instructions)
	parser.add_argument('-r', '--range', action='store', dest='range', help = "Nucleotides to add before and after gene (default = '0:0').")
	parser.add_argument('-i', '--in', action='store', dest='search_term_file', help = "Name and path of file with search terms (default = 'searchterms.txt').")
	parser.add_argument('-s', '--search', action='store', dest='search_term', help = "A specific search term fed by user")
	parser.add_argument('-o', '--out', action='store', dest='output', help = "Output folder name (default = 'regiongrab').")
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

	if args.range:
		add_to_gene_fiveprime = int(args.range.split(":")[0])
		add_to_gene_threeprime = int(args.range.split(":")[1])
		output_foldername_suffix = f"({add_to_gene_fiveprime}_{add_to_gene_threeprime})"

#now create directories with a main grab directory with a gbff and fna directory inside
	folder_list = create_directories(output_foldername, suffix=output_foldername_suffix)
	print(folder_list)
	full_gbff_directory = folder_list[0]
	gbff_directory = folder_list[1]
	fna_directory = folder_list[2]
	
	search_terms = []
	if user_defined_search_term:
		search_terms.append(user_defined_search_term.rstrip())

	elif search_filename:
		with open(search_filename , 'r') as f:
			for line in f:
				search_terms.append(line.rstrip()) 

	input_table = pd.read_csv('data/mmseq_output/cluster_metadata.tab', sep='\t')	
	output_table = input_table[input_table.isin(search_terms).any(axis=1)]
	output_table.to_csv(f'{output_foldername}/search_output_metadata.tab', sep='\t', index=False)

#now read the new metatdata table, line by line, and extract the information necessary 
#to feed to the NCBI entrez server.
	for index, row in output_table.iterrows():
		locus_id = row["LOCUS_ID"]
		nucleotide_accession = row["NCBI_NUCLEOTIDE"]
		gene_location = row["LOCATION"]
		gene_direction = row["DIRECTION"]
		annotation = row["NCBI_ANNOTATION"]

#check below to see if gene is truncated or contains a "<" or ">" symbol in its location
		location_start = (gene_location.split("..")[0])
		if "<" in location_start or ">" in location_start:
			location_start = int(location_start[1:])-1
		else: location_start = int(location_start)-1
		location_end = (gene_location.split("..")[1])
		if "<" in location_end or ">" in location_end:
			location_end = int(location_end[1:])
		else: location_end = int(location_end)

#here we add bases to retrieve to the 5-prime and 3-prime ends of the query gene
		if gene_direction == "F":
			location_start = location_start - add_to_gene_fiveprime
			location_end = location_end + add_to_gene_threeprime
		if gene_direction == "R":
			location_start = location_start - add_to_gene_threeprime
			location_end = location_end + add_to_gene_fiveprime

#if adding more sequence to the 5' end of the sequence results in a negative number (you
#have reached beyond the contig), the retrieval will fail.	So we set the sequence range
#to start at the first base
		if location_start <= 0:
			location_start = 1
		print(f"{nucleotide_accession} from {location_start} to {location_end} ({gene_direction})")


#To limit the number of times you have to re-download a given gbff file from NCBI if you're
#trying different parsings, we save the full length genome or contig in the full_gbff directory
#if the contig/genome was previously downloaded you just parse that file instead.

		filename = f"{full_gbff_directory}/{nucleotide_accession}.gbff"
		if not os.path.exists(filename):
			handle = Entrez.efetch(db="nucleotide", id=nucleotide_accession, rettype="gbwithparts", retmode="text")
			record = SeqIO.read(handle, "genbank")
			SeqIO.write(record, filename, "genbank")
			print(f"downloaded {nucleotide_accession} from NCBI and saved as {filename}.")
		else:
			record = SeqIO.read(filename, "genbank")
			print(f"{filename} exists and will be parsed.")

#with BioPython, when you modify a genbank record you can lose its metadata and annotations
#to prevent this I first save the metadata into variables a, b, c, and d...THEN modify
#the genbank record (reversing it or truncating it).  Then the metadata that
#we saved is put back into the modified genbank record.

		a=record.id
		b=record.name
		c=record.description.split(',')[0]
		d=record.annotations
		
#if adding more sequence to the 3' end of the sequence extends beyond the end of the contig
#the retrieval will still work but the user will think they have grabbed more sequence than
#they acxtually have.  So we set the sequence range to end at the last base of the contig
		if location_end > (len(record)):
			location_end = (len(record))

		region = record[location_start:location_end]

#here we put the metadata back into the sequence file.
		if gene_direction == "F":
			region.id = a
			region.name = b
			modified_description = c + f" REGION: {location_start}..{location_end}"
			region.description = modified_description
			region.annotations = d

		if gene_direction == "R":
			region = region.reverse_complement()

#for whatever reason when you reverse a sequence in Biopython it loses the fact that
#the molecule type is DNA.	We have to re-insert that information here.
			region.annotations["molecule_type"] = "DNA"
			region.id = a
			region.name = b
			modified_description = c + f" REGION: complement({location_start}..{location_end})"
			region.description = modified_description
			region.annotations = d

#Write out the genbank flat file
		gbff_filename = f"{gbff_directory}/{locus_id}_({annotation})_add_{add_to_gene_fiveprime}up-{add_to_gene_threeprime}dwn_({gene_direction}).gbff"
		with open(gbff_filename, "w") as output_handle:
			SeqIO.write(region, output_handle, "genbank")

#write out the fasta file
		fna_filename = f"{fna_directory}/{locus_id}_({annotation})_add_{add_to_gene_fiveprime}-{add_to_gene_threeprime}_({gene_direction}).fna"
		with open(fna_filename, "w") as output_handle2:
			SeqIO.write(region, output_handle2, "fasta")