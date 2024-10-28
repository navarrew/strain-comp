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
import os
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio import Entrez

'''
┌----------------------------------------┐
| define functions						 |
└----------------------------------------┘
'''

def create_subdirectory(base_name):
    # Get current working directory
    cwd = os.getcwd()
    
    # Start with the base name
    new_dir = base_name
    suffix = 1
    
    # Check if the directory already exists and increment the name if it does
    while os.path.exists(os.path.join(cwd, new_dir)):
        new_dir = f"{base_name}({suffix})"
        suffix += 1

    # Create the new directories for the region grabs
    os.makedirs(os.path.join(cwd, new_dir))
    new_dir_gbff = f"{new_dir}/gbff"
    os.makedirs(os.path.join(cwd, new_dir_gbff))
    new_dir_fna = f"{new_dir}/fna"
    os.makedirs(os.path.join(cwd, new_dir_fna))
    print(f"Directory '{new_dir}' created.")
    return new_dir_gbff, new_dir_fna


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
	Entrez.email = "navarrelab@gmail.com"  # Always tell NCBI who you are
	instructions = "Put genes you want to grab in a file (default name 'searchterms.txt')\n" \
			"Use the -r flag to specify amount to add to 5' and 3' ends separated by a colon."
	
	parser = argparse.ArgumentParser(add_help=True, description=instructions)
	parser.add_argument('-m', '--meta', action='store', dest='metadata_file', help = "Name and path of input metatdata file (default = data/mmseq_output/cluster_metadata.tab).")
	parser.add_argument('-i', '--in', action='store', dest='search_terms', help = "Name and path of file with search terms (default = 'searchterms.txt').")
	parser.add_argument('-r', '--range', action='store', dest='range', help = "Nucleotides to add before and after gene (default = '0:0').")
	parser.add_argument('-o', '--out', action='store', dest='output', help = "Output folder name (default = 'grab').")
	
	args = parser.parse_args()

	if args.metadata_file:
		metadata_filename = args.metadata_file

	if args.search_terms:
		search_filename = args.search_terms

	if args.output:
		output_foldername = args.output

	if args.range:
		add_to_gene_fiveprime = int(args.range.split(":")[0])
		add_to_gene_threeprime = int(args.range.split(":")[1])
		output_foldername = f"{output_foldername}_{add_to_gene_fiveprime}_to_{add_to_gene_threeprime}"

#now create directories with a main grab directory with a gbff and fna directory inside
	gbff_directory, fna_directory = create_subdirectory(output_foldername)

	search_terms = []
	with open(search_filename , 'r') as f:
		for line in f:
			search_terms.append(line.rstrip()) 

	input_table = pd.read_csv('data/mmseq_output/cluster_metadata.tab', sep='\t')	
	output_table = input_table[input_table.isin(search_terms).any(axis=1)]
	output_table.to_csv(f'{output_foldername}/search_output_metadata.tab', sep='\t')


#now read the new metatdata table, line by line, and extract the information necessary 
#to feed to the NCBI entrez server.
	for index, row in output_table.iterrows():
		locus_id = row["LOCUS_ID"]
		nucleotide_accession = row["NCBI_NUCLEOTIDE"]
		gene_location = row["LOCATION"]
		gene_direction = row["DIRECTION"]
		annotation = row["NCBI_ANNOTATION"]

#check below to see if gene is truncated.
		if pd.isna(row["NOTE"]): #if gene is not truncated then this should yield 'true'
			location_start = int(gene_location.split("..")[0])-1
			location_end = int(gene_location.split("..")[1])

#here we add bases to retrieve to the 5-prime and 3-prime ends of the query gene
			if gene_direction == "F":
				location_start = location_start - add_to_gene_fiveprime
				location_end = location_end + add_to_gene_threeprime
			if gene_direction == "R":
				location_start = location_start - add_to_gene_threeprime
				location_end = location_end + add_to_gene_fiveprime

#if adding more sequence to the 5' end of the sequence results in a negative number (you
#have reached beyond the contig), the retrieval will fail.  So we set the sequence range
#to start at the first base
			if location_start <= 0:
				location_start = 1
				
			handle = Entrez.efetch(db="nucleotide", id=nucleotide_accession, rettype="gbwithparts", retmode="text")
			record = SeqIO.read(handle, "genbank")
			
#with BioPython, when you modify a genbank record you can lose its metadata and annotations
#to prevent this I first save the metadata into variables a, b, c, and d...THEN modify
#the genbank record (reversing it or truncating it).  Then I add the metadata that
#I saved back into the modified genbank record.

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

			if gene_direction == "F":
				region.id = a
				region.name = b
				modified_description = c + f" REGION: {location_start}..{location_end}"
				region.description = modified_description
				region.annotations = d

			if gene_direction == "R":
				region = region.reverse_complement()
#for whatever reason when you reverse a sequence in Biopython it loses the fact that
#the molecule type is DNA.  We have to re-insert that information here.
				region.annotations["molecule_type"] = "DNA"
				region.id = a
				region.name = b
				modified_description = c + f" REGION: complement({location_start}..{location_end})"
				region.description = modified_description
				region.annotations = d

			print(f"{nucleotide_accession} from {location_start} to {location_end} ({gene_direction})")

#Write out the genbank flat file
			gbff_filename = f"{gbff_directory}/{locus_id}_({annotation})_add_{add_to_gene_fiveprime}up-{add_to_gene_threeprime}dwn_({gene_direction}).gbff"
			with open(gbff_filename, "w") as output_handle:
				SeqIO.write(region, output_handle, "genbank")

#write out the fasta file
			fna_filename = f"{fna_directory}/{locus_id}_({annotation})_add_{add_to_gene_fiveprime}-{add_to_gene_threeprime}_({gene_direction}).fna"
			with open(fna_filename, "w") as output_handle2:
				SeqIO.write(region, output_handle2, "fasta")