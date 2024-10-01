#!/usr/bin/env python

import os
import glob
import shutil 
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
from typing import Iterator, Tuple

	
def get_locus_tag(x): #this will give you the locus tag at the front of the locusid.  There are two general formats of locus tags.
	import re
	if x is None:
		return
	if "_" in x:
		locus_tag = str(x.split("_")[0]) #this works for locus tags formatted like STM14_#######
		return locus_tag
	else:
		locus = re.match(r"\D+", x) #this works for locus tags formatted like STM#### (no hyphen...letters before a string of #s)
		locus_tag = locus.group(0)
		return locus_tag

def line_format(sequence): #this function returns a FASTA seqeunce with 80 characters per line
		newline_list = []
		while sequence:
				newline = sequence[0:80]
				newline_list.append(newline)
				sequence = sequence[80:]
		output = "\n".join(newline_list)
		return output

def get_gbkey(description_line):
	descriptors = str(description_line[1:-1]) #remove the brackets at the end of the line
	descriptor_list = descriptors.split("] [") #further split the tags into individual items by splitting at the "] [" between each tag
	for item in descriptor_list:
		if "gbkey=" in str(item):
			descrip = (str(item)[6:])
			return descrip
	descrip = str("none")
	return descrip

def convert_fna_to_faa(filename):
	from Bio import SeqIO
	output_handle = open("faa/" + filename[4:].split(".")[0] + ".faa", "w")
	for seq_record in SeqIO.parse(filename, "fasta"):
		if get_gbkey(seq_record.description) == "CDS":
			output_handle.write(">" + seq_record.description[4:]+"\n")
			output_handle.write(line_format(str(seq_record.seq.translate()[:-1]))+"\n") #the [:-1] removes the * for the stop codon from the sequence record
	output_handle.close()
	print("Converted "+ filename + "to faa file")

def tree2list(directory: str) -> Iterator[Tuple[str, str, str]]:
	import os
	for i in os.scandir(directory):
		if i.is_dir():
			yield i.path
			yield from tree2list(i.path)
		else:
			yield i.path



'''
┌------------------------------------------------------------------------------------┐
| 1. Setting up directories, initializing report file, making copies of input files. |
└------------------------------------------------------------------------------------┘
'''

if __name__ == '__main__':

#set up directories we'll need
	
	directory_list = list(tree2list('.'))  #get current directories and put into list

	#stop the program if certain files and directories are not already in place
	if "./data" not in directory_list:
		print("\nCannot find a 'data' directory in the current directory.")
		print("Are you sure you're in the right directory?\n")
		print("Have you unpacked your ncbi_dataset file properly?\n")
		exit()

	if "./data/ncbi/data" not in directory_list:
		print("\nCannot find the folder with GCF files in it.")
		exit()

	if "./data/ncbi/master_table.tab" not in directory_list:
		print("\nCannot find the master_table.tab file in the data/ncbi directory.")
		print("Are you sure you're in the right directory?\n")
		print("Have you unpacked your ncbi_dataset file properly?\n")
		exit()

	#if not already in the current directory, make the directories 'report' and 'fna' and 'faa'
	if "./data/report" not in directory_list:
		os.mkdir("data/report")
	if "./data/fna" not in  directory_list:
		os.mkdir("data/fna")
	if "./data/fna" not in  directory_list:
		os.mkdir("data/faa")


	#get the names of all the fna files in the data/ncbi/data directory
	fna_sequence_file_list = glob.glob('data/ncbi/data/*/*.fna')

	#make the 'strainlist.txt' file by getting metadata in the ncbi masterLtable.tab file.
	strain_lookup_dictionary = {}
	fh = open('data/ncbi/master_table.tab','r')
	lines = fh.readlines()[1:] #start reading after the first line.
	for line in lines:
		splitline = line[:-1].strip().split('\t')
		assembly_id = splitline[0]
		assembly_index_long = assembly_id[4:]   #this number comes after "GCA_" or "GCF_"
		assembly_index = assembly_index_long.split(".")[0]
		species = splitline[1].strip()
		strain = splitline[2]
		bioproject_id = splitline[3]
		biosample_id = splitline[4]
		assembly_status = splitline[5]
		if assembly_status == 'Complete Genome':
			assembly_status = 'Complete'
		strain_info_line = (species + " " + strain + " [" + assembly_id + "; " + biosample_id + "; " + bioproject_id  + "; " + assembly_status + "]")
		strain_lookup_dictionary[assembly_index] = strain_info_line
	fh.close()

	#rename the fna files in the fna directory and build the strainlist.txt file
	fh3 = open('strainlist.txt','w')	
	temp_fna_filename_list = []

	for fna_filename in fna_sequence_file_list:
		filename_tag_1 = fna_filename.split("_")[1]
		filename_tag_2 = filename_tag_1.split(".")[0]	
		strain_information = strain_lookup_dictionary[filename_tag_2]
	
		fh2 = open(fna_filename,'r')
		line = fh2.readline()
		locus_tag_1 = line.split("locus_tag=")[1]
		locus_tag_2 = locus_tag_1.split("]")[0]
		locus_prefix = get_locus_tag(locus_tag_2)
		fh2.close()
	
	#now write out the line of the dictionary to the "strainist.txt" file
		strain_info_line = (locus_prefix + " | " + strain_information)
		fh3.write(strain_info_line + "\n")

	#now rename the fna file to the locus_prefix."temp_fna"
		temp_filenames_for_both_faa_and_fna_outputs = []
		new_temp_fna_filename_with_locus_id = str("data/fna/" + locus_prefix + "_temp.fna")
		new_faa_filename_with_locus_id = str("data/faa/" + locus_prefix + ".faa")
		shutil.copy(fna_filename, new_temp_fna_filename_with_locus_id)
		temporary_file_combine_list = []
		temporary_file_combine_list.append(new_temp_fna_filename_with_locus_id)
		temporary_file_combine_list.append(new_faa_filename_with_locus_id)
		temp_fna_filename_list.append(temporary_file_combine_list)
		temporary_file_combine_list = []

	#Now convert the temp_fna files to fna files with GC contents
	#also gather a list of potential pseudogenes
	output_handle_pseudogenes = open('data/potential_pseudogenes.fna', 'w')
	output_handle_pseudogenes.write("#The following CDS may be psuedogenes or cut off at the end of a contig:\n\n")

	for temp_fna_filename in temp_fna_filename_list:
		output_handle_fna = open(temp_fna_filename[0].split("_temp.")[0] + ".fna", "w")
		output_handle_faa = open(temp_fna_filename[1], "w")
		for record in SeqIO.parse(temp_fna_filename[0], "fasta"):
			old_description_line = record.description
			fragmented_description = old_description_line.split("[gbkey")
			GC_content = str(100*gc_fraction(record.seq))
			middle_fragment = "[pctGC="+ str('{:.4}'.format(GC_content)) + "] [gbkey"
			record.description = fragmented_description[0] + middle_fragment + fragmented_description[1]
			output_handle_fna.write(">" + record.description+"\n")
			output_handle_fna.write(line_format(str(record.seq)+"\n"))
			if len(record.seq) %3 != 0:
				output_handle_pseudogenes.write(record.description + "\n")
			#NOTE below we remove the >lcl from the sequence annotation
			output_handle_faa.write(">" + str(record.description)[4:]+"\n")
			output_handle_faa.write(line_format(str(record.seq.translate()[:-1]))+"\n")
		
		output_handle_faa.close()		
		output_handle_fna.close()
		
		remove_command = "rm " + temp_fna_filename[0]
		os.system(remove_command)
		
		concatenate_fna_files_command = "cat data/fna/*.fna > data/fna/all.fna"
		os.system(concatenate_fna_files_command)
