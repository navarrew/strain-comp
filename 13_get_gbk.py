#!/usr/bin/env python3
# NC_002695.1 - you can use this as a test accession

import argparse
import os
from Bio import SeqIO
from Bio import Entrez
Entrez.email = "william.navarre@utoronto.ca"  # Always tell NCBI who you are


# the following fxn will query the user for 
def retrieve_request_info():
	print ("\n\n\n#######gbk_segment.py######\n")
	print ("This script will retrieve a portion of a genbank seqeunce,")
	print ("either from a file in this directory or downloaded from NCBI\n")
	
	filename = input("NCBI accession or local filename: ")

	all_or_part = input("All of the file or just a region (A or R)? ")

	if (	all_or_part == "R" or 
			all_or_part =="Region" or
			all_or_part == "r" or
			all_or_part == "region" or
			all_or_part == "REGION"):
		print("\nYou have opted to grab only part of the genome.  Please select the region.\n")
		strt = int(input("start: "))
		end = int(input("end: "))
		grab = 1
	else: 
		grab = 0
		strt = 0
		end = 0

	direction = input("In fwd or rev direction (F or R)? ")
	
	output = input("Output GBFF or FASTA nucleotide format (G or F)? ")
	
	return filename, direction, output, grab, strt, end

instructions = 		"***GRAB_GBK.py***\
					\nThis will retrieve all or part of any NCBI genbank file via its \
					nucleotide accession number (online) or will extract a region from\
					any .gbk or .gbff file you name that is in your directory.\n"

# STEP 1 = the following argparse fxn reads the command line for specific inputs
parser = argparse.ArgumentParser(add_help=True, description=instructions)
parser.add_argument('-in', action='store', dest='infile', help = "name of the input file or accession number")
parser.add_argument('-direct', action='store', dest='direction', help = "f or r for 'forward' or 'reverse'")
parser.add_argument('-out', action='store', dest='output', help = "fna or gbk for 'fasta nucleotide' or 'genbank'")
parser.add_argument('-range', action='store', dest='range', help = "(optional) region in format (start..end)")

args = parser.parse_args()
inlist = []
region = str()

# STEP 2 - if there is nothing in the command line - run the fxn above to query the user
if args.infile == None:
	inlist = retrieve_request_info()
	filename = inlist[0]
	direction = inlist[1]
	output = inlist[2]
	grab = inlist[3]
	strt = inlist[4]
	end = inlist[5]

else: # STEP 3 - if there is something in the command line parse it.
	filename = args.infile
	if args.direction:
		direction = args.direction
	else: direction = "fwd"
	if args.output:
		output = args.output
	else: output = "gbk"
	grab = 0
	if args.range != None:
		range = args.range
		strt = int(range.split('..')[0])
		end = int(range.split('..')[1])
		grab = 1

#STEP 4 - if the sequence file is not in the directory - go to NCBI entrez and download it and save it as a gbff file

#if there's a suffix like gbk or gbff at the end of the filename, remove it.
if filename.split('.')[-1] == "gbk" or filename.split('.')[-1] == "gbff":
	filename = ".".join(filename.split('.')[0:-1])
	print(filename)
	
filename_appended_gbk = filename + ".gbk"
filename_appended_gbff = filename + ".gbff"

if os.path.isfile(filename_appended_gbff):
	filename_appended = filename_appended_gbff

if os.path.isfile(filename_appended_gbk):
	filename_appended = filename_appended_gbk

if not os.path.isfile(filename_appended_gbk) and not os.path.isfile(filename_appended_gbff):
    print("Downloading...")
    net_handle = Entrez.efetch(db="nucleotide", id=filename, rettype="gbwithparts", retmode="text")
    filename = filename + ".gbff"
    out_handle = open(filename, "w")
    out_handle.write(net_handle.read())
    out_handle.close()
    net_handle.close()
    filename_appended = filename_appended_gbff
    print("Saved the downloaded sequence file as: " + str(filename_appended))

if grab == 1 or (direction == "rev" or direction == "r" or direction == "REV" or direction == "Rev" or direction == "R"):	
	print("Parsing...")
	
	# STEP 5 -  now get the file (either from the directory or from NCBI) and read into a "record" Seq_record object
	record = SeqIO.read(filename_appended, "gb")
	
	#STEP 6 - save the header info.  NOTE that the reverse complement and segmentation of 
	# a gbk file, for whatever reason, loses its header info.
	# Here I put in extra steps to save the header info before reversal and then put it back
	# after the sequence was reversed.
	
	a=record.id
	b=record.name
	c=record.description.split(',')[0]
	d=record.annotations
	
	if (direction == "rev" or direction == "r" or direction == "REV" or direction == "Rev" or direction == "R"):
		reverse_index = 1 # a reverse_index of "1" indicates the sequence should be reverse-complemented.
	else: reverse_index = 0
		
	if grab == 1: # a grab variable of "1" indicates that you only want a portion of the sequence
		real_start_point = strt - 1 #we have to adjust the index range when slicing the sequence by one.
		region = record[real_start_point:end]
	
	else: 
		region = record
		region.desription = c
		
	if reverse_index == 1:
		region = region.reverse_complement()
	
	region.id = a
	region.name = b
	region.description = c
	region.annotations = d
	
	#now to make the appendage to the filename and modify the sequence description to 
	# inform the reader that you only have a segment of the genome, etc...
	appendage = ""
	if grab == 1:
		appendage = "_"+str(strt)+"_to_"+str(end)
		region.description = c + ", (bases: "+str(strt)+":"+str(end)+")"
	
	if reverse_index == 1:
		appendage = appendage + "(rev)"
		region.description = region.description + " (reverse orientation)"


	if output == "G" or  output == "g" or output == "gbk" or output == "gbff":
		SeqIO.write(region, (filename + appendage + ".gbff"), "genbank")
	
	if output == "F" or output == "f" or output == "fna" or output == "fasta":
		SeqIO.write(region, (filename + appendage + ".fna"), "fasta")
