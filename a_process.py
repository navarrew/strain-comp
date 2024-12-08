#!/home/navarrelab/anaconda3/bin/python
"""
Process genomic data from NCBI into a format that can be used for strain comparisons.

This module can take either zipped or unzipped genbank or cds_from_genomic.fna files 
and process them for further analysis by the strain-comp pipeline.

Detailed instructions are on GitHub.
"""

import os
import re
import glob
import shutil
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import warnings
from Bio import BiopythonWarning
from typing import Iterator


def tree2list(directory: str) -> Iterator[str]:
	for i in os.scandir(directory):
		if i.is_dir():
			yield i.path
			yield from tree2list(i.path)
		else:
			yield i.path


def directory_setter(directory_name):
	os.makedirs(directory_name, exist_ok=True)


def zipextractor(myzip, file_name, output_directory_name, tail):
	"""Extracts a file from a zip archive and renames it."""
	fileaccession = os.path.dirname(file_name).split('/')[-1]
	myzip.getinfo(file_name).filename = fileaccession + tail
	myzip.extract(file_name, path=output_directory_name)


def ncbi_unzip(zip_filename='ncbi_dataset.zip'):
	"""Unzips NCBI dataset and organizes files."""
	from zipfile import ZipFile
	with ZipFile(zip_filename) as myzip:
		list_of_files = myzip.namelist()
		jsonl_files = [x for x in list_of_files if 'assembly_data_report.jsonl' in x]
		if not jsonl_files:
			print('Unable to find an assembly_data_report.jsonl file in the zip file.')
			metatable_present = "N"
		elif len(jsonl_files) == 1:
			myzip.getinfo(jsonl_files[0]).filename = 'assembly_data_report.jsonl'
			myzip.extract(jsonl_files[0], path='data/input/')
			make_assembly_table_from_json()
			metatable_present = "Y"
		else:
			print('There appears to be more than one assembly_data_report.jsonl file in this zip file:')
			print('\n'.join(jsonl_files))
		file_types = {
			'cds_from_genomic.fna': ('./data/input/cds/', '.fna'),
			'genomic.fna': ('./data/input/genomic_fasta/', '_genomic.fna'),
			'genomic.gbff': ('./data/input/gbff/', '.gbff'),
			'genomic.gtf': ('./data/input/gtf/', '.gtf'),
			'genomic.gff': ('./data/input/gff/', '.gff'),
			'protein.faa': ('./data/input/protein_faa/', '.faa')
		}
		type_of_files_list = []
		for file_type, (output_dir, suffix) in file_types.items():
			files = [x for x in list_of_files if file_type in x]
			if files:
				type_of_files_list.append(file_type)
				directory_setter(output_dir)
				for file_name in files:
					if file_type == 'genomic.fna' and file_name.endswith("cds_from_genomic.fna"):
						continue
					zipextractor(myzip, file_name, output_dir, suffix)
	erase_if_empty('data/input/genomic_fasta')
	return(type_of_files_list, metatable_present)


def make_assembly_table_from_json():
	"""this function takes the json file from NCBI and makes a table used to make the strainlist"""
	with open('data/input/ncbi_strain_assembly_table.tab', 'w') as ncbi_strain_assembly_table:
		ncbi_strain_assembly_table.write("Accession\tSpecies\tStrain\tBioProject\tBioSample\tLevel\n")
	os.system(
		"dataformat tsv genome --inputfile data/input/assembly_data_report.jsonl "
			"--fields accession,ani-submitted-species,assminfo-biosample-strain,"
			"assminfo-bioproject,assminfo-biosample-accession,assminfo-level "
			"--elide-header >> data/input/ncbi_strain_assembly_table.tab"
	)


def line_format(sequence: str) -> str:
	"""Returns a FASTA sequence with 80 characters per line."""
	return "\n".join([sequence[i:i+80] for i in range(0, len(sequence), 80)])


def erase_if_empty(directory: str) -> None:
	"""
	Checks if the given directory is empty and deletes it if so.
	:param directory: Path to the directory to check and delete.
	"""
	if os.path.exists(directory) and os.path.isdir(directory):
		if not os.listdir(directory):  # Check if the directory is empty
			os.rmdir(directory)  # Remove the empty directory
			#print(f"Directory '{directory}' was empty and has been deleted.")
		else:
			pass #print(f"Directory '{directory}' is not empty.")
	else:
		pass #print(f"Directory '{directory}' does not exist.")


def get_locus_prefix(locus_id):
	"""Extracts the locus prefix from a locus ID."""
	if locus_id is None:
		return
	return str(locus_id.split("_")[0]) if "_" in locus_id else re.match(r"\D+", locus_id).group(0)


def convert_biopython_location_format_to_NCBI(input_string):
	"""Converts Biopython 'location' format to the format found in NCBI RefSeq FASTA headers."""
	join_pattern = r"join\{((?:\[\d+:\d+\]\([\+\-]\),?\s*)+)\}"
	join_match = re.match(join_pattern, input_string)
	
	if join_match:
		segments = join_match.group(1)
		segment_pattern = r"\[(\d+):(\d+)\]\(([\+\-])\)"
		matches = re.findall(segment_pattern, segments)
		
		if not matches:
			raise ValueError("No valid segments found inside join{}.")
		
		formatted_segments = []
		use_complement = False
		
		for start, end, symbol in matches:
			start_digit = int(start) + 1
			new_start = str(start_digit)
			
			if symbol == "+":
				formatted_segments.append(f"{new_start}..{end}")
			elif symbol == "-":
				formatted_segments.append(f"{new_start}..{end}")
				use_complement = True
		
		formatted_string = f"join({','.join(formatted_segments)})"
		
		if use_complement:
			return f"complement({formatted_string})"
		return formatted_string
	
	else:
		pattern = r"\[([<>]?\d+):([<>]?\d+)\]\(([\+\-])\)"
		match = re.match(pattern, input_string)
	
		if match:
			start, end, symbol = match.groups()
			start_symbol = start[0] if start[0] in '<>' else ''
			start_digit = int(start[1:] if start_symbol else start) + 1
			new_start = f"{start_symbol}{start_digit}"
			
			if symbol == '+':
				return f"{new_start}..{end}"
			elif symbol == '-':
				return f"complement({new_start}..{end})"
		else:
			raise ValueError(f"Input string {input_string} does not match the expected format.")


def convert_gbff_cds_to_fasta(seq_feature, seq_record, cds_index, locus_number):
	"""Converts CDS features to FASTA format."""
	cds_index += 1
	seq_feature_sequence = line_format(str(seq_feature.extract(seq_record).seq))
	translation = seq_feature.qualifiers.get('translation', [None])[0]
	if translation:
		seq_feature_protein_sequence = line_format(translation)
	else:
		seq_feature_protein_sequence = line_format(str(seq_feature.extract(seq_record).seq.translate()[:-1]))
	gene_name = seq_feature.qualifiers.get('gene', [None])[0]
	gene_name_tag = f"[gene={gene_name}] " if gene_name else ""
	codon_frame = seq_feature.qualifiers.get('codon_start', [None])[0]
	if codon_frame == "2":
		frame_tag = "[frame=2] "
	elif codon_frame == "3":
		frame_tag = "[frame=3] "
	else:
		frame_tag = ""
	gc_content = f"{gc_fraction(seq_feature.extract(seq_record).seq) * 100:.2f}"
	location = convert_biopython_location_format_to_NCBI(str(seq_feature.location))
	locus_tag = seq_feature.qualifiers.get("locus_tag", [None])[0]
	if not locus_tag:
		title = seq_record.annotations["accessions"][0]
		locus_tag = f"{title}_{str(locus_number).zfill(5)}"
		seq_feature.qualifiers["locus_tag"] = [locus_tag]
		locus_number += 5
	protein_id = seq_feature.qualifiers.get("protein_id", [None])[0]
	protein_description_tag = f"[protein_id={protein_id}] " if protein_id else ""
	product = seq_feature.qualifiers.get("product", [None])[0]
	if "pseudo" in seq_feature.qualifiers:
		protein_id = locus_tag
		protein_description_tag = "[pseudo=true] "

	with open(f"data/fna/{locus_tag.split('_')[0]}.fna", "a") as output_fna_handle:
		output_fna_handle.write(
			f">lcl|{seq_record.id}_cds_{protein_id}_{cds_index} {gene_name_tag}[locus_tag={locus_tag}] [protein={product}] "
			f"{frame_tag}{protein_description_tag}[location={location}] [pctGC={gc_content}] [gbkey=CDS]\n{seq_feature_sequence}\n"
		)
	with open(f"data/faa/{locus_tag.split('_')[0]}.faa", "a") as output_faa_handle:
		output_faa_handle.write(
			f">{seq_record.id}_cds_{protein_id}_{cds_index} {gene_name_tag}[locus_tag={locus_tag}] [protein={product}] "
			f"{frame_tag}{protein_description_tag} [location={location}] [pctGC={gc_content}] [gbkey=CDS]\n{seq_feature_protein_sequence}\n"
		)
	return cds_index, locus_number


def convert_rna_to_fasta(seq_feature, seq_record, rna_index):
#Converts RNA features to FASTA format.
	rna_index += 1
	seq_feature_sequence = line_format(str(seq_feature.extract(seq_record).seq))
	locus_tag = seq_feature.qualifiers.get("locus_tag", [None])[0]
	product = seq_feature.qualifiers.get("product", [None])[0]
	location = convert_biopython_location_format_to_NCBI(str(seq_feature.location))
	gc_content = f"{gc_fraction(seq_feature.extract(seq_record).seq) * 100:.2f}"
	
	with open(f"data/rna/{locus_tag.split('_')[0]}_rna.fna", "a") as output_rna_handle:
		output_rna_handle.write(
			f">lcl|{seq_record.id}_{seq_feature.type.lower()}_{rna_index} "
			f"[locus_tag={locus_tag}] [product={product}] [location={location}] "
			f"[pctGC={gc_content}] [gbkey={seq_feature.type}]\n{seq_feature_sequence}\n"
		)
	if "16S" in product and "<" not in location and ">" not in location:
		with open("data/rna/16SrRNA.fna", "a") as output_rrna_handle:
			output_rrna_handle.write(
				f">{locus_tag}:{seq_record.id}_{seq_feature.type.lower()}_{rna_index} "
				f"[locus_tag={locus_tag}] [product={product}] [location={location}] "
				f"[pctGC={gc_content}] [gbkey={seq_feature.type}]\n{seq_feature_sequence}\n"
			)
	return rna_index


def extract_CRISPR_info(seq_feature, seq_record):
#Extracts CRISPR information from a sequence feature and writes it out to 
#a file called "CRISPR_regions.txt"
	description = seq_record.description.split(",")[0].split(' NODE')[0]
	rpt_family = seq_feature.qualifiers.get('rpt_family', ['Unknown'])[0]
	array_seq = seq_feature.extract(seq_record).seq
	dir_repeat = seq_feature.qualifiers.get('rpt_unit_seq', [''])[0].upper()
	spacers = "\n".join(part for part in str(array_seq).split(dir_repeat) if part.strip()) if dir_repeat else "No direct repeat sequence found."
	
	with open("data/report/CRISPR_regions.txt", "a") as output_CRISPR_handle:
		output_CRISPR_handle.write(
			f">{description}\n{seq_record.id}; location: {seq_feature.location}; {rpt_family}\n"
			f"full repeat region: {array_seq}\ndirect repeat: {dir_repeat}\nSpacers:\n{spacers}\n\n"
		)


def convert_genbank(filename: str) -> None:
#Converts GenBank files to a set of FASTA (nucleotide cds, RNA, and protein).
	with open(filename, "r") as input_handle:
		cds_index = 0
		rna_index = 0
		locus_number = 5
		for seq_record in SeqIO.parse(input_handle, "genbank"):
			for seq_feature in seq_record.features:
				if seq_feature.type == "CDS":
					cds_index, locus_number = convert_gbff_cds_to_fasta(seq_feature, seq_record, cds_index, locus_number)
				elif seq_feature.type in ["tRNA", "rRNA", "ncRNA", "tmRNA"]:
					rna_index = convert_rna_to_fasta(seq_feature, seq_record, rna_index)
				elif seq_feature.type == "repeat_region" and "CRISPR" in seq_feature.qualifiers.get('rpt_family', []):
					extract_CRISPR_info(seq_feature, seq_record)


def make_strainlist_dictionary():
#this will read the strain_assembly_table.tab file and return a dictionary where
#the key is the GCF_ or GCA_ assembly ID and the returned value is the 
#strain metatdata including the species, strain, and assembly, BioSample, and BioProject IDs
#as well as the assembly state (contig, scaffold, or complete)
	strain_lookup_dictionary = {}
	with open('data/input/ncbi_strain_assembly_table.tab', 'r') as assembly_table:
		for line in assembly_table.readlines()[1:]:
			splitline = line.strip().split('\t')
			assembly_id = splitline[0]
			assembly_index = assembly_id[4:].split(".")[0]
			strain_info_line = f"{splitline[1].strip()} {splitline[2]} [{assembly_id}; {splitline[4]}; {splitline[3]}; {'Complete' if splitline[5] == 'Complete Genome' else splitline[5]}]"
			strain_lookup_dictionary[assembly_index] = strain_info_line
	return strain_lookup_dictionary


def extract_assemblyid_from_filename(filename):
#This function takes the filename of an NCBI file that contains the
#NCBI assembly ID (GCF or GCA number) and returns the assembly id.
	# Define the regex pattern
	pattern = r'.*/GC[FA]_(\d+)\.\d+\.(gbff|fna)$'
	
	# Search for the pattern in the filename
	match = re.search(pattern, filename)
	
	if match:
		# Extract and return the number
		return match.group(1)
	else:
		return None


def extract_assemblyid_from_gbff_metadata(ncbi_gbff_filename):
#if you can't get the assembly_id from the filename you can try to get it from
#metadata in the gbff file (the dbxrefs feature)
	with open(ncbi_gbff_filename, "r") as input_handle:
		for seq_record in SeqIO.parse(input_handle, "genbank"):
				for ref in seq_record.dbxrefs:
					if ref[:9] == "Assembly:":
						return ref[9:].split('.')[0].split("_")[1]
	return None


def extract_locusprefix_from_gbff(filename):
#this combs through a gbff file and finds a locus prefix.
#we assume that every gbff file is represented by a single locus prefix
#it exits on the first one it finds.  So if the gbff file has multiple 
#locus prefixes (e.g. multiple gbk files concatenated into a single file)
#this function should not be used.
	with open(filename, "r") as input_handle:
		for seq_record in SeqIO.parse(input_handle, "genbank"):
			for seq_feature in seq_record.features:
				locus_tag = seq_feature.qualifiers.get("locus_tag", [None])[0]
				if locus_tag:
					locus_prefix = get_locus_prefix(locus_tag)
					if locus_prefix:
						return locus_prefix
	return None


def extract_locusprefix_from_cds(filename):
#this combs through a cds fasta file and finds a locus profix.
#we assume that every file is represented by a single locus prefix
#it exits on the first one it finds.  So if the file has multiple 
#locus prefixes (e.g. multiple fasta files from different strains or 
#assemblies that are concatenated into a single file)
#this function should not be used.
	pattern = re.compile(r'\[locus_tag=(.*?)\]')
	with open(filename, 'r') as file:
		for line in file:
			if line.startswith('>'):
				match = pattern.search(line)
				if match:
					locus_prefix = get_locus_prefix(match.group(1))
					return locus_prefix
	return None
	
	
def move_unzipped_data(filelist):
#This function moves NCBI files to their correct location if they are in an 
#unzipped folder instead of in a zip archive.
	json_present = "N"
	metatable_present = "N"
	strainlist_present = "N"
	gbff_type = "N"
	fna_type = "N"
	data_types = []
	
	os.makedirs('data/input/gbff', exist_ok=True)
	os.makedirs('data/input/cds', exist_ok=True)

	for file_path in filelist:
		if file_path.endswith('.gbff'):
			if os.path.basename(file_path) != 'genomic.gbff':
				shutil.copy(file_path, os.path.join('data/input/gbff', os.path.basename(file_path)))
				gbff_type = "Y"
			else:
				fileaccession = os.path.dirname(file_path).split('/')[-1]
				shutil.copy(file_path, f'data/input/gbff/{fileaccession}.gbff')
				gbff_type = "Y"

		elif file_path.endswith('.fna'):
			if os.path.basename(file_path) != 'cds_from_genomic.fna':
				shutil.copy(file_path, os.path.join('data/input/cds', os.path.basename(file_path)))
				fna_type = "Y"
			else:
				fileaccession = os.path.dirname(file_path).split('/')[-1]
				shutil.copy(file_path, f'data/input/cds/{fileaccession}.fna')
				fna_type = "Y"

		elif os.path.basename(file_path) == 'assembly_data_report.jsonl':
			shutil.copy(file_path, os.path.join('data/input', os.path.basename(file_path)))
			json_present = "Y"

		elif os.path.basename(file_path) == 'ncbi_strain_assembly_table.tab':
			shutil.copy(file_path, os.path.join('data/input', os.path.basename(file_path)))
			metatable_present = "Y"

		elif os.path.basename(file_path) == 'strainlist.txt':
			shutil.copy(file_path, os.path.join('./', os.path.basename(file_path)))
			strainlist_present = "Y"

		if gbff_type == "Y":
			data_types.append("genomic.gbff")
		if fna_type == "Y":
			data_types.append("cds_from_genomic.fna")
	
	erase_if_empty('data/input/gbff')
	erase_if_empty('data/input/cds')		
	return data_types, json_present, metatable_present, strainlist_present


def convert_ncbi_cds(filepath, locus_prefix):
	fna_output_filename = f"data/fna/{locus_prefix}.fna"
	faa_output_filename = f"data/faa/{locus_prefix}.faa"
	with open(fna_output_filename, "w") as output_handle_fna, open(faa_output_filename, "w") as output_handle_faa:
		for record in SeqIO.parse(filepath, "fasta"):
				GC_content = float(100 * gc_fraction(record.seq))
				record.description = f"{record.description.split('[gbkey')[0]}[pctGC={GC_content:.4}] [gbkey{record.description.split('[gbkey')[1]}"
				output_handle_fna.write(f">{record.description}\n{line_format(str(record.seq))}\n")
				output_handle_faa.write(f">{record.description[4:]}\n{line_format(str(record.seq.translate()[:-1]))}\n")


def extract_strainlist_info_from_gbff (locustag, filename):
	with open(filename, "r") as input_handle:
		for seq_record in SeqIO.parse(input_handle, "genbank"):
			description = seq_record.description.split(',')[0]
			accession = ", ".join(seq_record.dbxrefs)
			strainlist_entry = f"{locustag} | [{accession}] | {description}"
			return strainlist_entry
	return None


def clear_directory(directory: str) -> None:
	"""
	Erases all contents of a specified directory without removing the directory itself.

	Parameters:
		directory (str): The path to the directory to clear.

	Returns:
		None
	"""
	if not os.path.isdir(directory):
		raise ValueError(f"The specified path '{directory}' is not a directory or does not exist.")

	for item in os.listdir(directory):
		item_path = os.path.join(directory, item)
		try:
			if os.path.isfile(item_path) or os.path.islink(item_path):
				os.unlink(item_path)  # Remove files or symlinks
			elif os.path.isdir(item_path):
				shutil.rmtree(item_path)  # Remove directories
		except Exception as e:
			print(f"Failed to delete '{item_path}'. Reason: {e}")



if __name__ == "__main__":

	import argparse

	# set up defaults for the variables

	strainlist = "strainlist.txt"
	input_zipfile = "./ncbi_dataset.zip"
	input_sequence_datatype = "gbff"
	zipped_data = "Y"
	input_directory = ""
	json_present = "N"
	metatable_present = "N"
	strainlist_present = "N"

#parse the command line for user input
	instructions_for_argparse = "This script processes NCBI data. Provide paths to zipped or unzipped data."

# Initialize ArgumentParser
	parser = argparse.ArgumentParser(add_help=True, description=instructions_for_argparse)

# Argument for zipped file
	parser.add_argument('-z', '--zip', dest='input_zipfile', help="Name and path of the zipped file with search terms (default = 'ncbi_dataset.zip').")

# Argument for unzipped data directory
	parser.add_argument('-i', '--input', dest='unzipped_input_directory', help="Path to the directory containing unzipped NCBI data.")

# Verbose flag
	parser.add_argument('-v', '--verbose', action='store_true', help="Turn on Biopython warnings.")

# Parse arguments
	args = parser.parse_args()

	if args.input_zipfile:
		input_zipfile = args.input_zipfile

	if args.unzipped_input_directory:
		zipped_data = 'N'
		input_directory = args.unzipped_input_directory

	"""if not verbose (default) then do warnings.simplefilter('ignore', BiopythonWarning)"""
	if args.verbose:
		print("BioPython warnings will be displayed.\n")
		verbose_setting = 'Y'
	else:
		warnings.simplefilter('ignore', BiopythonWarning)
		verbose_setting = 'N'
		
	#create the necessary directories in the project directory
	directory_setter('data')
	directory_setter('data/input')
	directory_setter('data/fna')
	directory_setter('data/faa')
	directory_setter('data/tables')
	directory_setter('data/report')
	

###Look for json file (assembly_data_report.jsonl)
###or assembly table (ncbi_strain_assembly_table.tab)
###or completed strainlist (strainlist.txt)

	#if not zipped then move files into correct directories.
	if zipped_data == 'N':
		directory_list = tree2list(input_directory)
		data_types, json_present, metatable_present, strainlist_present = move_unzipped_data(directory_list)

	#if zipped then unzip into correct directories and get back what types of data files you have.
	if zipped_data == 'Y':
		data_types, metatable_present = ncbi_unzip(input_zipfile)

#get both the assembly ID and locus prefix for each file
	if 'genomic.gbff' in data_types:
		ncbi_gbff_file_list = glob.glob('data/input/gbff/*.gbff')
		directory_setter('data/rna')
		clear_directory('data/rna')
		clear_directory('data/fna')
		clear_directory('data/faa')
		query_list = []
		for ncbi_gbff_filename in ncbi_gbff_file_list:
			search_terms = []
			assembly_id = extract_assemblyid_from_filename(ncbi_gbff_filename)
			
			#if the filename is not the assembly ID...
			if not assembly_id:
				assembly_id = extract_assemblyid_from_gbff_metadata(ncbi_gbff_filename)
				
			locus_prefix = extract_locusprefix_from_gbff(ncbi_gbff_filename)
			if not locus_prefix:
				locus_prefix = "UNKNOWN"
			search_terms = [assembly_id, locus_prefix, ncbi_gbff_filename]
			query_list.append(search_terms)
			
			#Now write out the faa, fna and rna files from the genbank file
			convert_genbank(ncbi_gbff_filename)

	elif 'cds_from_genomic.fna' in data_types:
		ncbi_cds_file_list = glob.glob('data/input/cds/*.fna')
		query_list = []
		for ncbi_cds_filename in ncbi_cds_file_list:
			search_terms = []
			assembly_id = extract_assemblyid_from_filename(ncbi_cds_filename)
			if not assembly_id:
				print(f"Can't retrieve an an assembly ID for {ncbi_cds_filename}")
				assembly_id = "None"
			locus_prefix = extract_locusprefix_from_cds(ncbi_cds_filename)
			if not locus_prefix:
				locus_prefix = "None"
			search_terms = [assembly_id, locus_prefix, ncbi_cds_filename]
			query_list.append(search_terms)
		
		#Now write out the faa, fna and rna files from the cds file by sending
		#the locus prefix and filepath of the cds file to the convert function
		#the locus prefix is used to rename the fna and faa files
		for item in query_list:
			convert_ncbi_cds(item[2], item[1])


	#now create a strainlist if one was not supplied by the user
	if strainlist_present == "N":
		strainlist_dictionary = {}

		#if theres a metadata table present make a dictionary from it
		if metatable_present == "Y":
			strainlist_dictionary = make_strainlist_dictionary()
		#if there's an NCBI JSON metafile then first make a metadata table then
		#convert that into a dictionary variable 
		elif json_present == "Y":
			make_assembly_table_from_json()
			strainlist_dictionary = make_strainlist_dictionary()

		#if there's a strainlist dictionary made from the previous steps then 
		#query the dictionary with the assembly ID
		if strainlist_dictionary:
			with open('strainlist.txt', 'w') as strainlist_output_handle:
				for item in query_list:
					strainlist_output_handle.write(f"{item[1]} | {strainlist_dictionary.get(item[0], 'Unknown')}\n")

		#if there's no user-provided strainlist.txt AND there's no metadata available
		#in either tab or json format then build a crude one directly from info
		#in the genbank files or the cds files.
		elif 'genomic.gbff' in data_types:
			with open('strainlist.txt', 'w') as strainlist_output_handle:
				for item in query_list:
					strainlist_entry = extract_strainlist_info_from_gbff (item[1], item[2])
					strainlist_output_handle.write(f"{strainlist_entry}\n")

		elif 'cds_from_genomic.fna' in data_types:
			with open('strainlist.txt', 'w') as strainlist_output_handle:
				for item in query_list:
					strainlist_output_handle.write(f"{item[1]} | {item[2].split('/')[-1]}\n")
					
		else:
			print('Unable to construct a strainlist.txt file!')
