#!/usr/bin/env python

from zipfile import ZipFile
import os
import shutil
import pathlib
from typing import Iterator, Tuple

# need to move this "assembly_data_report.jsonl"

# need to make a function that makes the right directories as needed.
def directory_setter(directory_name, directory_list2):
	if directory_name[0:-1] not in directory_list2:
		os.makedirs(directory_name)
	return()

def tree2list(directory: str) -> Iterator[Tuple[str, str, str]]:
	for i in os.scandir(directory):
		if i.is_dir():
			yield i.path
			yield from tree2list(i.path)
		else:
			yield i.path

def zipextractor (input_file_name, output_directory_name, tail):
	file_directory = os.path.dirname(input_file_name)
	fileaccession = file_directory.split('/')[-1]
	#rename the zip file in the zip'object' before extracting it.
	myzip.getinfo(input_file_name).filename = (fileaccession + tail) 
	myzip.extract(input_file_name, path=output_directory_name)
	
if __name__ == '__main__':
	
	#get current directories and put into list
	directory_list = list(tree2list('.')) 

	with ZipFile('ncbi_dataset.zip') as myzip:
		#get list of files within the zip folder
		list_of_files_inside_zip_folder = myzip.namelist()
		
		jsonl_files = [x for x in list_of_files_inside_zip_folder if 'assembly_data_report.jsonl' in x]
		if len(jsonl_files) == 0:
			print('Unable to find an assembly_data_report.jsonl file in the zip file.')
		if len(jsonl_files) == 1:
			myzip.getinfo(jsonl_files[0]).filename = 'assembly_data_report.jsonl' 
			myzip.extract(jsonl_files[0], path='data/ncbi/')
		elif len(jsonl_files) >= 2:
			print('There appears to be more than one assembly_data_report.jsonl file in this zip file:')
			print('\n'.join(jsonl_files))
			
		#now iterate through the ZipFile list and extract the filenames you want
		allfna_files_from_zip = [x for x in list_of_files_inside_zip_folder if 'genomic.fna' in x]

		genomic_files_from_zip = [x for x in allfna_files_from_zip if '/cds_from_genomic.fna' not in x]
		if len(genomic_files_from_zip) != 0:
			output_directory_name = './data/ncbi/genomic_fasta/'
			file_suffix = '_genomic.fna'
			directory_setter(output_directory_name, directory_list)
			for file_name in genomic_files_from_zip:
				zipextractor(file_name, output_directory_name, file_suffix)

		cdsfna_files_from_zip = [x for x in allfna_files_from_zip if '/cds_from_genomic.fna' in x]
		if len(cdsfna_files_from_zip) != 0:
			output_directory_name = './data/ncbi/cds_fasta/'
			file_suffix = '.fna'
			directory_setter(output_directory_name, directory_list)
			for file_name in cdsfna_files_from_zip:
				zipextractor(file_name, output_directory_name, file_suffix)

		gbff_files_from_zip = [x for x in list_of_files_inside_zip_folder if 'genomic.gbff' in x]
		if len(gbff_files_from_zip) != 0:
			output_directory_name = './data/ncbi/gbff/'
			file_suffix = '.gbff'
			directory_setter(output_directory_name, directory_list)
			for file_name in gbff_files_from_zip:
				zipextractor(file_name, output_directory_name, file_suffix)

		gtf_files_from_zip = [x for x in list_of_files_inside_zip_folder if 'genomic.gtf' in x]
		if len(gtf_files_from_zip) != 0:
			output_directory_name = './data/ncbi/gtf/'
			file_suffix = '.gtf'
			directory_setter(output_directory_name, directory_list)
			for file_name in gtf_files_from_zip:
				zipextractor(file_name, output_directory_name, file_suffix)
				
		gff_files_from_zip = [x for x in list_of_files_inside_zip_folder if 'genomic.gff' in x]
		if len(gff_files_from_zip) != 0:
			output_directory_name = './data/ncbi/gff/'
			file_suffix = '.gff'
			directory_setter(output_directory_name, directory_list)
			for file_name in gff_files_from_zip:
				zipextractor(file_name, output_directory_name, file_suffix)
								
		protein_files_from_zip = [x for x in list_of_files_inside_zip_folder if 'protein.faa' in x]
		if len(protein_files_from_zip) != 0:
			output_directory_name = './data/ncbi/protein_faa/'
			file_suffix = '_protein.faa'
			directory_setter(output_directory_name, directory_list)
			for file_name in protein_files_from_zip:
				zipextractor(file_name, output_directory_name, file_suffix)
				
	