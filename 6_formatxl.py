#!/usr/bin/env python

import os
import argparse
import pandas as pd
import xlsxwriter
from datetime import datetime, timedelta
from typing import Iterator, Tuple

def tree2list(directory: str) -> Iterator[Tuple[str, str, str]]:
	for i in os.scandir(directory):
		if i.is_dir():
			yield i.path
			yield from tree2list(i.path)
		else:
			yield i.path


if __name__ == '__main__':

#make a list of all the files in the 'tables' folder
	directory_list = list(tree2list('tables'))
	filelist = []
	for filename in directory_list:
		if 'tables/cluster_table' in filename:
			if '.tab' in filename:
				filelist.append(filename)
	filelist.sort()
	for i, x in enumerate(filelist):
		print(str(i+1) + '. ' + x)
	
	input_filepath_index = int(input('Which file to convert? '))
	input_filepath = filelist[input_filepath_index - 1]
	input_filename = input_filepath.split('/')[-1][0:-4]


	#read the file into a pandas dataframe
	timestart = datetime.now()
	print("Putting table into a pandas dataframe.  Started at: " + timestart.strftime("%Y-%m-%d %H:%M:%S"))
	input_table = pd.read_csv(input_filepath, sep='\t')

	timestart = datetime.now()
	print("Putting pandas dataframes into an excelwriter object.  Started at: " + timestart.strftime("%Y-%m-%d %H:%M:%S"))

	# Create a Pandas 'ExcelWriter' object using XlsxWriter as the engine.
	writer = pd.ExcelWriter(r"output/"+ input_filename + ".xlsx", engine="xlsxwriter")

	# Convert the pandas dataframe to sheets in the ExcelWriter object.
	input_table.to_excel(writer, sheet_name="cluster_table")

	workbook = writer.book
	worksheet1 = writer.sheets["cluster_table"]

# Define the header formats - start with a default set that will be applied to all.  Then modify for color.
	
	default_parameters = {"bold": True, "text_wrap": True, "valign": "bottom", "border": 1}
	
	#for the cluster header
	header_format_cluster = workbook.add_format(default_parameters)
	header_format_cluster.set_font_color('white')
	header_format_cluster.set_fg_color('black')

	#for the info like GC content and protein length
	header_format_lead = workbook.add_format(default_parameters)
	header_format_lead.set_fg_color('#5DEFF0')

	#for the info extracted from NCBI	
	header_format_ncbi = workbook.add_format(default_parameters)
	header_format_ncbi.set_fg_color('#FCD5B4')

	#for the info extracted from COG	
	header_format_cog = workbook.add_format(default_parameters)
	header_format_cog.set_fg_color('#FFBF00')


	#for the strain and hit count column headers
	header_format_counts = workbook.add_format(default_parameters)
	header_format_counts.set_fg_color('#00F900')

	#strain headers will be rotated so the locus prefix is visible.
	rotated_parameters = {"bold": True, "text_wrap": True, "valign": "bottom", "border": 1, "rotation": 90} 
	#for complete genomes use a bright green.  Get less bright with each lower assembly step. 
	header_format_complete = workbook.add_format(rotated_parameters)
	header_format_complete.set_fg_color('#00F900')
	header_format_chromosome = workbook.add_format(rotated_parameters)
	header_format_chromosome.set_fg_color('#00F900')	
	header_format_scaffold = workbook.add_format(rotated_parameters)
	header_format_scaffold.set_fg_color('#D4FB79')
	header_format_contig = workbook.add_format(rotated_parameters)
	header_format_contig.set_fg_color("#E3FBD6")

	timestart = datetime.now()
	print("Formatting columns in object.  Started at: " + timestart.strftime("%Y-%m-%d %H:%M:%S"))
	
	cluster_name_example = input_table.iloc[1, input_table.columns.get_loc('CLUSTER')]
	cluster_name_length = len(cluster_name_example)
	
	# Write the column headers with the defined format.
	for col_num, value in enumerate(input_table.columns.values):
		if " | " in value:
			worksheet1.set_column(col_num + 1, col_num + 1, 2)	
			if "; Complete]" in value:
				worksheet1.write(0, col_num + 1, value, header_format_complete)
			elif "; Chromosome]" in value:
				worksheet1.write(0, col_num + 1, value, header_format_chromosome)
			elif "; Scaffold]" in value:
				worksheet1.write(0, col_num + 1, value, header_format_scaffold)
			elif "; Contig]" in value:
				worksheet1.write(0, col_num + 1, value, header_format_contig)
			
		elif value[0:4] == "NCBI":
			worksheet1.write(0, col_num + 1, value, header_format_ncbi)
			if value == "NCBI annotations":
				worksheet1.set_column(col_num + 1, col_num + 1, 45)	

		elif value[0:4] == "COG_":
			worksheet1.write(0, col_num + 1, value, header_format_cog)
			if value == "NCBI annotations":
				worksheet1.set_column(col_num + 1, col_num + 1, 45)	
				
		elif value == 'total count' or value == 'strain count':
			worksheet1.write(0, col_num + 1, value, header_format_counts)
			worksheet1.set_column(col_num + 1, col_num + 1, 5)	
			
		else:
			worksheet1.write(0, col_num + 1, value, header_format_lead)
			worksheet1.set_column(col_num + 1, col_num + 1, 8)	

	worksheet1.write(0, 0, "original order", header_format_cluster)
	worksheet1.write(0, 1, "CLUSTER", header_format_cluster)
	worksheet1.set_column(0, 0, 6)
	worksheet1.set_column(1, 1, cluster_name_length+1)
	worksheet1.freeze_panes(1, 2)
	
	
	timestart = datetime.now()
	print("Writing out to excel file.  Started at: " + timestart.strftime("%Y-%m-%d %H:%M:%S"))

	# Close the Pandas Excel writer and output the Excel file.
	writer.close()
	
	timestart = datetime.now()
	print("Ended at: " + timestart.strftime("%Y-%m-%d %H:%M:%S"))
