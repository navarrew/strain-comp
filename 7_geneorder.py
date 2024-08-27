#!/usr/bin/env python

import os
import argparse
import pandas as pd

def get_locations(cell_item):
	first_cell_fragment = cell_item.split("]|")[1].rstrip()
	cell_fragments = first_cell_fragment.split(", ")
	position_list = []
	for fragment in cell_fragments:
		contig_info = fragment.split("|")[-1]
		positional_info = contig_info.split("_")
		positional_info.insert(0, positional_info[-1].zfill(5))
		position_list.append("_".join(positional_info).replace("_", ":", 1))
	position_list.sort()
	position_string = str(", ".join(map(str, position_list)))
	return position_string
	
def tree2list(directory: str) -> Iterator[Tuple[str, str, str]]:
	import os
    for i in os.scandir(directory):
        if i.is_dir():
            yield i.path[2:]
            yield from tree2list(i.path)
        else:
            yield i.path[2:]
            
if __name__ == '__main__':

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

	input_table = pd.read_csv(input_filepath, sep='\t', index_col='CLUSTER')

	end_of_annotations_index = input_table.columns.get_loc('strain count')
	number_of_columns = input_table.shape[1]+1 #had to add one because it doesn't count the index column I guess...or its indexed by 0.
	annotations_table = input_table.iloc[:,0:end_of_annotations_index+1]
	strain_data_table = input_table.iloc[:,(end_of_annotations_index+1):]	

	strain_data_table.to_csv('temp1.tab', sep='\t')

	with open('temp1.tab', 'r') as f:
		output_table = []
		for line in [line.rstrip() for line in f.readlines()]:
			line_list = []
			line_items = line.split('\t')
			line_list.append(line_items[0])
			for cell_item in line_items[1:]:
				if cell_item[0] == "[":
					processed_cell = get_locations(cell_item)
					line_list.append(processed_cell)
				else:
					line_list.append(cell_item)
			new_output_line = "\t".join(line_list)
			output_table.append(new_output_line)	
	f.close()

	header = output_table.pop(0).split('\t')
	strain_positions_table = pd.DataFrame([x.split('\t') for x in output_table], columns=header).set_index('CLUSTER')

	final_output_table = pd.concat([annotations_table, strain_positions_table], axis=1).reindex(annotations_table.index)
	final_output_table.reset_index(drop=False, inplace=True)
	final_output_table.to_csv('tab/cluster_table_geneordered.tab', sep='\t')

	del annotations_table
	del strain_data_table
	del strain_positions_table
	os.system('rm temp1.tab')
	
