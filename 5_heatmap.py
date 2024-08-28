#!/usr/bin/env python

from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import sys
import os
import argparse
from typing import Iterator, Tuple

def tree2list(directory: str) -> Iterator[Tuple[str, str, str]]:
	import os
	for i in os.scandir(directory):
		if i.is_dir():
			yield i.path
			yield from tree2list(i.path)
		else:
			yield i.path

if __name__ == '__main__':

#sometimes a heatmap needs a higher number of recursive calls to function.	We're setting a high limit here.
	sys.setrecursionlimit(100000) 

#Setting us up with directories we'll need if they don't exist already.
	directory_list = list(tree2list('.'))

	if "./tables" not in directory_list:
		print("\nCannot find a 'tables' directory in the current directory.")
		print("Are you sure you're in the right directory?\n")
		exit()
	if "./strainlist.txt" not in directory_list:
		print("\nCannot find strainlist.txt in the current directory.")
		print("Please make a copy of it and put it in the main project directory.\n")
		exit()
	if "./tables/cluster_hit_count_table.tab" not in directory_list:
		print("\nCannot find the 'tables/cluster_hit_count_table.tab file.")
		print("Are you sure you're in the right directory?\n")
		exit()
	#if not already in the tab directory make the directory 'tables/unsorted'
	if "./tables/archive" not in directory_list:
		os.mkdir("tables/archive")

#if not already in the current directory make the directories 'report' and 'heatmap'
	if "./data/report" not in directory_list:
		os.mkdir("data/report")
	if "./data/heatmap" not in  directory_list:
		os.mkdir("data/heatmap")
	if "./output" not in  directory_list:
		os.mkdir("output")

	instructions = "The heatmap.py script allows you to cluster strains together based on " \
				   "gene/presence absence (August 5, 2024)"

	'''
	┌----------------------------------------------------------------------------------┐
	| 1. Getting user input for hierarchical clustering and generation of the heatmap. |
	└----------------------------------------------------------------------------------┘
	'''
# STEP 1 = the following argparse fxn reads the command line for specific inputs
	parser = argparse.ArgumentParser(add_help=True, description=instructions)
	parser.add_argument('-m', '--method', action='store', dest='method', help = "Clustering method (default = average).")
	parser.add_argument('-c', '--color', action='store', dest='cmap', help = "Plot color (default = Blues).")
	parser.add_argument('-hm', '-hitmax', action='store', type=int, dest='vmax', help = "Max hits expected per cluster (default = 1).")
	parser.add_argument('-l', '--label', action='store', dest='xticklabels', help = "Label the axes? (Y/N)")
	parser.add_argument('-o', '--order', action='store', dest='reorder_opt', help = "Reorder the strainlist and table after clustering? (Y/N)")
	parser.add_argument('-r', '--range', action='store', dest='strainrange', help = "Min and max strain hits to use in clustering analysis (format: min#:max#)")

	args = parser.parse_args()

# STEP 2 - if there is nothing in the command line - use default values

	if args.method:
		method = args.method
	else: method = "average"

	if args.cmap:
		cmap = args.cmap
	else: cmap = "Blues"

	if args.vmax:
		vmax = args.vmax
	else: vmax = 1

	if args.xticklabels:
		xticklabels = args.xticklabels
	else: xticklabels = True

	if args.reorder_opt:
		reorder_opt = args.reorder_opt
	else: reorder_opt = "Y"

	if args.strainrange != None:
		strainrange = args.strainrange
		if ":" in strainrange:
			min_strain_hits = int(strainrange.split(':')[0])
			max_strain_hits = int(strainrange.split(':')[1])
		else: 
			print('range is in wrong format - should be min:max')
			strainrange = "WARN THE USER"
	else: strainrange = "WARN THE USER"


	# Write the date and time of execution to the report file,
	# the hierarchical clustering method,
	# whether strain list or tables were reordered,
	# the colormap used,
	# the maximum numebr of expected hits for a cluster in any given strain.
	with open('data/report/report.txt', 'a') as f:
		timestart = datetime.now()
		# MAKING A TIME STAMP FOR ALL SAVEFILES
		timestamp = str(timestart.strftime("%y%m%d%H%M%S"))
		f.write("\n\n#######################################"
				"\nHEATMAP AND HIERARCHICAL CLUSTERING: " + str(timestart.strftime("%Y-%m-%d %H:%M:%S")) +
				"\nFigure and unsorted data files given the timestamp: " + timestamp +
				"\n\nHierarchical clustering method: " + method +
				"\nMatplotlib color map used: " + cmap +
				"\n'vmax' parameter set to: " + str(vmax))
		if reorder_opt == 'Y':
			f.write("\nStrain list and tables were reordered based on hierarchical clustering.")
		else:
			f.write("\nStrain list and tables were not reordered based on hierarchical clustering.")

	'''
	┌---------------------------------------------------------------------┐
	| 2. Processing the data that will go into the clustermap.            |
	└---------------------------------------------------------------------┘
	'''

# importing cluster_hit_cout_table. tab data into a pandas dataframe called original_table

	original_table = pd.read_csv('tables/cluster_hit_count_table.tab', sep='\t', index_col='CLUSTER')

# grabbing the indicies (row index)
	original_table_clusters_list = original_table.index.tolist()
	original_table_clusters_set = set(original_table_clusters_list)
	
# we determine the maximum number of strains being analyzed by counting the columns with strain info
	strainmax = original_table.shape[1] - 2

# the settings below are the defaults...all strains = max.  Minimum strains = 1.
	if strainrange =="WARN THE USER":
		max_strain_hits = strainmax
		min_strain_hits = 1

# trimming the table to exclude data that won't be used in clustering.
	table2 = original_table.loc[(original_table['strain count'] <= max_strain_hits) & (original_table['strain count'] >= min_strain_hits)]
	table2 = table2.assign(diff=table2['total count'] - table2['strain count'])
	table2 = table2.loc[table2['diff'] == 0]
#we trimmed the table (table2) to create our final trimmed table, simply called 'table'
	table = table2.iloc[:,2:-1] #remove the columns for diff, strain counts, and total counts

# write the newly trimmed table out for the user to look at if they need to.
	table.to_csv('data/heatmap/trimmed_cluster_count_table_for_heatmap.tab', sep='\t')

#lets figure out which proteins were removed from the original protein clusters
	trimmed_table_clusters_list = table.index.tolist()
	trimmed_table_clusters_set = set(trimmed_table_clusters_list)
	proteins_excluded_in_clustering_set = original_table_clusters_set - trimmed_table_clusters_set
	excluded_list = list(proteins_excluded_in_clustering_set)

	excluded_table = original_table.reindex(index=excluded_list)
	excluded_table = excluded_table.sort_values(by='strain count', ascending=False)
	excluded_table.to_csv('data/heatmap/clusters_excluded_from_analysis_'+timestamp+'.tab', sep='\t')
	excluded_list = excluded_table.index.tolist()
	
	number_of_protein_clusters = table.shape[0]

	print("You will cluster a total of " + str(number_of_protein_clusters) +  \
			  " protein clusters across "+ str(strainmax) +" strains." + \
			  "\nTotal features = " + str(number_of_protein_clusters * strainmax))

	if strainrange =="WARN THE USER":
		print("\nVery large tables can take a long time to cluster and make a figure for.")
		print("If you want to reduce time you can use the -r option to reduce the size of the dataset.")


#writing out the number of strains and window to the report file
	with open('data/report/report.txt', 'a') as f:
		f.write("\nTotal number of strains in heatmap: " + str(strainmax) +
				"\nMax 'strain hits' included in clustering: " + str(max_strain_hits) +
				"\nMin 'strain hits' included in clustering: " + str(min_strain_hits) +
				"\nNumber of protein clusters analyzed: " + str(number_of_protein_clusters))

	'''
	┌---------------------------------------------------------------------┐
	| 3. Executing hierarchical clustering and generation of the heatmap. |
	└---------------------------------------------------------------------┘
	'''

	print("Generating heatmap.")
#its nice to know how long clustering is taking.  This will tell us.
	timestart = datetime.now()
	print("\nStarted clustering at: " + timestart.strftime("%Y-%m-%d %H:%M:%S"))

	# No whitespace in margins
	plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
	plt.margins(0, 0)

	# Set font scale
	sns.set(font_scale=0.05)

	# Enter user inputs into clustermap...a Seaborn object we're giving the name 'g'
	g = sns.clustermap(table,
					   method=method,
					   figsize=(50, 50),
					   cmap=cmap,
					   cbar=False,
					   vmin=0,
					   vmax=vmax,
					   xticklabels=True)

	# Set axis labels
	ax = g.ax_heatmap
	ax.set(xlabel='STRAINS',
		   ylabel='CLUSTERS')

#how much time did it take? 
	timeend = datetime.now()
	timedelta = timeend - timestart
	print("Finished clustering at: " + timeend.strftime("%Y-%m-%d %H:%M:%S"))
	hours, remainder = divmod(timedelta.total_seconds(), 3600)
	minutes, seconds = divmod(remainder, 60)
	print('Total time to cluster: ' + '{:02}:{:02}:{:02}'.format(int(hours), int(minutes), int(seconds)))

#Save the figure
	timestart = datetime.now()
	print("\nStarted making png figure at: " + timestart.strftime("%Y-%m-%d %H:%M:%S"))
	plt.savefig('data/heatmap/cluster_heatmap_'+timestamp+'.png', dpi=600, bbox_inches='tight', pad_inches=0)
	timeend = datetime.now()
	timedelta = timeend - timestart
	print("Finished making figure at: " + timeend.strftime("%Y-%m-%d %H:%M:%S"))
	hours, remainder = divmod(timedelta.total_seconds(), 3600)
	minutes, seconds = divmod(remainder, 60)
	print('Total time to make figure: ' + '{:02}:{:02}:{:02}'.format(int(hours), int(minutes), int(seconds)))
	
	cmd = 'cp data/heatmap/cluster_heatmap_'+timestamp+'.png output/heatmap.png'
	os.system(cmd)
	
	'''
	┌---------------------------------------------------------------------------------┐
	| 4. Generating reordered strain list and table based on hierarchical clustering. |
	└---------------------------------------------------------------------------------┘
	'''

	if reorder_opt == 'Y':

		print("\nReordering the strain list and table based on hierarchical clustering...\n")

		# Retrieve the new order of the columns and re-order the table according to the clustering
		col_order = g.dendrogram_col.reordered_ind
		table = table[table.columns[col_order]]

		# Retrieve the new order of the rows and re-order the table according to the clustering
		row_order = g.dendrogram_row.reordered_ind
		table = table.iloc[row_order, :]

		ordered_protein_clusters_list = table.index.tolist()
		ordered_protein_clusters_list.extend(excluded_list)
		
		ordered_index = []
		for proteincluster in ordered_protein_clusters_list:
			ordered_index.append(proteincluster.split(" - ")[0])
		
		#copy the old strainlist.txt to the 'tab/archive folder for archiving.'
		cmd = str('cp strainlist.txt tables/archive/strainlist_'+timestamp+'.txt')
		os.system(cmd)
		
		# Produce the new reordered strain list
		ordered_strains = list(table.columns)
		with open('strainlist.txt', 'w') as f:
			f.write("\n".join(ordered_strains))

		# Produce the reordered table of protein counts in each cluster for the heat map
		table.to_csv('data/heatmap/cluster_hit_count_table_sorted_'+timestamp+'.tab', sep='\t')

		# Produce a reordered cluster_table.tab file
		# rename the original cluster_table.tab as unsorted with a time stamp
		cmd = str('cp tables/cluster_table.tab tables/archive/cluster_table_'+timestamp+'.tab')
		os.system(cmd)
		
		#make a new set of headers for the sorted proteins
		#first you need to open the old cluster_table.tab file and grab its header line.
		#then find where the 'strain count' is in the header...and use that to mark where 
		#annotations end and strain info begins. 
		with open('tables/cluster_table.tab') as input_file:
			head = next(input_file)
		input_file.close()
		header_items = head.split('\t')
		index_of_straincount = header_items.index('strain count')
		annotation_header_items = header_items[0:index_of_straincount + 1]
		
		#now extend the 'annotation header' with the newly resorted strain list headers
		annotation_header_items.extend(ordered_strains)

		#use that list of headers to sort the table.
		reord_clus_tab = pd.read_csv('tables/archive/cluster_table_'+timestamp+'.tab', sep='\t')[annotation_header_items]
		reord_clus_tab.set_index('CLUSTER', inplace=True)
		reord_clus_tab = reord_clus_tab.reindex(index=ordered_index)
		reord_clus_tab.to_csv('tables/cluster_table.tab', sep='\t')


	elif reorder_opt == 'N':
		print("You opted out of reordering a strain list and table based on hierarchical clustering.\n")
	print("Finished!")
