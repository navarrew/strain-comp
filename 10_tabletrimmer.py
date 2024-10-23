#!/usr/bin/env python

import pandas as pd
import os

#this script allows you to take a strainlist that you modify by hand to include/exclude
#isolates from a larger list of strains...then use that new, shorter, list to make a new
#set of tab tables that can be used for heatmapping, etc.  It gets rid of proteins 
#specific to strains that are no longer present in the list.

#to start you need to:
#start a new main directory for this project 
#copy over your tab files from the original project to a tab directory within the new main dir
#copy over your original strainlist.txt (or strainlist_sorted.txt), parse it however 
#you want by removing or adding lines, then put it in the main directory.
#make sure it's named 'strainlist.txt'
#then go!  Your newly trimmed tab files will appear in a new 'trimmed' directory


if __name__ == '__main__':

# FIRST LETS SET UP OUR DIRECTORIES
# Importing the strainlist to rearrange the table.
    directory_list = []
    for entry in os.scandir('.'):
        if entry.is_dir():
            directory_list.append(entry.name)
        
    #if not already in the current directory make the directories reports and modified_tabs
    if "trimmed" not in directory_list:
        os.mkdir("trimmed")

    with open('strainlist.txt') as f:
        str_list = [line.rstrip() for line in f]
        # Count the total number of strains in new strainlist
        strain_max = len(str_list)

    print("You have a total number of " + str(strain_max))

    full_sorted_heads = ['CLUSTER', 'NCBI gene names', 'NCBI annotations', 'GCpct', 'GC spread', 'protein length', 'flags', 'total count', 'strain count'] + str_list
    count_sorted_heads = ['CLUSTER', 'total count', 'strain count'] + str_list

    reord_count_tab = pd.read_csv('tab/cluster_hit_count_table.tab', sep='\t')[count_sorted_heads]
    reord_count_tab.set_index('CLUSTER', inplace=True)
    #now counting total hits and strain counts
    reord_count_tab['total count'] = reord_count_tab[str_list].sum(axis=1)
    reord_count_tab['strain count'] = strain_max - (reord_count_tab[str_list].eq(0).sum(axis=1))
    total_counts = reord_count_tab['total count'].values.tolist()
    strain_counts = reord_count_tab['strain count'].values.tolist()
    reord_count_tab = reord_count_tab.loc[reord_count_tab['total count'] != 0]
    reord_count_tab.to_csv('trimmed/cluster_hit_count_table.tab', sep='\t')

    reord_full_tab = pd.read_csv('tab/cluster_table.tab', sep='\t')[full_sorted_heads]
    reord_full_tab.set_index('CLUSTER', inplace=True)
    reord_full_tab['total count'] = total_counts
    reord_full_tab['strain count'] = strain_counts
    reord_full_tab = reord_full_tab.loc[reord_full_tab['total count'] != 0]
    reord_full_tab.to_csv('trimmed/cluster_table.tab', sep='\t')
    
    os.rename('tab', 'old-tab')
    os.rename('trimmed', 'tab')
