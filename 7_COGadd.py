#!/usr/bin/env python
import subprocess
import shlex
import os
import argparse
from typing import Iterator, Tuple

# code modified from https://stackoverflow.com/questions/72938098/list-of-entries-files-and-folders-in-a-directory-tree-by-os-scandir-in-pytho
def tree2list(directory: str) -> Iterator[Tuple[str, str, str]]:

    for i in os.scandir(directory):
        if i.is_dir():
            yield i.path[2:]
            yield from tree2list(i.path)
        else:
            yield i.path[2:]
    
    
def make_COG_function_entry(COG_function_codes):
    COG_dict = dict(J="Translation, ribosomal structure and biogenesis",
                    A="RNA processing and modification",
                    K="Transcription",
                    L="Replication, recombination and repair",
                    B="Chromatin structure and dynamics",
                    D="Cell cycle control, cell division, chromosome partitioning",
                    Y="Nuclear structure",
                    V="Defense mechanisms",
                    T="Signal transduction mechanisms",
                    M="Cell wall/membrane/envelope biogenesis",
                    N="Cell motility",
                    Z="Cytoskeleton",
                    W="Extracellular structures",
                    U="Intracellular trafficking, secretion, and vesicular transport",
                    O="Posttranslational modification, protein turnover, chaperones",
                    X="Mobilome: prophages, transposons",
                    C="Energy production and conversion",
                    G="Carbohydrate transport and metabolism",
                    E="Amino acid transport and metabolism",
                    F="Nucleotide transport and metabolism",
                    H="Coenzyme transport and metabolism",
                    I="Lipid transport and metabolism",
                    P="Inorganic ion transport and metabolism",
                    Q="Secondary metabolites biosynthesis, transport and catabolism",
                    R="General function prediction only",
                    S="Function unknown")
    codes = list(COG_function_codes)
    function_list = []
    for code in codes:
        function = COG_dict[code]
        function_list.append(function)
    function_description = "; ".join(function_list)
    function_output = str(COG_function_codes + ": " + function_description)
    return(function_output)


if __name__ == '__main__':

    directory_list = list(tree2list("."))
    cwd = os.getcwd()
    if "Dropbox/bio" in cwd:
        Dropbox_directory_first_part = cwd.split('Dropbox/bio')[0]

#     dropbox_directory_list = list(tree2list("~/Dropbox/bio"))
#     print(dropbox_directory_list)


# STEP - get info from the USER

    instructions = "##############\n# COGadd.py  #\n##############\n" \
                    "\nThis program will run deepNOG on the cluster_representative_sequences.faa file" \
                    "\nin the main project directory.  This will output files into the COG directory.\n\n" \
                    "RULES:  You must...\n" \
                    "1. be in an environment where deepnog has been installed.\n" \
                    "2. be in a project directory with a tab folder that has cluster_table.tab in it.\n" \
                    "3. have the file cluster_representative_sequences.faa in the project folder.\n"\
                    "4. know where your cog-20.def.tab file is.  We have a default location set.\n" \
                    "....if you violate these rules you have to tell the script where these files are."

    parser = argparse.ArgumentParser(add_help=True, description=instructions,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-f', '--fasta', action='store', dest='infile', help = "name/location of input faa file (default = cluster_representative_sequences.faa).")
    parser.add_argument('-t', '--tab', action='store', dest='tabfolder', help = "location of cluster_table.tab (default = tab/cluster_table.tab)")
    parser.add_argument('-od', '--outdir', action='store', dest='outdir', help = "set name/location of output directory without leading or trailing '/' symbols or dots (default = data/COG).")
    parser.add_argument('-o', '--outfile', action='store', dest='outfile', help = "set name of the output file (default = cluster_table_with_COGs.tab)")
    parser.add_argument('-d', '--definitions', action='store', dest='def_loc', help = "location of the COG definitions table (default = (custom)/Dropbox/bio/tools/COG/cog-20.def.tab)")
    parser.add_argument('-co', '--confidence', action='store', dest='confidence', help = "Set confidence threshold for deepnog prediction.")

    args = parser.parse_args() #args is a 'parser object' that stores flags entered in on the command line

# if there is nothing in the command line we can set default values below

    if args.infile:
        infile = args.infile
    else: infile = "data/mmseq_output/cluster_representative_sequences.faa"

    if args.outdir:
        outdir = args.outdir
        if outdir[-1] == "/": #remove a trailing '/' from the directory entry if user added it
            outdir = outdir[0:-1] 
    else: outdir = "data/COG"

    if args.tabfolder:
        tabfolder = args.tabfolder
    else: tabfolder = "tables/cluster_table.tab"

    if args.outfile:
        outfile = args.outfile
    else: outfile = "default"

    if args.def_loc:
        def_loc = args.def_loc
    else: def_loc = Dropbox_directory_first_part + 'Dropbox/bio/tools/COG/cog-20.def.tab'

    from pathlib import Path
    my_file = Path(def_loc)
    if not my_file.is_file():
        print('\nCannot find your definitions file at:' + def_loc + '!')
        print('Please enter an absolute path to the COG definitions file.')
        print('If it helps, your current working directory is:')
        print(cwd+'\n')
        exit()

# NEXT STEP - setting us up with directories we'll need if they don't exist already.
    directory_list = list(tree2list("."))

#     #if not already in the current directory make the output directory

    if outdir not in directory_list:
        os.mkdir(outdir)

        
# NEXT STEP - Now we make the command to send to deepnog
    deepnog_command_option_list = []
    #setting the deepnog output filename:
    if args.confidence:
        confidence = args.confidence
        deepnog_output_filename = str( outdir + '/deepnog_output_' +str(confidence) + '.tab')
    else:
        confidence = "0"
        deepnog_output_filename = str( outdir + '/deepnog_output.tab')
    deepnog_command_option_list.append(str("--out " + deepnog_output_filename))
    if args.confidence:
        deepnog_command_option_list.append("-c " + str(args.confidence))
    deepnog_command_option_list.append(str(infile))
    deepnog_command = "deepnog infer -of tsv " + str(" ".join(deepnog_command_option_list))


    print("Sending the following command to deepnog:\n> " + deepnog_command)
    t = subprocess.run(shlex.split(deepnog_command))
#we're using the subprocess.run module to run deepnog a little bit 'outside' of the script
#shlex.split module makes the deepnog_command readable by subprocess.run() by splitting it into tokens

# NEXT STEP- take the output of deepnog and format it for our purposes with CLUSTER as its own column.
# Read in the COG tsv list
    with open(deepnog_output_filename) as f:
        cluster_info_list = [line.rstrip() for line in f]
        new_cluster_output_list = []
        for cluster_info in cluster_info_list:
            cluster_info_exploded = cluster_info.split('\t')
            cluster_name_and_fasta_header = cluster_info_exploded[0]
            cluster_name = "_".join(cluster_name_and_fasta_header.split('_')[0:2])
            cluster_fasta_header = "_".join(cluster_name_and_fasta_header.split('_')[2:])
            if len(cluster_info_exploded) == 3:
                cluster_cog_assignment = cluster_info_exploded[1]
                cluster_confidence = cluster_info_exploded[2]
                if cluster_confidence != 'confidence':
                    cluster_confidence = f'{float(cluster_info_exploded[2]):.3f}'
            else:
                cluster_cog_assignment = "*"
                cluster_confidence = "*"
                
            new_cluster_output_line = str(cluster_name + "\t" +
                                          cluster_fasta_header + "\t" +
                                          cluster_cog_assignment + "\t" +
                                          cluster_confidence)
            new_cluster_output_list.append(new_cluster_output_line)

#if the output list has the default header given by deepNOG
    old_header = new_cluster_output_list[0]
    old_header_start = old_header.split("\t")[0]
    if old_header_start == "sequence_id":
    #remove the header
        del new_cluster_output_list[0]
    #replace with a new header (at position zero in the list)
        new_cluster_output_list.insert(0, "CLUSTER\tFASTA header\tCOG_id\tCOG_deepnog_confidence (limit " + str(confidence) + ")\tCOG_functional_category\tCOG_definition\tCOG_gene_name\tCOG_pathway\tPDB")

#write it all out to a new half-finished temporary file
    output_temp_file = open(outdir + '/COG_annotations_tmp.tab', 'w')
    for cluster_line in new_cluster_output_list:
        output_temp_file.write(cluster_line + '\n')

    f.close()
    output_temp_file.close()

# NEXT STEP - read the COG definition file, make a dictionary out of it

    COG_definition_dictionary = {}
    with open(def_loc, 'r') as f:
        COGdef_lines = f.readlines()
        for line in COGdef_lines:
            line = line.rstrip()
            definition_list = line.split('\t')
            if len(definition_list) == 3:
                definition_list.extend(["*","*","*","*"])
            if len(definition_list) == 4:
                definition_list.extend(["*","*","*"])
            if len(definition_list) == 5:
                definition_list.extend(["*","*"])
            if len(definition_list) == 6:
                definition_list.extend(["*"])
            for i, n in enumerate(definition_list):
                if n == "":
                    definition_list[i] = "*"
            #remove pubmed IDs
            pubmed = definition_list.pop(5)
            #pop the COG_ID from the front of the line
            COGID = definition_list.pop(0)
            #pop the COG code (single letters) from the front of the list
            COG_function_codes = definition_list.pop(0)
            #use the make_COG_function_entry function to match the COG letter to its overall function
            COG_function_description = make_COG_function_entry(COG_function_codes)
            #reinsert the COG letter and its definition at position 0
            definition_list.insert(0, COG_function_description)
            definition_string = "\t".join(definition_list)
            COG_definition_dictionary[COGID] = definition_string
    f.close()

#NEXT STEP - append the definitions to the
    cluster_cog_definitions_file = open(outdir + '/COG_annotations.tab', 'w')
    cluster_cog_definitions_list = []
    with open(outdir + '/COG_annotations_tmp.tab', 'r') as g:
        definition_to_fill_lines = g.readlines()
        header_line = definition_to_fill_lines[0]
        cluster_cog_definitions_file.write(header_line)
        cluster_cog_definitions_list.append(header_line)
        for line in definition_to_fill_lines[1:]:
            line = line.rstrip()
            COG_hit_list = line.split('\t')
            COG_key_id = COG_hit_list[2]
            if COG_key_id in COG_definition_dictionary:
                outputline = line + "\t" + COG_definition_dictionary[COG_key_id] + "\n"
            else:
                outputline = line + ('\t*' * 5) + "\n"
            cluster_cog_definitions_file.write(outputline)
            cluster_cog_definitions_list.append(outputline)
    g.close()
    
    remove_command = 'rm ' + outdir + '/COG_annotations_tmp.tab'
    os.system(remove_command)
    
#now make a new dictionary with key = CLUSTER_XXX and value = COG definitions
    cluster_COG_definition_dictionary = {}
    for cog_definition_item in cluster_cog_definitions_list:
        cluster_id_key = cog_definition_item.split("\t")[0]
        cluster_cog_defline = '\t'.join(cog_definition_item.split("\t")[2:])
        cluster_COG_definition_dictionary[cluster_id_key] = cluster_cog_defline
        
#now add COG columns to the cluster_table.tab file
    if 'tables' in directory_list:
        cluster_table_with_COGs = open('tables/cluster_table_COG.tab', 'w')
    else:
         cluster_table_with_COGs = open('data/COG/cluster_table_COG.tab', 'w')
         print("WARNING: Could not find tables directory...saving output to data/COG/.")
    with open(tabfolder, 'r') as h:
        cluster_table_line = h.readlines()
        for line in cluster_table_line:
            cluster_id = line.split("\t")[0]
            front_of_line = '\t'.join(line.split("\t")[0:3])
            end_of_line = '\t'.join(line.split("\t")[3:])
            middle_of_line = cluster_COG_definition_dictionary[cluster_id]
            newline = front_of_line + "\t" + middle_of_line.rstrip() + "\t" + end_of_line
            cluster_table_with_COGs.write(newline)
    h.close()
    cluster_table_with_COGs.close()

#now add COG columns to the cluster_table_with_keggs.tab file
    if 'tables/cluster_table_with_keggs.tab' in directory_list:
        cluster_table_COGs_and_KEGG = open('tables/cluster_table_COG_KEGG.tab', 'w') 
        with open('tables/cluster_table_KEGG.tab', 'r') as k:
            kegg_cluster_table_line = k.readlines()
            for line in kegg_cluster_table_line:
                cluster_id = line.split("\t")[0]
                front_of_line = '\t'.join(line.split("\t")[0:5])
                end_of_line = '\t'.join(line.split("\t")[5:])
                middle_of_line = cluster_COG_definition_dictionary[cluster_id]
                newline = front_of_line + "\t" + middle_of_line.rstrip() + "\t" + end_of_line
                cluster_table_COGs_and_KEGG.write(newline)
        k.close()
        cluster_table_COGs_and_KEGG.close()
    
