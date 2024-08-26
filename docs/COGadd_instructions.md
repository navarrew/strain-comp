# Adding COG annotations to the cluster_table.tab file with deepnog and 8_COGadd.py #

[DEEPNOG](https://github.com/univieCUBE/deepnog) is a package that can rapidly assigns proteins, provided in FASTA format, to one or more putative COGs (Clusters of Orthologous Groups).

>Roman Feldbauer, Lukas Gosch, Lukas Lüftinger, Patrick Hyden, Arthur Flexer, Thomas Rattei, DeepNOG: Fast and accurate protein orthologous group assignment, Bioinformatics, 2020, btaa1051, https://doi.org/10.1093/bioinformatics/btaa1051

When combined with the NCBI annotations, this extra metadata is useful to reveal functional links between proteins and sometimes succeeds in assigning a function to a protein that NCBI missed.  E-values are provided for each assignment to give the user a sense of how confident the prediction is.

We developed the COGadd.py to run the Deepnog prediction and to import associated metadata from a table from NOG. 

Because DEEPNOG was developed by somebody else and the metadata is a fairly large table we suggest creating a separate conda environment for running deepnog.

## Setting up the environment for COGadd.py. ##

This script uses standard Python libraries.  Setting up a new environment should be simple.  Create it by typing:

`$ conda env create -n deepnog python=3.11.5 deepnog`


## Files you will need to provide. ##
The only file you need to provide that isn't already in place after running the first five scripts in the pipeline is a table of metadata for each COG ID.  This table is called **'cog-20.def.tab'**.  We provide a recent version of the table in this repository under the 'sample data' folder.  

The cog-20.def.tab file can also be found here [https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.def.tab](https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.def.tab)


## Input parameters. ##

`$ 8_COGadd.py -co 20 -d where/you/put/the/cog-20.def.tab_file`


- '-f', '--fasta'    name/location of input faa file (default = data/mmseq-output/cluster_representative_sequences.faa)
- '-t', '--tab'    location of cluster_table.tab (default = tables/cluster_table.tab
- '-od', '--outdir',    set name/location of output directory without leading or trailing '/' symbols or dots (default = data/COG)
- '-o', '--outfile',    set name of the output file (default = tables/cluster_table_with_COGs.tab)
- '-d', '--definitions'    provide location of the COG definitions table (default = project directory)
- '-co', '--confidence'    set the confidence threshold for deepnog prediction in percent (should enter at least 20)


## How the script operates. ##
The 8_COGadd.py script, when run from project directory, will invoke deepnog with the parameters you include.  

First the script will tell deepnog to grab protein sequence data from the mmseq2 output file called 'data/mmseq-output/cluster_representative_sequences.faa' that was created when you ran the 3_mmseqcluster.py script.  

After assigning a COG ID to each protein cluster (if possible) it will output this data to a table (**data/COG/deepnog_output.tab**). 

Then the script will grab metatdata for each COG from the **cog-20.def.tab** table and put it together with the deepnog_output table to make a new table called **data/COG/COG_annotations.tab**.

As a last step it will merge the **data/COG/COG_annotations.tab** file and the **tables/cluster_table.tab** file.  

A new file will be made called **tables/cluster_table_COG.tab**.  


## Understanding the output. ##
