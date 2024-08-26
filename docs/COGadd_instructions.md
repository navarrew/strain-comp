# Adding COG annotations to the cluster_table.tab file with deepnog and 8_COGadd.py #

[DEEPNOG](https://github.com/univieCUBE/deepnog) is a package that can rapidly assigns proteins, provided in FASTA format, to one or more putative COGs (Clusters of Orthologous Groups).

>Roman Feldbauer, Lukas Gosch, Lukas Lüftinger, Patrick Hyden, Arthur Flexer, Thomas Rattei, DeepNOG: Fast and accurate protein orthologous group assignment, Bioinformatics, 2020, btaa1051, https://doi.org/10.1093/bioinformatics/btaa1051

When combined with the NCBI annotations, this extra metadata is useful to reveal functional links between proteins and sometimes succeeds in assigning a function to a protein that NCBI missed.  E-values are provided for each assignment to give the user a sense of how confident the prediction is.

We developed the COGadd.py to run the Deepnog prediction and to import associated metadata from a table from NOG. 

Because DEEPNOG was developed by somebody else and the metadata is a fairly large table we suggest creating a separate conda environment for running deepnog.

## Setting up the environment for COGadd.py. ##


## Files you will need to provide. ##
The only file you need to provide that isn't already in place after running the first five scripts in the pipeline is a table of metadata for each COG ID.  This table is called **'cog-20.def.tab'**.  We provide a recent version of the table in this repository under the 'sample data' folder.  


The file can also be found here [https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.def.tab](https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.def.tab)


The 8_COGadd.py script, when run from project directory, will first grab protein sequence data from the mmseq2 output file called 'cluster_representative_sequences.faa'.  It can be found in the data/mmseq-output/ directory after you run the 3_mmseqcluster.py script.  
After assigning each protein a COG ID it will output this data to a table now located in the data/COG/ directory that it creates. 

Then the script will link metatdata for each COG from the **cog-20.def.tab** table.  

Finally it will use the Cluster ID/COG ID associations to weave the metadata for each COG into several new columns on the cluster_table.tab file.  

A new file will be made called **cluster_table_COG.tab**.  