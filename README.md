# strain-comp
## Scripts to compare bacterial genomes

This is a set of scripts developed in the Navarre lab to rapidly comapre and group strains by gene content.  The output is an easy to read .xlsx worksheet that can be parsed in any number of ways.  Strains can be grouped easily using hierarchical clustering methods and common gene sets can be identified. Data from KEGG and COG/Deepnog/EGGnog can be added easily at later steps.  This type of analysis has been successfuly performed across over a thousand strains in a single run of the pipeline in less than an hour.  Analyzing a few hundred strains can be completed in less than a few minutes.

## Dependencies ##

All but one of these scripts are coded in Python, which also must be installed.  The first script is a (bash or zsh) shell script.
  1. **python** >= 3.11 (see conda/anaconda below)
  2. **ncbi-datasets-cli** [https://github.com/ncbi/datasets] - _conda install -c conda-forge ncbi-datasets-cli_
  3. **biopython** [https://github.com/biopython/biopython] - _conda install -c conda-forge biopython_
  4. **mmseqs2** [https://github.com/soedinglab/MMseqs2] - _conda install -c conda-forge -c bioconda mmseqs2_
  5. **pandas** [https://github.com/pandas-dev/pandas] - _conda install -c conda-forge pandas_
  6. **seaborn** [https://github.com/mwaskom/seaborn] - _conda install -c conda-forge seaborn_
  7. **xlsxwriter** [https://github.com/jmcnamara/XlsxWriter] - _conda install -c conda-forge xlsxwriter_
  8. numpy [https://github.com/numpy/numpy] *probably was installed already when installing biopython
  9. scipy [https://github.com/scipy/scipy] *probably was installed when installing seaborn
  10. matplotlib [https://github.com/matplotlib/matplotlib] *probably was installed when installing seaborn
  11. (optional) fastcluster [https://github.com/fastcluster/fastcluster?tab=readme-ov-file]  _fastcluster may need to be constructed from source.  It otherwise is pinned to python version 3.5, which isn't compatible with most other parts of this pipeline.  It isn't essential but makes hierarchical clustering much faster in heatmap.py if the dataset is very large._

      
### Using conda/anaconda [https://anaconda.org] ###
We have been using anaconda to manage our python environments.  It's worked well across both mac and linux platforms and many of the issues we had in earlier days have been resolved.

After installing conda you should create a new environment for this pipeline.  We have posted the **'environment.yml'** file that includes all packages we have installed via conda and that can be used to recreate a functional environment.


`conda create -n whatever_you_want_to_name_your_environment -c conda-forge python=3.11 biopython mmseqs2 pandas seaborn xlsxwriter ncbi-datasets-cli`

or use the 'environment.yml' file provided here.  
`conda env create -f environment.yml`



## Getting started analyzing genomes - the basic steps involved ##

1.	Make a new directory (we'll call the project directory) for your project.
2.	Download data using the ncbi-datasets command line interface.
3.	Put the data in your project directory.  Unzip it.
4.	Rename the datafiles using the 1_rename.sh script.
5.	Extract metadata, translate nucleotide to protein seqeuences with 2_process_ncbi.py script.
6.	Cluster similar proteins and name the clusters with 3_mmseqcluster.py
7.	Make a tab delimited file that provides metadata for all protein clusters and compares them across all strains with the 4_maketable.py script.
8.	Use hierarchical clustering to group strains and protein clusters together with the 5_heatmap.py script
9.	Get pretty readable output with the 6_formatxl.py script.


## Getting genomic data using the NCBI datasets command line interface
When ncbi-datasets-cli is installed you can type commands directly into the terminal to download genomic data and its associated metadata.

### Downloading genomes using a taxonomic name

If you want all the Gardnerella RefSeq annotated ORFs (cds) from genomes without MAGs and atypical genomes type:

  `$ datasets download genome taxon gardnerella --include cds --assembly-source RefSeq --exclude-atypical --annotated --exclude-multi-isolate --mag exclude --filename ncbi_dataset.zip`

### Downloading genomes using accession numbers

If you know the accessions of a set of genomes.

  `$ datasets download genome accession GCF_01234567.1,GCF_022662295.1 --include cds --assembly-source RefSeq --filename ncbi_dataset.zip`

### Downloading several genomes from a file of accession numbers ###

You can put a lot of accession numbers into a single text file (one accession per line) and feed them into the download program using the --inputfile flag.  The file does not need to be named accessions.txt.

 `$ datasets download genome accession --inputfile accessions.txt --include cds --assembly-source RefSeq --filename ncbi_dataset.zip`

### Downloading other types of genomic data ###
If you want more than just the open reading frames use the --include tag.

To get the 'genbank flat file' (gbff) format:
 `$ datasets download genome accession --inputfile accessions.txt --include gbff --assembly-source RefSeq --filename ncbi_dataset.zip`

To get the 'genbank flat file' (gbff) and cds formats in the same package:
 `$ datasets download genome accession --inputfile accessions.txt --include cds,gbff --assembly-source RefSeq --filename ncbi_dataset.zip`

To get the 'genbank flat file' (gbff) and cds formats in the same package:
 `$ datasets download genome accession --inputfile accessions.txt **--include genome** --assembly-source RefSeq --filename ncbi_dataset.zip`



## Step 1 - prepare raw NCBI data files for use with 1_rename.sh 

Make a directory for your project and put the zip file of genomic data inside of it.  
In the terminal go to the project directory, activate the conda environment you have created for this pipeline, then type:

`unzip ncbi_dataset.zip`

then...

`1_rename.sh`


## Step 2 - generate nucleotide and protein files for each genome and to consolidate the metadata for all strains into a single file (strainlist.txt) with 2_process_ncbi.py

`2_process_ncbi.py`

## Step 3 - create a set of clustered proteins across all strains with _3_mmseqcluster.py_  
This script takes the protein files from all the strains and concatenates them.  It uses mmseq2 clusters the proteins together by relatedness (similar proteins from each strain are grouped together into a single cluster).  Each protein cluster group is given a unique identifier with a prefix you can define using the -n flag (default = CLUSTER).  The output files are stored in a folder data/mmseq_output and are used by 4_maketable.py.


`3_mmseqcluster.py -n PREFIX_FOR_CLUSTERS -p PERCENT_IDENTITY (default = 80)`
example:
`3_mmseqcluster.py -n LACTOBACILLUS -p 85`



## Step 4 - _4_maketable.py_ to generate a table of strains vs. clusters.
This script combines output from the mmseq clustering with the straintable.txt file made by the process_ncbi.py script to make a tab file that compares each protein cluster across all strains in the analysis.  The first few columns of the output table have metadata about the protein cluster (ie. what the protein does, how big it is, its GC content) the remaining columns are for each strain.  If a protein is present in the strain the cell will give basic information about the protein.  If the protein is absent in the strain the table cell will be filled by an asterisk.  

## Step 5 - _5_heatmap.py_ - use hierarchical clustering to group strains and protein clusters together by similarity of presence/absence. 
This script uses hierarchical clustering by the 'seaborn' library to group strains together by relatedness based on protein presence/absence.  It will also group proteins together based on their similarily in distribution across the strains. It will produce a graphic png file and will re-sort the strainlist.txt and tab file to put related strains and proteins next to each other.

## Step 6 - _6_formatxl.py_ -

X

## Step 7 - _7_geneorder.py_ -

