# strain-comp
## Scripts to compare bacterial genomes

This is a set of scripts developed in the Navarre lab to rapidly comapre and group strains by gene content.  The output is an easy to read .xlsx worksheet that can be parsed in any number of ways.  Strains can be grouped easily using hierarchical clustering methods and common gene sets can be identified. Data from KEGG and COG/Deepnog/EGGnog can be added easily at later steps.  This type of analysis has been successfuly performed across over a thousand strains in a single run of the pipeline in less than an hour.  Analyzing a few hundred strains can be completed in less than a few minutes.

## Getting started with installing dependencies ##

Key to having these scripts work are the following packages/software - all but one use Python, which also must be installed.
  1. python >= 3.11 (see conda/anaconda below)
  2. ncbi-datasets-cli [https://github.com/ncbi/datasets] - _conda install -c conda-forge ncbi-datasets-cli_
  3. biopython [https://github.com/biopython/biopython] - _conda install -c conda-forge biopython_
  4. mmseqs2 [https://github.com/soedinglab/MMseqs2] - _conda install -c conda-forge -c bioconda mmseqs2_
  5. pandas [https://github.com/pandas-dev/pandas] - _conda install -c conda-forge pandas_
  6. seaborn [https://github.com/mwaskom/seaborn] - _conda install -c conda-forge seaborn_
  7. xlsxwriter [https://github.com/jmcnamara/XlsxWriter]
  8. numpy [https://github.com/numpy/numpy] *seems to be installed already when installing biopython
  9. scipy [https://github.com/scipy/scipy] *seems to be installed when installing seaborn
  10. matplotlib [https://github.com/matplotlib/matplotlib] *seems to be installed when installing seaborn
  11. fastcluster [https://github.com/fastcluster/fastcluster?tab=readme-ov-file]
      
### Using conda/anaconda [https://anaconda.org] ###
We have been using anaconda to manage our python environments.  It's worked well across both mac and linux platforms and many of the issues we had in earlier days have been resolved.

We have posted the **'environment.yml'** file that includes all packages we have installed.
This file includes some packages that we tried but found weren't necessary or useful.
You can use our environment.yml file or just download all the packages above to get a fresh start.


## Getting started analyzing genomes ##

1.	Make a new directory (we'll call the project directory) for your project.

2.	Download the necessary genomic data from NCBI.

5.	Unzip the downloaded genomic data.  
  a.	You should end up with a folder called **ncbi_dataset**.
  b.	Inside the ncbi_dataset folder should be a folder called data.
  c.	Inside the data folder should be several folders with genome data for each strain.
    1.	The genome datafiles in each will be named ‘cds_from_genomic.fna’
    2.	There should also be a file ‘assembly_data_report.jsonl’ in the data folder.
    3.	Put the ncbi_dataset folder in your project directory.


## Getting genomic data using the NCBI datasets command line interface
In any command line environment where ncbi-datasets-cli is installed you can type commands directly into the terminal to download genomic data and its associated metadata.

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

## Step 1 - _1_rename.sh_ to prepare raw NCBI data files for use

X

## Step 2 - _2_process_ncbi.py_ to generate nucleotide and protein files for each genome and to consolidate the metadata for all strains into a single file (strainlist.txt)

X

## Step 3 - _3_mmseq.py_ to generate clusters of similar proteins across all strains. 

X

## Step 4 - _4_maketable.py_ to generate a table of strains vs. clusters.

X

## Step 5 - _5_heatmap.py_ - use hierarchical clustering to group strains and protein clusters together by similarity of presence/absence. 

X

## Step 6 - _6_geneorder.py_ -

X

## Step 7 - _7_formatxl.py_ -

