# strain-comp
Scripts to compare bacterial genomes

This is a set of scripts that can be run in sequential order to process NCBI 'cds' files to compare multiple strains to each other by gene content.


## Getting started ##

1.	Make a new folder (directory) for your project. 
  a.	Give the project folder a descriptive name like wn_iners_analysis_2024.
  b.	If it’s easier in the moment you could call it ‘iners’ and then go back and change the name to be a more descriptive one when you’re done.  But it’s good to let people know what you did and what is in the folder.

2.	Download the necessary data from NCBI.

3.	Unzip the downloaded genomic data.  
   a.	You should end up with a folder called **ncbi_dataset**.  
   b.	Inside the ncbi_dataset folder should be a folder called data.
   c.	Inside the data folder should be several folders with genome data for each strain.
1.	The genome datafiles should be named ‘cds_from_genomic.fna’
iv.	There should also be a file ‘assembly_data_report.jsonl’ in the data folder.
d.	Put the ncbi_dataset folder in your project folder.  Then get started with the rename.sh script.


## Getting genomic data using the NCBI datasets command line interface
In any command line environment where ncbi-datasets is installed you can type commands directly into the terminal to download genomic data and its associated metadata.

### Downloading genomes using a taxonomic name

If you want all the Gardnerella RefSeq annotated ORFs (cds) from genomes without MAGs and atypical genomes type:

  `$ datasets download genome taxon gardnerella --include cds --assembly-source RefSeq --exclude-atypical --annotated --exclude-multi-isolate --mag exclude --filename ncbi_dataset.zip`


### Download genomes using an accession number

If you know the accessions of a set of genomes.

  `$ datasets download genome accession GCF_01234567.1,GCF_022662295.1 --include cds --assembly-source RefSeq --filename ncbi_dataset.zip`

### Download several genomes from a file of accession numbers ###

You can put a lot of accession numbers into a single text file (one accession per line) and feed them into the download program using the --inputfile flag.  The file does not need to be named accessions.txt.

 `$ datasets download genome accession --inputfile accessions.txt --include cds --assembly-source RefSeq --filename ncbi_dataset.zip`
