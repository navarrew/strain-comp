# strain-comp
Scripts to compare bacterial genomes

This is a set of scripts that can be run in sequential order to process NCBI 'cds' files to compare multiple strains to each other by gene content.

## Getting started with installing dependencies ##

Key to having these scripts work are the following packages/software - all but one use Python, which also must be installed.
  1. python (see conda/anaconda below)
  2. biopython [https://github.com/biopython/biopython]
  3. pandas [https://github.com/pandas-dev/pandas]
  4. seaborn [https://github.com/mwaskom/seaborn]
  5. numpy [https://github.com/numpy/numpy]
  6. scipy [https://github.com/scipy/scipy]
  7. fastcluster [https://github.com/fastcluster/fastcluster?tab=readme-ov-file]
  8. xlsxwriter [https://github.com/jmcnamara/XlsxWriter]
  9. ncbi-datasets-cli [https://github.com/ncbi/datasets]
  10. mmseqs2 [https://github.com/soedinglab/MMseqs2]

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

 
