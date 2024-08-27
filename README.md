# strain-comp
## Compare and group bacterial genomes by gene/protein content

![Navarre lab logo2](/docs/logo.png)

This is a set of scripts developed in the Navarre lab to rapidly comapre and group strains by gene content.  The output is an easy to read .xlsx worksheet that can be parsed in any number of ways.  Strains can be grouped easily using hierarchical clustering methods and common gene sets can be identified. Data from KEGG and COG/Deepnog/EGGnog can be added easily at later steps.  This type of analysis has been successfuly performed across over a thousand strains in a single run of the pipeline in less than an hour.  Analyzing a few hundred strains can be completed in less than a few minutes.


## Dependencies ##

All but one of these scripts are coded in Python.  The other script is a (bash or zsh) shell script.

These scripts depend on several other packages (each of which also depend on packages).  In total over 90 different packages will be installed to get these scripts to work.  Fortunately most of these are installed automatically by just intalling a few important top-level packages like biopython, pandas, and seaborn.
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


```
conda create -n whatever_you_want_to_name_your_environment -c conda-forge python=3.11.5 biopython mmseqs2 pandas seaborn xlsxwriter ncbi-datasets-cli
```

Or use the 'environment.yml' file provided in this repository. The environment will be named 'navpipe' unless you modify the environment.yml file in a text editor before using it to create the environment (_change the first line: name: navpipe to name: something_else_).


```
conda env create -f environment.yml
```

## Additional steps - make these scripts findable and executable. ##
Put the scripts into a directory that is searchable via your $PATH variable.  
Check their permissions to see if they are 'executable'.  If not you should make them executable with the following commands (from within the directory where the scripts are kept).


`chmod +x *.py`  #this will make all the python scripts executable


`chmod +x 1_rename.sh` #this will make the 1_rename.sh script executable

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

![strain-comp workflow](/docs/pipeline.png)

## Getting genomic data using the NCBI datasets command line interface
When ncbi-datasets-cli is installed you can type commands directly into the terminal to download genomic data and its associated metadata.

### Downloading genomes using a taxonomic name

If you want all the Gardnerella RefSeq annotated ORFs (cds) from genomes without MAGs and atypical genomes type:

```
datasets download genome taxon gardnerella --include cds --assembly-source RefSeq --exclude-atypical --annotated --exclude-multi-isolate --mag exclude --filename ncbi_dataset.zip
```

### Downloading genomes using accession numbers

If you know the accessions of a set of genomes.

```
datasets download genome accession GCF_01234567.1,GCF_022662295.1 --include cds --assembly-source RefSeq --filename ncbi_dataset.zip
```

### Downloading several genomes from a file of accession numbers ###

You can put a lot of accession numbers into a single text file (one accession per line) and feed them into the download program using the --inputfile flag.  The file does not need to be named accessions.txt.

```
datasets download genome accession --inputfile accessions.txt --include cds --assembly-source RefSeq --filename ncbi_dataset.zip
```


## Step 1 - prepare raw NCBI data files for use with 1_rename.sh 

Make a directory for your project and put the zip file of genomic data inside of it.  
In the terminal go to the project directory, activate the conda environment you have created for this pipeline, then type:

```
unzip ncbi_dataset.zip
```

then...

```
1_rename.sh
```


## Step 2 - generate nucleotide and protein files for each genome and to consolidate the metadata for all strains into a single file (strainlist.txt) with 2_process_ncbi.py
This script takes the poorly named NCBI cds files and renames them by the more readable locus tags for each strain, then it adds the GC% content for each nucleotide sequence and saves then in a folder 'fna'.  Then it translates all of the sequences into proteins (faa format) and saves them in the 'faa' folder.  It also makes the 'strainlist.txt' file that has the imnportant metadata for each strain. 
To execute the script simply type:

```
2_process_ncbi.py
```

## Step 3 - create a set of clustered proteins across all strains with _3_mmseqcluster.py_  
This script takes the protein files from all the strains and concatenates them.  It uses mmseq2 clusters the proteins together by relatedness (similar proteins from each strain are grouped together into a single cluster).  Each protein cluster group is given a unique identifier with a prefix you can define using the -n flag (default = CLUSTER).  


```
3_mmseqcluster.py -n PREFIX_FOR_CLUSTERS -p PERCENT_IDENTITY (default = 80)
```


example:


```
3_mmseqcluster.py -n LACTOBACILLUS -p 85
```

The output files are stored in a folder data/mmseq_output and are used by 4_maketable.py.


## Step 4 - generate a table of strains vs. protein clusters with 4_maketable.py.
This script combines output from the mmseq clustering with the straintable.txt file made by the process_ncbi.py script to make a tab file that compares each protein cluster across all strains in the analysis.  The first few columns of the output table have metadata about the protein cluster (ie. what the protein does, how big it is, its GC content) the remaining columns are for each strain.  If a protein is present in the strain the cell will give basic information about the protein.  If the protein is absent in the strain the table cell will be filled by an asterisk.  
Simply enter the following on the command prompt.

```
4_maketable.py
```

The output files are found in a folder called 'tables'.

> _You could stop here and use the cluster_table.tab output for most purposes._

## Step 5 - use hierarchical clustering to group strains and protein clusters together by presence/absence with _5_heatmap.py_. 
This script uses hierarchical clustering by the 'seaborn' library to group strains together by relatedness based on protein presence/absence.  It will also group proteins together based on their similarily in distribution across the strains. It will produce a graphic png file and will re-sort the strainlist.txt and tab file to put related strains and proteins next to each other.

To run the script with default settings type:

```
5_heatmap.py
```
You can alter the clustering methods and other parameters with the following flags:

```
-m', '--method', clustering method (average, ward, others, default=average)
-r', '--range', Min and max strain hits to use in clustering analysis (format: min#:max#, default = all)
-c', '--color', heatmap color (Blues, Reds, ...other matplotlib colors, default=Blues)
-hm', '-hitmax', Max hits expected per cluster (int, default = 1)
-l', '--label', Label the axes? (Y/N, default = Y)
-o', '--order', Reorder the strainlist and table after clustering? (Y/N, default = Y)
```

> _The 'cluster_table.tab' file will be rearranged with closely related strains next to each other in the table.  The original cluster_table.tab file is archived in a folder (tables/archived) with a datestamp.
> The remaining scripts below involve making the data in the cluster_table.tab pretty (formatxl), sortable by position (geneorder), or adding additional data from other annotation pipelines to it (COGadd and KEGGadd)._

**Below is an example of a clustered 'heatmap' produced by the heatmap.py script.  Strains are clustered along the x-axis and genes are clustered along the y-axis.**

![heatmap created by heatmap.py](/sample_data/heatmap.jpg)

## Step 6 - make the data tables easy to read in Excel format with _6_formatxl.py_.
This script formats the cluster_table.tab file in the 'tables' directory as an Excel file.  This both cuts the size of the table by approximately 50% but also makes it easy to read.  The xlsx file will appear in a folder called 'output'. To run the script simply type:

```
6_formatxl.py
```


## Step 7 - make a table that's easy to sort by gene order in the genome with _7_geneorder.py_.
The cluster_table.tab file created by 4_maketable.py serves as the primary 'database' created by this pipeline.  Each cell holds data for a given protein (row) in a given strain (column).  The data in these cells can be easily rearranged so that the 'gene position' appears up front in each cell.  This makes it easy to sort genes by their position in a genome by simply selecting that column and asking Excel to sort it in ascending or descending order. 

To execute the script simply type:

```
7_geneorder.py
```

The output will appear in the `tables` directory as `cluster_table_geneordered.tab`

# Adding data from other annotation tools like kofamscan (KEGG) and deepnog (COG)
We have developed scripts to take output from kofamscan and deepnog searches and incorporate them into the cluster_table.tab file (which can also be parsed by formatxl.py).

### Incorporating kofamscan and KEGG (Kyoto Encyclopedia of Genes and Genomes) annotations ###

Because kofamscan is a separate program maintained by KEGG, and because it has a large database file that comes with it, we recommend doing the kofamscan search in a separate conda environment specifically set up for kofamscan.  Then take the resulting output file (in detail-tsv format, _see instructions below_) and then use KEGGadd.py to incorporate that data into the cluster_table.tab file. 

[Instructions here](/docs/KEGGadd_instructions.md#adding-kegg-annotations-to-the-cluster_tabletab-file-with-kofamscan-and-keggaddpy)

### Incorporating deepnog and COG (Clusters of Orthologous Groups) annotations ###

Because deepnog is a separate program and because it has a large database file required to use it, we recommend doing the deepnog search in a separate conda environment specifically set up for deepnog and our COGadd.py script to incorporate COG annotation data into the cluster_table.tab file. 

[Instructions here](/docs/COGadd_instructions.md#adding-cog-annotations-to-the-cluster_tabletab-file-with-deepnog-and-cogaddpy)
