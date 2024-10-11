# Understanding file inputs and outputs from the strain-comp pipeline

![Navarre lab logo2](/docs/logo.png)


Author: William Navarre


# Repository Contents

<table>
<tr><th>File/Directory</th><th>Description</th></tr>
<tr><td valign="top">/docs</td><td valign="top">Directory containing markdown files (like this one) and images that are found within the markdown docs. </td><tr>
<tr><td valign="top">/sample_data</td><td valign="top">Directory containing sample data for you to work with when testing the directory.</td><tr>
<tr><td valign="top">README.md</td><td valign="top">The main markdown file explaining the purpose and contents of the repository, listing of links to specific content.</td><tr>
<tr><td valign="top">1_unzip_ncbi.py</td><td valign="top">A script that renames data files downloaded from NCBI into something useable, generates a strain table, and establishes the folders needed for the next step. </td> <tr>
<tr><td valign="top">2_process_ncbi.py</td><td valign="top">Python script. This script takes the NCBI cds files and converts them to protein FASTA files (suffix .faa).  It also makes the final strainlist with the metadata used for the column headers.</td><tr>
<tr><td valign="top">3_mmseqcluster.py</td><td valign="top">Python script.  Invokes the mmseq2 program based on user input to cluster the proteins of the strains based on similarity.  It generates a 'cluster_metadata.tab' file and a 'reprsentative proteins.faa' file.</td><tr>
<tr><td valign="top">4_maketable.py</td><td valign="top">Python script. </td><tr>
<tr><td valign="top">5_heatmap.py</td><td valign="top">Python script. </td><tr>
<tr><td valign="top">6_formatxl.py</td><td valign="top">Python script. </td><tr>
<tr><td valign="top">7_geneorder.py</td><td valign="top">Python script to rearrange the data in the cluster_table cells to put genomic gene order up front.  This allows fast sorting of genes in a genome by their location and makes it easy to spot genomic islands. </td><tr>
<tr><td valign="top">8_COGadd.py</td><td valign="top">Python script to incorporate deepnog COG annotations and metatdata into the cluster_table.tab/xlsx files. </td><tr>
<tr><td valign="top">9_KEGGadd.py</td><td valign="top">Python script to incorporate KEGG annotations from their program 'kofamscan' into the cluster_table.tab/xlsx files. </td><tr>
<tr><td valign="top">10_tabletrimmer.py</td><td valign="top">Python script to remove strains from the analysis without having to start the whole process from step 1. </td><tr>
<tr><td valign="top">11_sequence_grabber.py</td><td valign="top">Python script to incorporate retrieve protein and/or nucleotide sequences based on their CLUSTER ID, their NCBI WP number, or their locus id. </td><tr>
</table>


# The files you need to input

##NCBI 'cds_from_genomic' files has the nucleotide sequences for all the genes in each of the strains you are analyzing. 

.

##NCBI's assembly_data_report.jsonl has most of the metadata for all the assemblies you downloaded.



# The most important output files and how to read them

## cluster_table.tab (cluster_table.xlsx)

.

### variant - cluster_table_COG

.

### variant - cluster_table_KEGG

.

## heatmap.png

x

## strainlist.txt


