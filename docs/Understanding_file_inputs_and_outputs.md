# Understanding file inputs and outputs from the strain-comp pipeline

![Navarre lab logo2](/docs/logo.png)


Author: William Navarre


# Repository Contents

<table>
<tr><th>File/Directory</th><th>Description</th></tr>
<tr><td valign="top">/docs</td><td valign="top">Directory containing markdown files (like this one) and images that are found within the markdown docs. </td><tr>
<tr><td valign="top">/sample_data</td><td valign="top">Directory containing sample data for you to work with when testing the directory.</td><tr>
<tr><td valign="top">README.md</td><td valign="top">The main markdown file explaining the purpose and contents of the repository, listing of links to specific content.</td><tr>
<tr><td valign="top">1_rename.sh</td><td valign="top">A shell (bash or zsh) script that renames data files downloaded from NCBI into something useable. </td><tr>
<tr><td valign="top">2_process_ncbi.py</td><td valign="top">Python script. </td><tr>
<tr><td valign="top">3_mmseqcluster.py</td><td valign="top">Python script. </td><tr>
<tr><td valign="top">4_maketable.py</td><td valign="top">Python script. </td><tr>
<tr><td valign="top">5_heatmap.py</td><td valign="top">Python script. </td><tr>
<tr><td valign="top">6_formatxl.py</td><td valign="top">Python script. </td><tr>
<tr><td valign="top">7_geneorder.py</td><td valign="top">Python script to rearrange the data in the cluster_table cells to put genomic gene order up front.  This allows fast sorting of genes in a genome by their location and makes it easy to spot genomic islands. </td><tr>
<tr><td valign="top">8_COGadd.py</td><td valign="top">Python script to incorporate deepnog COG annotations and metatdata into the cluster_table.tab/xlsx files. </td><tr>
<tr><td valign="top">9_KEGGadd.py</td><td valign="top">Python script to incorporate KEGG annotations from their program 'kofamscan' into the cluster_table.tab/xlsx files. </td><tr>
</table>



# The files you need to input

##NCBI 'cds_from_genomic' files

##NCBI's .json file with metadata

# The most important output files

## cluster_table.tab and/or cluster_table.xlsx