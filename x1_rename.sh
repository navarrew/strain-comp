#!/bin/bash

shopt -s nullglob #this suppresses error messages for 'file not found', etc.

echo "This will rename all files in the ncbi_datset folder with their GCF number"

#lets get rid of the README.md file and move the assembly report to the top level directory
rm README.md
rm ncbi_dataset/data/dataset_catalog.json
mv ncbi_dataset/data/assembly_data_report.jsonl .
wait

#lets go through each folder and rename the various files to their assembly IDs.
for file in ncbi_dataset/data/*/cds_from_genomic.fna
do
directory_name=$(dirname $file)
accession=$(basename $directory_name)
mv "${file}" "${directory_name}/${accession}.fna"
done
wait

for file in ncbi_dataset/data/*/genomic.gbff
do
directory_name=$(dirname $file)
accession=$(basename $directory_name)
mv "${file}" "${directory_name}/${accession}.gbff"
done
wait

for file in ncbi_dataset/data/*/genomic.gff
do
directory_name=$(dirname $file)
accession=$(basename $directory_name)
mv "${file}" "${directory_name}/${accession}.gff"
done
wait

for file in ncbi_dataset/data/*/genomic.gtf
do
directory_name=$(dirname $file)
accession=$(basename $directory_name)
mv "${file}" "${directory_name}/${accession}.gtf"
done
wait

for file in ncbi_dataset/data/*/protein.faa
do
directory_name=$(dirname $file)
accession=$(basename $directory_name)
mv "${file}" "${directory_name}/${accession}.faa"
done
wait

for file in ncbi_dataset/data/*/GC*_genomic.fna
do
directory_name=$(dirname $file)
accession=$(basename $directory_name)
mv "${file}" "${directory_name}/${accession}_genomic.fna"
done
wait

for file in ncbi_dataset/data/*/sequence_report.jsonl
do
directory_name=$(dirname $file)
accession=$(basename $directory_name)
mv "${file}" "${directory_name}/${accession}_report.jsonl"
done
wait

mkdir data
mkdir data/ncbi

#move the data folder to the top level ('.') directory.
mv ncbi_dataset/data/ data/ncbi/
mv data/ncbi/cds_fasta data/ncbi/data
#remname the cds_fasta folder 


#delete the ncbi_dataset folder which should be empty anyway at this point.
rm -rf ncbi_dataset
wait



#now we make the master table for the strains using NCBIs dataformat program and the assembly_data_report.jsonl file.
#first we'll make a header for the master table
echo "Accession	Species	Strain	BioProject	BioSample	Level" > data/ncbi/master_table.tab
#then we append to the master table all the data we want from the assembly_data_report.jsonl file
mv assembly_data_report.jsonl data/ncbi/
dataformat tsv genome --inputfile data/ncbi/assembly_data_report.jsonl --fields accession,ani-submitted-species,assminfo-biosample-strain,assminfo-bioproject,assminfo-biosample-accession,assminfo-level --elide-header >> data/ncbi/master_table.tab
