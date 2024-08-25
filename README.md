# strain-comp
Scripts to compare bacterial genomes

## Getting genomic data using the NCBI datasets command line interface
In any command line environment where ncbi-datasets is installed you can type commands directly into the terminal to download genomic data and its associated metadata.

### Downloading genomes using a taxonomic name

If you want all the Gardnerella RefSeq annotated ORFs (cds) from genomes without MAGs and atypical genomes type:

  **$ datasets download genome taxon gardnerella --include cds --assembly-source RefSeq --exclude-atypical --annotated --exclude-multi-isolate --mag exclude --filename ncbi_dataset.zip**


### Download genomes using an accession number

If you know the accessions of a set of genome and want both the RefSeq ORFs and the full gbff file.

  **$ datasets download genome accession GCF_01234567.1,GCF_022662295.1 **--include cds,gbff** --assembly-source RefSeq --filename ncbi_dataset.zip**

### Download several genomes from a file of accession numbers ###

You can put a lot of accession numbers into a single text file (one accession per line) and feed them into the download program using the --inputfile flag.  The file does not need to be named accessions.txt.

$ datasets download genome accession --inputfile accessions.txt --include cds,gbff --assembly-source RefSeq --filename ncbi_dataset.zip
