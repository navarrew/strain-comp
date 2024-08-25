# strain-comp
Scripts to compare bacterial genomes

# Getting genomic data
In any command line environment where ncbi-datasets is installed you can type commands directly into the terminal to download genomic data and its associated metadata.

$ conda activate ncbi (or use the bi environment)


Downloading genomes using a taxonomic name

If you want all the Gardnerella RefSeq annotated ORFs (cds) from genomes without MAGs and atypical genomes type:

$ datasets download genome taxon gardnerella --include cds --assembly-source RefSeq --exclude-atypical --annotated --exclude-multi-isolate --mag exclude --filename gardnerella_dataset.zip


Download a genome using an accession number

If you know the accession of a particular genome and want both the RefSeq ORFs and the full gbff file.

$ datasets download genome accession GCF_01234567.1 --include cds,gbff --assembly-source RefSeq --filename strainXYZ.zip

Or if you want a few genomes, separate the accession numbers with commas

$ datasets download genome accession GCF_01234567.1,GCF_022662295.1 --include cds,gbff --assembly-source RefSeq --filename ncbi_dataset.zip

