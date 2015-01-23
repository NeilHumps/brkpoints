# Breakage Point Finder

##Overview
This is a "real" Readme file.

## Installation
1. Make a virtual environment (checkout virtualenvwrapper)
mkvirtualenv brkpoints

2. ```python setup.py develop```

### Dependencies
The following dependencies are required:
1. Biopython
Is there a specific version?
Developed using Biopython version 1.64

2. Blastn
Describe where to get this and how it should be set up
?```sudo apt-get install ncbi-blast+```
Download and install the blast commandline software using the NCBI FTP service, details:
http://www.ncbi.nlm.nih.gov/books/NBK1763/
by default, the software will install in /usr/local/ncbi/blast/
This was developed using blast version 2.2.28

3. A local copy of the Human genome
Obtain the human_genomic Blast database from:
ftp://ftp.ncbi.nlm.nih.gov/blast/db
human_genomic.*tar.gz : Human RefSeq (NC_######) chromosome records with 
                        gap adjusted concatenated NT_ contigs
in the ncbi/blast/ directory, make a directory called 'db'
download them into the 'db' directory and unzip

## Usage

Document how to use the script here.
Fragile_finder is a python module to deduce the locations of fragile site/region origins of translocation events. These translocation events are obtained from http://www.unav.es/genetica/TICdb/ as a tab-delimited file.
Both sides of the translocation junction are mapped to the human genome using BLAST, to create a list of fragile sites. Loci that appear to originate from the same fragile region, i.e. in close proximity on the genome, are merged into a single list entry, therefore removing duplicate entries. The list of fragile sites is output as a csv file, which documents fragile site name, details of original translocation event, suffixed with _a (represents left-hand-side of translocation) or _b (represents right-hand-side of translocation), chromosome number, start and end coordinates and strand (if strand is different in duplicates, strand=1).

Fragile_finder_test is used to determine if the fragile sites identified colocalise with the genes specified in the translocation data. 100bp of the gene sequence was obtained from GenBank using the HGNC gene name. This sequence is subjected to the same BLAST search as the fragile sites to obtain a genomic coordinate to compare to the fragile site list. The output states the % of sites that colocalise to the named gene, and a csv file of sites that dont colocalise, with comments. Colocalising sites are also listed below non-colocalising sites.

Run the test suite with

```py.test tests```
