# Ploidetect - Estimation of tumour purity and aneuploidy from whole-genome sequence data
By Luka Culibrk at Canada's Michael Smith Genome Sciences Centre

## Installation

Ploidetect is an R package. Install Ploidetect by running ```devtools::install_github("lculibrk/Ploidetect-package")```

## Input data
The input to Ploidetect is a tab-delimited file containing the following columns, where each column represents a genomic bin (similar to a bed file):
* Tumour read counts in the genomic bin
* Germline read counts in the genomic bin
* Window id (formatted as chromosome#_start where chromosome# is the chromosome (without "chr") and start is the start position of the genomic bin)
* Window size (Size of the genomic bin in bp)
* Tumour allele frequencies of heterozygous germline SNPs
* GC content of the genomic bin

## Running Ploidetect
Ploidetect is intended to be as simple as possible to run. After loading Ploidetect with ```library(Ploidetect)```, simply load the input data file, for instance:

```dat <- read.table("example.txt", sep = "\t", header = T, stringsAsFactors = F```)

and then run Ploidetect:

output <- ```ploidetect(dat)```

It's really that simple!

...To be updated with detailed instructions on interpretation of results, manually forcing models, etc...
