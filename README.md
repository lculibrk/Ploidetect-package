# Ploidetect - Estimation of tumour purity and aneuploidy from whole-genome sequence data
By Luka Culibrk at Canada's Michael Smith Genome Sciences Centre

## Installation

Ploidetect is an R package. Install Ploidetect by running ```devtools::install_github("lculibrk/Ploidetect-package")``` in R.

## Input data
The input to Ploidetect is a tab-delimited file containing the following columns, where each column represents a genomic bin (similar to a bed file):
* Tumour read counts in the genomic bin
* Germline read counts in the genomic bin
* Tumour allele frequencies of heterozygous germline SNPs
* Window id (formatted as chromosome#_start where chromosome# is the chromosome (without "chr") and start is the start position of the genomic bin)
* Window size (Size of the genomic bin in bp)
* GC content of the genomic bin

Ploidetect requires you to supply the column indices for each of these variables, and by default they are assumed to be in the above order

## Running Ploidetect
Ploidetect is intended to be as simple as possible to run. After loading Ploidetect with ```library(Ploidetect)```, simply load the input data file, for instance:

```dat <- read.table("example.txt", sep = "\t", header = T, stringsAsFactors = F)```

and then run Ploidetect:

output <- ```ploidetect(dat)```

The outputs can be accessed from the output list:

```output$TC_calls``` will give you the purity calls

```output$plots``` will print the model plots, in the same order as ```output$TC_calls```

```output$CN_calls``` will give you the segmented copy number calls

If you believe Ploidetect made an error in automatic model estimation,you can supply it with the information of the model you believe to be correct using the ```cndiff```, ```lowest```, and ```comp``` parameters for ```ploidetect()```.Please see the [demo](demo/Demo.md) (Currently under construction) for more information on this. 
