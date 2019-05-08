# Ploidetect - Estimation of tumour purity and aneuploidy from whole-genome sequence data
By Luka Culibrk at Canada's Michael Smith Genome Sciences Centre

## Installation

Ploidetect is an R package. Install Ploidetect by running ```devtools::install_github("lculibrk/Ploidetect-package")``` in R.

## Input data
The input to Ploidetect is a tab-delimited .bed file containing the following columns in the following order, where each column represents a genomic bin:
* chromosome without the "chr" ie. 1, 2, X
* Start coordinate (0-based)
* End coordinate 
* Raw Tumour read depth (e.g. calculated by bedtools coverage or a similar tool such as hts-nim-tools) 
* Raw normal read depth)
* Heterozygous germline SNP allele frequencies in tumour (e.g. calculated from an mpileup or similar)
* GC content of the genomic bin (e.g. calculated using bedtools nuc)

Column names are optional

## Running Ploidetect
Ploidetect is intended to be as simple as possible to run. After loading Ploidetect with ```library(Ploidetect)```, simply load the input data file, for instance:

```dat <- read_ploidetect("example.txt")```

and then run Ploidetect:

output <- ```ploidetect(dat)```

The outputs can be accessed from the output list:

```output$TC_calls``` will give you the purity calls

```output$plots``` will print the model plots, in the same order as ```output$TC_calls```

```output$CN_calls``` will give you the segmented copy number calls

If you believe Ploidetect made an error in automatic model estimation,you can supply it with the information of the model you believe to be correct using the ```cndiff```, ```lowest```, and ```comp``` parameters for ```ploidetect()```.Please see the [demo](demo/Demo.md) (Currently under construction) for more information on this. 
