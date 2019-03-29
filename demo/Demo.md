Ploidetect Example Case
================
Luka Culibrk
2/11/2019

Data import
-----------

First we import the data, formatted as discussed in the README.md. This data is a metastatic tumor biopsy sequenced as part of the Personalized Oncogenomics Project.

``` r
dat <- read.table("examplePOG.txt", stringsAsFactors = F, header = T)
str(dat)
```

    ## 'data.frame':    126254 obs. of  6 variables:
    ##  $ V4    : int  14904 13843 14757 14472 14044 13728 14379 14222 14155 14054 ...
    ##  $ V5    : int  1000040 1000014 1000008 1000032 1000016 1000018 1000007 1000022 1000016 1000002 ...
    ##  $ V6    : chr  "." "." "." "." ...
    ##  $ window: chr  "1_0" "1_126352" "1_267214" "1_570656" ...
    ##  $ size  : int  126352 140862 303442 142132 28522 30467 26608 23202 23910 28241 ...
    ##  $ V7    : num  0.379 0.282 0.295 0.415 0.413 ...

Running Ploidetect
------------------

Now we run Ploidetect by calling `ploidetect()`, and supplying the indices for the columns in our data.

``` r
result <- ploidetect(all_data = dat, normal = 2, tumour = 1, avg_allele_freq = 3, window_id = 4, window_size = 5, GC = 6, verbose = F, CNA_call = F)
```

Let's get a look at how the results look.

``` r
purity_calls <- result$TC_calls
plots <- result$plots
knitr::kable(purity_calls)
```

|     |  reads\_per\_copy|  zero\_copy\_depth|  Ploidy|  tumour\_purity|  lowest\_peak\_CN|  maf\_error|  CN\_diff|  Comparator|  model\_error|  avg\_ploidy|  unmatchedpct|
|-----|-----------------:|------------------:|-------:|---------------:|-----------------:|-----------:|---------:|-----------:|-------------:|------------:|-------------:|
| 1   |          10498.76|         32943.4237|       2|       0.3892687|                 1|   0.0318992|         2|           3|      3.626962|         2.87|     0.0000000|
| 2   |           7243.14|         33012.1694|       3|       0.3049843|                 2|   0.0313687|         4|           5|    106.748614|         4.15|     0.0000000|
| 21  |          18581.03|          -286.1983|       3|       1.0077611|                 2|   0.0660855|         1|           3|    433.156900|         3.33|     0.5788985|

``` r
plots[[3]]
```

![](Demo_files/figure-markdown_github/unnamed-chunk-3-1.png)

TC\_calls object
----------------

ploidetect() returns a list of three objects; TC\_calls, or "tumor content calls" describes the models that Ploidetect has used to estimate tumour purity using a variety of different assumptions. It is ordered by the strength of the model, so the first model is most likely to describe the true tumour purity. In this case the purity is predicted to be 39%.

There are 11 columns in the TC\_calls object. The columns correspond to the following values:

-   reads\_per\_copy: The number of reads in a genomic bin that are expected to come from a single copy of the genome

-   zero\_copy\_depth: The number of reads in a genomic bin which is homozygously deleted in the tumor sample. This is the number of reads coming from the contaminating normal.

-   Ploidy: the estimated most common copy number state in the genome.

-   tumour\_purity: the tumour purity

-   lowest\_peak\_CN: used in Ploidetect's modeling. Corresponds to the copy number of the lowest well-represented copy number state in the genome

-   maf\_error: The median difference between expected and observed germline-heterozygous SNP allele frequencies

-   CN\_diff: Used in Ploidetect's modeling. See below for more information.

-   Comparator: Used in Ploidetect's modeling. See below for more information.

-   Model\_error: The error computed by Ploidetect for this model. Used to rank models.

-   Avg\_ploidy: The mean copy number of the genome.

-   Unmatchedpct: An estimate of the proportion of the genome that was excluded from the model fit. See below for more information.

A brief lesson in NGS copy number variation data
------------------------------------------------

Copy number variation is a common type of mutation in most cancers. To explain the methodology, let's use a simple toy genome to go over a few concepts that Ploidetect uses. First, we generate a small toy genome with a decent amount of chromosomal instability.

``` r
set.seed(42069)
expand.segments <- function(copynumbers){
  out_cns <- c()
  out_segs <- c()
  for(segment in 1:length(copynumbers)){
    new_size <- round(rnorm(n = 1, mean = 20, sd = 4), digits = 0)
    out_cns <- c(out_cns, rep(copynumbers[segment], times = max(1, new_size)))
    out_segs <- c(out_segs, rep(segment, times = max(1, new_size)))
  }
  out <- data.frame("copynumber" = out_cns, "segment" = out_segs)
  return(out)
}

genome_positions <- 1:10
genome_copynumber <- round(rnorm(n = length(genome_positions), mean = 2, sd = 1), digits = 0)
genome_obj <- expand.segments(genome_copynumber)
genome_obj$start <- 1:nrow(genome_obj)
genome_obj$end <- genome_obj$start + 1

ggplot(genome_obj, aes(x = start, y = copynumber)) + geom_point() + theme_bw()
```

![](Demo_files/figure-markdown_github/unnamed-chunk-4-1.png)

Okay. So we have a genome. Now let's simulate depth based on the copy numbers we've generated. We'll have the counts follow a normal distribution for simplicity, and assume that the depth is 40x (ie. two-copy regions get 40 aligned reads, one-copy regions get 20, and so forth)

``` r
set.seed(42069)
genome_obj$counts <- NA
for(segment in genome_obj$segment){
  current_segment <- genome_obj[genome_obj$segment == segment,]
  genome_obj$counts[genome_obj$segment == segment] <- rnorm(n = nrow(current_segment), mean = 20 * current_segment$copynumber[1])
}
ggplot(genome_obj, aes(x = start, y = counts, color = factor(copynumber))) + geom_point() + theme_bw() + scale_color_viridis(discrete = T)
```

![](Demo_files/figure-markdown_github/unnamed-chunk-5-1.png)

In this example, the take-home message is that the difference in depth between a one-copy region and a two-copy region is the same as the difference in depth between a two and three-copy region - in this case, the difference is 20. Here's another way to look at the copy number landscape of this genome:

``` r
ggplot(genome_obj, aes(x = counts)) + geom_density() + scale_x_continuous(limits = c(0, 100)) + theme_bw()
```

![](Demo_files/figure-markdown_github/unnamed-chunk-6-1.png)

The density plot communicates the relative quantity and approximate depth of the copy number variants in this genome. In this example genome, we have assumed that there is zero normal contamination. In a real tumour biopsy, this is almost never the case.

Consider that the total reads that align to a given locus can be represented as such:

![\\sum reads = reads\_t + reads\_n](https://latex.codecogs.com/png.latex?%5Csum%20reads%20%3D%20reads_t%20%2B%20reads_n "\sum reads = reads_t + reads_n")

Where ![reads\_t](https://latex.codecogs.com/png.latex?reads_t "reads_t") is the number tumour reads and ![reads\_n](https://latex.codecogs.com/png.latex?reads_n "reads_n") is the amount of normal reads.

In the case of 50% tumour purity, the read depth of a locus with equal copy number in both tumour and normal should be equal:

![\\sum reads = reads\_t + reads\_n = 2reads\_t = 2reads\_n](https://latex.codecogs.com/png.latex?%5Csum%20reads%20%3D%20reads_t%20%2B%20reads_n%20%3D%202reads_t%20%3D%202reads_n "\sum reads = reads_t + reads_n = 2reads_t = 2reads_n")

Or where the copy number is inequal, it is a weighted average (![C\_n](https://latex.codecogs.com/png.latex?C_n "C_n") for normal copy number and ![C\_t](https://latex.codecogs.com/png.latex?C_t "C_t") for tumour copy number):

![\\sum reads = \\frac{C\_treads\_t + C\_nreads\_n}{2}](https://latex.codecogs.com/png.latex?%5Csum%20reads%20%3D%20%5Cfrac%7BC_treads_t%20%2B%20C_nreads_n%7D%7B2%7D "\sum reads = \frac{C_treads_t + C_nreads_n}{2}")

Since (nearly) the entire normal genome should be diploid, we can simplify the above by removing ![C\_n](https://latex.codecogs.com/png.latex?C_n "C_n"):

![\\sum reads = \\frac{C\_treads\_t + 2reads\_n}{2} = \\frac{C\_treads\_t}{2} + reads\_n](https://latex.codecogs.com/png.latex?%5Csum%20reads%20%3D%20%5Cfrac%7BC_treads_t%20%2B%202reads_n%7D%7B2%7D%20%3D%20%5Cfrac%7BC_treads_t%7D%7B2%7D%20%2B%20reads_n "\sum reads = \frac{C_treads_t + 2reads_n}{2} = \frac{C_treads_t}{2} + reads_n")

The conclusion of this (long-winded) string of equations is that for any given locus, the number of reads that originate from contaminating normal is a constant, ![reads\_n](https://latex.codecogs.com/png.latex?reads_n "reads_n"). Since ![reads\_n](https://latex.codecogs.com/png.latex?reads_n "reads_n") corresponds to the depth of a two-copy region, we can get the degree of normal contamination (![\\alpha](https://latex.codecogs.com/png.latex?%5Calpha "\alpha")) from the following equation:

![\\alpha = \\frac{reads\_n}{\\sum reads}](https://latex.codecogs.com/png.latex?%5Calpha%20%3D%20%5Cfrac%7Breads_n%7D%7B%5Csum%20reads%7D "\alpha = \frac{reads_n}{\sum reads}")

Where ![\\sum reads](https://latex.codecogs.com/png.latex?%5Csum%20reads "\sum reads") was measured at copy number two.

Going back to the toy genome, let's add some contaminating counts - we'll aim for 50% tumour purity here.

``` r
set.seed(42069)
genome_obj$counts <- genome_obj$counts/2 + rnorm(nrow(genome_obj), mean = 20)
ggplot(genome_obj, aes(x = start, y = counts, color = factor(copynumber))) + geom_point() + theme_bw() + scale_color_viridis(discrete = T)
```

![](Demo_files/figure-markdown_github/unnamed-chunk-7-1.png)

``` r
ggplot(genome_obj, aes(x = counts)) + geom_density() + scale_x_continuous(limits = c(0, 100), breaks = seq(from = 0, to = 100, by = 10)) + theme_bw()
```

![](Demo_files/figure-markdown_github/unnamed-chunk-7-2.png)

First, how do we determine ![reads\_n](https://latex.codecogs.com/png.latex?reads_n "reads_n")?

Going back to ![\\sum reads = reads\_t + reads\_n](https://latex.codecogs.com/png.latex?%5Csum%20reads%20%3D%20reads_t%20%2B%20reads_n "\sum reads = reads_t + reads_n"), we just have to find where ![reads\_t = 0](https://latex.codecogs.com/png.latex?reads_t%20%3D%200 "reads_t = 0"), which occurs at loci that are homozygously deleted. In our toy case (and in many tumour biopsies), there are no homozygously deleted regions. However, we can predict where it will be. From the above plot, we can see that the difference in depth between each copy number is about 10 counts. This information can be used to determine the depth of a zero-copy region using regression:

``` r
fit <- lm(formula = counts~copynumber, data = genome_obj)
predict(object = fit, data.frame("copynumber" = 0))
```

    ##        1 
    ## 20.04273

Now ![\\sum reads](https://latex.codecogs.com/png.latex?%5Csum%20reads "\sum reads") for two-copy can be challenging to find, since we need to identify which loci are two-copy. In this example, we have it in `genome_obj`

``` r
genome_obj %>% group_by(copynumber) %>% dplyr::summarise("mean_cov" = mean(counts))
```

    ## # A tibble: 4 x 2
    ##   copynumber mean_cov
    ##        <dbl>    <dbl>
    ## 1          1     29.8
    ## 2          2     40.1
    ## 3          3     50.2
    ## 4          4     59.7

Now we just plug these numbers in:

``` r
20/40.1
```

    ## [1] 0.4987531

Tumour purity is 0.4998187 for our toy example.

The above is the methodology that Ploidetect uses to determine tumour purity. Let's do the estimate for the cancer case ourselves.

``` r
result$plots[[3]]
```

![](Demo_files/figure-markdown_github/unnamed-chunk-11-1.png)

CN=2 occurs at approximately 55000 depth, and CN=1 occurs at about 45000. CN=0 must occur at about 35000. As an aside, Ploidetect knows the copy number state of the peaks in the histogram by fitting the SNP allele frequencies indicted by "MAF = x" in the plot.

``` r
1 - 35000/55000
```

    ## [1] 0.3636364

Eyeballing it gives us 36% purity, and Ploidetect predicted 39%.

However, determining the depth of a homozygous deletion may be challenging the case of a messy genome. Notably, in cases of subclonal copy number variation, not all peaks in the density histogram correspond to integer copy number states. In these cases, it is challenging to estimate the depth difference between different integer copy number states. To illustrate this, here's another metastatic tumour, however this one contains subclonal copy number variation, and we'll go through Ploidetect's process to explain how it manages to determine the purity of this biopsy.

``` r
clonalcase <- read.table("clonal_example.txt", sep = "\t", header = T, stringsAsFactors = F)
str(clonalcase)
```

    ## 'data.frame':    152304 obs. of  6 variables:
    ##  $ V4    : int  14334 16122 26613 16676 16590 16570 16715 16624 16854 17290 ...
    ##  $ V5    : int  1000002 1000034 1000005 1000021 1000004 1000017 1000019 1000008 1000020 1000011 ...
    ##  $ V6    : chr  "." "." "." "." ...
    ##  $ window: chr  "1_0" "1_126577" "1_536971" "1_600543" ...
    ##  $ size  : int  126577 410394 63572 113805 23282 25173 23165 22789 18340 25249 ...
    ##  $ V7    : num  0.379 0.275 0.439 0.426 0.401 ...

Here we're going to go through Ploidetect's internals to demonstrate what it's doing, step by step. First let's set the variables that would otherwise be handled by the parameters of ploidetect()

``` r
all_data <- clonalcase
normal = 2
tumour = 1
avg_allele_freq =3
window_id = 4
window_size = 5
GC = 6
verbose = F
rerun = F
lowest = NA
cndiff = NA
comp = NA
nomaf = F
```

Note that you won't be able to call these functions yourself as they aren't exported by the namespace of the package.

First, we call ploidetect\_preprocess, which performs basic preprocessing of the data. First, Ploidetect converts combines sequential bins into ~100kb sized new bins for more accurate modeling. It corrects for lingering GC biases using a loess fit and also corrects the read depth using the normal read depth, since using constant depth bins can fail at the terminal end of the chromsome due to being truncated by the end of the chromosome. We get three objects, one of which is important for our sake.

``` r
preprocess_output <- ploidetect_preprocess(all_data = all_data)
str(preprocess_output)
```

    ## List of 4
    ##  $ x           :'data.frame':    149837 obs. of  8 variables:
    ##   ..$ y_raw          : num [1:149837] 26614 16590 16570 16715 16624 ...
    ##   ..$ x_raw          : int [1:149837] 1000005 1000004 1000017 1000019 1000008 1000020 1000011 1000032 1000021 1000006 ...
    ##   ..$ maf            : num [1:149837] NaN NaN 0.385 0.514 0.425 ...
    ##   ..$ window         : chr [1:149837] "1_536971" "1_714348" "1_737630" "1_762803" ...
    ##   ..$ size           : int [1:149837] 63572 23282 25173 23165 22789 18340 25249 29427 27180 27812 ...
    ##   ..$ GC             : num [1:149837] 0.439 0.401 0.428 0.465 0.46 ...
    ##   ..$ residual       : num [1:149837] 8417 -1429 -1573 -1595 -1668 ...
    ##   ..$ normalized_size: int [1:149837] 63572 23282 25173 23165 22789 18340 25249 29427 27180 27812 ...
    ##  $ maxpeak     : num 18840
    ##  $ highoutliers:'data.frame':    479 obs. of  8 variables:
    ##   ..$ merge : int [1:479] 6216 6217 6218 6219 6220 6221 6222 6223 6224 6225 ...
    ##   ..$ chr   : chr [1:479] "1" "1" "1" "1" ...
    ##   ..$ tumour: num [1:479] 45520 46647 47303 47287 46170 ...
    ##   ..$ normal: int [1:479] 1008473 1037335 1033976 1028034 1045385 1074295 1025875 1004450 1111839 1028179 ...
    ##   ..$ maf   : num [1:479] NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN ...
    ##   ..$ wind  : chr [1:479] "1_121484470" "1_121484480" "1_121484490" "1_121484500" ...
    ##   ..$ size  : int [1:479] 10 10 10 10 10 10 9 8 8 7 ...
    ##   ..$ gc    : num [1:479] 0.6 0.2 0.3 0.4 0.4 ...
    ##  $ merged      : num 5

``` r
filtered <- preprocess_output$x
maxpeak <- preprocess_output$maxpeak
highoutliers <- preprocess_output$highoutliers
```

filtered contains the corrected and filtered read depths for our data. Let's compare before and after (and some filtering to visualize the difference more easily).

``` r
filtered %>% filter (y_raw < 200000) %>% ggplot(aes(x = size)) + geom_point(aes(y = y_raw), alpha = 0.05, size = 0.1)
```

![](Demo_files/figure-markdown_github/unnamed-chunk-17-1.png)

``` r
filtered %>% filter (y_raw < 200000) %>% ggplot(aes(x = size)) + geom_point(aes(y = residual), alpha = 0.05, size = 0.1)
```

![](Demo_files/figure-markdown_github/unnamed-chunk-17-2.png)

Ploidetect uses a kernel density estimation to find peaks in the depth histogram. We use a simple heuristic to get a decent bandwidth for peak calling.

``` r
bw = maxpeak/40
```

It's not much. In some cases, this step can be extremely helpful in reducing the mapping biases in the data. Ploidetect does a quick check to see if this step manages to actually reduce the variability. The idea behind this is that the worse data will have points far more spread out, resuting in a higher average normalized density

``` r
## Compute density for mapping-corrected data
den <- density(filtered$residual, n = nrow(filtered), bw = bw)
## Normalize the density to 0->1 range
den$y <- (den$y - min(den$y))/(max(den$y) - min(den$y))
## Sanity check - see if ploidetect_preprocess resulted in a worse separation of peaks than the initial data
den2 <- density(filtered$y_raw, n = nrow(filtered), bw = bw)
den2$y <- (den2$y - min(den2$y))/(max(den2$y) - min(den2$y))
if(mean(den2$y) < mean(den$y)){
  filtered$residual <- filtered$y_raw - maxpeak
}
```

If you were wondering, highoutliers contains exceptionally high depth outliers that are either caused by sequencing error or other causes, and we exclude these bins to limit their interference on the modeling. maxpeak is the most common read depth that we observe in the data.

Next, we need to call peaks in the density histogram. The histogram, for reference:

``` r
plot(density(filtered$residual, bw = bw, n = nrow(filtered$residual)))
```

![](Demo_files/figure-markdown_github/unnamed-chunk-20-1.png)

Ploidetect calls peaks using a simple derivative approach. It takes the second derivative of the kernel density estimation and detects the point where the sign goes negative.

``` r
library(kedd)
plot(dkde(filtered$residual, deriv.order = 2, h = bw))
```

![](Demo_files/figure-markdown_github/unnamed-chunk-21-1.png)

Of course, since we used a simple heuristic for bandwidth, we should probably throw another heuristic in for good measure. What if we detect only one peak in the data? Either there is no copy number variation, or Ploidetect needs to reduce the bandwidth.

``` r
allPeaks <- peakcaller(filtered[findInterval(filtered$residual, vec = quantile(filtered$residual, probs = c(0, 0.999))) == 1,], bw)
knitr::kable(allPeaks)
```

|       start|         end|     height|         pos|     trough|  ratiotrough|  troughdiff|       diff|
|-----------:|-----------:|----------:|-----------:|----------:|------------:|-----------:|----------:|
|    237.2232|   1405.7505|  1.0000000|    772.4014|  0.3646409|    0.3646409|   0.6353591|      0.000|
|  -3502.0643|  -2391.9633|  0.5719458|  -2952.7582|  0.2293659|    0.4010272|   0.3425800|   3725.160|
|  -1398.7151|   -756.0251|  0.3646409|   -756.1020|  0.3646409|    1.0000000|   0.0000000|   1528.503|
|   3918.0842|   5028.1852|  0.1345760|   3976.4056|  0.1344579|    0.9991227|   0.0001181|   3204.004|
|  11455.0854|  13091.0236|  0.0497407|  12198.2201|  0.0061578|    0.1237973|   0.0435829|  11425.819|
|   8066.3562|   8650.6198|  0.0061578|   8650.5778|  0.0061578|    1.0000000|   0.0000000|   7878.176|
|  -7358.2044|  -6306.5298|  0.0059025|  -6833.0282|  0.0009978|    0.1690457|   0.0049047|   7605.430|

So we've called more than one peak, we can continue. If Ploidetect called only one then it would try again by halving the bandwidth and giving a warning

``` r
if(nrow(allPeaks) == 1){
  warning("Zero peaks detected. Attempting peak calling with lower bandwidth")
  bw = bw/2
  allPeaks <- peakcaller(filtered[findInterval(filtered$residual, vec = quantile(filtered$residual, probs = c(0, 0.999))) == 1,], bw)
}
```

Next we do a bit of housekeeping and center the positions and read depths of the peaks around zero

``` r
filtered$residual <- filtered$residual - allPeaks$pos[1]
allPeaks$start <- allPeaks$start - allPeaks$pos[1]
allPeaks$end <- allPeaks$end - allPeaks$pos[1]
allPeaks$pos <- allPeaks$pos - allPeaks$pos[1]
```

next we need to do a bit of work with the allele frequency data. Specifically, we want to determine the SNP allele frequency that best represents each peak. This is done by ploidetect\_processallpeaks. We also map the data to their respective peaks. This doesn't map them all, because peak boundaries end up being fairly conservative due to the kernel density approach.

``` r
output <- ploidetect_processallpeaks(filtered, allPeaks)
filtered <- output$filtered
allPeaks <- output$allPeaks
filtered %>% ggplot(aes(x = size, y = residual, color = peak)) + geom_point(size = 0.1, alpha = 0.1) + scale_color_viridis(discrete = F)
```

![](Demo_files/figure-markdown_github/unnamed-chunk-25-1.png)

This is where Ploidetect presents you with its first plot:

``` r
filteredforplot <- filtered[findInterval(filtered$residual, vec = quantile(filtered$residual, probs = c(0.01, 0.999))) <= 1,]
filteredforplot$residual <- filteredforplot$residual + maxpeak
plot <- ggplot(data = filteredforplot, mapping = aes_string(x = "size", y = "residual", color = "mafflipped")) + geom_point(size = 0.1, alpha = 0.1) +
  #xlab("Window size") + 
  xlab("Window size of constant-coverage bins") + 
  #ylab("Reads mapping to bins in Somatic") + 
  ylab("Normalized somatic read counts") + 
  #ggtitle(paste0("Coverage vs normalized bin size (Surrogate for Mappability + GC bias)")) + 
  ggtitle(paste0("Coverage Plot for Filtered and Normalised Data")) + 
  scale_colour_viridis(option = "plasma", name = "Major Allele\n Frequency") +
  #scale_x_continuous(limits = quantile(filtered$size, probs = c(0.05, 0.99))) +
  theme_bw(base_size = 12)
plot
```

![](Demo_files/figure-markdown_github/unnamed-chunk-26-1.png)

And now it generates the density plot, which also illustrates the peaks that have been called in the density histogram

``` r
den <- density(filteredforplot$residual, bw = bw)
dendf <- data.frame(x = den$x, y = den$y)
# Normalize the density to 0->1 range
dendf$y <- (dendf$y - min(dendf$y))/(max(dendf$y) - min(dendf$y))
plot <- ggplot(data = dendf, mapping = aes(x = x, y = y)) + geom_line() + 
  xlab("Normalized somatic read counts") + 
  ylab("Normalized Density") + 
  ggtitle(paste0("Kernel Density Estimate of Count Data")) + 
  geom_vline(data = allPeaks, aes(xintercept = pos + maxpeak), linetype = 2) + 
  geom_text(data = allPeaks, aes(x = pos + maxpeak, y = allPeaks$height + 0.05, label = paste0("MAF = ", round(mainmaf, digits = 3)))) +
  theme_bw()
plot
```

![](Demo_files/figure-markdown_github/unnamed-chunk-27-1.png)

I've already told you that this case has subclonal copy number variation. I'm going to spoil it right away and tell you that the peaks immediately adjacent to the tallest one are the result of subclonal copy number events. The reasons for this will become more clear in a bit.

If you recall with the toy case, we used the plot information to estimate the number of reads that align to a single copy of the tumour genome, which we will refer to as the number of reads per copy. This is done by Ploidetect by selecting a set of peaks that it reasonably believes to correspond to integer copy number variants, and using those peaks to determine the reads per copy. Here's the table of peaks again:

``` r
allPeaks %>% arrange(pos) %>% knitr::kable()
```

|       start|         end|     height|        pos|     trough|  ratiotrough|  troughdiff|       diff|  npeak|    mainmaf|
|-----------:|-----------:|----------:|----------:|----------:|------------:|-----------:|----------:|------:|----------:|
|  -8130.6058|  -7078.9312|  0.0059025|  -7605.430|  0.0009978|    0.1690457|   0.0049047|   7605.430|      1|  0.5191255|
|  -4274.4656|  -3164.3647|  0.5719458|  -3725.160|  0.2293659|    0.4010272|   0.3425800|   3725.160|      2|  0.6257158|
|  -2171.1164|  -1528.4264|  0.3646409|  -1528.503|  0.3646409|    1.0000000|   0.0000000|   1528.503|      3|  0.5583460|
|   -535.1782|    633.3491|  1.0000000|      0.000|  0.3646409|    0.3646409|   0.6353591|      0.000|      4|  0.5146988|
|   3145.6829|   4255.7838|  0.1345760|   3204.004|  0.1344579|    0.9991227|   0.0001181|   3204.004|      5|  0.5837145|
|   7293.9548|   7878.2185|  0.0061578|   7878.176|  0.0061578|    1.0000000|   0.0000000|   7878.176|      6|  0.7259628|
|  10682.6840|  12318.6223|  0.0497407|  11425.819|  0.0061578|    0.1237973|   0.0435829|  11425.819|      7|  0.6796519|

So one example model might select peaks 1, 2, 4, 6, 7, and 9 (peaks are numbered by depth). Then, we can compute the tumour purity using this information. To generate these models, Ploidetect needs to first perform a rough estimation step. Iteratively, Ploidetect selects the tallest peak (peak 4 in this case) and another peak in the dataset to serve as a comparator. For this case, let's select peak 2.

``` r
allPeaks_plot <- allPeaks
allPeaks_plot$selected <- F
allPeaks_plot$selected[allPeaks$npeak %in% c(2, 4)] <- T
plot <- ggplot(data = dendf, mapping = aes(x = x, y = y)) + geom_line() + 
  xlab("Normalized somatic read counts") + 
  ylab("Normalized Density") + 
  ggtitle(paste0("Peaks 1 and 2")) + 
  geom_vline(data = allPeaks_plot, aes(xintercept = pos + maxpeak, color = selected), linetype = 2) + 
  geom_text(data = allPeaks_plot, aes(x = pos + maxpeak, y = allPeaks$height + 0.05, label = paste0("MAF = ", round(mainmaf, digits = 3)))) +
  theme_bw() + 
  scale_color_viridis(discrete = T)
plot
```

![](Demo_files/figure-markdown_github/unnamed-chunk-29-1.png)

The difference in coverage between peaks 1 and 2:

``` r
allPeaks$diff[which(allPeaks$npeak == 2)]
```

    ## [1] 3725.16

If we move left and right in steps of 3725.1595823 reads, we would notice that we land very close to peak 1, peak 6, peak 7, and peak 9:

``` r
diffvalue <- allPeaks$diff[which(allPeaks$npeak == 2)]
positions <- seq(from = -2*diffvalue, by = diffvalue, length.out = 6) + maxpeak

allPeaks_plot <- allPeaks
allPeaks_plot$selected <- F
allPeaks_plot$selected[allPeaks$npeak %in% c(2, 4)] <- T
plot <- ggplot(data = dendf, mapping = aes(x = x, y = y)) + geom_line() + 
  xlab("Normalized somatic read counts") + 
  ylab("Normalized Density") + 
  ggtitle(paste0("Extrapolating model for peaks 1 and 2 with a diff of one")) + 
  geom_vline(data = data.frame(positions = positions), aes(xintercept = positions), linetype = 2) + 
  geom_text(data = allPeaks_plot, aes(x = pos + maxpeak, y = allPeaks$height + 0.05, label = paste0("MAF = ", round(mainmaf, digits = 3)))) +
  theme_bw() + 
  scale_color_viridis(discrete = T)
plot
```

![](Demo_files/figure-markdown_github/unnamed-chunk-31-1.png)

This model was generated by selecting peak 2 as the comparator, and assuming a copy number difference between peak 1 and 2 of one. This model has matched with peaks 1, 2, 4, 6, 7, and 9.

Alternatively, we might assume that peak 2 is two copies away from peak 1, in which case we would move left and right in steps of 9340 reads.

``` r
diffvalue <- allPeaks$diff[which(allPeaks$npeak == 2)]/2
positions <- seq(from = -4*diffvalue, by = diffvalue, length.out = 12) + maxpeak

allPeaks_plot <- allPeaks
allPeaks_plot$selected <- F
allPeaks_plot$selected[allPeaks$npeak %in% c(2, 4)] <- T
plot <- ggplot(data = dendf, mapping = aes(x = x, y = y)) + geom_line() + 
  xlab("Normalized somatic read counts") + 
  ylab("Normalized Density") + 
  ggtitle(paste0("Extrapolating model for peaks 1 and 2 with a diff of two")) + 
  geom_vline(data = data.frame(positions = positions), aes(xintercept = positions), linetype = 2) + 
  geom_text(data = allPeaks_plot, aes(x = pos + maxpeak, y = allPeaks$height + 0.05, label = paste0("MAF = ", round(mainmaf, digits = 3)))) +
  theme_bw() + 
  scale_color_viridis(discrete = T)
plot
```

![](Demo_files/figure-markdown_github/unnamed-chunk-32-1.png)

This model was generated again by selecting peak 2 as the comparator, and assuming a copy number difference of two. While this model matched with all of the peaks, it predicted that there should be peaks in locations we don't find any, for example between peaks 1 and 2 or 6 and 7.

ploidetect\_roughmodels will do this for each peak, as well as a variety of copy number differences

``` r
xdists <- ploidetect_roughmodels(allPeaks = allPeaks, maxpeak, verbose = verbose, rerun = rerun)
knitr::kable(xdists)
```

| unmatched  | predunmatched                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |  reads\_per\_copy|  Copy\_number\_difference\_between\_peaks|  Comparator\_peak\_rank|
|:-----------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------:|-----------------------------------------:|-----------------------:|
|            | 13135.7231768991\_20741.1526938378\_24543.8674523071\_28346.5822107764                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |         1713.6151|                                         4|                       7|
|            | 13252.055941143\_20702.3751057564\_24427.5346880632\_28152.6942703699                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |         1772.2979|                                         2|                       2|
|            | 14033.788937233\_12431.786811443\_20441.7974403931\_23645.8016919732\_25247.8038177632\_28451.8080693432                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |         1700.8908|                                         2|                       4|
|            | 14112.8894220116\_12537.2541244811\_20415.4306121336\_23566.7012071946\_25142.3365047251\_28293.6070997861                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |         1711.8964|                                         5|                       6|
|            | 14254.2851603703\_12725.7817756261\_20368.2986993473\_23425.3054688359\_24953.8088535801\_28010.8156230686                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |         1666.3049|                                         1|                       3|
|            | 14276.5376044399\_12755.4517010522\_20360.8812179908\_23403.0530247663\_24924.138928154\_27966.3107349295\_29487.3966383172                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |         1711.6442|                                         5|                       7|
|            | 16356.3555930653\_13872.9158715275\_12631.1960107586\_20081.515175372\_21323.2350361409\_23806.6747576787\_25048.3946184476\_27531.8343399854\_28773.5542007543                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |         1289.8250|                                         3|                       2|
|            | 16703.7924802164\_14567.7896458297\_13499.7882286363\_12431.786811443\_19907.7967317964\_20975.7981489898\_23111.8009833765\_24179.8024005698\_25247.8038177632\_27383.8066521499\_28451.8080693432\_29519.8094865366                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |         1278.9220|                                         3|                       4|
|            | 17908.5054190264\_16045.9256278731\_14183.3458367197\_13252.055941143\_12320.7660455663\_19771.0852101798\_20702.3751057564\_22564.9548969098\_23496.2447924865\_24427.5346880632\_25358.8245836398\_27221.4043747932\_28152.6942703699\_29083.9841659465                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |          907.2244|                                         4|                       2|
|            | 18038.7942517081\_16436.792125918\_15635.791063023\_14033.788937233\_13232.787874338\_12431.786811443\_19640.7963774981\_20441.7974403931\_21242.7985032881\_22844.8006290782\_23645.8016919732\_24446.8027548682\_25247.8038177632\_26048.8048806582\_27650.8070064482\_28451.8080693432\_29252.8091322383                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |          765.5157|                                         4|                       4|
|            | 18075.543622231\_16547.0402374867\_15782.7885451146\_14254.2851603703\_13490.0334679982\_12725.7817756261\_11961.5300832539\_19604.0470069752\_20368.2986993473\_21132.5503917195\_22661.0537764637\_23425.3054688359\_24189.557161208\_24953.8088535801\_25718.0605459522\_27246.5639306965\_28010.8156230686\_28775.0673154408\_29539.3190078129                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |          760.2121|                                         2|                       3|
|            | 18094.7633981417\_16604.6995652191\_15859.6676487577\_14369.603815835\_13624.5718993737\_12879.5399829124\_12134.508066451\_19584.8272310644\_20329.8591475258\_21074.8910639871\_22564.9548969098\_23309.9868133711\_24055.0187298325\_24800.0506462938\_25545.0825627552\_26290.1144792165\_27780.1783121392\_28525.2102286005\_29270.2421450619                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |          754.2556|                                         5|                       2|
|            | 18198.9944642871\_16917.3927636551\_16276.591913339\_15635.791063023\_14354.189362391\_13713.388512075\_13072.587661759\_12431.786811443\_11790.985961127\_19480.5961649191\_20121.3970152351\_20762.1978655511\_21402.9987158671\_22684.6004164992\_23325.4012668152\_23966.2021171312\_24607.0029674472\_25247.8038177632\_25888.6046680792\_27170.2063687112\_27811.0072190272\_28451.8080693432\_29092.6089196593\_29733.4097699753                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |          723.4216|                                         5|                       4|
|            | 18330.294186355\_17820.7930581069\_16801.7908016107\_16292.2896733627\_15782.7885451146\_14763.7862886184\_14254.2851603703\_13744.7840321222\_13235.2829038741\_12725.7817756261\_12216.280647378\_11706.7795191299\_19349.2964428512\_19858.7975710993\_20368.2986993473\_20877.7998275954\_21387.3009558435\_22406.3032123397\_22915.8043405878\_23425.3054688359\_23934.8065970839\_24444.307725332\_24953.8088535801\_25463.3099818282\_25972.8111100763\_26991.8133665725\_27501.3144948205\_28010.8156230686\_28520.3167513167\_29029.8178795648\_29539.3190078129                                                                                                                                                                                                                                                                                                                                                                                                                                                               |          520.7406|                                         3|                       3|
|            | 18457.669468417\_18075.543622231\_17693.4177760449\_16929.1660836728\_16547.0402374867\_16164.9143913006\_15782.7885451146\_15400.6626989285\_14636.4110065564\_14254.2851603703\_13872.1593141843\_13490.0334679982\_13107.9076218121\_12725.7817756261\_12343.65592944\_11961.5300832539\_11579.4042370679\_19221.9211607892\_19604.0470069752\_19986.1728531613\_20368.2986993473\_20750.4245455334\_21132.5503917195\_21514.6762379055\_22278.9279302777\_22661.0537764637\_23043.1796226498\_23425.3054688359\_23807.4313150219\_24189.557161208\_24571.6830073941\_24953.8088535801\_25335.9346997662\_25718.0605459522\_26100.1863921383\_26482.3122383244\_27246.5639306965\_27628.6897768826\_28010.8156230686\_28392.9414692547\_28775.0673154408\_29157.1931616268\_29539.3190078129\_29921.444853999                                                                                                                                                                                                                        |          380.1060|                                         4|                       3|
|            | 18534.0946376542\_18228.3939607054\_17922.6932837565\_17616.9926068077\_17005.59125291\_16699.8905759611\_16394.1898990123\_16088.4892220634\_15782.7885451146\_15477.0878681657\_14865.686514268\_14559.9858373192\_14254.2851603703\_13948.5844834215\_13642.8838064726\_13337.1831295238\_13031.4824525749\_12725.7817756261\_12420.0810986772\_12114.3804217284\_11808.6797447795\_11502.9790678307\_19145.4959915519\_19451.1966685008\_19756.8973454496\_20062.5980223985\_20368.2986993473\_20673.9993762962\_20979.700053245\_21285.4007301939\_21591.1014071428\_22202.5027610405\_22508.2034379893\_22813.9041149382\_23119.604791887\_23425.3054688359\_23731.0061457847\_24036.7068227336\_24342.4074996824\_24648.1081766313\_24953.8088535801\_25259.509530529\_25565.2102074778\_25870.9108844267\_26176.6115613755\_26482.3122383244\_27093.7135922221\_27399.4142691709\_27705.1149461198\_28010.8156230686\_28316.5163000175\_28622.2169769663\_28927.9176539152\_29233.618330864\_29539.3190078129\_29845.0196847617 |          308.9638|                                         5|                       3|
| 2\_3\_5    |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |         5932.5924|                                         2|                       5|
| 2\_3\_5\_6 |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |        11211.8404|                                         1|                       5|
| 2\_3\_5\_7 |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |         7744.3688|                                         1|                       7|
| 3          |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |         3638.6490|                                         1|                       2|
| 3          | 13126.8859432517\_24552.7046859545                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |         3526.0733|                                         4|                       5|
| 3          | 13587.6776561681\_24091.9129730381                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |         3572.3631|                                         3|                       6|
| 3          | 13769.5089699773\_23910.0816592289\_28980.3680038546                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |         3571.9820|                                         3|                       7|
| 3          | 23410.1228116842\_27980.4503087654                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |         3491.3371|                                         5|                       5|
| 3          | 28451.8080693432                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |         3523.9135|                                         1|                       4|
| 7          | 12931.1629488637\_20809.3394365162\_24748.4276803425\_28687.5159241688                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |         1717.7592|                                         4|                       6|

You can ignore the "mean\_err" variable here. Each row of the xdists object describes a model, such as the ones above, that Ploidetect has fit to the data. The unmatched variable is a list of the peaks that were not matched by the model, and the predunmatched variable lists the positions of peaks that were predicted to occur but were not detected. This information is used as input to modelbuilder\_iterative(), which is a fairly verbose function that takes this information and selects a model that balances fitting a maximum number of peaks, predicting few imaginary peaks and minimizing the difference between predicted and observed allele frequencies using a dirty heuristic algorithm.

``` r
allPeaks$pos <- allPeaks$pos + maxpeak
TC_calls <- list()

for(i in 1:nrow(xdists)){
  modelbuilder_output <- modelbuilder_iterative(xdists[i,], allPeaks = allPeaks, lowest = lowest, filtered = filtered, strict = T, get_homd = F, mode = "TC", nomaf = nomaf, rerun = rerun, maxpeak = maxpeak, bw = bw)
  TC_calls <- c(TC_calls, list(modelbuilder_output$out))
  plots <- c(plots, list(modelbuilder_output$outplot))
}

TC_calls <- do.call(rbind.data.frame, TC_calls) %>% arrange(model_error)
knitr::kable(TC_calls)
```

|  reads\_per\_copy|  zero\_copy\_depth|  Ploidy|  tumour\_purity|  lowest\_peak\_CN|  maf\_error|  CN\_diff|  Comparator|  model\_error|  avg\_ploidy|  unmatchedpct|
|-----------------:|------------------:|-------:|---------------:|-----------------:|-----------:|---------:|-----------:|-------------:|------------:|-------------:|
|          3695.062|          11413.127|       2|       0.3930238|                 0|   0.1728321|         1|           2|      4.327698|         1.84|     0.1939974|
|          3207.020|          12219.616|       2|       0.3442180|                 0|   0.1752919|         1|           4|     20.097236|         1.87|     0.1939974|
|          3143.662|           6050.033|       4|       0.5096167|                 1|   0.0365015|         4|           5|     31.294220|         3.87|     0.1939974|
|          3143.662|           6050.033|       4|       0.5096167|                 1|   0.0365015|         3|           6|     31.294220|         3.87|     0.1939974|
|          2747.146|          12990.657|       2|       0.2972306|                 0|   0.1818888|         5|           5|     58.403267|         1.90|     0.1939974|
|          5833.416|          12997.097|       1|       0.4730322|                 0|   0.1051991|         2|           5|           Inf|         1.09|     0.9710928|
|         10989.455|           7889.823|       1|       0.7358501|                 0|         Inf|         1|           5|           Inf|         1.04|     0.9716307|
|          7744.655|          11096.765|       1|       0.5826099|                 0|   0.1454233|         1|           7|           Inf|         1.00|     0.9982511|

And we're back to the familiar output data! As an aside, we have essentially mirrored Ploidetect's default behavior:

``` r
subclonalcalls <- ploidetect(clonalcase)
subclonalcalls$TC_calls %>% knitr::kable()
```

|     |  reads\_per\_copy|  zero\_copy\_depth|  Ploidy|  tumour\_purity|  lowest\_peak\_CN|  maf\_error|  CN\_diff|  Comparator|  model\_error|  avg\_ploidy|  unmatchedpct|
|-----|-----------------:|------------------:|-------:|---------------:|-----------------:|-----------:|---------:|-----------:|-------------:|------------:|-------------:|
| 06  |         18792.711|           56603.12|       2|       0.3990445|                 0|   0.1671647|         2|           6|      26.46738|         1.84|     0.2217781|
| 012 |          9152.032|           57494.25|       4|       0.2414838|                 0|   0.4508861|         3|           8|      30.33578|         3.54|     0.0241590|
| 0   |          9286.932|           57046.53|       4|       0.2456198|                 0|   0.4511516|         1|           4|      31.39913|         3.69|     0.0000000|
| 01  |          9286.932|           57046.53|       4|       0.2456198|                 0|   0.4511516|         5|           7|      31.39913|         3.69|     0.0000000|
| 02  |          9286.932|           57046.53|       4|       0.2456198|                 0|   0.4511516|         4|           8|      31.39913|         3.69|     0.0000000|
| 03  |          9286.932|           57046.53|       4|       0.2456198|                 0|   0.4511516|         2|           2|      31.39913|         3.69|     0.0000000|
| 04  |          9286.932|           57046.53|       4|       0.2456198|                 0|   0.4511516|         1|           3|      31.39913|         3.69|     0.0000000|
| 08  |         18569.862|           56945.13|       2|       0.3947471|                 0|   0.1669210|         1|           2|      32.13835|         1.83|     0.2230469|
| 1   |          9312.878|           56927.71|       4|       0.2465242|                 1|   0.0257199|         5|           6|      33.09912|         3.69|     0.0000000|
| 15  |          9191.438|           57331.35|       4|       0.2427929|                 1|   0.0256816|         4|           7|      35.11734|         3.54|     0.0241590|
| 09  |         17804.138|           58191.69|       2|       0.3796193|                 0|   0.1666533|         3|           7|      41.16327|         1.76|     0.2459371|
| 010 |         17804.138|           58191.69|       2|       0.3796193|                 0|   0.1666533|         1|           5|      41.16327|         1.76|     0.2459371|
| 011 |         17804.138|           58191.69|       2|       0.3796193|                 0|   0.1666533|         2|           8|      41.16327|         1.76|     0.2459371|
| 11  |          8582.935|           42258.50|       6|       0.2888692|                 1|   0.0240996|         2|           5|      79.17911|         5.72|     0.0000000|
| 07  |         20955.395|           53192.70|       2|       0.4406861|                 0|   0.1749926|         3|           6|     102.32783|         1.80|     0.2230469|
| 14  |         15855.421|           29793.76|       4|       0.5155849|                 1|   0.0337180|         4|           6|     151.53217|         3.88|     0.2217781|
| 12  |          8038.968|           45183.86|       6|       0.2624464|                 1|   0.0243461|         5|           8|     166.63327|         5.74|     0.0000000|
| 05  |         54974.059|           39417.02|       1|       0.7361030|                 0|         Inf|         1|           6|           Inf|         1.04|     0.9747319|
| 13  |         44282.827|            5708.16|       2|       0.9394513|                 1|   0.2307235|         1|           7|           Inf|         2.01|     0.9976222|

So now, an interesting thought arises; what if we decided that Ploidetect's top model (Ploidy = 2, purity = 0.3990445) was wrong? What if, based on our own inspection of the data and what we have learned so far, we decided that Ploidetect made a mistake and the correct model should have been the second on the list? Or what if the model looks correct but the ploidy is wrong? This is why Ploidetect provides a ton of plots for you - it helps with manual review.

The first two plots that Ploidetect returns will be a scatter plot of the coverage and a plot of all of the called peaks:

``` r
subclonalcalls$plots[1:2]
```

    ## [[1]]

![](Demo_files/figure-markdown_github/unnamed-chunk-36-1.png)

    ## 
    ## [[2]]

![](Demo_files/figure-markdown_github/unnamed-chunk-36-2.png)

The next plots will illustrate the models that Ploidetect used for each of the purity predictions

``` r
subclonalcalls$plots[c(3, 4)]
```

    ## [[1]]

![](Demo_files/figure-markdown_github/unnamed-chunk-37-1.png)

    ## 
    ## [[2]]

![](Demo_files/figure-markdown_github/unnamed-chunk-37-2.png)

So let's say that the other model looks more correct for whatever reason. We can force Ploidetect to pick this model, using the information in the table or the plot (note the CN\_diff and comparator values):

``` r
newsubclonalcalls <- ploidetect(clonalcase, cndiff = 3, comp = 8, lowest = 0)
newsubclonalcalls$TC_calls %>% knitr::kable()
```

|     |  reads\_per\_copy|  zero\_copy\_depth|  Ploidy|  tumour\_purity|  lowest\_peak\_CN|  maf\_error|  CN\_diff|  Comparator|  model\_error|  avg\_ploidy|  unmatchedpct|
|-----|-----------------:|------------------:|-------:|---------------:|-----------------:|-----------:|---------:|-----------:|-------------:|------------:|-------------:|
| 0   |          9152.032|           57494.25|       4|       0.2414838|                 0|   0.4508861|         3|           8|      30.33578|         3.54|      0.024159|

``` r
newsubclonalcalls$plots[3]
```

    ## [[1]]

![](Demo_files/figure-markdown_github/unnamed-chunk-38-1.png)

These features (and indeed, this tutorial) is there for the inevitable cases where Ploidetect gets it wrong. It isn't helpful here, but Ploidetect includes a helpful function for testing allele frequency values at specific copy numbers and tumour purities. Here we'll demonstrate why things can get a bit unclear sometimes by comparing predicted allele frequency values at the top two tumour purities. First, the predicteds for Peak \#2 under the two models:

``` r
testMAF(CN = 1, tp = 0.39)
```

    ##        1 
    ## 0.621118

``` r
testMAF(CN = 2, tp = 0.24)
```

    ##    1    2 
    ## 0.50 0.62

We notice that a diploid genome would predict an allele frequency of 0.62 for peak 2, however a tetraploid model would also predict a similar allele frequency. Which one is correct?

This is an example of the [identifiability problem](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4103592/). Essentially both assumptions of purity and ploidy explain the data equally well. But which fit looks better? The conservative model that skips the smaller peaks immediately adjacent to peak \#1 or the greedy model which fits essentially everything? In practice we tend to see that integer copy number gains and losses relative to tumour ploidy have a steadily declining quantity - that is, you would usually expect the \#2 peak by height to represent a copy number difference of one. Smaller peaks falling in between, as seen in this case, generally correspond to subclones.

Now that you are familiar with the rationale behind Ploidetect, there are two more cases to go over to further illustrate the manual review process.

``` r
very_crazy_case <- read.table("horrible_no_good_very_bad_data.bed", header = T, stringsAsFactors = F)
str(very_crazy_case)
```

    ## 'data.frame':    145014 obs. of  6 variables:
    ##  $ V4    : int  17791 17481 15982 16539 16496 16035 16115 16688 15947 16147 ...
    ##  $ V5    : int  1000006 1000187 1000002 1000012 1000003 1000031 1000018 1000019 1000005 1000012 ...
    ##  $ V6    : chr  "." "." "." "." ...
    ##  $ window: chr  "1_0" "1_251674" "1_569481" "1_693498" ...
    ##  $ size  : int  251674 317807 124017 36723 27929 24812 23377 19706 23959 27395 ...
    ##  $ V7    : num  0.326 0.298 0.41 0.432 0.421 ...

``` r
output <- ploidetect(very_crazy_case)
```

First, the purity calls:

``` r
output$TC_calls %>% knitr::kable()
```

|        |   reads\_per\_copy|  zero\_copy\_depth|  Ploidy|  tumour\_purity|  lowest\_peak\_CN|  maf\_error|  CN\_diff|  Comparator|  model\_error|  avg\_ploidy|  unmatchedpct|
|--------|------------------:|------------------:|-------:|---------------:|-----------------:|-----------:|---------:|-----------:|-------------:|------------:|-------------:|
| 212    |          12780.104|          22242.505|       4|       0.5347020|                 2|   0.0343097|         5|           9|      42.46588|         5.28|     0.0681942|
| 214    |          12419.224|          23648.108|       4|       0.5122749|                 2|   0.0321826|         3|           6|      43.55131|         5.22|     0.1301886|
| 210    |          10783.954|          30362.612|       4|       0.4153224|                 2|   0.0383839|         3|           7|      47.16489|         5.50|     0.0681942|
| 211    |          10783.954|          30362.612|       4|       0.4153224|                 2|   0.0383839|         4|           8|      47.16489|         5.50|     0.0681942|
| 112    |          21193.422|          30980.631|       2|       0.5777334|                 1|   0.0427840|         1|          12|      49.11205|         2.71|     0.6836550|
| 113    |          21193.422|          30980.631|       2|       0.5777334|                 1|   0.0427840|         3|          11|      49.11205|         2.71|     0.6836550|
| 27     |          13182.326|          20579.537|       4|       0.5616169|                 2|   0.0379463|         5|           8|      64.64948|         5.27|     0.1038306|
| 11     |          22168.291|          28973.707|       2|       0.6047798|                 1|   0.0388941|         3|           9|      70.36832|         2.66|     0.7272209|
| 12     |          22168.291|          28973.707|       2|       0.6047798|                 1|   0.0388941|         3|          10|      70.36832|         2.66|     0.7272209|
| 216    |          10356.196|          32156.624|       4|       0.3917681|                 2|   0.0427555|         5|          11|     100.93556|         5.45|     0.1246034|
| 2      |           8089.872|          31526.581|       5|       0.3391530|                 2|   0.0367286|         3|          12|     187.91119|         7.39|     0.0202601|
| 21     |           7677.958|          34012.680|       5|       0.3110462|                 2|   0.0403044|         5|           6|     198.29951|         7.67|     0.0010285|
| 13     |          21095.248|          31229.938|       2|       0.5746424|                 1|   0.0431699|         2|           6|     199.80839|         2.73|     0.6826265|
| 14     |          21095.248|          31229.938|       2|       0.5746424|                 1|   0.0431699|         4|           9|     199.80839|         2.73|     0.6826265|
| 15     |          21095.248|          31229.938|       2|       0.5746424|                 1|   0.0431699|         1|           3|     199.80839|         2.73|     0.6826265|
| 16     |          21095.248|          31229.938|       2|       0.5746424|                 1|   0.0431699|         5|          10|     199.80839|         2.73|     0.6826265|
| 17     |          21095.248|          31229.938|       2|       0.5746424|                 1|   0.0431699|         3|           5|     199.80839|         2.73|     0.6826265|
| 18     |          21095.248|          31229.938|       2|       0.5746424|                 1|   0.0431699|         3|           8|     199.80839|         2.73|     0.6826265|
| 110    |          21095.248|          31229.938|       2|       0.5746424|                 1|   0.0431699|         4|          11|     199.80839|         2.73|     0.6826265|
| 24     |          11213.070|          28556.149|       4|       0.4398810|                 2|   0.0345789|         4|           7|     212.86976|         5.51|     0.1028021|
| 25     |          11213.070|          28556.149|       4|       0.4398810|                 2|   0.0345789|         1|           4|     212.86976|         5.51|     0.1028021|
| 213    |          12517.469|          23168.487|       4|       0.5193601|                 2|   0.0327334|         4|           5|     220.43327|         5.25|     0.1291601|
| 1      |          22330.046|          28572.855|       2|       0.6098361|                 1|   0.0383666|         2|           5|     223.64591|         2.69|     0.7261924|
| 23     |          11154.662|          28851.950|       4|       0.4360588|                 2|   0.0350963|         5|           5|     240.94864|         5.52|     0.1028021|
| 26     |          11154.662|          28851.950|       4|       0.4360588|                 2|   0.0350963|         2|          12|     240.94864|         5.52|     0.1028021|
| 19     |          21455.141|          30323.021|       2|       0.5859395|                 1|   0.0418044|         4|          10|     242.36634|         2.73|     0.6826265|
| 111    |          21455.141|          30323.021|       2|       0.5859395|                 1|   0.0418044|         2|           7|     242.36634|         2.73|     0.6826265|
| 28     |          10853.018|          30238.232|       4|       0.4178716|                 2|   0.0379424|         2|           3|     242.83079|         5.50|     0.1051463|
| 215    |          12271.869|          24312.278|       4|       0.5023688|                 2|   0.0316262|         1|           2|     243.08693|         5.27|     0.1235749|
| 22     |          11088.591|          29190.314|       4|       0.4317357|                 2|   0.0357226|         4|           6|     303.25275|         5.52|     0.1028021|
| 29     |           9357.819|          27159.610|       5|       0.4079681|                 2|   0.0341634|         5|           7|     439.77360|         6.71|     0.1051463|
| 114    |          48807.625|          25024.887|       1|       0.7959486|                 1|   0.1107488|         1|           5|    2117.33144|         1.27|     0.8753966|
| 115    |          86543.791|         -13033.579|       1|       1.0814324|                 1|         Inf|         1|          10|           Inf|         1.07|     0.9853251|
| 116    |          76514.017|          -3125.201|       1|       1.0208482|                 1|         Inf|         1|          11|           Inf|         1.07|     0.9820841|
| 117    |          83339.716|          -9868.285|       1|       1.0629310|                 1|         Inf|         1|           9|           Inf|         1.08|     0.9797399|
| 0      |          43035.108|          31135.536|       1|       0.7343515|                 0|         Inf|         2|          10|           Inf|         1.22|     0.9338975|
| 01     |          43035.108|          31135.536|       1|       0.7343515|                 0|         Inf|         1|           7|           Inf|         1.22|     0.9338975|
| 118    |          38596.738|          -3186.629|       2|       1.0430586|                 1|   0.3426988|         2|           8|           Inf|         2.24|     0.8926759|
| 02     |          40851.202|          33073.690|       1|       0.7118416|                 0|         Inf|         2|          11|           Inf|         1.26|     0.8903316|
| 03     |          40851.202|          33073.690|       1|       0.7118416|                 0|         Inf|         4|          13|           Inf|         1.26|     0.8903316|
| 04     |          40851.202|          33073.690|       1|       0.7118416|                 0|         Inf|         1|           6|           Inf|         1.26|     0.8903316|
| 05     |          40851.202|          33073.690|       1|       0.7118416|                 0|         Inf|         2|           9|           Inf|         1.26|     0.8903316|
| 06     |          40851.202|          33073.690|       1|       0.7118416|                 0|         Inf|         5|          13|           Inf|         1.26|     0.8903316|
| Lots o |  f models. Let's s|      ee the plots:|        |                |                  |            |          |            |              |             |              |

``` r
output$plots[1:3]
```

    ## [[1]]

![](Demo_files/figure-markdown_github/unnamed-chunk-42-1.png)

    ## 
    ## [[2]]

![](Demo_files/figure-markdown_github/unnamed-chunk-42-2.png)

    ## 
    ## [[3]]

![](Demo_files/figure-markdown_github/unnamed-chunk-42-3.png)

...yuck. This is an example of when things have gone quite wrong at some stage. Most likely this is the result of sequencing error, a high degree of subclonal copy number variation, or both.

Basically cases that look like this can't be trusted to be called correctly. It might be, but it's quite likely that it isn't. As a start, we can eliminate a bunch of possibilities based on the allele frequencies. The top peak has allele frequency of 0.825 - quite high. The top model claims that this is four-copy.

``` r
testMAF(4, 0.53)
```

    ##         2         3         4 
    ## 0.5000000 0.6732026 0.8464052

So this is possible, but it would require four-copy homozygosity. This seems fairly extreme. We can next test the two-copy, 57% or 60% model:

``` r
testMAF(2, 0.57)
```

    ##     1     2 
    ## 0.500 0.785

``` r
testMAF(2, 0.6)
```

    ##   1   2 
    ## 0.5 0.8

Also not the best, although 60% is fairly close. Let's check out how that fit looks:

``` r
output$plots[8]
```

    ## [[1]]

![](Demo_files/figure-markdown_github/unnamed-chunk-45-1.png)

This model assumes that CN = 1 has allele frequency 0.752 and CN = 3 has allele frequency of 0.63. Let's test:

``` r
testMAF(CN = 1, 0.60)
```

    ##         1 
    ## 0.7142857

``` r
testMAF(CN = 3, 0.60)
```

    ##         2         3 
    ## 0.6153846 0.8461538

It looks like the 60% purity model predicts allele frequencies slightly lower than what we observe - indicating that this tumour is slightly above 60% purity. However, this trend is consistent. The model is likely correct, but inaccuracies in fitting the model resulted in an inaccurate estimation of purity.

``` r
testMAF(CN = 1, 0.67)
```

    ##         1 
    ## 0.7518797

The error in the estimate is about 7%. In practice this is probably the worst case example for Ploidetect's accuracy, although 7% isn't likely to severely affect downstream variant calling.

One last example, just to drill in these points a bit further.
