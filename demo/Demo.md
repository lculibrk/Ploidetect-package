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

First, we call ploidetect\_preprocess, which performs basic preprocessing of the data. First, Ploidetect converts combines sequential bins into ~100kb sized new bins for more accurate modeling. It corrects for lingering GC biases using a loess fit and also corrects the read depth using the normal read depth, since using constant depth bins can fail at the terminal end of the chromsome due to being truncated by the end of the chromosome. We get four objects, one of which is important for our sake.

``` r
preprocess_output <- ploidetect_preprocess(all_data = all_data)
str(preprocess_output)
```

    ## List of 4
    ##  $ x           :'data.frame':    30000 obs. of  8 variables:
    ##   ..$ y_raw          : num [1:30000] 84054 88603 84898 87641 86285 ...
    ##   ..$ x_raw          : int [1:30000] 5000075 5000083 5000093 5000118 5000150 5000062 5000085 5000119 5000071 5000112 ...
    ##   ..$ maf            : num [1:30000] 0.432 NaN NaN 0.607 0.382 ...
    ##   ..$ window         : chr [1:30000] "1_737630" "1_852346" "1_994324" "1_1130229" ...
    ##   ..$ size           : int [1:30000] 114716 141978 135905 143522 133304 133542 158099 106973 125470 128645 ...
    ##   ..$ GC             : num [1:30000] 0.469 0.626 0.59 0.627 0.575 ...
    ##   ..$ residual       : num [1:30000] -7529 -3982 -7554 -4946 -6093 ...
    ##   ..$ normalized_size: int [1:30000] 114716 141978 135905 143522 133304 133542 158099 106973 125470 128645 ...
    ##  $ maxpeak     : num 94203
    ##  $ highoutliers:'data.frame':    107 obs. of  8 variables:
    ##   ..$ merge : num [1:107] 1243 1244 1245 1246 1247 ...
    ##   ..$ chr   : chr [1:107] "1" "1" "1" "1" ...
    ##   ..$ tumour: num [1:107] 232944 250380 274998 302823 323815 ...
    ##   ..$ normal: int [1:107] 5153203 5244638 5287191 5364167 5255397 5427384 5216989 5522157 5375805 5595164 ...
    ##   ..$ maf   : num [1:107] NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN ...
    ##   ..$ wind  : chr [1:107] "1_121484470" "1_121484520" "1_121484562" "1_121484594" ...
    ##   ..$ size  : int [1:107] 50 42 32 28 26 32 37 39 31 28 ...
    ##   ..$ gc    : num [1:107] 0.38 0.469 0.286 0.493 0.267 ...
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

It's not much. In some cases, this step can be extremely helpful in reducing the mapping biases in the data. Ploidetect does a quick check to see if this step manages to actually reduce the variability. The idea behind this is that the worse data will have points far more spread out, resuting in a higher average normalized density

Ploidetect uses a kernel density estimation to find peaks in the depth histogram. We use a simple heuristic to get a decent bandwidth for peak calling.

``` r
bw = maxpeak/40
```

Now we perform the density estimation that the rest of the analysis is based on. Afterward, we do a quick check to ensure that the correction step actually helped in correcting biases in the depth - if it didn't, we revert back to the raw depth.

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

|       start|         end|     height|           pos|     trough|  ratiotrough|  troughdiff|       diff|
|-----------:|-----------:|----------:|-------------:|----------:|------------:|-----------:|----------:|
|   -2529.101|    2802.118|  1.0000000|      49.21603|  0.2992419|    0.2992419|   0.7007581|      0.000|
|  -21048.070|  -15997.442|  0.5620062|  -18607.51988|  0.1824274|    0.3246003|   0.3795789|  18656.736|
|  -11227.405|   -7579.729|  0.3217848|   -8797.92659|  0.2992419|    0.9299443|   0.0225429|   8847.143|
|    9536.289|   11219.831|  0.1603730|    9539.53482|  0.1604204|    1.0002956|  -0.0000474|   9490.319|
|   15709.279|   20759.907|  0.1241323|   16850.45812|  0.1179252|    0.9499961|   0.0062071|  16801.242|
|   53588.990|   60884.341|  0.0510317|   57051.28192|  0.0117187|    0.2296357|   0.0393130|  57002.066|
|   44610.095|   47416.000|  0.0117187|   47412.89380|  0.0117187|    1.0000000|   0.0000000|  47363.678|
|  -40128.221|  -35358.183|  0.0055449|  -37801.00712|  0.0004143|    0.0747157|   0.0051306|  37850.223|

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

|       start|         end|     height|         pos|     trough|  ratiotrough|  troughdiff|       diff|  npeak|    mainmaf|
|-----------:|-----------:|----------:|-----------:|----------:|------------:|-----------:|----------:|------:|----------:|
|  -40177.437|  -35407.399|  0.0055449|  -37850.223|  0.0004143|    0.0747157|   0.0051306|  37850.223|      1|  0.5424196|
|  -21097.286|  -16046.658|  0.5620062|  -18656.736|  0.1824274|    0.3246003|   0.3795789|  18656.736|      2|  0.6181312|
|  -11276.621|   -7628.945|  0.3217848|   -8847.143|  0.2992419|    0.9299443|   0.0225429|   8847.143|      3|  0.5591691|
|   -2578.317|    2752.902|  1.0000000|       0.000|  0.2992419|    0.2992419|   0.7007581|      0.000|      4|  0.5265049|
|    9487.073|   11170.615|  0.1603730|    9490.319|  0.1604204|    1.0002956|  -0.0000474|   9490.319|      5|  0.5468677|
|   15660.063|   20710.691|  0.1241323|   16801.242|  0.1179252|    0.9499961|   0.0062071|  16801.242|      6|  0.5727502|
|   44560.879|   47366.784|  0.0117187|   47363.678|  0.0117187|    1.0000000|   0.0000000|  47363.678|      7|  0.7438366|
|   53539.773|   60835.125|  0.0510317|   57002.066|  0.0117187|    0.2296357|   0.0393130|  57002.066|      8|  0.6797716|

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

    ## [1] 18656.74

If we move left and right in steps of 1.865673610^{4} reads, we would notice that we land very close to peak 1, peak 6, peak 7, and peak 9:

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

| unmatched     | predunmatched                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |  reads\_per\_copy|  Copy\_number\_difference\_between\_peaks|  Comparator\_peak\_rank|
|:--------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------:|-----------------------------------------:|-----------------------:|
|               | 117004.087557932\_128404.500735733                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |          9229.667|                                         5|                       6|
|               | 65732.3048410155\_122674.217563646\_132164.536350751                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |          9180.292|                                         1|                       4|
|               | 65785.0545407585\_122621.467863903\_132094.20341776                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |          9143.207|                                         5|                       7|
|               | 65815.5938406098\_122590.928564052\_132053.484351292                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |          9141.524|                                         4|                       8|
|               | 66218.1573386489\_122188.365066013\_131516.733020573                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |          9198.917|                                         2|                       2|
|               | 67661.8333316166\_120744.689073045\_129591.831696616                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |          9079.596|                                         1|                       3|
|               | 69001.3980750914\_60600.7770326784\_119405.12432957\_127805.745371983\_136206.366414396                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |          9072.627|                                         2|                       5|
|               | 71493.127312954\_63923.0826831617\_116913.395091707\_124483.4397215\_132053.484351292\_147193.573610876                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |          9133.076|                                         5|                       8|
|               | 81765.4372629165\_69327.6133235024\_63108.7013537953\_100422.173172038\_119078.909081159\_125297.821050866\_131516.733020573\_137735.64499028                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |          5704.451|                                         3|                       2|
|               | 88602.8471740553\_71801.6050892291\_66201.1910609537\_60600.7770326784\_99803.6752306061\_116604.917315432\_122205.331343708\_127805.745371983\_133406.159400258\_144606.987456809                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |          4812.490|                                         3|                       5|
|               | 89458.1018087782\_79967.7830216731\_70477.464234568\_65732.3048410155\_60987.145447463\_98948.4205958832\_108438.739382988\_117929.058170093\_122674.217563646\_127419.376957198\_132164.536350751\_136909.695744304\_146400.014531409                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |          4590.146|                                         2|                       4|
|               | 89539.0772250504\_80210.7092704898\_70882.3413159292\_66218.1573386489\_61553.9733613686\_98867.445179611\_108195.813134172\_117524.181088732\_122188.365066013\_126852.549043293\_131516.733020573\_136180.916997853\_145509.284952414                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |          4599.459|                                         4|                       2|
|               | 89779.689890545\_80932.5472669736\_72085.4046434022\_67661.8333316166\_63238.2620198309\_58814.6907080452\_98626.8325141164\_107473.975137688\_116321.117761259\_120744.689073045\_125168.260384831\_129591.831696616\_134015.403008402\_138438.974320188\_147286.116943759                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |          4524.594|                                         2|                       3|
|               | 90002.9506811242\_81602.3296387111\_73201.708596298\_69001.3980750914\_64801.0875538849\_60600.7770326784\_98403.5717235372\_106804.19276595\_115204.813808363\_119405.12432957\_123605.434850777\_127805.745371983\_132006.05589319\_136206.366414396\_144606.987456809\_148807.297978016                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |          4518.902|                                         4|                       5|
|               | 90471.9140205065\_83009.219656858\_79277.8724750337\_71815.1781113852\_68083.830929561\_64352.4837477367\_60621.1365659125\_97934.6083841549\_101665.955565979\_109128.649929628\_116591.344293276\_120322.6914751\_124054.038656925\_127785.385838749\_131516.733020573\_135248.080202397\_138979.427384222\_146442.12174787                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |          3773.748|                                         5|                       2|
|               | 90843.0127853655\_87482.7643684002\_80762.2675344698\_77402.0191175045\_70681.5222835741\_67321.2738666088\_63961.0254496436\_60600.7770326784\_97563.5096192959\_100923.758036261\_107644.254870192\_114364.751704122\_117725.000121087\_121085.248538053\_124445.496955018\_127805.745371983\_131165.993788948\_134526.242205914\_137886.490622879\_144606.987456809\_147967.235873774                                                                                                                                                                                                                                                                                                                                                                                                                                                               |          3135.317|                                         5|                       5|
|               | 91039.821606629\_87876.3820109273\_81549.5028195239\_78386.0632238222\_72059.1840324189\_68895.7444367172\_65732.3048410155\_62568.8652453138\_59405.4256496121\_97366.7007980324\_100530.140393734\_106857.019585137\_113183.898776541\_116347.338372243\_119510.777967944\_122674.217563646\_125837.657159348\_129001.096755049\_132164.536350751\_135327.975946453\_138491.415542154\_144818.294733558\_147981.734329259                                                                                                                                                                                                                                                                                                                                                                                                                            |          3109.862|                                         3|                       4|
|               | 91254.2136611402\_88305.1661199498\_82407.0710375689\_79458.0234963784\_73559.9284139975\_70610.880872807\_67661.8333316166\_64712.7857904261\_61763.7382492356\_58814.6907080452\_97152.3087435212\_100101.356284712\_105999.451367093\_108948.498908283\_114846.593990664\_117795.641531854\_120744.689073045\_123693.736614235\_126642.784155426\_129591.831696616\_132540.879237807\_135489.926778997\_138438.974320188\_144337.069402569\_147286.116943759                                                                                                                                                                                                                                                                                                                                                                                        |          3018.727|                                         3|                       3|
|               | 91830.6815055544\_89458.1018087782\_87085.5221120019\_82340.3627184494\_79967.7830216731\_77595.2033248968\_72850.0439313443\_70477.464234568\_68104.8845377918\_65732.3048410155\_63359.7251442392\_60987.145447463\_58614.5657506867\_96575.840899107\_98948.4205958832\_101321.00029266\_106066.159686212\_108438.739382988\_113183.898776541\_115556.478473317\_117929.058170093\_120301.63786687\_122674.217563646\_125046.797260422\_127419.376957198\_129791.956653975\_132164.536350751\_134537.116047527\_136909.695744304\_139282.27544108\_144027.434834632\_146400.014531409\_148772.594228185                                                                                                                                                                                                                                             |          2321.732|                                         4|                       4|
|               | 91991.4755464379\_89779.689890545\_87567.9042346522\_83144.3329228665\_80932.5472669736\_78720.7616110808\_74297.1902992951\_72085.4046434022\_69873.6189875094\_67661.8333316166\_65450.0476757237\_63238.2620198309\_61026.476363938\_58814.6907080452\_96415.0468582235\_98626.8325141164\_100838.618170009\_105262.189481795\_107473.975137688\_109685.760793581\_114109.332105366\_116321.117761259\_118532.903417152\_120744.689073045\_122956.474728938\_125168.260384831\_127380.046040723\_129591.831696616\_131803.617352509\_134015.403008402\_136227.188664295\_138438.974320188\_142862.545631973\_145074.331287866\_147286.116943759\_149497.902599652                                                                                                                                                                                   |          2289.343|                                         4|                       3|
|               | 92305.1974449097\_90407.1336874887\_88509.0699300677\_86611.0061726467\_82814.8786578046\_80916.8149003836\_79018.7511429626\_77120.6873855416\_73324.5598706996\_71426.4961132785\_69528.4323558575\_67630.3685984365\_65732.3048410155\_63834.2410835945\_61936.1773261735\_60038.1135687525\_58140.0498113314\_96101.3249597517\_97999.3887171727\_99897.4524745937\_101795.516232015\_105591.643746857\_107489.707504278\_109387.771261699\_113183.898776541\_115081.962533962\_116980.026291383\_118878.090048804\_120776.153806225\_122674.217563646\_124572.281321067\_126470.345078488\_128368.408835909\_130266.47259333\_132164.536350751\_134062.600108172\_135960.663865593\_137858.727623014\_139756.791380435\_143552.918895277\_145450.982652698\_147349.046410119\_149247.11016754                                                     |          1852.647|                                         5|                       4|
|               | 92433.8326776164\_90664.4041529021\_88894.9756281879\_87125.5471034736\_83586.690054045\_81817.2615293308\_80047.8330046165\_78278.4044799022\_76508.9759551879\_72970.1189057594\_71200.6903810451\_69431.2618563308\_67661.8333316166\_65892.4048069023\_64122.976282188\_62353.5477574737\_60584.1192327595\_58814.6907080452\_95972.689727045\_97742.1182517593\_99511.5467764735\_101280.975301188\_104819.832350616\_106589.260875331\_108358.689400045\_111897.546449473\_113666.974974188\_115436.403498902\_117205.832023616\_118975.260548331\_120744.689073045\_122514.117597759\_124283.546122473\_126052.974647188\_127822.403171902\_129591.831696616\_131361.260221331\_133130.688746045\_134900.117270759\_136669.545795473\_138438.974320188\_140208.402844902\_143747.25989433\_145516.688419045\_147286.116943759\_149055.545468473 |          1827.913|                                         5|                       3|
| 2\_3\_5\_6\_7 |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |         56015.204|                                         1|                       6|
| 2\_3\_5\_6\_8 |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |         44140.844|                                         1|                       7|
| 3\_5          |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |         19091.233|                                         2|                       6|
| 3\_5          | 65702.2282578287\_122704.294146833                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |         17772.563|                                         4|                       6|
| 3\_5\_7       |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |         18433.697|                                         3|                       6|
| 3\_5\_7       | 131516.733020573                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |         18494.796|                                         1|                       2|
| 3\_5\_8       | 125779.046381855                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |         18242.713|                                         3|                       7|
| 3\_5\_8       | 127805.745371983                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |         18061.251|                                         1|                       5|
| 3\_5\_8       | 132053.484351292                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |         18288.627|                                         2|                       8|
| 8             | 117885.100086974\_129726.019529296                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |          9151.164|                                         4|                       7|
| 8             | 68969.7791030232\_119436.743301638\_132053.484351292                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |          9140.834|                                         3|                       8|

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
|         18792.711|           56603.12|       2|       0.3990445|                 0|   0.1671647|         2|           6|      26.46738|         1.84|     0.2217781|
|          9152.032|           57494.25|       4|       0.2414838|                 0|   0.4508861|         3|           8|      30.33578|         3.54|     0.0241590|
|          9286.932|           57046.53|       4|       0.2456198|                 0|   0.4511516|         1|           4|      31.39913|         3.69|     0.0000000|
|          9286.932|           57046.53|       4|       0.2456198|                 0|   0.4511516|         5|           7|      31.39913|         3.69|     0.0000000|
|          9286.932|           57046.53|       4|       0.2456198|                 0|   0.4511516|         4|           8|      31.39913|         3.69|     0.0000000|
|          9286.932|           57046.53|       4|       0.2456198|                 0|   0.4511516|         2|           2|      31.39913|         3.69|     0.0000000|
|          9286.932|           57046.53|       4|       0.2456198|                 0|   0.4511516|         1|           3|      31.39913|         3.69|     0.0000000|
|         18569.862|           56945.13|       2|       0.3947471|                 0|   0.1669210|         1|           2|      32.13835|         1.83|     0.2230469|
|          9312.878|           56927.71|       4|       0.2465242|                 1|   0.0257199|         5|           6|      33.09912|         3.69|     0.0000000|
|          9191.438|           57331.35|       4|       0.2427929|                 1|   0.0256816|         4|           7|      35.11734|         3.54|     0.0241590|
|         17804.138|           58191.69|       2|       0.3796193|                 0|   0.1666533|         3|           7|      41.16327|         1.76|     0.2459371|
|         17804.138|           58191.69|       2|       0.3796193|                 0|   0.1666533|         1|           5|      41.16327|         1.76|     0.2459371|
|         17804.138|           58191.69|       2|       0.3796193|                 0|   0.1666533|         2|           8|      41.16327|         1.76|     0.2459371|
|          8582.935|           42258.50|       6|       0.2888692|                 1|   0.0240996|         2|           5|      79.17911|         5.72|     0.0000000|
|         20955.395|           53192.70|       2|       0.4406861|                 0|   0.1749926|         3|           6|     102.32783|         1.80|     0.2230469|
|         15855.421|           29793.76|       4|       0.5155849|                 1|   0.0337180|         4|           6|     151.53217|         3.88|     0.2217781|
|          8038.968|           45183.86|       6|       0.2624464|                 1|   0.0243461|         5|           8|     166.63327|         5.74|     0.0000000|
|         54974.059|           39417.02|       1|       0.7361030|                 0|         Inf|         1|           6|           Inf|         1.04|     0.9747319|
|         44282.827|            5708.16|       2|       0.9394513|                 1|   0.2307235|         1|           7|           Inf|         2.01|     0.9976222|

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
