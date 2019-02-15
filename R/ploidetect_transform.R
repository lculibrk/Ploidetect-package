ploidetect_transform <- function(x, bw, normal = 2, tumour = 1, avg_allele_freq = 3, window_id = 4, window_size = 5, GC = 6, verbose = F){
  ## First linearTransform with slope 0
  filtered <- linearTransform(x, tumour = tumour, normal = normal, avg_allele_freq = avg_allele_freq, window_id = window_id, slope = 0, verbose = verbose)
## Generate list of peaks called with transformation slope = 0 and heavy outlier filtering
  allPeaks <- peakcaller(filtered[findInterval(filtered$residual, vec = quantile(filtered$residual, probs = c(0.01, 0.99))) == 1,], bw)
  ## If we call zero peaks (due to poor bw, for example)
  if(nrow(allPeaks) == 1){
    warning("Zero peaks detected. Attempting peak calling with lower bandwidth")
    bw = bw/2
    allPeaks <- peakcaller(filtered[findInterval(filtered$residual, vec = quantile(filtered$residual, probs = c(0.01, 0.99))) == 1,], bw)
  }
  ## Perform bisection search for optimal slope
  ## Because of the way this is set up, we can input a heavily quantile-filtered "filtered" data.frame into optimizesmoothing
  ## which linearTransform will use to transform the unfiltered data, thus preventing the outliers from interfering with the transformation
  if(verbose){
    print("Optimizing slope for linear transformation")
  }
  filtered <- linearTransform(x = x, 
                              tumour = tumour, 
                              normal = normal, 
                              avg_allele_freq = avg_allele_freq, 
                              window_id = window_id, 
                              slope = optimizesmoothing(x = x, tumour = tumour, normal = normal, avg_allele_freq = avg_allele_freq, window_id = window_id, filtered = filtered[findInterval(filtered$residual, vec = quantile(filtered$residual, probs = c(0.001, 0.999))) == 1,], allPeaks = allPeaks, verbose = verbose, bw = bw),
                              verbose = verbose)
  ## Perform peak calls with much more permissive quantile filtering
  allPeaks <- peakcaller(filtered[findInterval(filtered$residual, vec = quantile(filtered$residual, probs = c(0.001, 0.999))) == 1,], bw)
  allPeaks$meansig <- NULL
  ## Center the residual and peak data about the tallest peak
  filtered$residual <- filtered$residual - allPeaks$pos[1]
  allPeaks$pos <- allPeaks$pos - allPeaks$pos[1]
  output <- list("filtered" = filtered, "allPeaks" = allPeaks)
  return(output)
}