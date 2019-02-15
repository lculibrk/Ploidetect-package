optimizesmoothing <- function(x, tumour = tumour, normal = normal, avg_allele_freq = avg_allele_freq, window_id = window_id, filtered, allPeaks, verbose = F, bw = bw){
  if(verbose){
    print("Beginning linear transformation step")
  }
  return(0)
  if(is.null(bw)){
    stop("BW is not set!")
  }
  ## Record number of peaks called
  npeaks <- nrow(allPeaks)
  
  ## Record the peak positions
  origpeakpos <- allPeaks$pos
  
  ## Record the max peak height under slope = 0
  zeroslope <- mean(allPeaks$meansig)
  
  ## Should be renamed; "troughs" here actually refers to the max peak height. Record all in a vector
  troughs <- c(zeroslope)
  ## Record all slopes used in a vector
  slopes <- c(0)
  ## Check if we should be going into negatives or positives for this by checking whether we have a better result at -1 or 1
  allPeaks <- peakcaller(linearTransform(x, tumour = tumour, normal = normal, avg_allele_freq = avg_allele_freq, window_id = window_id, slope = 1, verbose = F), bw = bw)
  slopes <- c(slopes, 1)
  troughs <- c(troughs, allPeaks$meansig[1])
  allPeaks <- peakcaller(linearTransform(x, tumour = tumour, normal = normal, avg_allele_freq = avg_allele_freq, window_id = window_id, slope = -1, verbose = F), bw = bw)
  slopes <- c(slopes, -1)
  troughs <- c(troughs, allPeaks$meansig[1])
  
  if(which.min(troughs[2:3]) == 1){
    slope = 1
  }else{
    slope = -1
  }
  
  if(verbose){
    if(slope > 0){
      print("Data appears to have a positive trend between germline mappability and read depth in tumour")
    }else{
      print("Data appears to have a negative trend between germline mappability and read depth in tumour")
    }
  }
  
  slopes <- c(0, slope)
  troughs <- c(zeroslope, troughs[which.min(troughs[2:3]) + 1])
  ## Start iterating at +-0.5
  slope <- slope/2
  for(i in 3:10){
    print(paste0("Testing slope offset=",slope))
    bestsofar <- slopes[which.min(troughs)]
    slopes[i] <- slope
    transformedData <- linearTransform(x, tumour = tumour, normal = normal, avg_allele_freq = avg_allele_freq, window_id = window_id, slope = slope, verbose = F)
    allPeaks <- peakcaller(transformedData, bw)
    troughs[i] <- allPeaks$meansig[1]
    slope <- mean(c(bestsofar, slope))
  }
  
  best = slopes[which.min(troughs)]
  
  if(verbose){
    print(paste0("Slope offset has been set to ", best))
  }
  return(best)
}