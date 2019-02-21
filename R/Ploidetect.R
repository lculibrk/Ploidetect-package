#' @export
ploidetect <- function(all_data, normal = 2, tumour = 1, avg_allele_freq = 3, window_id = 4, window_size=5, GC = 6, limited = F, top = Inf, plots = F, verbose = F, nomaf = F, bw = 800, lowest = NA, runCNAs = F, comp=NA, cndiff=NA, segmentation_threshold = 0.75, CNA_call = F, debugPlots = F){
  if(!is.numeric(bw)){
    stop("Bandwidth must be numeric!")
  }
  
  ## Run ploidetect_preprocess
  
  output <- ploidetect_preprocess(all_data = all_data, verbose = verbose, debugPlots = debugPlots, tumour = tumour, normal = normal, avg_allele_freq = avg_allele_freq, window_id = window_id, window_size = window_size, GC = GC)
  
  ## Unpack the output list from ploidetect_preprocess
  maxpeak <- output$maxpeak
  x <- output$x
  highoutliers <- output$highoutliers
  
  bw = maxpeak/25
  
  ## Run ploidetect_transform
  output <- ploidetect_transform(x, bw, verbose = verbose, tumour = tumour, normal = normal, avg_allele_freq = avg_allele_freq, window_id = window_id, window_size = window_size)
  
  ## Unpack the output list from ploidetect_transform
  allPeaks <- output$allPeaks
  filtered <- output$filtered 
  ## Run processallpeaks
  output <- ploidetect_processallpeaks(filtered, allPeaks)
  
  ## Unpack the output list from processallpeaks
  filtered <- output$filtered
  allPeaks <- output$allPeaks
  nomaf <- output$nomaf
  ## Generate coverage plots for interpretation
  
  filteredforplot <- filtered[findInterval(filtered$residual, vec = quantile(filtered$residual, probs = c(0.001, 0.999))) <= 1,]
  filteredforplot$residual <- filteredforplot$residual + maxpeak
  plot <- ggplot(data = filteredforplot, mapping = aes_string(x = "size", y = "residual", color = "mafflipped")) + geom_point(size = 0.1, alpha = 0.1) +
    #xlab("Window size") + 
    xlab("Normalized window size of constant-coverage bins") + 
    #ylab("Reads mapping to bins in Somatic") + 
    ylab("Residuals of reads mapping to constant-coverage bins in Somatic") + 
    #ggtitle(paste0("Coverage vs normalized bin size (Surrogate for Mappability + GC bias)")) + 
    ggtitle(paste0("Coverage Plot for Filtered and Normalised Data")) + 
    scale_colour_viridis(option = "plasma", name = "Major Allele\n Frequency") +
    #scale_x_continuous(limits = quantile(filtered$size, probs = c(0.05, 0.99))) +
    theme_bw(base_size = 12)
  print(plot)
  den <- density(filteredforplot$residual)
  # Normalize the density to 0->1 range
  den$y <- (den$y - min(den$y))/(max(den$y) - min(den$y))
  plot <- ggplot(mapping = aes(x = den$x, y = den$y)) + geom_line() + 
    xlab("Residuals of Reads Mapping to 50kb windows in Somatic") + 
    ylab("Normalized Density") + 
    ggtitle(paste0("Kernel Density Estimate of Residual Data")) + 
    geom_vline(aes(xintercept = allPeaks$pos + maxpeak), linetype = 2) + 
    geom_text(mapping = aes(x = allPeaks$pos + maxpeak, y = allPeaks$height + 0.05, label = paste0("MAF = ", round(allPeaks$mainmaf, digits = 3)))) +
    theme_bw()
  print(plot)
  
  
  rerun = F
  if(!all(is.na(c(comp, cndiff)))){
    rerun=T
  }
  
  xdists <- ploidetect_roughmodels(allPeaks = allPeaks, maxpeak, verbose = verbose, rerun = rerun)
  
  ## Can't really do much else if there's only one peak in the data, so we return a message explaining this and exit
  if(nrow(allPeaks) == 1){
    return(xdists)
  }
  
  ## Normalize allPeaks to maxPeak
  allPeaks$pos <- allPeaks$pos + maxpeak
  
  ## If this is a "rerun" case, ie. if the user is specifying the model for Ploidetect to use
  
  
  if(rerun){
    xdists <- xdists %>% filter(Comparator_peak_rank == comp & Copy_number_difference_between_peaks == cndiff)
  }
  if(nrow(xdists) == 0 & rerun){
    stop("You are trying to force models which do not exist! Check the model information you entered and try again")
  }
  
  TC_calls <- list()
  plots <- list()
  for(i in 1:nrow(xdists)){
    modelbuilder_output <- modelbuilder_iterative(xdists[i,], allPeaks = allPeaks, lowest = NA, filtered = filtered, strict = T, get_homd = F, mode = "TC", nomaf = nomaf, rerun = rerun, maxpeak = maxpeak, bw = bw)
    TC_calls <- c(TC_calls, list(modelbuilder_output$out))
    plots <- c(plots, list(modelbuilder_output$outplot))
  }

  #TC_calls <- lapply(xdists, function(x) modelbuilder_iterative(xdists = x, allPeaks = allPeaks, lowest = NA, filtered = filtered, strict = T, get_homd = F, mode = "TC", nomaf = nomaf, rerun = rerun, maxpeak = maxpeak, bw = bw))
  
  TC_calls <- do.call(rbind.data.frame, TC_calls)
  
  plots <- plyr::compact(plots)
  
  ordering <- order(TC_calls$newerr)
  
  TC_calls <- TC_calls[ordering,]
  
  plots <- plots[ordering]
  
  #TC_calls <- do.call(rbind.data.frame, TC_calls) %>% arrange(newerr)
  
  if(!CNA_call){
    return(list("TC_calls" = TC_calls, "plots" = plots, "CN_calls" = NULL))
  }
  
  best_model = xdists[which(xdists$Comparator_peak_rank == TC_calls$Comparator[1] & xdists$Copy_number_difference_between_peaks == TC_calls$CN_diff[1]),]
  
  output <- modelbuilder_iterative(xdists = best_model[1,], allPeaks = allPeaks, lowest = TC_calls$lowest_peak_CN[1], filtered = filtered, strict = F, get_homd = F, mode = "CNA", nomaf = nomaf, rerun = rerun, maxpeak = maxpeak, bw = bw)
  
  predictedpositions <- output$predictedpositions
  
  matchedPeaks = output$matchedPeaks
  
  depthdiff <- output$depthdiff
  
  CN_calls <- ploidetect_segmentator(filtered, matchedPeaks, maxpeak, predictedpositions, highoutliers, depthdiff, avg_allele_freq = avg_allele_freq, window_size = window_size, window_id = window_id, tumour = tumour, segmentation_threshold = segmentation_threshold)
  
  return(list("TC_calls" = TC_calls, "CN_calls" = CN_calls))
}