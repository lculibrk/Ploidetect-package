#' @export
ploidetect <- function(all_data, normal = 2, tumour = 1, avg_allele_freq = 3, window_id = 4, window_size=5, GC = 6, plots = F, verbose = F, nomaf = F, lowest = NA, runCNAs = F, comp=NA, cndiff=NA, segmentation_threshold = 0.75, CNA_call = F, debugPlots = F){

  plots <- list()
  ## Run ploidetect_preprocess

  output <- ploidetect_preprocess(all_data = all_data, verbose = verbose, debugPlots = debugPlots, tumour = tumour, normal = normal, avg_allele_freq = avg_allele_freq, window_id = window_id, window_size = window_size, GC = GC)
  
  ## Unpack the output list from ploidetect_preprocess
  maxpeak <- output$maxpeak
  filtered <- output$x
  highoutliers <- output$highoutliers
  
  bw = maxpeak/40
  
  ## Perform peak calling with a heavily filtered "filtered" object to see how many we call:
  allPeaks <- peakcaller(filtered[findInterval(filtered$residual, vec = quantile(filtered$residual, probs = c(0.01, 0.99))) == 1,], bw)
  
  ## If we call zero peaks (due to poor bw, for example), try with half bandwidth and throw a warning
  if(nrow(allPeaks) == 1){
    warning("Zero peaks detected. Attempting peak calling with lower bandwidth")
    bw = bw/2
    allPeaks <- peakcaller(filtered[findInterval(filtered$residual, vec = quantile(filtered$residual, probs = c(0.01, 0.99))) == 1,], bw)
  }
  
  ## Now we peak call with a much more permissive quantile filter
  allPeaks <- peakcaller(filtered[findInterval(filtered$residual, vec = quantile(filtered$residual, probs = c(0.001, 0.999))) == 1,], bw)
  
  ## Center the residual and peak data about the tallest peak
  filtered$residual <- filtered$residual - allPeaks$pos[1]
  allPeaks$start <- allPeaks$start - allPeaks$pos[1]
  allPeaks$end <- allPeaks$end - allPeaks$pos[1]
  allPeaks$pos <- allPeaks$pos - allPeaks$pos[1]
  
  ## Run processallpeaks
  output <- ploidetect_processallpeaks(filtered, allPeaks)
  
  ## Unpack the output list from processallpeaks
  filtered <- output$filtered
  allPeaks <- output$allPeaks %>% filter(pos + maxpeak > 0)
  nomaf <- output$nomaf
  
  ## Generate coverage plots for interpretation
  
  filteredforplot <- filtered[findInterval(filtered$residual, vec = quantile(filtered$residual, probs = c(0.01, 0.99))) <= 1,]
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
  plots <- c(plots, list(plot))
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
  plots <- c(plots, list(plot))
  
  
  rerun = F
  if(!all(is.na(c(comp, cndiff)))){
    rerun=T
  }
  
  xdists <- ploidetect_roughmodels(allPeaks = allPeaks, maxpeak, verbose = verbose, rerun = rerun)
  
  ## Can't really do much else if there's only one peak in the data, so we return a message explaining this and exit
  if(nrow(allPeaks) == 1){
    return(list("TC_calls" = xdists, "plots" = NA, "CN_calls" = NA))
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

  for(i in 1:nrow(xdists)){
    modelbuilder_output <- modelbuilder_iterative(xdists[i,], allPeaks = allPeaks, lowest = lowest, filtered = filtered, strict = T, get_homd = F, mode = "TC", nomaf = nomaf, rerun = rerun, maxpeak = maxpeak, bw = bw)
    TC_calls <- c(TC_calls, list(modelbuilder_output$out))
    plots <- c(plots, list(modelbuilder_output$outplot))

  }

  #TC_calls <- lapply(xdists, function(x) modelbuilder_iterative(xdists = x, allPeaks = allPeaks, lowest = NA, filtered = filtered, strict = T, get_homd = F, mode = "TC", nomaf = nomaf, rerun = rerun, maxpeak = maxpeak, bw = bw))

  TC_calls <- do.call(rbind.data.frame, TC_calls)
  
  if(nrow(TC_calls) == 0){
    TC_calls <- list()
    plots <- list()
    for(i in 1:nrow(xdists)){
      modelbuilder_output <- modelbuilder_iterative(xdists[i,], allPeaks = allPeaks, lowest = lowest, filtered = filtered, strict = F, get_homd = F, mode = "TC", nomaf = nomaf, rerun = rerun, maxpeak = maxpeak, bw = bw)

      TC_calls <- c(TC_calls, list(modelbuilder_output$out))

      plots <- c(plots, list(modelbuilder_output$outplot))

    }
    TC_calls <- do.call(rbind.data.frame, TC_calls)
  }
  
  if(length(TC_calls) == 1){
    return(list("TC_calls" = TC_calls, "plots" = plots, "CN_calls" = NULL))
  }
  
  plots <- plyr::compact(plots)
  
  ordering <- order(TC_calls$model_error)
  
  TC_calls <- TC_calls[ordering,]
  
  plots <- plots[c(1, 2, ordering + 2)]
  
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
  
  result_segments <- CN_calls %>% group_by(chr, segment) %>% dplyr::summarise(pos = min(pos), end = max(end), residual = median(residual))
  
  segmentation_plot <- ggplot(CN_calls, aes(x = pos, y = raw_residual + maxpeak)) + geom_point(size = 0.5) + scale_color_viridis() + geom_segment(data = result_segments, mapping = aes(y = residual, yend = residual, x = pos, xend = end), size = 1, color = "red") + facet_wrap(~chr, ncol = 3) + ggtitle("Segmentation results")

  CNA_plot <- ggplot(CN_calls, aes(x = pos, y = raw_residual + maxpeak, color = CN)) + geom_point(size = 0.5) + scale_color_viridis() + facet_wrap(~chr, ncol = 3) + ggtitle("Copy number profile")

  plots <- c(plots, list(segmentation_plot), list(CNA_plot))
  
  return(list("TC_calls" = TC_calls, "plots" = plots, "CN_calls" = CN_calls))
}