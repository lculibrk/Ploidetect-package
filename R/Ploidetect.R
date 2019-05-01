#' Runs Ploidetect on read count and allele frequency data
#' @param all_data input data, formatted as a data.frame where each row is a genomic bin. Columns correspond to
#' somatic read counts, normal read counts, allele frequencies, window id (in the format chr_pos), window size, and GC content
#' @param normal column index of normal read depth
#' @param tumour column index of tumour read depth
#' @param avg_allele_freq column index of snp allele frequencies in the tumour
#' @param window_id column index of window id (formatted as chr_pos)
#' @param window_size column index for size of the genomic bin in bp
#' @param GC column index of GC content as a percentage 
#' @param plots Logical. should plots be output?
#' @param verbose Logical. Print verbose messages?
#' @param lowest Integer from 0-2. Optional. Forces the copy number identity of the lowest fit peak in the read depth KDE
#' @param comp Integer. Forces selection of a specific comparator peak in the read depth kernel density estimate (KDE) for tumour content modeling
#' @param cndiff Integer. Forces the copy number difference between tallest peak in read depth KDE and the peak selected in comp
#' @param segmentation_threshold Float between 0 and 1. Threshold for segmentation of adjacent bins during copy number calling.
#' Default is 0.75. Smaller values will produce more segments, but may result in over-segmentation of copy number data.
#' @param CNA_call Logical. Should copy number calling be performed?
#' @param debugplots Logical. Should plots of data normalization be output?
#' @return A named list, containing three elements:
#' TC_calls: A data.frame containing all models considered for modeling tumour purity and ploidy
#' plots: a list of plots. Plot 1 is always a scatter plot where the x-axis is window size, y-axis is corrected somatic read counts, and
#' points are colored by SNP allele frequency. Plot 2 shows the KDE of the data, with all peaks shown. Next, the models in TC_calls are plotted,
#' one plot per model. Finally, if CNA_call was TRUE, copy number data is plotted, one plot per chromosome.
#' CNA_calls: a data.frame containing segmented copy number and LOH calls
#' @examples 
#' ## Run Ploidetect without specifying a model
#' ploidetect(all_data)
#' ## Run Ploidetect and force a model which fits the second most common copy number as being a true integer copy number 
#' ## with copy number equal to Ploidy +- 1, and the lowest common copy number being homozygous deletion
#' ploidetect(all_data, comp = 2, cndiff = 1, lowest = 0)
#' ## Run Ploidetect with the above model and also output copy number data
#' ploidetect(all_data, comp = 2, cndiff = 1, lowest = 0, CNA_call = T)
#' ## Run Ploidetect with a custom segmentation_threshold 
#' ploidetect(all_data, comp = 2, cndiff = 1, lowest = 0, CNA_call = T, segmentation_threshold = 0.5)
#' @export
ploidetect <- function(all_data, normal = 2, tumour = 1, avg_allele_freq = 3, window_id = 4, window_size=5, GC = 6, plots = F, verbose = F, lowest = NA, comp=NA, cndiff=NA, segmentation_threshold = 0.75, CNA_call = F, debugPlots = F){
  ## Initialize plots object
  plots <- list()
  
  ## Simplify_size
  simplify_size = 100000
  
  ## Load centromeres
  centromeres <- centromeres
  
  ## Run ploidetect_preprocess

  output <- ploidetect_preprocess(all_data = all_data, verbose = verbose, debugPlots = debugPlots, tumour = tumour, normal = normal, avg_allele_freq = avg_allele_freq, window_id = window_id, window_size = window_size, GC = GC, simplify = T, simplify_size = simplify_size)

  ## Unpack the output list from ploidetect_preprocess
  maxpeak <- output$maxpeak
  filtered <- output$x
  highoutliers <- output$highoutliers
  unaltered <- output$merged
  
  bw = maxpeak/40
  
  den <- density(filtered$residual, n = nrow(filtered), bw = bw)
  ## Normalize the density to 0->1 range
  den$y <- (den$y - min(den$y))/(max(den$y) - min(den$y))
  
  ## Sanity check - see if ploidetect_preprocess resulted in a worse separation of peaks than the initial data
  
  den2 <- density(filtered$y_raw, n = nrow(filtered), bw = bw)
  den2$y <- (den2$y - min(den2$y))/(max(den2$y) - min(den2$y))
  
  #den3 <- density(filtered$fan_correction, n = nrow(filtered), bw = bw)
  #den3$y <- (den3$y - min(den3$y))/(max(den3$y) - min(den3$y))
  
  ## If ploidetect_preprocess made things worse, then revert to raw input data
  
  decision <- which.min(c(mean(den$y), mean(den2$y)))#, mean(den3$y)))
  
  if(decision != 1){
    if(decision == 2){
      filtered$residual <- filtered$y_raw - maxpeak
    }
    #if(decision == 3){
    #  filtered$residual <- filtered$fan_correction
    #}
    if(verbose){
    }
  }
  
  
  
  ## Perform peak calling
  allPeaks <- peakcaller(filtered[findInterval(filtered$residual, vec = quantile(filtered$residual, probs = c(0, 0.999))) == 1,], bw, maxpeak = maxpeak)
  
  ## If we call zero peaks (due to poor bw, for example), try with half bandwidth and throw a warning
  if(nrow(allPeaks) == 1){
    warning("Zero peaks detected. Attempting peak calling with lower bandwidth")
    bw = bw/2
    allPeaks <- peakcaller(filtered[findInterval(filtered$residual, vec = quantile(filtered$residual, probs = c(0, 0.999))) == 1,], bw, maxpeak = maxpeak)
  }

  
  ## Run processallpeaks
  output <- ploidetect_processallpeaks(filtered, allPeaks)
  
  ## Unpack the output list from processallpeaks
  filtered <- output$filtered
  allPeaks <- output$allPeaks %>% filter(pos + maxpeak > 0)
  nomaf <- output$nomaf
  
  ## Generate coverage plots for interpretation
  
  filteredforplot <- filtered %>% filter(residual < max(allPeaks$pos) + maxpeak)
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
  den <- density(filteredforplot$residual, n = nrow(filteredforplot), bw = bw)
  dendf <- data.frame(x = den$x, y = den$y)
  # Normalize the density to 0->1 range and plot
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
    return(list("TC_calls" = xdists, "plots" = plots, "CN_calls" = NA))
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
    modelbuilder_output <- modelbuilder_iterative(xdists[i,], allPeaks = allPeaks, lowest = lowest, filtered = filtered, strict = T, get_homd = F, mode = "TC", nomaf = nomaf, rerun = rerun, maxpeak = maxpeak, bw = bw, verbose = verbose)
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
  
  if(!CNA_call){
    return(list("TC_calls" = TC_calls, "plots" = plots, "CN_calls" = NULL))
  }
  
  best_model = xdists[which(xdists$Comparator_peak_rank == TC_calls$Comparator[1] & xdists$Copy_number_difference_between_peaks == TC_calls$CN_diff[1]),]
  
  output <- modelbuilder_iterative(xdists = best_model[1,], allPeaks = allPeaks, lowest = TC_calls$lowest_peak_CN[1], filtered = filtered, strict = F, get_homd = F, mode = "CNA", nomaf = nomaf, rerun = rerun, maxpeak = maxpeak, bw = bw)
  
  predictedpositions <- output$predictedpositions
  
  matchedPeaks = output$matchedPeaks
  
  depthdiff <- output$depthdiff
  
  CN_calls <- ploidetect_segmentator(filtered, matchedPeaks, maxpeak, predictedpositions, highoutliers, depthdiff, avg_allele_freq = avg_allele_freq, window_size = window_size, window_id = window_id, tumour = tumour, segmentation_threshold = segmentation_threshold, verbose = verbose, GC = GC, all_data = all_data)
  
  #CN_calls %>% filter(chr == 20) %>%  ggplot(aes(x = pos,  y = raw_residual, color = segment)) + geom_point() + scale_color_viridis()
  
  new_CN_calls <- ploidetect_fineCNAs(all_data = all_data, CNAout = CN_calls, depthdiff = depthdiff, maxpeak = maxpeak, TC = TC_calls$tumour_purity[1], ploidy = TC_calls$Ploidy[1], verbose = verbose, tumour = tumour, normal = normal, avg_allele_freq = avg_allele_freq, window_id = window_id, window_size = window_size, GC = GC, decision = decision, simpsize = simplify_size, unaltered = unaltered)
  
  
  new_CN_calls <- do.call(rbind.data.frame, new_CN_calls)
  
  
  #new_CN_calls %>% ungroup %>%  filter(chr == 11, CN < 10) %>%  ggplot(aes(x = pos,  y = raw_residual, color = CN)) + geom_point() + scale_color_viridis()


  iters <- 2
  
  target_iters <- ceiling(log2(simplify_size/median(all_data$size)))
  #new_CN_calls %>% ungroup %>%  filter(chr == "17", CN < 10) %>%  ggplot(aes(x = pos, y = raw_residual, color = CN)) + geom_point() + scale_color_viridis()
  while(iters < (target_iters)){
    initial_points <- nrow(new_CN_calls)
    new_CN_calls <- ploidetect_fineCNAs(all_data = all_data, CNAout = new_CN_calls, depthdiff = depthdiff, maxpeak = maxpeak, TC = TC_calls$tumour_purity[1], ploidy = TC_calls$Ploidy[1], verbose = verbose, tumour = tumour, normal = normal, avg_allele_freq = avg_allele_freq, window_id = window_id, window_size = window_size, GC = GC, decision = decision, simpsize = simplify_size/(2^(iters-1)), unaltered = unaltered)
    new_CN_calls <- do.call(rbind.data.frame, new_CN_calls)
    print(paste0("Completed iteration ", iters))
    iters <- iters+1
    if(nrow(new_CN_calls) == initial_points){
      iters <- Inf
    }
  }
  if(verbose){
    print("Calling LOH")
  }
  CN_calls <- ploidetect_loh(purity = TC_calls$tumour_purity[1], CNA_object = new_CN_calls)
  if(verbose){
    print("Finished calling LOH")
  }
  CN_calls$state <- factor(CN_calls$state)
  
  CN_palette <- c("0" = "#cc0000", 
                  "1" = "#000066", 
                  "2" = "#26d953", 
                  "3" = "#609f70", 
                  "4" = "#cccc00", 
                  "5" = "#80804d",
                  "6" = "#cc6600", 
                  "7" = "#856647", 
                  "8" = "#cc0000"
                  )
  
  CN_calls <- split(CN_calls, f = CN_calls$chr)
  
  CNA_plot <- lapply(CN_calls, function(x){
    chr = x$chr[1]
    x %>% filter(end < centromeres$pos[which(centromeres$chr %in% chr)[1]] | pos > centromeres$end[which(centromeres$chr %in% chr)[2]]) %>% ggplot(aes(x = pos, y = log(raw_residual + maxpeak), color = as.character(state))) + 
      geom_point(size = 0.5) + 
      scale_color_manual(name = "State",
                         values = CN_palette, 
                         labels = c("0" = "HOMD", 
                                    "1" = "CN = 1", 
                                    "2" = "CN = 2 HET", 
                                    "3" = "CN = 2 HOM", 
                                    "4" = "CN = 3 HET", 
                                    "5" = "CN = 3 HOM", 
                                    "6" = "CN = 4 HET", 
                                    "7" = "CN = 4 HOM", 
                                    "8" = "CN = 5+")) + 
      ylab("log(Read Depth)") + 
      xlab("position") + 
      ggtitle(paste0("Chromosome ", chr, " copy number profile")) + 
      theme_bw()
  })
  
  vaf_plot <- lapply(CN_calls, function(x){
    chr = x$chr[1]
    x %>% filter(end < centromeres$pos[which(centromeres$chr %in% chr)][1] | pos > centromeres$end[which(centromeres$chr %in% chr)][2]) %>% filter(!is.na(mafflipped)) %>% ggplot(aes(x = pos, y = mafflipped, color = as.character(state))) + 
      geom_point(size = 0.5) + 
      scale_color_manual(name = "State",
                         values = CN_palette, 
                         labels = c("0" = "HOMD", 
                                    "1" = "CN = 1", 
                                    "2" = "CN = 2 HET", 
                                    "3" = "CN = 2 HOM", 
                                    "4" = "CN = 3 HET", 
                                    "5" = "CN = 3 HOM", 
                                    "6" = "CN = 4 HET", 
                                    "7" = "CN = 4 HOM", 
                                    "8" = "CN = 5+")) + 
      ylab("Major allele frequency") + 
      xlab("position") + 
      ggtitle(paste0("Chromosome ", chr, " allele frequency profile")) + 
      theme_bw()
  })
  
  cna_plots <- list()
  
  for(i in 1:length(CNA_plot)){
    cna_plots[i] <- list(plot_grid(CNA_plot[[i]], vaf_plot[[i]], align = "v", axis = "l", ncol = 1))
  }
  
  CN_calls <- do.call(rbind.data.frame, CN_calls)
  
  plots <- c(plots, cna_plots)
  
  return(list("TC_calls" = TC_calls, "plots" = plots, "CN_calls" = CN_calls))
}
