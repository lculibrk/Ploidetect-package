ploidetect_segmentator <- function(filtered, matchedPeaks, maxpeak, predictedpositions, highoutliers, depthdiff, verbose = F, avg_allele_freq = avg_allele_freq, window_size = window_size, window_id = window_id, tumour = tumour, GC = GC, segmentation_threshold = segmentation_threshold, all_data = all_data){
  if(verbose){
    print("Performing segmentation and copy number variant calling")
  }
  # Start with all NA values for CN
  filtered$CN <- NA
  # Get chr and "pos" from the "window" field

  #filtered$chr <- gsub(pattern = "_.*", replacement = "", filtered$wind)
  #filtered$pos <- as.numeric(gsub(pattern = ".*_", replacement = "", filtered$wind))
  #for(window in 1:nrow(filtered)){
  #  if(!is.na(filtered$peak[window]) & filtered$peak[window] %in% matchedPeaks$npeak){
  #    filtered$CN[window] <- matchedPeaks$CN[which(matchedPeaks$npeak == filtered$peak[window])]
  #  }
  #}
  # Filter for points with maf values
  filterednomafna <- filtered %>% filter(!is.na(mafflipped))
  # Get end value for intervals
  #filtered$end <- filtered$pos + filtered$size
  # Extract relevant data for CNA calling
  data <- filtered[,c("chr", "pos", "end", "mafflipped", "residual", "CN")]
  data$residual <- data$residual + maxpeak
  # Extract windows that are extreme outliers in coverage and include them in CNA calling (since they were filtered out in data pre-processing for TC/ploidy modeling)
  #highoutliers <- highoutliers %>% filter(V5 > 1000)
  if(nrow(highoutliers) > 0){
    # Populate chr, pos, end fields
    #highoutliers$chr <- gsub(pattern = "_.*", replacement = "", highoutliers$wind)
    #highoutliers$pos <- as.numeric(gsub(pattern = ".*_", replacement = "", highoutliers$wind))
    #highoutliers$end <- highoutliers$pos + highoutliers$size
    # Flip allele frequencies as we've done for most of the data
    highoutliers$mafflipped <- abs(highoutliers$maf - 0.5)+0.5
    # Populate the requisite fields to merge this with "data"
    highoutliers$residual <- highoutliers[,tumour] + maxpeak
    highoutliers$CN <- NA
    highoutliers <- highoutliers[,c("chr", "pos", "end", "mafflipped", "residual", "CN")]
    data <- rbind.data.frame(data, highoutliers)

  }
  ## Split data by chromosome
  datsplit <- split(data, data$chr)
  ## Generate model for regression-based CNA calling
  df_train <- data.frame("CN" = names(predictedpositions), "median_segment" = predictedpositions, stringsAsFactors = T)
  train <- lm(CN ~ median_segment, data = df_train)
  ## Run segmentation by compression on all chromosomes
  if(verbose){
    print("Performing segmentation of copy number data")
  }
  compressedalldat <- lapply(datsplit, runiterativecompression, x = depthdiff, segmentation_threshold = segmentation_threshold, verbose = verbose)
  
  ## Generate a "median_segment" column for median coverage per segment
  compressedalldat <- lapply(compressedalldat, function(x){
    x <- x %>% group_by(segment) %>% dplyr::mutate("median_segment" = median(residual), "median_maf" = median(mafflipped, na.rm = T))
  })
  if(verbose){
    print("Calling copy numbers for each segment")
  }
  compressedalldat <- lapply(compressedalldat, function(x){
    x$median_segment <- x$median_segment
    x$CN <- round(predict(train, x), digits = 1)
    return(x)
  })
  ## Breakpoint calling from segmentation data
  if(verbose){
    print("Calling breakpoints between segments")
  }
  compressedalldat <- lapply(compressedalldat, callbreakpoints, predictedpositions = predictedpositions, maxpeak = maxpeak)
  
  CNAout <- do.call(rbind.data.frame, compressedalldat)
  names(CNAout)[which(names(CNAout) == "residual")] <- "raw_residual"
  names(CNAout)[which(names(CNAout) == "median_segment")] <- "residual"
  CNAout <- CNAout %>% arrange(chr, pos)
  if(verbose){
    print("Copy number analysis complete!")
  }
  return(CNAout)
}
