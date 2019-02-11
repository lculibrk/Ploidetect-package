ploidetect_segmentator <- function(filtered, matchedPeaks, maxpeak, predictedpositions, highoutliers, depthdiff, verbose = F, avg_allele_freq = avg_allele_freq, window_size = window_size, window_id = window_id, tumour = tumour, segmentation_threshold = segmentation_threshold){
  if(verbose){
    print("Performing segmentation and copy number variant calling")
  }
  # Start with all NA values for CN
  filtered$CN <- NA
  # Get chr and "pos" from the "window" field
  filtered$chr <- gsub(pattern = "_.*", replacement = "", filtered$wind)
  filtered$pos <- as.numeric(gsub(pattern = ".*_", replacement = "", filtered$wind))
  for(window in 1:nrow(filtered)){
    if(!is.na(filtered$peak[window]) & filtered$peak[window] %in% matchedPeaks$npeak){
      filtered$CN[window] <- matchedPeaks$CN[which(matchedPeaks$npeak == filtered$peak[window])]
    }
  }
  # Filter for points with maf values
  filterednomafna <- filtered %>% filter(!is.na(mafflipped))
  # Get end value for intervals
  filtered$end <- filtered$pos + filtered$size
  # Extract relevant data for CNA calling
  data <- filtered[,c("chr", "pos", "end", "mafflipped", "residual", "CN")]
  # Extract windows that are extreme outliers in coverage and include them in CNA calling (since they were filtered out in data pre-processing for TC/ploidy modeling)
  #highoutliers <- highoutliers %>% filter(V5 > 1000)
  if(nrow(highoutliers) > 0){
    # Populate chr, pos, end fields
    highoutliers$chr <- gsub(pattern = "_.*", replacement = "", highoutliers[,window_id])
    highoutliers$pos <- as.numeric(gsub(pattern = ".*_", replacement = "", highoutliers[,window_id]))
    highoutliers$end <- highoutliers$pos + highoutliers[,window_size]
    # Flip allele frequencies as we've done for most of the data
    highoutliers[,avg_allele_freq][highoutliers[,avg_allele_freq] == "."] <- NA
    highoutliers[,avg_allele_freq] <- as.numeric(highoutliers[,avg_allele_freq])
    highoutliers$mafflipped <- abs(highoutliers[,avg_allele_freq] - 0.5)+0.5
    # Populate the requisite fields to merge this with "data"
    highoutliers$residual <- highoutliers[,tumour]
    highoutliers$CN <- NA
    highoutliers <- highoutliers[,c("chr", "pos", "end", "mafflipped", "residual", "CN")]
    data <- rbind.data.frame(data, highoutliers)
  }
  ## Split data by chromosome
  datsplit <- split(data, data$chr)
  ## Generate model for regression-based CNA calling
  df_train <- data.frame("CN" = as.numeric(names(predictedpositions)), "median_segment" = predictedpositions, stringsAsFactors = F)
  train <- lm(CN ~ median_segment, data = df_train)
  ## Run segmentation by compression on all chromosomes
  if(verbose){
    print("Performing segmentation of copy number data")
  }
  compressedalldat <- lapply(datsplit, runiterativecompression, x = depthdiff, segmentation_threshold = segmentation_threshold)
  ## Change generate a "median_segment" column for median coverage per segment
  compressedalldat <- lapply(compressedalldat, function(x){
    x <- x %>% group_by(segment) %>% mutate("median_segment" = median(residual))
  })
  if(verbose){
    print("Calling copy numbers for each segment")
  }
  compressedalldat <- lapply(compressedalldat, function(x){
    x$median_segment <- x$median_segment + maxpeak
    x$CN <- round(predict(train, x), digits = 0)
    return(x)
  })
  ## Breakpoint calling from segmentation data
  if(verbose){
    print("Calling breakpoints between segments")
  }
  compressedalldat <- lapply(compressedalldat, function(t){
    t$breakpoint <- F
    for(window in 1:(nrow(t)-1)){
      if(t$segment[window] != t$segment[window+1]){
        if(t$CN[window] == t$CN[window+1]){
          next
        }
        diffs <- abs(c(min(abs(t$residual[window] - predictedpositions + maxpeak)), min(abs(t$residual[window+1] - predictedpositions + maxpeak))))
        decision <- which.max(diffs) - 1
        t$breakpoint[window + decision] <- T
      }
    }
    return(t)
  })
  compressedalldat <- do.call(rbind.data.frame, compressedalldat)
  CNAout <- compressedalldat
  names(CNAout)[which(names(CNAout) == "residual")] <- "raw_residual"
  names(CNAout)[which(names(CNAout) == "median_segment")] <- "residual"
  CNAout <- CNAout %>% arrange(chr, pos)
  if(verbose){
    print("Copy number analysis complete!")
  }
  return(CNAout)
}
