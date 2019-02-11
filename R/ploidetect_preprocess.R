ploidetect_preprocess <- function(all_data, normal = 2, tumour = 1, avg_allele_freq = 3, window_id = 4, window_size = 5, GC = 6, debugPlots = F, verbose = F){
  if(verbose){
    print("Ploidetect - Detection of tumour purity and aneuploidy from whole-genome sequence data")
    print("Thank you for using Ploidetect! Please remember to cite this tool if used in your research ^_^")
    print("Beginning data pre-processing steps")
  }
  # Load data
  x <- as.data.frame(all_data)
  
  
  # Test if data is configured and input properly
  if(!all(is.numeric(x[,normal]))){
    stop("At least one element in normals column is not numeric. Did you not specify the indices correctly/coerce the values to numeric?")
  }
  if(!all(is.numeric(x[,tumour]))){
    stop("At least one element in tumour column is not numeric. Did you not specify the indices correctly/coerce the values to numeric?")
  }
  if(!all(is.numeric(x[,window_size]))){
    stop("At least one element in window_size column is not numeric. Did you not specify the indices correctly/coerce the values to numeric?")
  }
  if(!grepl("[0-9]*_[0-9]*", x[1,window_id])){
    stop("Values in window_id column do not appear to match expected format: chrnumber_start")
  }
  
  
  ## The very first thing we do is measure the read depth at the highest density of read coverage
  maxpeak <- density(x[,tumour], n = nrow(x))$x[which.max(density(x[,tumour], n = nrow(x))$y)]
  
  ## Filter for chr1-22 and X for this stage of modeling
  x <- x[grepl(pattern = paste0(c(paste0("^", 1:22), "^X"), "_", collapse = "|"), x = x[,window_id]),]
  
  
  ## Set row names to window_ids
  row.names(x) <- as.data.frame(x)[,window_id]
  
  ## Perform basic pre-filtering, find the windows within the 90th percentile of tumour read counts
  rangedf <- x[findInterval(x[,tumour], 
                            quantile(x[,tumour], 
                                     probs=c(0,0.90))) == 1,]
  
  ## Obtain the range of values observed within the 90th percentile of tumour read counts
  range <- range(rangedf[,tumour])[2] - range(rangedf[,tumour])[1]
  
  ## Maximum allowable read depth for retention is the 90th percentile + the above computed range.
  max <- rangedf$V1[2] + range
  min <- 0
  
  ## Set outliers aside for later steps (these will be CNA called later, but are exempt from preprocessing steps)
  highoutliers <- x[findInterval(as.data.frame(x)[,tumour], c(min, max)) > 1,]
  
  ## Filter data for everything within the read depth range
  x <- x[findInterval(as.data.frame(x)[,tumour], c(min, max)) == 1,]
  
  # Extract X chromosome regions
  
  chrX <- x[grep("X_", x[,window_id]),]
  isMale=F
  if(which.min(abs((median(chrX[,window_size]) - median(x[,window_size])) - c(0, median(x[,window_size])))) == 2){
    x$V5[grep("X_", x[,window_id])] <- x$V5[grep("X_", x[,window_id])]/2
    isMale=T
  }
  if(verbose){
    if(isMale){
      print("Automated sex detection infers MALE sex")
    }else{
      print("Automated sex detection infers FEMALE sex")
    }
  }
  
  ## This is a very broad-strokes filtering step for window size. Basically removing extreme outliers w.r.t germline mappability, as we don't want to use these in modeling
  x <- x[(x[,window_size] > (median(x[,window_size])/10)) & (x$V5 < (median(x[,window_size])*5)),]
  
  
  x <- x[,c(tumour, normal, window_id, avg_allele_freq, window_size, GC)]
  names(x) <- c("y_raw", "x_raw", "window", "maf", "size", "GC")

  print(str(x))

  if(debugPlots){
    rawPlot <- x %>% ggplot(aes(x = size, y = y_raw)) + geom_point(size = 0.1, alpha = 0.1) + xlab("Window size") + ylab("Tumour Read counts") + ggtitle("Raw counts by window size") + theme_minimal()
    print(rawPlot)
  }
  
  ## Set names of dataframe
  
  
  if(debugPlots){
    GCplot <- ggplot(x, aes(x = GC * 100, y = y_raw)) + geom_point(size = 0.3, alpha = 0.2) + theme_bw() + xlab("GC Content %") + ylab("Tumour Read Counts") + ggtitle("Read count and GC content relationship")
    print(GCplot)  
  }
  
  ## Normalize for GC content
  GCnorm <- loessFit(y = x$y_raw, x = x$GC, span = 0.1)
  
  ## Residuals of this model can be thought of as the normalized read counts, so scale them back to something resembling read counts
  GCnorm$residuals <- GCnorm$residuals + median(x$y_raw)
  
  if(debugPlots){
    GCnormplot <- ggplot(x, aes(y = GCnorm$residuals, x = size)) + 
      geom_point(size = 0.3, alpha = 0.1) + 
      theme_bw() + 
      xlab("Window size") + 
      ylab("Normalized tumour read counts") + 
      ggtitle("Non-linear normalization of read counts by GC content")
    print(GCnormplot)
  }
  
  ## Perform a linear model fit so as to get a linear transformation of the data w.r.t. the fitted values of the loess fit
  
  #GCnorm2 <- rlm(x$y_raw ~ GCnorm$fitted)
  
  ## Normalized final values: Linear-normalized read depth divided by expected response of the linear model, then scaled to a practical scale
  ## This scales the values such that they are linear (ie. removes any x~y relationship and ensures the variation is only in the y-axis)
  
  
  x$normalized_tumour <- GCnorm$residuals
  
  ## We use window_size as a variable to represent the germline mappability in later steps
  x$normalized_size <- x$size
  
  if(debugPlots){
    GCnorm2plot <- ggplot(x, aes(y = normalized_tumour, x = normalized_size)) + 
      geom_point(size = 0.3, alpha = 0.1) + 
      theme_bw() + 
      xlab("Window size") + 
      ylab("Normalized tumour read counts") + 
      ggtitle("Linear normalization of read counts by GC content")
    print(GCnorm2plot)
  }
  output <- list("x" = x, "maxpeak" = maxpeak, "highoutliers" = highoutliers)
  return(output)
}