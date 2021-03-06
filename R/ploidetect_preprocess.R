ploidetect_preprocess <- function(all_data, normal = 2, tumour = 1, avg_allele_freq = 3, window_id = 4, window_size = 5, GC = 6, debugPlots = F, verbose = F, simplify = T, simplify_size = 500000){
  if(verbose & simplify){
    print("Ploidetect - Detection of tumour purity and aneuploidy from whole-genome sequence data")
    print("Thank you for using Ploidetect! Please remember to cite this tool if used in your research ^_^")
    print("Beginning data pre-processing steps")
  }
  # Load data
  x <- as.data.frame(all_data)
  
  # Process centromere data
  centromeres_preprocess <- centromeres %>% group_by(chr) %>% dplyr::summarise(pos = first(pos), end = last(end))
  #if(any("." %in% x$maf)){
  #  x$maf[which(x$maf == ".")] <- NA
  #  x$maf <- as.numeric(x$maf)
  #}
  # Test if data is configured and input properly
  if(any(grepl("chr", x$chr))){
    stop("Expected numeric chromosomes (1, 2, 3, ... X), not chr1, chr2, etc")
  }
  if(!all(is.numeric(x$normal))){
    stop("At least one element in normals column is not numeric.")
  }
  if(!all(is.numeric(x$tumour))){
    stop("At least one element in tumour column is not numeric.")
  }
  #if(!all(is.numeric(x$maf) | is.na(x$maf))){
  #  stop("VAF column must contain only numeric or NA values")
  #}
  if(!all(is.numeric(x$gc))){
    stop("At least one element in GC-content column is not numeric")
  }
  

  
  ## Step 1: Merge data such that window size is about 100Kb
  
  ## Filter for chr1-22 and X
  
  #x <- x[grepl(pattern = paste0(c(paste0("^", 1:22), "^X"), "_", collapse = "|"), x = x[,window_id]),]
  
  if(verbose){
    print("Filtering for chromosomes 1-22 and X")
  }
  x <- x %>% filter(chr %in% paste0(c(1:22, "X")))
  x <- split(x, x$chr)
  centromeres_split <- split(centromeres_preprocess, centromeres_preprocess$chr)
  
  if(verbose){
    print("Filtering out centromeric loci")
  }
  
  x <- lapply(x, function(k){
    chr <- k$chr[1]
    centro_start <- centromeres_split[[chr]]$pos %>% unlist()
    centro_end <- centromeres_split[[chr]]$end %>% unlist()
    #print(str(k))
    k <- k %>% filter(end < centro_start | pos > centro_end)
    #print(str(k))
    return(k)
  })
  if(verbose){
    print("Completed centromere filtering")
  }
  x <- do.call(rbind.data.frame, x) %>% arrange(chr, pos)
  
  x$window_size <- x$end - x$pos
  
  
  #mean_size <- mean(x[,window_size])
  mean_size <- mean(x$window_size)
  
  closest <- round(simplify_size/mean_size, digits = 0)
  if(simplify){
    # Compute the closest integer multiple of the window size to 100kb
    if(closest < 1){
      closest = 1
    }
    # Create a merge vector
    x$merge <- floor(seq(from = 0, by = 1/closest, length.out = nrow(x)))
  }else{x$merge <- 1:nrow(x)}
  
  # Process the allele frequency column into a numeric vector
  
  #x[,avg_allele_freq][which(x[,avg_allele_freq] == ".")] <- NA
  x$maf[which(x$maf == ".")] <- NA
  
  
  #x[,avg_allele_freq] <- as.numeric(x[,avg_allele_freq])
  #x$maf <- as.numeric(x$maf)
  
  #x$chr <- gsub("_.*", "", x[,window_id])
  
  # Sanitize the column names
  
  #x <- x[,c(tumour, normal, avg_allele_freq, window_id, window_size, GC, 7, 8)]
  
  #names(x) <- c("tumour", "normal", "maf", "wind", "size", "gc", "merge", "chr")
  #x <- x %>% group_by(merge, chr) %>% summarise(tumour = sum(tumour), normal = sum(normal), maf = merge_mafs(maf, na.rm = T), wind = dplyr::first(wind), size = sum(size), gc = mean(gc))
  x <- x %>% group_by(merge, chr) %>% dplyr::summarise(pos = first(pos), end = last(end), tumour = sum(tumour), normal = sum(normal), maf = merge_mafs(maf, na.rm = T, exp = T), gc = mean(gc), window_size = sum(window_size))
  
  ## Measure the read depth at the highest density of read coverage
  maxpeak <- density(x$tumour, n = nrow(x))$x[which.max(density(x$tumour, n = nrow(x))$y)]
  


  ## Remove tibble formatting
  x <- as.data.frame(x)
  
  ## Set row names to window_ids
  #row.names(x) <- x$wind

  ## Get median normal coverage
  median_normal <- median(x$normal)
  
  norm_factor <- x$normal/median_normal
  
  x$tumour <- x$tumour/norm_factor
  

  ## Perform basic pre-filtering, find the windows within the 90th percentile of tumour read counts
  rangedf <- x[findInterval(x$tumour, 
                            quantile(x$tumour, 
                                     probs=c(0,0.90))) == 1,]
  
  ## Obtain the range of values observed within the 90th percentile of tumour read counts
  range <- range(rangedf$tumour)[2] - range(rangedf$tumour)[1]
  
  ## Maximum allowable read depth for retention is the 90th percentile + the above computed range.
  max <- range(rangedf$tumour)[2] + range
  min <- 0
  
  ## Set outliers aside for later steps (these will be CNA called later, but are exempt from TC/Ploidy analysis)
  highoutliers <- x[findInterval(x$tumour, c(min, max)) > 1,]
  
  ## Filter data for everything within the read depth range
  x <- x[findInterval(x$tumour, c(min, max)) == 1,]
  
  # Extract X chromosome regions
  
  chrX <- x[x$chr == "X",]
  isMale=F
  if(which.min(abs((median(chrX$window_size) - median(x$window_size)) - c(0, median(x$window_size)))) == 2){
    x$normal[x$chr == "X"] <- x$normal[x$chr == "X"]/2
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
  x <- x[(x$window_size > (median(x$window_size)/10)) & (x$window_size < (median(x$window_size)*5)),]
  
  
  #x <- x[,c(tumour, normal, window_id, avg_allele_freq, window_size, GC)]
  #names(x) <- c("y_raw", "x_raw", "window", "maf", "size", "GC")
  
  #x <- x[,3:8]
  
  #names(x) <- c("y_raw", "x_raw", "maf", "window", "size", "GC")
  x <- x %>% dplyr::rename("y_raw" = "tumour", "x_raw" = "normal")
  
  if(debugPlots){
    rawPlot <- x %>% ggplot(aes(x = window_size, y = y_raw)) + geom_point(size = 0.1, alpha = 0.1) + xlab("Window size") + ylab("Tumour Read counts") + ggtitle("Raw counts by window size") + theme_minimal() + 
      theme(
        plot.title = element_text(size = 20),
        plot.caption = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)
      )
    print(rawPlot)
  }
  
  ## Set names of dataframe
  
  
  if(debugPlots){
    GCplot <- ggplot(x, aes(x = gc * 100, y = y_raw)) + geom_point(size = 0.3, alpha = 0.2) + theme_bw() + xlab("GC Content %") + ylab("Tumour Read Counts") + ggtitle("Read count and GC content relationship") + 
      theme(
        plot.title = element_text(size = 20),
        plot.caption = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)
      )
    print(GCplot)  
  }
  
  ## Normalize for GC content
  GCnorm <- loessFit(y = x$y_raw, x = x$gc, span = 0.75)
  
  ## Residuals of this model can be thought of as the normalized read counts, so scale them back to something resembling read counts
  GCnorm$residuals <- GCnorm$residuals
  
  if(debugPlots){
    GCnormplot <- ggplot(x, aes(y = GCnorm$residuals, x = window_size)) + 
      geom_point(size = 0.3, alpha = 0.1) + 
      theme_bw() + 
      xlab("Window size") + 
      ylab("Normalized tumour read counts") + 
      ggtitle("Non-linear normalization of read counts by GC content") + 
      theme(
        plot.title = element_text(size = 20),
        plot.caption = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)
      )
    print(GCnormplot)
  }
  
  
  ## Normalized final values: Linear-normalized read depth divided by expected response of the linear model, then scaled to a practical scale
  ## This scales the values such that they are linear (ie. removes any x~y relationship and ensures the variation is only in the y-axis)
  
  x$residual <- GCnorm$residuals
  
  ## We use window_size as a variable to represent the germline mappability in later steps
  x$normalized_size <- x$window_size
  
  ## Experimental scaling:
  
  #x$fan_correction <- ((x$residual/x$size) * median(x$size)) + median(x$size)
  
  
  if(debugPlots){
    GCnorm2plot <- ggplot(x, aes(y = residual, x = normalized_size)) + 
      geom_point(size = 0.3, alpha = 0.1) + 
      theme_bw() + 
      xlab("Window size") + 
      ylab("Normalized tumour read counts") + 
      ggtitle("Linear normalization of read counts by GC content")
    print(GCnorm2plot)
  }
  output <- list("x" = x, "maxpeak" = maxpeak, "highoutliers" = highoutliers, "merged" = closest)
  return(output)
}
