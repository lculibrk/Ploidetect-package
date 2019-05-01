ploidetect_fineCNAs <- function(all_data, CNAout, TC, ploidy, depthdiff = depthdiff, maxpeak = maxpeak, verbose = verbose, tumour = tumour, normal = normal, avg_allele_freq = avg_allele_freq, window_id = window_id, window_size = window_size, GC = GC, decision = decision, simpsize = simpsize, unaltered = unaltered){
  ## Generate processed data.frame for unmerged data
  unmerged_data <- ploidetect_preprocess(all_data = all_data, verbose = verbose, debugPlots = F, tumour = tumour, normal = normal, avg_allele_freq = avg_allele_freq, window_id = window_id, window_size = window_size, GC = GC, simplify = T, simplify_size = simpsize/2)

  den <- density(unmerged_data$x$residual, n = nrow(unmerged_data$x$residual))
  offset <- den$x[which.max(den$y)]
  unmerged_data$x$residual <- unmerged_data$x$residual - offset
  
  ## Unpack maxpeak
  unmerged_maxpeak <- unmerged_data$maxpeak
  
  ## Compute reads-per-copy and HOMD location for unmerged data based on merged data
  unmerged_diff <- depthdiff/(unaltered/unmerged_data$merged)
  unmerged_normalreads <- unmerged_maxpeak - ploidy*unmerged_diff
  
  predictedpositions <- seq(from = unmerged_normalreads, by = unmerged_diff, length.out = 11)
  names(predictedpositions) <- 0:10
  
  df.train <- data.frame("CN" = 0:10, "median_segment" = predictedpositions)
  model <- lm(CN ~ median_segment, data = df.train)
  #print(predictedpositions)
  
  ## Continue unpacking data
  unmerged_highoutliers <- unmerged_data$highoutliers
  unmerged_data <- unmerged_data$x
  if(decision == 2){
    unmerged_data$residual <- unmerged_data$y_raw - unmerged_maxpeak
  }
  
  ## Sanitise highoutliers and merge with rest of data
  unmerged_highoutliers <- unmerged_highoutliers[,c("tumour", "normal", "maf", "wind", "size", "gc", "tumour", "size")]
  names(unmerged_highoutliers) <- names(unmerged_data)
  unmerged_data$residual <- unmerged_data$residual + unmerged_maxpeak
  unmerged_data <- rbind.data.frame(unmerged_data, unmerged_highoutliers)
  
  ## Compute the positional columns (chr, pos, end) for each window
  unmerged_data$chr <- gsub("_.*", "", unmerged_data$window)
  unmerged_data$pos <- as.numeric(gsub(".*_", "", unmerged_data$window))
  unmerged_data$end <- unmerged_data$pos + unmerged_data$size
  unmerged_data <- unmerged_data %>% arrange(chr, pos)
  
  

  
  
  ## Split both merged and unmerged data by chr
  unmerged_data <- split(unmerged_data, f = unmerged_data$chr)
  merged_data <- split(CNAout, f = CNAout$chr)
  
  ## First map old CNs to new higher-res data
  for(chr in names(unmerged_data)){
    unmerged_chr <- unmerged_data[[chr]] %>% arrange(pos)
    merged_chr <- merged_data[[chr]]
    merged_segments <- merged_chr %>% group_by(chr, segment) %>% summarise(pos = dplyr::first(pos), end = last(end), CN = mean(CN))
    merged_segments$pos[1] <- 0
    unmerged_chr$segment <- findInterval(unmerged_chr$pos, merged_segments$pos)
    unmerged_chr <- left_join(unmerged_chr, merged_segments[,c("segment", "CN")], by = "segment")
    unmerged_data[[chr]] <- unmerged_chr
  }
  
  #unmerged_data$`11` %>% ggplot(aes(x = pos,  y = residual, color = segment)) + geom_point() + scale_color_viridis()

  ## Compute the standard deviation of read depth in the 50% longest segments
  grouped_data <- do.call(rbind.data.frame, unmerged_data)
  sd <- grouped_data %>% group_by(chr, segment) %>% summarise("sd" = sd(residual), "mean_residual" = mean(residual), "length" = n()) %>% ungroup %>% arrange(desc(length)) %>% slice(1:(n()/2)) %>%  summarise("medsd" = median(sd, na.rm = T)) %>% unlist
  
  #test_data$new_CN <- round(predict(model, data.frame("median_segment" = test_data$residual)), 0)
  #test_data$flagged <- test_data$CN != test_data$new_CN
  #test_data %>% filter(chr == 1, residual < 1e+06) %>% ggplot(aes(x = pos, y = residual, color = flagged)) + geom_point() + scale_color_viridis(discrete = T)
  
  #unmerged_data$`1` %>% ggplot(aes(x = pos, y = residual, color = CN)) + geom_point() + scale_color_viridis()
  
  for(chr in names(unmerged_data)){
    unmerged_chr <- unmerged_data[[chr]] %>% arrange(pos)
    merged_chr <- merged_data[[chr]]
    merged_segments <- merged_chr %>% group_by(chr, segment) %>% summarise(pos = dplyr::first(pos), end = last(end), CN = mean(CN))
    merged_segments$pos[1] <- 0
    #unmerged_chr$segment <- findInterval(unmerged_chr$pos, merged_segments$pos)
    #unmerged_chr <- left_join(unmerged_chr, merged_segments[,c("segment", "CN")], by = "segment")
    unmerged_chr$mafflipped <- abs(unmerged_chr$maf - 0.5) + 0.5
    #unmerged_chr$new_CN <- round(predict(model, data.frame("median_segment" = unmerged_chr$residual)), 0)
    unmerged_chr <- unmerged_chr %>% group_by(segment) %>% mutate("mean_residual" = mean(residual), "z" = (residual - mean(residual))/sd)
    #unmerged_chr %>% ggplot(aes(x = pos, y = residual, color = abs(z) > 3)) + geom_point() + scale_color_viridis(discrete = T)
    unmerged_chr <- unmerged_chr %>% group_by(segment) %>% mutate("median_segment" = median(residual), "median_maf" = median(mafflipped, na.rm = T))
    unmerged_chr <- unmerged_chr[,c("chr", "pos", "end", "mafflipped", "residual", "CN", "segment", "median_segment", "median_maf", "z")] %>% arrange(pos)
    unmerged_chr$flagged <- F
    #unmerged_chr$flagged[which(abs(unmerged_chr$residual - unmerged_chr$median_segment) > unmerged_diff * 0.5)] <- T
    unmerged_chr$flagged[which(abs(unmerged_chr$z) > 3)] <- T
    #unmerged_chr %>% ggplot(aes(x = pos, y = residual, color = flagged)) + geom_point() + scale_color_viridis(discrete = T)
    whichflagged <- which(unmerged_chr$flagged)
    for(flagged in whichflagged){
      if(!any(c(flagged - 1, flagged + 1) %in% whichflagged)){
        unmerged_chr$flagged[flagged] <- F
      }
    }
    #unmerged_chr %>% filter(CN < 10) %>%  ggplot(aes(x = pos, y = residual, color = CN)) + geom_point() + scale_color_viridis(discrete = F)
    #unmerged_chr$flagged <- unmerged_chr$CN != unmerged_chr$new_CN
    unmerged_chr$old_segment <- F
    #merged_chr[which(merged_chr$breakpoint),]
    
    ## Identify intervals of old breakpoints
    #breakpoints <- c(merged_chr$pos[which(merged_chr$breakpoint)], merged_chr$end[which(merged_chr$breakpoint)]) %>% sort()
    
    breakpoints <- merged_chr %>% group_by(segment) %>% summarise("5-term_pos" = first(pos) - 1, "5-term_end" = first(end) + 1, "3-term_pos" = last(pos) - 1, "3-term_end" = last(end) + 1)
    
    breakpoints <- as.matrix(breakpoints[,2:5]) %>% t() %>% as.vector() %>% sort()
    
    ## All windows that fall within old breakpoints need to be broken into length=1 segments
    ## First we generate intervals based on the old breakpoints
    unmerged_chr$atomize <- findInterval(unmerged_chr$pos, breakpoints)
    intervals <- unique(unmerged_chr$atomize)
    ## Find odd-numbered intervals, which denote the points which actually fell within the range of old breakpoints
    relevant_intervals <- intervals[which(intervals %% 2 == 1)]
    
    breakpoints <- c(unmerged_chr$pos[which(unmerged_chr$atomize %in% relevant_intervals)], unmerged_chr$end[which(unmerged_chr$atomize %in% relevant_intervals)]) %>% sort()
    
    #breakpoints <- c(breakpoints, breakpoints + 1)
    new_segment_interval <- unmerged_chr$pos[which(unmerged_chr$flagged)]

    new_segment_interval <- sort(c(new_segment_interval, new_segment_interval + 1))
    new_segment_interval <- sort(c(breakpoints, new_segment_interval))
    unmerged_chr$segment <- findInterval(unmerged_chr$pos, new_segment_interval, rightmost.closed = F) + 1
    #unmerged_chr %>% filter(CN < 10) %>% ggplot(aes(x = pos, y = residual, color = CN)) + geom_point() + scale_color_viridis(discrete = F)
    unmerged_to_compress <- unmerged_chr %>% mutate("npoints" = 1) %>% group_by(segment) %>% arrange(pos) %>% summarise("chr" = dplyr::first(chr), "pos" = dplyr::first(pos), "end" = last(end), "npoints" = sum(npoints), "residual" = sum(residual)) %>% arrange(pos) %>% mutate("mean_residual" = residual/npoints)
    #unmerged_to_compress %>% ggplot(aes(x = pos, xend = end, y = residual/npoints, yend = residual/npoints)) + geom_segment() + geom_point(data = unmerged_chr_filt, mapping = aes(x = pos, y = residual, color = segment), inherit.aes = F) + scale_color_viridis()
    new_segments <- runiterativecompression(t = unmerged_to_compress, x = unmerged_diff, segmentation_threshold = 0.5, verbose = verbose)
    #new_segments <- compressdata(t = unmerged_to_compress, x = unmerged_diff, segmentation_threshold = 0.25) %>% mutate("segment" = 1:n())
    #new_segments %>% filter(CN < 10) %>% ggplot(aes(x = pos, y = residual/npoints, color = segment)) + geom_point() + scale_color_viridis() # + geom_point(data = unmerged_chr_filt, mapping = aes(x = pos, y = residual, color = segment), inherit.aes = F) + scale_color_viridis()
    
    new_segments <- new_segments %>% group_by(segment) %>% summarise("chr" = dplyr::first(chr), "pos" = dplyr::first(pos), "end" = last(end), "npoints" = sum(npoints), "residual" = sum(residual), "len" = n())

    #long_segment <- 0
    #new_chr = F
    #if(nrow(new_segments) > 1){
    #  for(segment in 1:nrow(new_segments)){
    #    ## If there's a new chromosome, reset long_segment
    #    if(segment > 1){
    #      if(new_segments$chr[segment] != new_segments$chr[segment-1]){
    #        long_segment <- 0
    #      }
    #    }
    #    ### Lookahead to get new long_segment
    #    ## If long_segment is zero
    #    if(long_segment == 0){
    #      ## Record current chromosome
    #      current_chr = new_segments$chr[segment]
    #      ## Lookahead
    #      for(sub_segment in segment:nrow(new_segments)){
    #        ## If we get to a new chromosome, break
    #        if(current_chr != new_segments$chr[sub_segment]){
    #          new_chr = T
    #          break
    #        }
    #        ## When long_segment is found, record it
    #        if(new_segments$npoints[sub_segment] > 1){
    #          long_segment <- new_segments$segment[sub_segment]
    #          break
    #        }
    #      }
    #      if(new_chr){
    #        new_chr = F
    #        next
    #      }
    #    }
    #    ## If the current segment is "long" ie. supported by at least two adjacent points, make this the new "long_segment"
    #    if(new_segments$npoints[segment] > 1){
    #      long_segment <- new_segments$segment[segment]
    #    }
    #    ## If we haven't found a new long_segment yet, skip#

        ## If the current segment is length one, merge to long_segment
    #    if(new_segments$npoints[segment] == 1){
    #      new_segments$segment[segment] <- long_segment
    #    }
    #  }
    #}
    #new_segments <- new_segments %>% ungroup %>% group_by(segment) %>% summarise("chr" = dplyr::first(chr), "pos" = dplyr::first(pos), "end" = last(end), "npoints" = sum(npoints), "residual" = sum(residual), "len" = n())
    #print(new_segments)
    print("Generating segments based on compressiond ata")
    unmerged_chr$segment <- findInterval(unmerged_chr$pos, new_segments$pos, rightmost.closed = F)
    
    centromere_start <- centromeres$pos[centromeres$chr == chr][1]
    centromere_end <- centromeres$end[centromeres$chr == chr][2]
    
    unmerged_chr <- unmerged_chr %>% group_by(segment) %>% mutate("median_segment" = median(residual), "median_maf" = median(mafflipped, na.rm = T))
    print("calling copy number")
    unmerged_chr$CN <- round(predict(model, unmerged_chr), digits = 0)
    print("calling breakpoints")
    unmerged_chr <- callbreakpoints(unmerged_chr, predictedpositions = predictedpositions, maxpeak = unmerged_maxpeak)
    #unmerged_chr %>% filter(CN < 10) %>%  ggplot(aes(x = pos, y = residual, color = CN)) + geom_point() + scale_color_viridis(discrete = F)
    names(unmerged_chr)[which(names(unmerged_chr) == "residual")] <- "raw_residual"
    names(unmerged_chr)[which(names(unmerged_chr) == "residual")] <- "residual"
    unmerged_data[chr] <- list(unmerged_chr)
  }
  return(unmerged_data)
}
