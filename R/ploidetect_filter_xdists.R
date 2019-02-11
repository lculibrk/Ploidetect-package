ploidetect_filter_xdists <- function(xdists = xdists, allPeaks = allPeaks, lowest, filtered = filtered, strict = T, mode = "TC"){
  if(!is.na(lowest)){
    lowestrange <- lowest
  }else{
    if(allPeaks$npeak[which.max(allPeaks$height)] == 1){
      lowestrange <- 1:2
    }else{
      if(nomaf){
        lowestrange <- 1
      }else{
        lowestrange <- 0:2
      }
    }
  }
  out <- xdists[1,]
  for(i in 1:nrow(xdists)){
    # Ensure no carry-over from previous loops
    matchedPeaks = NULL; x = NULL; homd = NULL; lowestpeakmaf = NULL;
    # If we have models to compare
    # This step generates "matchedPeaks", a data.frame containing the peaks used in model computation
    if(nrow(xdists) > 1 | mode == "CNA" | rerun){
      # Filter out the biggest peak
      nonbiggest <- allPeaks[2:nrow(allPeaks),]
      # Down-weight the peaks with trough > peak
      #nonbiggest$height[which(nonbiggest$ratiotrough >= 1)] <- nonbiggest$height[which(nonbiggest$ratiotrough >= 1)]/2
      if(any(nonbiggest$trough > nonbiggest$height)){
        nonbiggest$trough[which(nonbiggest$trough > nonbiggest$height)] <- nonbiggest$height[which(nonbiggest$trough > nonbiggest$height)]
      }
      # Scale so that we at most penalize shoulders by 50%
      nonbiggest$ratiotrough <- abs((nonbiggest$trough/nonbiggest$height)/2 - 0.5) + 0.5
      #nonbiggest$troughdiff <- nonbiggest$height - nonbiggest$trough
      # Scale by the trough-ratio
      nonbiggest$height <- nonbiggest$height * nonbiggest$ratiotrough
      # Softmax the peak weights
      nonbiggest$height <- softmax(log(nonbiggest$height, base = 2))
      # Make dfs for matched and unmatched peaks
      matchedPeaks <- allPeaks[!(allPeaks$npeak %in% unlist(strsplit(xdists$unmatched[i], split = "_"))),]
      unmatchedPeaks <- allPeaks[allPeaks$npeak %in% unlist(strsplit(xdists$unmatched[i], split = "_")),]
      # Sum the weights of all unmatched peaks for peak-skipping evaluation
      unmatchederror <- sum(nonbiggest$height[nonbiggest$npeak %in% unmatchedPeaks$npeak])
      matchederror <- 1-unmatchederror
      # If we only had one model then just assume it's perfect
    }else{matchedPeaks <- allPeaks; matchederror = 1; unmatchederror = 0; unmatchedPeaks <- allPeaks[0,]}
    # Add predicted, unmatched peaks to matchedPeaks
    matchedPeaks <- plyr::rbind.fill(matchedPeaks, data.frame("pos" = as.numeric(unlist(strsplit(as.character(xdists$predunmatched[i]), split = "_"))))) %>% arrange(pos)
    # Prune off predicteds outside the range of real peaks' locations
    while(is.na(matchedPeaks$start[nrow(matchedPeaks)])){
      matchedPeaks <- matchedPeaks[-nrow(matchedPeaks),]
    }
    while(is.na(matchedPeaks$start[1])){
      matchedPeaks <- matchedPeaks[-1,]
    }
    # Some metrics for contiguity of the peak topology
    
    naindex <- c()
    peakindex <- c()
    
    # Generate a vector of NA ("ghost") peaks and true peaks
    
    for(k in 1:nrow(matchedPeaks)){
      if(is.na(matchedPeaks$start[k])){
        naindex <- c(naindex, k)
      }
      if(!is.na(matchedPeaks$start[k])){
        peakindex <- c(peakindex, k)
      }
    }
    
    # Test whether real peaks are sequential or separated by ghost peaks
    
    sequential <- c()
    if(nrow(matchedPeaks) > 2){
      for(k in 2:length(peakindex)){
        sequential[k-1] <- peakindex[k] == (peakindex[k-1] + 1)
      }
    }else{sequential <- T}
    # Hackiest part of the program
    # If we have more ghost peaks than half the true peaks (in #peaks > 3 cases), skip this model
    # If peak count is NOT 3 and less than half of the peaks are sequential, skip this model
    if(strict & !rerun){
      # If over half the peaks are ghost peaks
      if(length(naindex) > nrow(matchedPeaks)/2){
        next
      }
      # If few peaks are sequential in >3-peak models
      if(sum(as.numeric(sequential)) < length(peakindex)/2){
        if(nrow(allPeaks) > 3){
          next
        }
      }
    }
    
    #if(((length(naindex) > nrow(matchedPeaks)/2) | ((sum(as.numeric(sequential)) < length(peakindex)/2) & nrow(allPeaks) > 3)) & strict & !rerun){
    #  next
    #}
    # If we have three peaks, filter out if have zero sequentials or if we have more than one ghost IF strict mode is enabled
    if(nrow(allPeaks) == 3){
      if(strict){
        if((sum(as.numeric(sequential)) < 1) | (length(naindex) > 1)){
          next
        }
      }
    }
    out <- rbind.data.frame(out, xdists[i,])
  }
  return(list("xdists" = out[-1,], "matchedPeaks" = matchedPeaks))
}
