modelbuilder_iterative <- function(xdists = xdists, allPeaks = allPeaks, lowest = NA, filtered = filtered, strict = T, get_homd = F, mode = "TC", nomaf = nomaf, rerun = rerun, maxpeak = maxpeak,bw = bw){
  out <- list()
  outplots <- list()
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
  xdists <- as.list(xdists)
  # Ensure no carry-over from previous loops
  matchedPeaks = NULL; x = NULL; homd = NULL; lowestpeakmaf = NULL;
  # If we have models to compare
  # This step generates "matchedPeaks", a data.frame containing the peaks to be used in model computation
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
  matchedPeaks <- allPeaks[!(allPeaks$npeak %in% unlist(strsplit(xdists$unmatched, split = "_"))),]
  unmatchedPeaks <- allPeaks[allPeaks$npeak %in% unlist(strsplit(xdists$unmatched, split = "_")),]
  # Sum the weights of all unmatched peaks for peak-skipping evaluation
  unmatchederror <- sum(nonbiggest$height[nonbiggest$npeak %in% unmatchedPeaks$npeak])
  matchederror <- 1-unmatchederror
  # If we only had one model then just assume it's perfect
  # Add predicted, unmatched peaks to matchedPeaks
  matchedPeaks <- plyr::rbind.fill(matchedPeaks, data.frame("pos" = as.numeric(unlist(strsplit(as.character(xdists$predunmatched), split = "_"))))) %>% arrange(pos)
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
  # If we have more ghost peaks than  half the true peaks (in #peaks > 3 cases), filter out model
  # If peak count is NOT 3 and less than half of the peaks are sequential, filter out
  if(strict & !rerun){
    # If over half the peaks are ghost peaks
    if(length(naindex) > nrow(matchedPeaks)/2){
      return()
    }
    # If few peaks are sequential in >3-peak models
    if(sum(as.numeric(sequential)) < length(peakindex)/2){
      if(nrow(allPeaks) > 3){
        return()
      }
    }
  }
  
  #if(((length(naindex) > nrow(matchedPeaks)/2) | ((sum(as.numeric(sequential)) < length(peakindex)/2) & nrow(allPeaks) > 3)) & strict & !rerun){
  #  next
  #}
  # If we have three peaks, filter out if have zero sequentials or if we have more than one ghost IF strict mode is enabled
  if(nrow(allPeaks) == 3 & strict){
    if((sum(as.numeric(sequential)) < 1) | (length(naindex) > 1)){
        return()
    }
  }
  # Number of peaks
  npeaks <- nrow(matchedPeaks)
  
  peakslist <- rep(list(matchedPeaks), times = length(lowestrange))
  names(peakslist) <- lowestrange
  mafdev <- list()
  mafdeveven <- list()
  for(lowestpeak in lowestrange){
    peakslist[paste0(lowestpeak)] <- peakslist[paste0(lowestpeak)] %>% 
      as.data.frame() %>%
      arrange_at(4) %>% 
      mutate("CN" = seq(from = lowestpeak, length.out = nrow(as.data.frame(peakslist[paste0(lowestpeak)])))) %>% 
      filter(CN >= 0) %>% 
      list()
  }
  for(lowestpeak in lowestrange){
    ## Get matchedPeaks for this particular assumption of lowestpeak
    matchedPeaks <- as.data.frame(peakslist[paste0(lowestpeak)])
    ## Fix column names of matchedPeaks
    names(matchedPeaks)[1:10] <- substring(names(matchedPeaks)[1:10], first = 7)
    names(matchedPeaks)[11] <- substring(names(matchedPeaks)[11], first = 4)
    # Fit a weighted linear model to the matchedPeaks model to find slope (reads per copy) and intercept (HOMD coverage)
    CNs <- seq(from = lowestpeak, length.out = npeaks, by = 1)
    model <- lm(matchedPeaks$pos~CNs, weights = matchedPeaks$height)
    predictedpositions <- predict(model, data.frame("CNs" = seq(from = 0, to = 100)))
    #predictedpositions <- predict(model, data.frame("CNs" = seq(from = 0, to = max(matchedPeaks$CN, 2))))
    names(predictedpositions) <- seq(from = 0, to = 100)
    #model <- lm(y_raw ~ CN, data = filteredonlypeaks)
    # Degrees of freedom
    df <- df.residual(model)
    # For three-peak models (these are really annoying and hard) we allow the possibility of mapping to only two of three peaks
    # However the last peak is used in error estimation, unlike normally
    threepkerr <- 0
    if(nrow(allPeaks) == 3 & mode == "TC"){
      if(nrow(unmatchedPeaks) == 1){
        threepkerr <- min(abs(unmatchedPeaks$pos - model$fitted.values)) * unmatchedPeaks$height
        df <- 1
        warning("Only three peaks exist - allowing possibility of one of three peaks being subclonal. Dubious results possible")
      }
    }
    fitpeaks <- model$fitted.values
    
    # Summarise the model
    model <- summary(model)
    x <- model$coefficients[2,1]
    # HOMD from intercept
    homd <- model$coefficients[1,1]
    # Purity & impurity
    purity <- 1-(homd/(2*x + homd))
    impurity <- 1-purity
    
    if(mode == "CNA"){
      return(list("predictedpositions" = predictedpositions, "matchedPeaks" = matchedPeaks, "depthdiff" = x))
    }
    
    # Error function: Mean Standard Error multiplied by 1 + the number of ghost peaks, all divided by the logistic-transformed proportion of the data that was matched
    newerrors <- (sqrt(sum(model$residuals^2)/df) *  (length(naindex) + 1))/(matchederror)
    # For three-peak models, add the error computed a few lines above
    if(nrow(allPeaks) == 3){
      newerrors <- newerrors + threepkerr
    }
    # If zero degrees of freedom (ie. two peak model), error is zero (otherwise gets coerced to Inf)
    if(df == 0){
      newerrors <- 0
    }
    if(nomaf == F){
      # Initialization
      filterednomafna <- filtered %>% filter(!is.na(mafflipped))
      matchedPeaks$expected <- 0
      matchedPeaks$mafdeviation <- 0
      error <- c()
      toterr <- c()
      # Looping over all peaks
      for(j in 1:nrow(matchedPeaks)){
        if(!(matchedPeaks$CN[j] %in% 0:4)){
          next
        }
        # If it's not a ghost peak
        if(!is.na(matchedPeaks$mainmaf[j])){
          error <- c()
          err <- c()
          mafs <- filterednomafna$maf[filterednomafna$peak %in% matchedPeaks$npeak[j]]
          # CN for testing
          CN <- matchedPeaks$CN[j]
          # For below, we essentially make a vector of all possible MAF values at each CN state given the TC
          # "err" is a vector of the closest predicted MAF to the real MAF
          # We do this for CN 0-4, as beyond this the possibilities get numerous and doesn't really help
          if(CN == 0){
            err <- c()
            matchedPeaks$expected[j] = 0.5
            err <- abs(mafs - matchedPeaks$expected[j])
            #error <- c(err, error)
          }
          if(CN == 1){
            err <- c()
            comp <- c((0*purity + 1*impurity)/(CN*purity+2*impurity), 
                      (1*purity + 1*impurity)/(CN*purity+2*impurity))
            for(g in 1:length(mafs)){
              err[g] <- min(abs(mafs[g] - comp))
            }
          }
          if(CN == 2){
            err <- c()
            comp <- c(0.5,
                      (0*purity + 1*impurity)/(CN*purity+2*impurity),
                      (2*purity + 1*impurity)/(CN*purity+2*impurity))
            for(g in 1:length(mafs)){
              err[g] <- min(abs(mafs[g] - comp))
            }
          }
          if(CN == 3){
            err <- c()
            comp <- c((0*purity + 1*impurity)/(CN*purity+2*impurity),
                      (1*purity + 1*impurity)/(CN*purity+2*impurity),
                      (2*purity + 1*impurity)/(CN*purity+2*impurity),
                      (3*purity + 1*impurity)/(CN*purity+2*impurity))
            for(g in 1:length(mafs)){
              err[g] <- min(abs(mafs[g] - comp))
            }
          }
          if(CN == 4){
            err <- c()
            comp <- c(0.5,
                      (1*purity + 1*impurity)/(CN*purity + 2*impurity),
                      (3*purity + 1*impurity)/(CN*purity + 2*impurity),
                      (0*purity + 1*impurity)/(CN*purity + 2*impurity),
                      (4*purity + 1*impurity)/(CN*purity + 2*impurity))
            for(g in 1:length(mafs)){
              err[g] <- min(abs(mafs[g] - comp))
            }
          }
          error <- c(error, err)
          matchedPeaks$mafdeviation[j] <- mean(error)
          toterr <- c(toterr, error)
        }
      }
    }
    # Get ploidy
    ploidy <- matchedPeaks$CN[which.max(matchedPeaks$height)]
    # Set height of ghosts to zero
    matchedPeaks$height[is.na(matchedPeaks$height)] <- 0
    # Get average ploidy
    average_ploidy <- round(weighted.mean(matchedPeaks$CN, w = matchedPeaks$height), digits = 2)
    #
    # Mafdev is the median difference between predicted and observed MAFs
    mafdev[paste0(lowestpeak)] <- mean(toterr)
    lowestsize <- matchedPeaks$height[1]
    mafdeveven[paste0(lowestpeak)] <- mean(toterr)*(1+(max(ploidy-2, 0)*(1-lowestsize)))
    if(any(matchedPeaks$mafdeviation > 0.15)){
      newerrors <- Inf
      mafdev[paste0(lowestpeak)] <- Inf
    }
    # If tc<1.1 or < 0 then set mafdev to Inf (this is an impossible model) (1.1 to allow slight inaccuracies)
    if(!is.na(purity)){
      if(purity > 1.1 | purity < 0){
        mafdev[paste0(lowestpeak)] <- Inf
        #if(rerun){
        #  stop("You are trying to force a model which would result in TC > 110% or lower than 0%")
        #}
      }
      
    }
    #if we have a case of ridiculous by-chance fit with peak skipping (< 10^-8 works well), filter out
    if((unmatchederror > 0 & newerrors < 10^-8 & newerrors > 0)){
      x <-  NA
    }
    # If no maf in data
    if(nomaf){
      mafdev[paste0(lowestpeak)] <- NA
      lowestpeakmaf <- NA
    }
    # Arrange
    matchedPeaks <- matchedPeaks %>% arrange(pos)
    # If the biggest peak is somehow a HOMD, correct that
    if(nrow(matchedPeaks) > 0){
      if(ploidy == 0){
        mafdev[paste0(lowestpeak)] <- Inf
      }
    }
    if(strict & unmatchederror > 0.9){
      newerrors <- Inf
    }
    out[paste0(lowestpeak)] <- list(data.frame("reads_per_copy" = x, 
                                               "zero_copy_depth" = homd, 
                                               "Ploidy" = ploidy, 
                                               "tumour_purity" = purity, 
                                               "lowest_peak_CN" = lowestpeak, 
                                               "maf_error" = as.numeric(mafdev[paste0(lowestpeak)]), 
                                               "CN_diff" = xdists$Copy_number_difference_between_peaks, 
                                               "Comparator" = xdists$Comparator_peak_rank,
                                               "model_error" = newerrors * as.numeric(mafdev[paste0(lowestpeak)]),
                                               "avg_ploidy" = average_ploidy,
                                               "unmatchedpct" = unmatchederror))
  }
  out <- do.call(rbind.data.frame, out)
  if(nrow(out) > 1){
    out <- out[-which.max(out$maf_error),]
  }
  if(all(out$Ploidy %% 2 == 0)){
    for(level in out$lowest_peak_CN){
      out$maf_error[which(out$lowest_peak_CN == level)] <- as.numeric(mafdeveven[paste0(level)])
    }
  }
  out <- out[which.min(out$maf_error),]
  
  fitpeaks <- fitpeaks[fitpeaks > 0]
  color_frame <- data.frame("positions" = fitpeaks, "col" = "#ED553B", stringsAsFactors = F)
  matchedPeaks$CN <- seq(from = out$lowest_peak_CN[1], length.out = nrow(matchedPeaks))
  matchedPeaks <- matchedPeaks %>% arrange(pos)
  for(i in 1:nrow(matchedPeaks)){
    if(matchedPeaks$CN[i] == 0){
      color_frame$col[i] <- "#000000"
    }else if(matchedPeaks$CN[i] == 1){
      color_frame$col[i] <- "#00245e"
    }else if(matchedPeaks$CN[i] == 2){
      color_frame$col[i] <- "#25ba00"
    }else if(matchedPeaks$CN[i] == 3){
      color_frame$col[i] <- "#F6D55C"
    }else{
      color_frame$col[i] <- "#ED553B"
    }
  }
  
  newden <- filtered$residual[which(filtered$residual < max(matchedPeaks$pos) + x)] %>% density()
  plot <- ggplot(mapping = aes(x = newden$x + maxpeak, y = (newden$y - min(newden$y))/(max(newden$y) - min(newden$y)))) + 
    geom_line() + 
    xlab("Read Counts") + 
    ylab("Normalized Density") + 
    ggtitle(paste0("Model for \"CN_diff\" = ", xdists$Copy_number_difference_between_peaks, ",\"comparator\" = ", xdists$Comparator_peak_rank)) + 
    geom_vline(data = color_frame, aes(xintercept = positions, color = factor(1:nrow(color_frame))), linetype = 2) + 
    geom_text(mapping = aes(x = allPeaks$pos, y = allPeaks$height + 0.05, label = paste0("MAF = ", round(allPeaks$mainmaf, digits = 3)))) +
    theme_bw(base_size = 12) + 
    scale_color_manual(values = color_frame$col, labels = paste0("CN = ", matchedPeaks$CN), name = "Absolute copy number") #+ theme(legend.position = "none")
  
  return(list("out" = out, "outplot" = plot))
}
