ploidetect_roughmodels <- function(allPeaks, maxpeak, verbose = F, rerun = rerun){
  if(verbose){
    print("Generating list of peak models")
  }
  ## Check if we have enough peaks to model with this approach
  if(nrow(allPeaks) >= 3){
    ## Arrange by height
    allPeaks <- allPeaks %>% arrange(desc(height))
    # Blank data frame to rbind to later
    models <- data.frame("mean_err" = 0, 
                         "reads_per_copy" = 0, 
                         "Copy_number_difference_between_peaks" = 0, 
                         "Comparator_peak_rank" = 0, 
                         "unmatched" = 0,
                         "predunmatched" = 0)
    # For all peaks after peak intensity #1:
    for(l in 2:(nrow(allPeaks))){
      # For the possibility of having up to five levels of CNA between l and biggestpeak (this is slightly overkill):
      for(j in 1:5){
        if(verbose){
          print(paste0("Generating model for CN difference between fitted peaks = ", j,  " with secondary peak as #", l, " in terms of intensity"))
        }
        # Experimental
        weights <- allPeaks$height[l]
        # Compute the estimate of reads-per-copy
        x = abs((allPeaks$pos[1] - allPeaks$pos[l])/j)
        # Simulate peak topology given x for positive and negative peaks
        pos <- seq(from = allPeaks$pos[1], to = max(allPeaks$pos) + x, by = x)
        neg <- seq(from = allPeaks$pos[1], to = min(allPeaks$pos) - x, by = -x)
        # All predicted peaks positions
        allpredpeaks <- c(neg[-1], pos)
        #Extract separate lists of + and - peaks
        posPeaks <- allPeaks %>% filter(allPeaks$pos > 0)
        negPeaks <- allPeaks %>% filter(allPeaks$pos < 0)
        # Set initialized values
        t <- 0
        matchedPeaks <- c()
        pos.errors <- c()
        # For the positive peaks
        if(nrow(posPeaks) > 0){
          for(k in 2:length(pos)){
            pos <- seq(from = allPeaks$pos[1], to = max(allPeaks$pos) + x, by = x)
            if(k > length(pos)){
              break
            }
            # If we have positive peaks...
            # If we haven't matched this peak yet in the previous iteration
            if(t != which.min(abs(pos[k] - posPeaks$pos))){
              # If the peak is actually somewhere in the vicinity of where we expect it
              if(min(abs(pos[k] - posPeaks$pos)) <= x/2){
                # Which peak is closest to predicted?
                t <- which.min(abs(pos[k] - posPeaks$pos))
                x <- weighted.mean(c(x, posPeaks$pos[t]/(k-1)), w = c(weights, posPeaks$height[t]))
                weights <- weights + posPeaks$height[t]
                pos <- seq(from = allPeaks$pos[1], to = max(allPeaks$pos) + x, by = x)
                if(verbose){
                  print(paste0("matched positive peak indices in allPeaks is ", posPeaks$npeak[t]))
                }
                # Compute the difference * height of peak
                pos.errors <- c(pos.errors, posPeaks$height[t] * abs(pos[k] - posPeaks$pos)[t])
                # Add the peak to the index
                matchedPeaks <- c(matchedPeaks, posPeaks$npeak[t])
              }
            }
          }
        }
        # Same as before, but with negatives
        t <- 0
        neg.errors <- c()
        neg <- seq(from = allPeaks$pos[1], to = min(allPeaks$pos) - x, by = -x)
        if(nrow(negPeaks) > 0){
          for(k in 2:length(neg)){
            neg <- seq(from = allPeaks$pos[1], to = min(allPeaks$pos) - x, by = -x)
            if(k > length(neg)){
              break
            }
            if(t != which.min(abs(neg[k] - negPeaks$pos))){
              if(min(abs(neg[k] - negPeaks$pos)) <= x/2){
                t <- which.min(abs(neg[k] - negPeaks$pos))
                x <- weighted.mean(c(x, abs(negPeaks$pos[t])/(k-1)), w = c(weights, negPeaks$height[t]))

                weights <- weights + negPeaks$height[t]
                neg.errors <- c(neg.errors, negPeaks$height[t] * abs(neg[k] - negPeaks$pos)[t])
                matchedPeaks <- c(matchedPeaks, negPeaks$npeak[t])
              }
            }
          }
        }
        # Add the biggest peak to matchedPeaks
        matchedPeaks <- c(matchedPeaks, allPeaks$npeak[1])
        # Make list of all unmatched peaks
        unmatchedPeaks <- allPeaks$npeak[!(allPeaks$npeak %in% matchedPeaks)]
        # all unmatched peaks' data
        unmatcheddat <- allPeaks[allPeaks$npeak %in% unmatchedPeaks,]
        # Arrange
        allPeaks <- allPeaks %>% arrange(pos)
        # List of predicted peaks that were NOT matched by the model
        prednotmatched <- allpredpeaks
        predictednotmatched <- c()
        for(i in 1:nrow(allPeaks)){
          predictednotmatched <- c(predictednotmatched, which.min(abs(allPeaks$pos[i] - prednotmatched)))
        }
        prednotmatched <- allpredpeaks[-predictednotmatched]
        # Remove predicted but not matched peaks that are above and below the real peaks (ie. removing peaks predicting CNAs one level beyond what exist in the data)
        ind <- c()
        if(length(prednotmatched) > 0){
          for(i in 1:length(prednotmatched)){
            if(prednotmatched[i] < min(allPeaks$pos) | prednotmatched[i] > max(allPeaks$pos)){
              ind <- c(ind, i)
            }
          }
        }
        if(length(ind) > 0){
          prednotmatched <- prednotmatched[-ind]
        }
        # Normalize to maxpeak's value
        prednotmatched <- prednotmatched + maxpeak
        # Arrange
        allPeaks <- allPeaks %>% arrange(desc(height))
        # If we have data to compute error at all
        if(length(c(pos.errors, neg.errors)) > 1 | nrow(allPeaks) == 3 | rerun){
          # Deprecated error function
          error <- sqrt(mean(c(pos.errors, neg.errors)^2))/allPeaks$height[l]
          # Put model as a row of a df
          model <- data.frame("mean_err" = error, 
                              #"mean_err" = sum(unmatchedIntensities) * unmatchedPeaks * j * mean(c(pos.errors, neg.errors)), 
                              "reads_per_copy" = x, 
                              "Copy_number_difference_between_peaks" = j, 
                              "Comparator_peak_rank" = l,
                              "unmatched" = paste(unmatchedPeaks, collapse = "_"), stringsAsFactors = F,
                              "predunmatched" = paste(prednotmatched, collapse = "_"))
          # Bind the models
          models <- rbind(models, model)
        }
      }
    }
    # Remove the dummy models row from the beginning
    models <- models[-1,]
    # xdists, a data.frame of all the models
    xdists <- models %>% arrange(mean_err)
  }
  if(nrow(allPeaks) == 2){
    x <- max(allPeaks$diff)
    xdists <- data.frame("mean_err" = c(0, 0), 
                         "reads_per_copy" = c(x, x/2), 
                         "Copy_number_difference_between_peaks" = c(1, 2), 
                         "Comparator_peak_rank" = c(2, 2), 
                         "unmatched" = as.character(c("", "")),
                         "predunmatched" = c("", as.character(allPeaks$pos[which.max(allPeaks$height)] + x/2 + maxpeak)), 
                         stringsAsFactors = F)
  }
  if(nrow(allPeaks) == 1){
    xdists <- "No Copy number alteration detected. TC cannot be determined by this method"
    return(xdists)
  }
  xdists <- xdists %>% group_by(unmatched, predunmatched) %>% dplyr::summarise(reads_per_copy = mean(reads_per_copy), Copy_number_difference_between_peaks = first(Copy_number_difference_between_peaks), Comparator_peak_rank = first(Comparator_peak_rank))
  return(xdists)
}
