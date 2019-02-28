peakcaller <- function(filtered, bw = bw, verbose = F){
  if(verbose){
    print("Performing Kernel Density Estimation and Peak Calling")
  }
  # Filter out all regions without a MAF call
  filterednomafna <- filtered %>% filter(!is.na(maf))
  # Density estimate
  den <- density(filtered$residual, n = nrow(filtered), bw = bw)
  ## Normalize the density to 0->1 range
  den$y <- (den$y - min(den$y))/(max(den$y) - min(den$y))
  # Take second derivative
  ddx <- dkde(filtered$residual, deriv.order = 2, h = bw)
  dx <- dkde(filtered$residual, deriv.order = 1, h = bw)
  dx$est.fx <- abs(dx$est.fx)
  # Call peaks, essentially scan along ddx and find where it goes from negative to positive
  peaks <- list()
  for(i in 1:length(ddx$eval.points)){
    if(i > 1){
      if(ddx$est.fx[i] < 0 & ddx$est.fx[i-1] >= 0){
        peakstart <- ddx$eval.points[i]
      }
      if(ddx$est.fx[i] > 0 & ddx$est.fx[i-1] <= 0){
        peakend <- ddx$eval.points[i]
        peak <- c("start" = peakstart, "end" = peakend)
        peaks <- c(peaks, list(peak))
      }
    }
  }
  ## Get peak height
  newpeaks <- lapply(peaks, function(x){
    if(length(den$y[findInterval(x = den$x, vec = unlist(x)) == 1]) > 0){
      c(x, 
        height = max(den$y[findInterval(x = den$x, vec = unlist(x)) == 1]), 
        pos = (den$x[(den$y == max(den$y[findInterval(x = den$x, vec = x[1:2]) == 1])) & (den$x >= x[1]) & (den$x) <= x[2]]))
    }
  })
  # Merge into a df, arrange
  allPeaks <- as.data.frame(do.call(rbind, newpeaks))
  allPeaks <- allPeaks %>% arrange(pos)
  # Call the "troughs" (max(lowest region before the peak, lowest region after the peak))
  if(nrow(allPeaks) > 1){
    val <- c()
    for(i in (nrow(allPeaks)):1){
      # This is for the final iteration, otherwise indexes into the negatives. Probably hacky
      val2 <- 0
      tryCatch(
        val2 <- min(den$y[findInterval(x = den$x, vec = allPeaks$pos[i:(i+1)]) == 1]), error = function(e) 0)
      val3 <- min(den$y[findInterval(x = den$x, vec = allPeaks$pos[(i-1):i]) == 1])
      val4 <- max(val2, val3)
      val <- c(val,val4)
    }
  }else {
    val <- 0
  }
  # Merge troughs into new df
  allPeaks <- cbind.data.frame(allPeaks, "trough" = rev(val))
  # Calculate some more stuff from the troughs (one normalizes to height and one normalizes by absolute)
  allPeaks$ratiotrough <- allPeaks$trough/allPeaks$height
  allPeaks$troughdiff <- allPeaks$height - allPeaks$trough
  # Calculate the distance between peak and tallest peak
  allPeaks$diff <- abs(allPeaks$pos - allPeaks$pos[which(allPeaks$height == max(allPeaks$height))])
  # Arrange by height and filter for peaks representing density < 0.5%
  allPeaks <- allPeaks %>% arrange(desc(height)) %>% filter(height > 0.005)
  
  if(nrow(allPeaks %>% filter(ratiotrough < 1)) == 1){
    allPeaks <- allPeaks %>% filter(ratiotrough < 1)
  }
  return(allPeaks)
}
