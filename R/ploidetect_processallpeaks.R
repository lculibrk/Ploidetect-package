ploidetect_processallpeaks <- function(filtered, allPeaks, verbose = F){
  ## Assign IDs to each peak
  allPeaks <- allPeaks %>% arrange(pos) %>% mutate(npeak = seq(1, nrow(allPeaks))) %>% arrange(desc(height))
  
  ## Map peaks to data
  filtered$peak <- NA
  for(i in 1:nrow(allPeaks)){
    interval <- unlist(c(allPeaks[i,c("start", "end")]))
    filtered$peak[findInterval(x = filtered$residual, vec = interval) == 1] <- allPeaks$npeak[i]
  }
  ## Filter data for only data with allelic frequency values
  filterednomafna <- filtered %>% filter(!is.na(mafflipped))
  ## Check to see if we have allelic frequency data
  nomaf <- F
  if(nrow(filterednomafna) == 0){
    nomaf <- TRUE
    warning("No allelic frequency data detected. Either a bug or you never provided any in the first place, in which case ignore this")
  }
  ## Compute the median allele frequency per peak
  if(verbose){
    print("Computing allelic frequencies per peak")
  }
  allPeaks$mainmaf <- NA
  for(i in 1:nrow(allPeaks)){
    if(length(na.omit(filterednomafna$maf[filterednomafna$peak == allPeaks$npeak[i]])) > 1){
      md <- density(filterednomafna$mafflipped[filterednomafna$peak == allPeaks$npeak[i]], na.rm = T)
      allPeaks$mainmaf[i] <- md$x[which.max(md$y)]
    }else{
      allPeaks$mainmaf[i] <- NA
    }
  }
  output <- list("filtered" = filtered, "allPeaks" = allPeaks, "nomaf" = nomaf)
  return(output)
}