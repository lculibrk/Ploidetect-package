ploidetect_processallpeaks <- function(filtered, allPeaks, verbose = F){
  ## Center the residual and peak data about the tallest peak
  filtered$residual <- filtered$residual - allPeaks$pos[1]
  allPeaks$start <- allPeaks$start - allPeaks$pos[1]
  allPeaks$end <- allPeaks$end - allPeaks$pos[1]
  allPeaks$pos <- allPeaks$pos - allPeaks$pos[1]
  ## Assign IDs to each peak
  allPeaks <- allPeaks %>% arrange(pos) %>% mutate(npeak = seq(1, nrow(allPeaks))) %>% arrange(desc(height))
  #filtered$mafflipped <- abs(filtered$maf - 0.5)+0.5
  ## Map peaks to data
  filtered$peak <- NA
  for(i in 1:nrow(allPeaks)){
    interval <- unlist(c(allPeaks[i,c("start", "end")]))
    filtered$peak[findInterval(x = filtered$residual, vec = interval) == 1] <- allPeaks$npeak[i]
  }
  ## Filter data for only data with allelic frequency values
  filterednomafna <- filtered %>% filter(!is.na(maf))
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
    peakdata <- filterednomafna[filterednomafna$peak == allPeaks$npeak[i],]
    peakdata <- peakdata %>% filter(!is.na(maf))
    if(nrow(peakdata) > 1){
      peakdata <- filterednomafna[filterednomafna$peak == allPeaks$npeak[i],]
      peakdata <- peakdata %>% filter(!is.na(residual))
      position <- allPeaks$pos[i]
      peakdata$dev <- peakdata$residual - position
      peakdata <- peakdata %>% arrange(dev)
      peakdata <- peakdata[1:max(c(round(nrow(peakdata)/10, digits = 0), 2)),]
      md <- density(unmerge_mafs(peakdata$maf, flip = T), na.rm = T)
      allPeaks$mainmaf[i] <- md$x[which.max(md$y)]
    }else{
      allPeaks$mainmaf[i] <- NA
    }
  }
  filtered_maf <- filtered %>% filter(!is.na(maf))
  filtered_nomaf <- filtered %>% filter(is.na(maf))
  filtered_nomaf$mafflipped <- NA
  mafs <- unmerge_mafs(filtered$maf)
  mafs <- lapply(mafs, function(x){paste0(abs(as.numeric(x) - 0.5) + 0.5, collapse = ";")}) %>% unlist()
  filtered_maf$mafflipped <- mafs
  filtered <- rbind.data.frame(filtered_maf, filtered_nomaf) %>% arrange(chr, pos)
  output <- list("filtered" = filtered, "allPeaks" = allPeaks, "nomaf" = nomaf)
  return(output)
}
