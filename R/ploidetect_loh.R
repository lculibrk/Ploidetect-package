ploidetect_loh <- function(purity, CNA_object){
  CNA_object <- CNA_object %>% arrange(CN)
  CNA_object <- CNA_object %>% group_by(chr, segment) %>% dplyr::mutate("median_maf" = median(mafflipped, na.rm = T))
  CNA_object$LOH <- F
  CNA_object_mafs <- CNA_object[which(!is.na(CNA_object$median_maf)),]
  CNA_object_nomafs <- CNA_object[which(is.na(CNA_object$median_maf)),]
  allCNs <- unique(CNA_object_mafs$CN)
  allCNs <- allCNs[allCNs > 1]
  allCNs <- allCNs[allCNs <= 8]
  maf_possibilities <- list()
  for(copynumber in allCNs){
    mafs <- testMAF(copynumber, tp = purity)
    LOH_maf <- mafs[length(mafs)]
    almost_LOH_maf <- mafs[length(mafs) - 1]
    maf_possibilities[paste0(copynumber)] <- list(c(LOH_maf, almost_LOH_maf))
  }
  for(row in 1:nrow(CNA_object_mafs)){
    CN <- paste0(CNA_object_mafs$CN[row])
    if(CN == "1"){
      CNA_object_mafs$LOH[row] <- T
    }
    if(CN %in% names(maf_possibilities)){
      candidates <- unlist(maf_possibilities[CN])
      decision <- which.min(abs(CNA_object_mafs$median_maf[row] - candidates))
      if(decision == 1){
        CNA_object_mafs$LOH[row] <- T
      }else{
        CNA_object_mafs$LOH[row] <- F
      }
    }
  }
  CNA_object <- rbind.data.frame(CNA_object_mafs, CNA_object_nomafs) %>% arrange(chr, pos)
  return(CNA_object)
}