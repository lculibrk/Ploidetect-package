merge_mafs <- function(inputmafs, na.rm = T){
  if(na.rm){
    inputmafs <- inputmafs[which(!is.na(inputmafs))]
  }
  ## Compute magnitude about 0.5
  mafs_about_zero <- inputmafs - 0.5
  negatives <- which(mafs_about_zero < 0)
  positives <- which(mafs_about_zero > 0)
  ## Compute absolute magnitude about 0.5 and add 0.5
  flipped_mafs <- abs(mafs_about_zero)
  if(length(positives) > length(negatives)){
    out <- mean(0.5 + flipped_mafs)
  }else{
    out <- mean(0.5 - flipped_mafs)
  }
  return(out)
}
