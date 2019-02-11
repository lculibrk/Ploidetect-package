runiterativecompression <- function(t, x, segmentation_threshold = segmentation_threshold){
  converged <- F
  compress <- t
  if(nrow(t) == 1){
    converged = T
  }
  while(!converged){
    windows <- nrow(compress)
    compress <- compressdata(compress, x, segmentation_threshold)
    if(nrow(compress) == windows){
      converged <- T
    }
  }
  t$segment <- findInterval(t$pos, vec = compress$pos)
  return(t)
}