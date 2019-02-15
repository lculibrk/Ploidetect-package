#' @export
testMAF <- function(CN, tp){
  if(CN < 0 | CN > 8){
    stop("Please enter a CN between 0 and 8")
  }
  np <- 1-tp
  halfcn <- ceiling(CN/2)
  major_allele_possibilities = seq(from = halfcn, to = CN, by = 1)
  output <- c()
  for(al in major_allele_possibilities){
    majoraf <- ((al * tp) + (1 * np))/((CN * tp) + (2 * np))
    output[as.character(al)] <- majoraf
  }
  return(output)
}
