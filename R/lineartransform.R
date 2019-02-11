linearTransform <- function(x, tumour, normal, avg_allele_freq, window_id, slope = 0, verbose = verbose){
  ## Fit a linear model with specified slope
  if(verbose){
    print(paste0("Performing a linear transformation using slope = ", slope))
  }
  xfit <- lm(formula = normalized_tumour ~ offset(normalized_size*slope), data = x)
  
  ## Regenerate input data.frame
  fitted_model <- cbind.data.frame(x[,tumour], x[, normal], xfit$residuals, x$normalized_size, x$GC, x$size) %>% `row.names<-`(row.names(x))
  names(fitted_model) <- c("y_raw", "x_raw", "residual", "norm_size", "GC", "size")
  # Make window column from row.names
  fitted_model$wind <- row.names(fitted_model)
  filtered <- fitted_model
  # Add the minor allelic fraction data to the df and convert it from the weird BEDTOOLS nomenclature (e.g. having . instead of NA or NULL)
  filtered <- left_join(x = filtered, y = data.frame("maf" = x[,avg_allele_freq], "wind" = x[,window_id]), by = "wind", stringsAsFactors = F)
  filtered$maf[filtered$maf == "."] <- NA
  filtered$maf <- as.numeric(filtered$maf)
  ## Flip the minor allelic fractions to a 0.5-1 range
  filtered$mafflipped <- abs(filtered$maf - 0.5)+0.5
  return(filtered)
}