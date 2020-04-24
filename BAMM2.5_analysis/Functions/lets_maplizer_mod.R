lets.maplizer.mod<- function(x, y, z, func = mean, ras = FALSE) {
  
  # Change factor to numbers
  if (is.factor(y)) {
    y <- as.numeric(levels(y))[y]
  }
  
  # Get the matrix without coordinates
  p <- x[[1]][, -(1:2)]
  
  for(i in 1:ncol(p)) {
    pos <- x[[3]][i] == z
    if (length(pos) > 0) {
      p[, i] <- p[, i] * y[pos][1]## This bit was modifying. The value of the species is a unique value. without the [1] is a redundancy. 
      pos2 <- p[, i] == 0
      p[pos2, i] <- NA
    } else {
      p[, i] <- NA
    }
  }
  
  resum <- apply(p, 1, func, na.rm=T)
  
  # Matrix of result 
  resultado <- cbind(x[[1]][, 1:2], resum)
  resu2 <- na.omit(resultado)
  
  # Name change
  name <- paste("Variable", as.character(substitute(func)),
                sep = "_")
  colnames(resultado)[3] <- name 
  
  # Return result with or without the raster
  if (ras) {
    r <- rasterize(resu2[, 1:2], x[[2]], resu2[, 3])
    return(list(Matrix = resultado, Raster = r))
  } else {
    return(resultado)
  }
}