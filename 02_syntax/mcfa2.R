alpha <- function(dat){
  covar <- dat[lower.tri(dat)] #get unique covariances
  n <- ncol(dat) #number of items in the scale
  al <- ((sum(covar) / length(covar)) * n^2 / sum(dat))
  cat("alpha:", return(al), "\n")   
}

mcfa.input <- function(gp, dat){
  dat1 <- dat[complete.cases(dat), ]
  g <- dat1[ ,gp] #grouping
  freq <- data.frame(table(g))
  gn <- grep(gp, names(dat1)) #which column number is the grouping var
  dat2 <- dat1[ ,-gn] #raw only
  G <- length(table(g))
  n <- nrow(dat2)
  k <- ncol(dat2)
  scaling <- (n^2 - sum(freq$Freq^2)) / (n*(G - 1))
  varn <- names(dat1[ ,-gn])
  ms <- matrix(0, n, k)
  for (i in 1:k){
    ms[,i] <- ave(dat2[ ,i], g)
  }   
  cs <- dat2 - ms #deviation matrix, centered scores
  colnames(ms) <- colnames(cs) <- varn
  b.cov <- (cov(ms) * (n - 1)) / (G - 1) #group level cov matrix
  w.cov <- (cov(cs) * (n - 1)) / (n - G) #individual level cov matrix
  pb.cov <- (b.cov - w.cov)/scaling #estimate of pure/adjusted between cov matrix
  w.cor <- cov2cor(w.cov) #individual level cor matrix
  b.cor <- cov2cor(b.cov) #group level cor matrix
  pb.cor <- cov2cor(pb.cov) #estimate of pure between cor matrix
  icc <- round(diag(pb.cov) / (diag(w.cov) + diag(pb.cov)), 3) #iccs
  return(list(b.cov = b.cov, pw.cov = w.cov, ab.cov = pb.cov, pw.cor = w.cor,
              b.cor = b.cor, ab.cor = pb.cor,
              n = n, G = G, c. = scaling, sqc = sqrt(scaling),
              icc = icc, dfw = n - G, dfb = G, 
              pw.data = data.frame(cs),
              b.data = data.frame(ms)
  ) )
}
