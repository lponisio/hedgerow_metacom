
## add transparency to named colors
add.alpha <- function(col, alpha=0.2){
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
        rgb(x[1], x[2], x[3],
            alpha=alpha))
}


pdf.f <- function(f, file, ...) {
    ## writes plot to a pdf
    cat(sprintf('Writing %s\n', file))
    pdf(file, ...)
    on.exit(dev.off())
    f()
}


## return sorted unique values
id <- function(x) as.character(unique(sort(x)))

## standardize a vector
standardize <- function(x)
  (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)

## expit function
expit <- function(x) 1 / (1 + exp(-x))

## inverse logit
inv.logit <- function(a){
  exp(a)/(exp(a) + 1)
}


## function to clean up white-space in a column of data (replaces all
## instances of white-space with " " and empty cells with ""
fix.white.space <- function(d) {
  d <- as.character(d)
  remove.first <- function(s) substr(s, 2, nchar(s))
  d <- gsub("      ", " ", d, fixed=TRUE)
  d <- gsub("     ", " ", d, fixed=TRUE)
  d <- gsub("    ", " ", d, fixed=TRUE)
  d <- gsub("   ", " ", d, fixed=TRUE)
  d <- gsub("  ", " ", d, fixed=TRUE)

  tmp <- strsplit(as.character(d), " ")
  d <- sapply(tmp, function(x) paste(x, collapse=" "))

  first <- substr(d, 1, 1)
  d[first==" "] <- remove.first(d[first==" "])
  d
}

## load and return loaded object
load.local <- function(file) {
 v <- load(file)
 stopifnot(length(v) == 1)
 get(v)
}

## function to make pollinator visitation matrices
make.mats <- function(pollinator.id,
                      null.mat,
                      pollinator,
                      var1,
                      var2) {
  make.mat <- function(P) {
    var1 <- var1[pollinator==P]
    var2 <- var2[pollinator==P]
    m <- tapply(rep(1, length(var1)),
                list(sites=var1, dates=var2), sum)
    null.mat[rownames(m), colnames(m)][!is.na(m)] <- m[!is.na(m)]
    null.mat
  }
  mats <- lapply(pollinator.id, function(x) make.mat(x))
  names(mats) <- pollinator.id
  mats
}



## convert occurrence data into a site by species matrix
samp2site.spp <- function(site, yr, abund) {
  x <- tapply(abund, list(site=site,yr=yr), sum)
  x[is.na(x)] <- 0
  return(x)
}
