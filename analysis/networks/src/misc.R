
species.lev <- function(x){
    sl <- specieslevel(x)
    sl$'higher level'$tot.int <- colSums(x)
    sl$'lower level'$tot.int <- rowSums(x)
    return(sl)
}


checkDirExists <- function(save.dir){
    if(!dir.exists(save.dir)) {
        cat(paste("Needed dir",
                  save.dir,
                  "does not exist. OK to create? (Type 'yes' if ok.)"))
        okToMakeDir <- readlines()
        if(!identical(okToMakeDir, "yes"))
            stop("Stopping because permission to make save.dir denied.")
        dir.create(save.dir, showWarnings = FALSE)
    }
}

## standard error function
se <- function(x) sqrt(var(x)/(length(x)-1))

## write to a pdf
pdf.f <- function(f, file, ...) {
  cat(sprintf('Writing %s\n', file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}

## convert occurrence data into a site by species matrix
samp2site.spp <- function(site,spp,abund) {
  x <- tapply(abund, list(site=site,spp=spp), sum)
  x[is.na(x)] <- 0

  return(x)
}

## convert a site by species matrix into occurrences
comm.mat2sample <-  function (z) {
  temp <- data.frame(expand.grid(dimnames(z))[1:2],
                     as.vector(as.matrix(z)))
  temp <- temp[sort.list(temp[, 1]), ]
  data.frame(Site = temp[, 1], Abund = temp[, 3],
             GenSp = temp[, 2])
}

## add transparency to named colors
add.alpha <- function(col, alpha=0.2){
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
        rgb(x[1], x[2], x[3],
            alpha=alpha))
}

## inverse logit
inv.logit <- function(a){
  exp(a)/(exp(a) + 1)
}

## standardize a vector
standardize <- function(x)
  (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
