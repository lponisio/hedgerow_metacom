
calc.metric <- function(dat.web) {
  ## calculates modularity
  calc.mod <- function(dat.web){
    ## converts a p-a matrix to a graph for modularity computation
    mut.adj <- function(x) {
      nr <- dim(x)[1]
      nc <- dim(x)[2]
      to.fill <- matrix(0, ncol=nc + nr, nrow=nc + nr)
      to.fill[1:nr,(nr+1):(nc+nr)] <- x
      adj.mat <- graph.adjacency(to.fill, mode= "upper", weighted=TRUE)
      return(adj.mat)
    }
    graph <- mut.adj(dat.web)
    weights <- as.vector(dat.web)
    weights <- weights[weights != 0]

    ## if matrix is binary, modularity calculate is not affected by
    ## the weights
    greedy <- modularity(graph,
                         membership(fastgreedy.community(graph,
                                                         weights=weights)),
                         weights=weights)

    random.walk <-  modularity(graph,
                               membership(walktrap.community(graph,
                                                             weights=weights)),
                               weights=weights)
    dendro <-  modularity(graph,
                          membership(edge.betweenness.community(graph,
                                                                weights=
                                                                weights)),
                          weights=weights)
    return(c(greedy, random.walk, dendro))
  }

  dat.web <- as.matrix(empty(dat.web))
  ## matrix of all the same number
  if(min(dat.web) == max(dat.web)){
    return(c(mets=rep(0, 1),
             mod.met=rep(0,3)))

  }else{
    ## wbinary=TRUE ensures binary matrices are calcualted like
    ## unweighted
    nodf <- nestednodf(dat.web,
                       weighted=TRUE,
                       wbinary=TRUE)$statistic["NODF"]
    mets <-  c(nodf,
               networklevel(dat.web, index="H2",
                            H2_integer=TRUE))
  }
  mod.met <- calc.mod(dat.web)
  return(c(mets, mod.met= mod.met))

}



##  function that computes summary statistics on simulated null matrices
##  (nulls simulated from web N times)
cor.metrics <- function (true.stat, null.stat, N) {
  ## calculate pvalues
  pvals <- function(stats, nnull){
    colSums(stats >= stats[rep(1, nrow(stats)),])/(nnull + 1)
  }
  ## calculate zvalues two different ways
  zvals <-function(stats){
    z.sd <- (stats[1,] -
             apply(stats, 2, mean, na.rm = TRUE))/
               apply(stats, 2, sd, na.rm = TRUE)
    z.sd[is.infinite(z.sd)] <- NA
    return(z.sd)
  }
  out.mets <- rbind(true.stat, null.stat)
  ## compute z scores
  zvalues <- zvals(out.mets)
  ## compute p-values
  pvalues <- pvals(out.mets, N)
  out <- c(true.stat, zvalues, pvalues)
  names(out) <- c("NODF", "H2",
                      "modularityG", "modularityR","modularityD",
                      "zNODF", "zH2", "zmodG", "zmodR", "zmodD",
                      "pNODF", "pH2", "pmodG", "pmodR", "pmodD")
  return(out)
}


prep.dat <- function(cor.stats, spec.dat){
  dats <- do.call(rbind, cor.stats)
  out <- data.frame(dats)
  out$Site <- sapply(strsplit(names(cor.stats), "\\."),
                 function(x) x[1])
  out$Year <-  sapply(strsplit(names(cor.stats), "\\."),
                 function(x) x[2])
  out$SiteStatus <- spec.dat$SiteStatus[match(paste(out$Site, out$Year),
                                              paste(spec.dat$Site,
                                                    spec.dat$Year))]
  out$ypr <- spec.dat$ypr[match(paste(out$Site, out$Year),
                                paste(spec.dat$Site,
                                      spec.dat$Year))]
  rownames(out) <- NULL
  return(out)
}



