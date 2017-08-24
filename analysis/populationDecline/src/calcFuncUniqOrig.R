
calcFuncUniqOrig <- function(traits, traits.2.keep,
                             bee.nonbee,
                             abund.col,
                             weights){
  these.traits <- traits[traits$bee.syr == bee.nonbee,
                         colnames(traits) %in% traits.2.keep]

  rownames(these.traits) <- these.traits$GenusSpecies
  these.traits$GenusSpecies <- NULL

  abund <- traits[traits$bee.syr == bee.nonbee, abund.col]
  abund[abund == 0 | is.na(abund)] <- 1
  names(abund) <- rownames(these.traits)

  coords <- dbFD(these.traits, abund, w=weights,
                 corr="cailliez", print.pco=T)$x.axes
  centr <- apply(coords, 2, weighted.mean, w = abund)

  ## add centroid coords as last row in dataframe
  coords2 <- data.frame(rbind(coords, centr))
  rownames(coords2)[dim(coords)[1]+1] <- "centr"

  ## create a matrix of distances between all species coordinates and
  ## centroid
  dists_centr <- as.matrix(dist(coords2, diag=TRUE, upper=TRUE))
  for (i in 1:dim(dists_centr)[1]) {
    dists_centr[i,i] <- NA
  }
  ## Originality: Distance to centroid of the species present
  originality <- dists_centr[, dim(dists_centr)[1]]
  originality <- originality[-(length(originality))]

  ## Uniqueness: nearest neighbour among present species
  ## the the minimum disance for each row
  uniq <- apply(dists_centr[-(dim(dists_centr)[1]),
                            -(dim(dists_centr)[1])],
                1, min, na.rm=TRUE)
  names(uniq) <- names(originality)
  out <- cbind(scale(uniq),
           scale(originality))
  colnames(out) <- c("uniq", "originality")
  return(out)
}
