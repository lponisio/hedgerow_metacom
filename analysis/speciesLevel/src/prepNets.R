
## the purpose of this function is to break up data with many
## sites/years and prepare it for network analysis.

dropNet <- function(z){
  z[!sapply(z, FUN=function(q){
    any(dim(q) < 3)
  })]
}

breakNet <- function(spec.dat, site, year){
  ## puts data together in a list and removes empty matrices
  sites <- split(spec.dat, spec.dat[,site])
  networks <- lapply(sites, function(x){
    lapply(split(x, f=x[,year]), as.matrix)
  })
  ## formats data matrices appropriate for network analysis
  comms <- rapply(networks, function(y){
    samp2site.spp(site=y[,"PlantGenusSpecies"],
                  spp=y[,"GenusSpecies"],
                  abund=rep(1, nrow(y)))
  }, how="replace")
  adj.mat <- unlist(lapply(comms, dropNet), recursive=FALSE)
  if(year == "assem"){
    names(adj.mat) <- sub("\\.", "_", names(adj.mat),
                          perl=TRUE)
  }
  return(adj.mat)
}


getSpecies <- function(networks, FUN){
  species.site <- lapply(networks, FUN)
  site.plant <- rep(names(species.site), lapply(species.site, length))
  species <- data.frame(species=do.call(c, species.site),
                        siteStatus=site.plant,
                        site= sapply(strsplit(site.plant, "_"), function(x)
                          x[1]),
                        status= sapply(strsplit(site.plant, "_"), function(x)
                          x[2]))
  return(species)
}


## saves each element of a list with corresponding name
saveDats <- function(x, y, f.path){
  mapply(function(a, b)
         write.csv(a, file=file.path(f.path,
                        sprintf("%s.csv", b))),
         a=x, b=y,
         SIMPLIFY=FALSE)
}


## calculates the degree of species in a network
getDegree <- function(x, MARGIN){
  apply(x, MARGIN, function(y){
    length(y[y != 0])/length(y)
  })
}

## calculate various stats
calcStats <- function(x){
  means=mean(x)
  medians=median(x)
  mins <- min(x)
  maxs <- max(x)
  sds <- sd(x)
  return(c(mean=means,
           median=medians,
           min=mins,
           max=maxs,
           sd=sds))
}


##  species differences between early and late stages
getColExt <- function(dats){
  out <- lapply(unique(dats$site), function(x){
    this.site <- dats[dats$site == x,]
    early <- this.site$species[this.site$status == "early"]
    late <- this.site$species[this.site$status == "late"]
    colonists <- as.character(late[!late %in% early])
    ext <- as.character(early[!early %in% late])
    species <- c(colonists, ext)
    return(data.frame(species=species,
                      class=c(rep("colonist", length(colonists)),
                        rep("extinction", length(ext))),
                      site=rep(x, length(species))))
  })
  out <- do.call(rbind, out)
  return(out)
}


## number of species that interact
getCon <- function(x, INDEX){
  apply(x, INDEX, function(y) sum(y > 0))
}


## extreact specialization scores from specieslevel function and
## return data frame
getSpec <- function(species.lev, names.net, seps="_"){
  n.pp <- sapply(species.lev, nrow)
  pp <- c(unlist(sapply(species.lev, rownames)))
  names(pp) <- NULL
  all.pp <- do.call(rbind, species.lev)
  rownames(all.pp) <- NULL
  try(all.pp$GenusSpecies <- pp)
  all.pp$speciesType <- c(rep("pollinator", n.pp[1]),
                          rep("plant", n.pp[2]))
  all.pp$Site <- strsplit(names.net, seps)[[1]][1]
  all.pp$assem <- strsplit(names.net, seps)[[1]][2]
  return(all.pp)
}

