## setwd('~/Dropbox/hedgerow_metacom')
setwd('analysis/networks')
source('src/initialize.R')
Sys.setenv('R_MAX_VSIZE'=64000000000)
library(parallel)
options(ncores=10)

species.lev <- function(x){
    present.x <- bipartite::empty(x)
    sl <- specieslevel(present.x,
                       index=c("degree", "betweenness", "closeness"))
    if(ncol(x) != ncol(present.x)){
        not.present <- colnames(x)[!colnames(x) %in%
                                   rownames(sl$"higher level")]
        all.na.mat <- matrix(NA, ncol=ncol(sl$'higher level'),
                             nrow=length(not.present))
        rownames(all.na.mat) <- not.present
        colnames(all.na.mat) <- colnames(sl$'higher level')
        sl$'higher level' <- rbind(sl$'higher level', all.na.mat)
        sl$'higher level' <-
            sl$'higher level'[order(rownames(sl$'higher level')),]
    }

    if(nrow(x) != nrow(present.x)){
        not.present <- rownames(x)[!rownames(x) %in%
                                   rownames(sl$"lower level")]
        all.na.mat <- matrix(NA, ncol=ncol(sl$'lower level'),
                             nrow=length(not.present))
        rownames(all.na.mat) <- not.present
        colnames(all.na.mat) <- colnames(sl$'lower level')
        sl$'lower level' <- rbind(sl$'lower level', all.na.mat)
        sl$'lower level' <-
            sl$'lower level'[order(rownames(sl$'lower level')),]
    }
    return(sl)
}

calcMetrics <- function(j, z.array){
    this.year <- z.array[,j,]
    specs.year <- species.lev(this.year)
    ll <- specs.year$`lower level`
    hl <- specs.year$`higher level`
    return(list(ll=ll, hl=hl))
}

## load the GIANT model file that tracks all of the latent states
load('~/Dropbox/hedgerow_metacom_saved/occupancy/runs/Z_all_2500_350.Rdata')

zs <- grepl("Z", colnames(ms.ms.nimble[[1]]))
samples.z <- lapply(ms.ms.nimble, function(x) x[,zs])

these.samples <- samples.z[[1]]
these.samples <- these.samples[sample(1:nrow(these.samples), size=100, replace=FALSE),]

## specs.year.sp <- vector(mode="list",
##                         length=nrow(these.samples))
## specs.year.site <- vector(mode="list",
##                           length=nrow(these.samples))

calcYearNet <- function(i, these.samples, model.input){
    print(i)
    z.array <- array(these.samples[i,],
                     dim=dim(model.input$data$Z),
                     dimnames=dimnames(model.input$data$Z))
    metrics <- mclapply(1:dim(z.array)[2], calcMetrics, z.array)
    specs.year.site <- lapply(metrics, function(x) x[["ll"]])
    specs.year.sp <- lapply(metrics, function(x) x[["hl"]])
    return(list(site=specs.year.site, sp=specs.year.sp))
}

year.nets <- mclapply(1:nrow(these.samples), calcYearNet,
                      these.samples, model.input)


## for(i in 1:nrow(these.samples)){
##     print(i)
##     z.array <- array(these.samples[i,],
##                      dim=dim(model.input$data$Z),
##                      dimnames=dimnames(model.input$data$Z))
##     metrics <- lapply(1:dim(z.array)[2], calcMetrics, z.array)
##     specs.year.site[[i]] <- lapply(metrics, function(x) x[["ll"]])
##     specs.year.sp[[i]] <- lapply(metrics, function(x) x[["hl"]])

## }


getPosteriorNet <- function(j, nets, model.input, site.sp){
    print(j)
    this.year <- lapply(nets, function(x) as.matrix(x[[j]]))
    this.year.array <- simplify2array(this.year)

    year.mean <- try(apply(this.year.array, c(1,2), mean,
                           na.rm=TRUE), silent=TRUE)
    if(inherits(year.mean, "try-error")) browser()
    year.sd <- apply(this.year.array, c(1,2), sd,
                          na.rm=TRUE)
    colnames(year.sd) <- paste0("sd.", colnames(year.sd))
    year <- as.data.frame(cbind(year.mean, year.sd))
    year$Year <- dimnames(model.input$data$Z)[[2]][j]
    if(site.sp == "site"){
        year$Site <- rownames(year)
    } else if(site.sp == "sp"){
        year$GenusSpecies <- rownames(year)
    }
    rownames(year) <- NULL
    return(year)
}

year.site <- lapply(1:dim(model.input$data$Z)[2], getPosteriorNet,
                               mclapply(year.nets, function(x) x$site),
                    model.input, "site")

year.sp <- lapply(1:dim(model.input$data$Z)[2], getPosteriorNet,
                               mclapply(year.nets, function(x) x$sp),
                    model.input, "sp")





all.years.site <- all.years.sp <-
    vector(mode="list", length=length(specs.year.site[[1]]))

for(j in 1:length(specs.year.site[[1]])){
    browser()
    this.year <- lapply(specs.year.site, function(x) as.matrix(x[[j]]))
    this.year.array <- simplify2array(this.year)
    year.mean.site <- apply(this.year.array, c(1,2), mean,
                                 na.rm=TRUE)
    year.sd.site <- apply(this.year.array, c(1,2), sd,
                          na.rm=TRUE)
    colnames(year.sd.site) <- paste0("sd.", colnames(year.sd.site))
    year.site <- as.data.frame(cbind(year.mean.site, year.sd.site))
    year.site$Year <- dimnames(model.input$data$Z)[[2]][j]
    year.site$Site <- rownames(year.site)
    rownames(year.site) <- NULL
    all.years.site[[j]] <- year.site

    this.year.sp <- lapply(specs.year.sp, function(x) as.matrix(x[[j]]))
    this.year.sp.array <- simplify2array(this.year.sp)
    year.mean.sp <- apply(this.year.sp.array, c(1,2), mean,
                               na.rm=TRUE)
    year.sd.sp <- apply(this.year.sp.array, c(1,2), sd, na.rm=TRUE)

    colnames(year.sd.sp) <- paste0("sd.", colnames(year.sd.sp))
    year.sp <- as.data.frame(cbind(year.mean.sp, year.sd.sp))
    year.sp$Year <- dimnames(model.input$data$Z)[[2]][j]
    year.sp$GenusSpecies <- rownames(year.sp)
    rownames(year.sp) <- NULL
    all.years.sp[[j]] <- year.sp
}


fin.years.sp <- do.call(rbind,  all.years.sp)
fin.years.sp <-  merge(fin.years.sp,
                       traits[,c("GenusSpecies", "r.degree",
                                 "MeanITD")])

fin.years.site <- do.call(rbind,  all.years.site)
fin.years.site <- merge(fin.years.site, by.site)

## **********************************************************
## pollinators
## **********************************************************

## anything outputted by specieslevel
ys <- c("k", "betweenness")

xvar.species.form <- paste0("scale(", xvar.species, ")")

formulas.species <-lapply(ys, function(x) {
    as.formula(paste(x, "~",
                     paste(paste(xvar.species.form , collapse="+"),
                           "(1|GenusSpecies)",
                           sep="+")))
})

## within a year, across sites
mod.years.pol <- list()
mod.years.pol[[1]] <- lmer(formulas.species[[1]], data=fin.years.sp,
                           weights=(1/fin.years.sp$sd.k))

fin.years.sp$inv.sd.betweeness <- 1/fin.years.sp$sd.betweenness
fin.years.sp$inv.sd.betweeness[
                 !is.finite(fin.years.sp$inv.sd.betweeness)] <- NA
fin.years.sp$inv.sd.betweeness <- log(fin.years.sp$inv.sd.betweeness)

mod.years.pol[[2]] <- lmer(formulas.species[[2]], data=fin.years.sp,
                           weights=(fin.years.sp$inv.sd.betweeness))

## within a site across years
mod.sites.pol <- lapply(formulas.species, function(x){
    lmer(x, data=specs.site.pol)
})


names(mod.years.pol) <- names(mod.sites.pol) <- ys

## floral degree positivly related to all measures of network
## importance in space/time
print("metacommunity pollinator spatial network")
print(lapply(mod.years.pol, summary))
print("metacommunity pollinator temporal network")
print(lapply(mod.sites.pol, summary))

## **********************************************************
## sites
## **********************************************************
xvar.site.form <- paste0("scale(", xvar.site, ")")
formulas.site <-lapply(ys, function(x) {
    as.formula(paste(x, "~",
                     paste(paste(xvar.site.form, collapse="+"),
                           "(1|Site)",
                           sep="+")))
})

## within a year across sites
mod.years.site <- lapply(formulas.site, function(x){
    lmer(x, data=fin.years.site)
})

mod.years.site[[1]] <- lmer(formulas.site[[1]], data=fin.years.site,
                           weights=(1/fin.years.site$sd.k))

fin.years.site$inv.sd.betweeness <- 1/fin.years.site$sd.betweenness
fin.years.site$inv.sd.betweeness[
                 !is.finite(fin.years.site$inv.sd.betweeness)] <- NA
fin.years.site$inv.sd.betweeness <- log(fin.years.site$inv.sd.betweeness)

mod.years.site[[2]] <- lmer(formulas.site[[2]], data=fin.years.site,
                           weights=(fin.years.site$inv.sd.betweeness))


## name them the same as the pollinators
names(mod.years.site)  <- ys

print("metacommunity patch network")
print(lapply(mod.years.site, summary))

save(mod.years.pol, mod.sites.pol, mod.years.site,
     file=sprintf('saved/mods/zmets_drop_li_ht%s.Rdata', drop.li.ht))
