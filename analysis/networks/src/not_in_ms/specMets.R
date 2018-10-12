rm(list=ls())
setwd('~/Dropbox/hedgerow_metacom/analysis/speciesLevel')
source('src/initialize.R')

## **********************************************************
## species level metrics
## **********************************************************
spec.metrics <- calcSpec(nets, spec, spec.metric = "degree", 3)
save(spec.metrics, file=file.path(save.path, 'spec.metrics.Rdata'))

## linear models
load(file=file.path(save.path, 'spec.metrics.Rdata'))

## SiteStatus or ypr
xvar <- "s.Year"

## anything outputted by specieslevel
ys <- c("proportional.generality", "d", "degree",
        "tot.int",
        "partner.diversity")

formulas <-lapply(ys, function(x) {
  as.formula(paste(x, "~",
                   paste(paste(xvar, "overall.spec", sep="*"),
                             ## paste0("I(", xvar, "^2)*overall.spec"),
                         "(1|Site)",
                          "(1|GenusSpecies)",
                         sep="+")))
})

## formulas <-lapply(ys, function(x) {
##   as.formula(paste(x, "~",
##                    paste(xvar,
##                          "(1|Site)",
##                           "(1|GenusSpecies)",
##                          sep="+")))
## })

mod.pols.cont <- lapply(formulas, function(x){
  lmer(x,
       data=spec.metrics[spec.metrics$speciesType == "pollinator" &
         spec.metrics$SiteStatus == "control",])
})


mod.pols.mat <- lapply(formulas, function(x){
  lmer(x,
       data=spec.metrics[spec.metrics$speciesType == "pollinator" &
         spec.metrics$SiteStatus == "mature",])
})

mod.plants.cont <- lapply(formulas, function(x){
  lmer(x,
       data=spec.metrics[spec.metrics$speciesType == "plant" &
         spec.metrics$SiteStatus == "control",])
})

mod.plants.mat <- lapply(formulas, function(x){
  lmer(x,
       data=spec.metrics[spec.metrics$speciesType == "plant" &
         spec.metrics$SiteStatus == "mature",])
})

names(mod.pols) <- names(mod.plants) <- ys

lapply(mod.pols.cont, summary)
lapply(mod.pols.mat, summary)

lapply(mod.plants.cont, summary)
lapply(mod.plants.mat, summary)


## [[1]]
## [1] -456.7538

## [[2]]
## [1] -220.4532

## [[3]]
## [1] 2951.248

## [[4]]
## [1] 2951.248

## [[5]]
## [1] 1165.832


save(mod.pols, mod.plants, ys, file=file.path(save.path,
            sprintf('mods/spec.metrics_%s.Rdata', xvar)))

## **********************************************************
## degree distributions (abundance distributions)
## **********************************************************

baci.sites <- c("Barger", "Butler", "MullerB", "Sperandio", "Hrdy")
spec.metrics <- spec.metrics[spec.metrics$Site %in% baci.sites,]

layout(matrix(1:6, nrow=2))
cols <- rainbow(length(unique(spec.metrics$ypr)))
lapply(unique(spec.metrics$Site), function(x){
  print(x)
  this.spec.metrics <- spec.metrics[spec.metrics$Site == x, c("degree", "ypr")]
  plot(NA, ylim=c(0,0.8), xlim=c(0,25),
       ylab="Frequency",
       xlab="Abundance",
       main=x)
  for(i in 1:length(unique(this.spec.metrics$ypr))){
    this.ypr <- unique(this.spec.metrics$ypr)[i]
    print(this.ypr)
    points(density(this.spec.metrics$degree[this.spec.metrics$ypr == this.ypr]),
           col=cols[i], type="l", lwd=2)
  }
})

plot(NA, ylim=c(0,1), xlim=c(0,1), xaxt="n", yaxt="n", ylab="", xlab="")
legend("center", col=cols, lwd="2",
       legend=sort(unique(spec.metrics$ypr)),
       bty="n")
