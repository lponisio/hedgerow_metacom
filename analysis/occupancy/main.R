## ************************************************************
## setwd('~/Dropbox/hedgerow_metacom')
rm(list=ls())
setwd('analysis/occupancy')
library('abind')
library('nimble')
library('R2jags')
source('src/misc.R')
source('src/plotting.R')
source('src/prep.R')
source('src/initialize.R')
w.ypr <- FALSE

## ************************************************************
## prep data
## ************************************************************

model.input <- prepOccModelInput(nzero=0,
                    threshold=5,
                    save.dir=save.dir,
                    spec,
                    sr.sched,
                    all.traits,
                    col.name.trait1 = "r.degree",
                    col.name.trait2 = "BodyLength",
                    HRarea=sum.dist.area, ##spstats
                    natural.mat=nat.area.sum, ## natural
                    kinda.natural.mat=NULL, ## kinda natural
                    natural.decay="350",
                    veg=by.site,
                    w.ypr=w.ypr,
                    load.inits=FALSE)

scale <- 1e3
burnin <- 1e1*scale
niter <- (1e3)*scale
nthin <- 2
nchain <- 3

model.input$data$kinda.natural <- NULL

source('src/complete_noRain.R')

input1 <- c(code=ms.ms.occ,
            model.input)

## *****************************************************************
## run in nimble using mcmc suite
## *****************************************************************
ms.ms.nimble <- compareMCMCs_withMonitors(input1,
                                          MCMCs=c('nimble'),
                                          niter=niter,
                                          burnin = burnin,
                                          thin=nthin,
                                          summary=FALSE,
                                          check=FALSE,
                                          monitors=model.input$monitors)

save(ms.ms.nimble, file=file.path(save.dir,
                                  'runs/nimble_bees.Rdata'))

## *****************************************************************
## cross level sampler
## *****************************************************************
## source('../../../occupancy/analysis/all/samplers/sampler_crossLevel_new.R')

## MCMCdefs.crosslevel <- list('nimbleCrosslevel' = quote({
##     customSpec <- configureMCMC(Rmodel)
##     customSpec$removeSamplers('Z')
##     ## find node names of each species for random effects
##     base.names <- c(p.0,
##       p.day.1,
##       p.day.2,
##       phi.0,
##       phi.hr.area,
##       phi.nat.area,
##       phi.kinda.nat.area,
##       phi.fra,
##       phi.k,
##       phi.B,
##       phi.nat.area.fra,
##       phi.hr.area.fra,
##       phi.nat.area.k,
##       phi.hr.area.k,
##       gam.0,
##       gam.hr.area,
##       gam.nat.area,
##       gam.kinda.nat.area,
##       gam.fra,
##       gam.k,
##       gam.B,
##       gam.hr.area.fra,
##       gam.nat.area.fra,
##       gam.hr.area.k,
##       gam.nat.area.k)

##     customSpec$removeSamplers(base.names)
##     customSpec$addSampler(target = base.names,
##                           type ='sampler_crossLevelBinary',
##                           print=FALSE)
##     customSpec
## }))


## ms.ms.crosslevel <- compareMCMCs_withMonitors(input1,
##                                               MCMCs=c('nimbleCrosslevel'),
##                                               MCMCdefs = MCMCdefs.crosslevel,
##                                               niter= niter,
##                                               burnin = burnin,
##                                               thin=nthin,
##                                               summary=FALSE,
##                                               check=FALSE,
##                                               monitors=model.input$monitors)

## save(ms.ms.crosslevel, file=file.path(save.dir,
##                                       'runs/crosslevel.Rdata'))

## ## *****************************************************************
## ## run in jags as a check
## ## *****************************************************************
## source('src/complete_jags.R')

## ms.ms.jags <- ms.ms(d=model.input)




## ## *****************************************************************
## ## not using mcmc suite
## ## *****************************************************************

## ms.ms.model <- nimbleModel(code=ms.ms.occ,
##                            constants=model.input$constants,
##                            data=model.input$data,
##                            inits=model.input$inits,
##                            check=FALSE,
##                            calculate=FALSE)

## ## configure and build mcmc
## mcmc.spec <- configureMCMC(ms.ms.model,
##                            print=FALSE,
##                            monitors = model.input$monitors)
## mcmc <- buildMCMC(mcmc.spec)

## ## compile model in C++
## C.model <- compileNimble(ms.ms.model)
## C.mcmc <- compileNimble(mcmc, project = ms.ms.model)

## ## run model
## C.mcmc$run(niter)

## samples <- as.matrix(C.mcmc$mvSamples)

## means <- apply(samples, 2, mean)
