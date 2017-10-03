## ************************************************************
rm(list=ls())
setwd('~/Dropbox/hedgerow_metacom/analysis/occupancy')
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
                    natural.mat=nat.area.sum,
                    natural.decay="350",
                    veg=by.site,
                    w.ypr=w.ypr,
                    load.inits=FALSE)

## buffer "d2000"

scale <- 1e2
burnin <- 1e1*scale
niter <- (1e3)*scale
nthin <- 5

source('src/complete_noRain.R')

input1 <- c(code=ms.ms.occ,
            model.input)

## *****************************************************************
## run in nimble
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
                                  'runs/nimble_bees_noRain_short.Rdata'))

## *****************************************************************
## cross level sampler
## *****************************************************************
source('../../../occupancy/analysis/all/samplers/sampler_crossLevel_new.R')

MCMCdefs.crosslevel <- list('nimbleCrosslevel' = quote({
    customSpec <- configureMCMC(Rmodel)
    customSpec$removeSamplers('Z')
    ## find node names of each species for random effects
    base.names <- c('p.0',
                    'p.day.1',
                    'p.day.2',
                    'phi.0',
                    'phi.hr.area',
                    'phi.nat.area',
                    'phi.fra',
                    'phi.traits1',
                    'phi.traits2',
                    'phi.traits1.fra',
                    'phi.hr.area.fra',
                    'phi.nat.area.fra',
                    'gam.0',
                    'gam.hr.area',
                    'gam.nat.area',
                    'gam.fra',
                    'gam.traits1',
                    'gam.traits2',
                    'gam.traits1.fra',
                    'gam.hr.area.fra',
                    'gam.nat.area.fra',
                    )
    customSpec$removeSamplers(base.names)
    customSpec$addSampler(target = base.names,
                          type ='sampler_crossLevelBinary',
                          print=FALSE)
    customSpec
}))


ms.ms.crosslevel <- compareMCMCs_withMonitors(input1,
                                              MCMCs=c('nimbleCrosslevel'),
                                              MCMCdefs = MCMCdefs.crosslevel,
                                              niter= niter,
                                              burnin = burnin,
                                              thin=nthin,
                                              summary=FALSE,
                                              check=FALSE,
                                              monitors=model.input$monitors)

save(ms.ms.crosslevel, file=file.path(save.dir,
                                      'runs/crosslevel.Rdata'))

## *****************************************************************
## not using mcmc suite
## *****************************************************************


ms.ms.model <- nimbleModel(code=ms.ms.occ,
                           constants=model.input$constants,
                           data=model.input$data,
                           inits=model.input$inits,
                           check=FALSE,
                           calculate=FALSE)

## configure and build mcmc
mcmc.spec <- configureMCMC(ms.ms.model,
                           print=FALSE,
                           monitors = model.input$monitors)
mcmc <- buildMCMC(mcmc.spec)

## compile model in C++
C.model <- compileNimble(ms.ms.model)
C.mcmc <- compileNimble(mcmc, project = ms.ms.model)

## run model
C.mcmc$run(niter)

samples <- as.matrix(C.mcmc$mvSamples)

means <- apply(samples, 2, mean)


## *****************************************************************
## run in jags as a check
## *****************************************************************
source('src/complete_jags.R')

ms.ms.jags <- ms.ms(d=model.input)
