## ************************************************************
## setwd('~/Dropbox/hedgerow_metacom')
rm(list=ls())
setwd('analysis/occupancy')
source('src/initialize.R')
w.ypr <- FALSE
## allInt, "no_noncrop"
include.int <- "allInt"
## 350, 1000, 2500
natural.decay <- "350"
filtering <- FALSE


## spec <- spec[spec$SiteStatus == "mature" | spec$SiteStatus ==
##              "maturing",]
## spec <- spec[spec$SiteStatus == "control",]

## sr.sched <- sr.sched[sr.sched$Site %in% unique(spec$Site),]

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
                    natural.decay=natural.decay,
                    veg=by.site,
                    w.ypr=w.ypr,
                    load.inits=FALSE,
                    model.type=include.int,
                    col.name.div.type = "Div") ## div.visits, Div

scale <- 1e2
burnin <- 1e1*scale
niter <- (1e3)*scale
nthin <- 2
nchain <- 3


source(sprintf('src/models/complete_%s.R', include.int))

if(filtering){
    source('src/dynamicOcc.R')
    model.input$data$Z <- NULL
    model.input$inits$Z <- NULL
    ## We do not want any X element equal to NA or they will not be
    ## considered data and will be sampled.
    model.input$data$X[ is.na(model.input$data$X) ] <- -1000
    source(sprintf('src/models/complete_%s_filter_cppp.R', include.int))
}

## ## ****************************************************************
## ## not using mcmc suite
## ##
## *****************************************************************

ms.ms.model <- nimbleModel(code=ms.ms.occ,
                           constants=model.input$constants,
                           data=model.input$data,
                           inits=model.input$inits,
                           check=FALSE,
                           calculate=FALSE)
C.model <- compileNimble(ms.ms.model)

## configure and build mcmc
mcmc.spec <- configureMCMC(ms.ms.model,
                           print=FALSE,
                           monitors = model.input$monitors)
mcmc <- buildMCMC(mcmc.spec)
C.mcmc <- compileNimble(mcmc, project = ms.ms.model)

## run model
ms.ms.nimble <- runMCMC(C.mcmc, niter=niter,
                        nchains=nchain,
                        nburnin=burnin,
                        WAIC=TRUE)


save(ms.ms.nimble, file=file.path(save.dir,
                                  sprintf('runs/nimble_bees_%s_%s.Rdata',
                                          natural.decay, include.int)))



## ****************************************************************
## posterior probabilities
##
## *****************************************************************

load(file=file.path(save.dir,
                    sprintf('runs/nimble_bees_%s_%s.Rdata',
                            natural.decay, include.int)))


if(is.list(ms.ms.nimble$samples)){
    samples.4.cppp <- do.call(rbind, ms.ms.nimble$samples)
} else{
    samples.4.cppp <- ms.ms.nimble$samples
}


ms.ms.model <- nimbleModel(code=ms.ms.occ,
                           constants=model.input$constants,
                           data=model.input$data,
                           inits=model.input$inits,
                           check=FALSE,
                           calculate=FALSE)


samples.4.cppp <- samples.4.cppp[, colnames(samples.4.cppp) %in%
                                   ms.ms.model$getNodeNames(includeData=FALSE,
                                                            stochOnly=TRUE)]


## H param > 0
h1 <- apply(samples.4.cppp, 2, function(x) sum(x > 0)/length(x))
## H param == 0
h0 <- apply(samples.4.cppp, 2, function(x) dnorm(0, mean(x), sd(x))/length(x))
## H param < 0
h2 <- apply(samples.4.cppp, 2, function(x) sum(x < 0)/length(x))

posterior.probs <- cbind(h1,h0,h2)

save(posterior.probs, file=file.path(save.dir,
                                  sprintf('runs/post_probs_nimble_bees_%s_%s.Rdata',
                                          natural.decay, include.int)))


## ****************************************************************
## runn cppp on model
##
## *****************************************************************
source(sprintf('src/models/complete_%s_cppp.R', include.int))
source("src/cppp.R")

likeDiscFuncGenerator <- nimbleFunction(
    setup = function(model, ...){},
    run = function(){
        output <- -calculate(model)
        returnType(double(0))
        return(output)
    },
    contains = discrepancyFunction_BASE
)


mcmcCreator <- function(model){
  mcmc.spec <- configureMCMC(model,
                             print=FALSE,
                             monitors = colnames(samples.4.cppp))
  mcmc <- buildMCMC(mcmc.spec)
}


model.cppp <- runCPPP(ms.ms.model,
                      dataNames= "X",
                      discrepancyFunction= likeDiscFuncGenerator,
                      mcmcCreator = mcmcCreator,
                      origMCMCOutput= samples.4.cppp,
                      cpppControl = list(nBootstrapSDReps = 100),
                      mcmcControl = list(nMCMCiters = 10000,
                                         burnInProp= 0.1),
                      nCores = 4)

## ability to deal with a multi chain input
## catch unpassed arguments eariler, origMCMCOutput, discFunction,
## dataNames
## breaks if you monitor things other than parameters, why?

## ## *****************************************************************
## ## run in nimble using mcmc suite
## ## *****************************************************************
## input1 <- c(code=ms.ms.occ,
##             model.input)

## ms.ms.nimble <- compareMCMCs_withMonitors(input1,
##                                           MCMCs=c('nimble'),
##                                           niter=niter,
##                                           burnin = burnin,
##                                           thin=nthin,
##                                           summary=FALSE,
##                                           check=FALSE,
##                                           monitors=model.input$monitors)

## save(ms.ms.nimble, file=file.path(save.dir,
##                                   sprintf('runs/nimble_bees_%s.Rdata',
##                                           natural.decay)))


## ## *****************************************************************
## ## run in jags as a check
## ## *****************************************************************
## source('src/complete_jags.R')

## ms.ms.jags <- ms.ms(d=model.input)


