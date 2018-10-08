## ************************************************************
## setwd('~/Dropbox/hedgerow_metacom')
rm(list=ls())
setwd('analysis/occupancy')
args <- commandArgs(trailingOnly=TRUE)
source('src/initialize.R')
source("src/cppp.R")
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
                      cpppControl = list(nBootstrapSDReps = scale/10),
                      mcmcControl = list(nMCMCiters = scale*10,
                                         burnInProp= 0.1),
                      nCores = ncores)

## ability to deal with a multi chain input
## catch unpassed arguments eariler, origMCMCOutput, discFunction,
## dataNames
## breaks if you monitor things other than parameters, why?


