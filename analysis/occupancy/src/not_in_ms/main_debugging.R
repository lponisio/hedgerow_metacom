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
filtering <- TRUE

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
                    col.name.div.type = "div.visits") ## div.visits, Div

source(sprintf('src/models/complete_%s.R', include.int))

if(filtering){
    source('src/dynamicOcc.R')
    model.input$data$Z <- NULL
    model.input$inits$Z <- NULL
    ## We do not want any X element equal to NA or they will not be considered data and
    ## will be sampled.
    model.input$data$X[ is.na(model.input$data$X) ] <- -1000
    source(sprintf('src/models/complete_%s_filter.R', include.int))
}

## ## ****************************************************************
## ## not using mcmc suite
## ##
## *****************************************************************
## until next release so can use waic
install_github("nimble-dev/nimble",
               ref = "devel",
               subdir = "packages/nimble")

ms.ms.model <- nimbleModel(code=ms.ms.occ,
                           constants=model.input$constants,
                           data=model.input$data,
                           inits=model.input$inits,
                           check=FALSE,
                           calculate=FALSE)
C.model <- compileNimble(ms.ms.model)

## This steps through the calculations in order
for(node in C.model$getMaps('nodeNamesLHSall')) {
    writeLines(node)
    C.model$calculate(node)
}
## It crashes on X[2, 1:10, 1:7, 1]
## Now I've reloaded and will look at X[2, 1:10, 1:7, 1]
C.model$X[2, 1:10, 1:7, 1] ## This does have missing data from row 1.

C.model$calculate() ## NA!
dim(C.model$logProb_X)
C.model$logProb_X[, 1, 1, ] ## All NA
XnodeNames <- C.model$expandNodeNames('X')
C.model$getLogProb(XnodeNames[1])
dim(C.model$phi)
C.model$X[1, , , 1] ## This has ragged entries, so nrep may be important
C.model$phi[, 1, ]
C.model$gam[, 1, ]
## Focus on site 1, sp 1:
C.model$nrep[1, , 1] ## Hmmm, I wonder if 0 trips at error
C.model$psi[1, 1, 1] ## This has ragged entries too.
C.model$phi[1, , 1]
C.model$gam[1, , 1]
C.model$p[1, , , 1]
C.model$calculate("psi[1,1,1]")
C.model$calculate("psi.1[1,1]")
nonXnodes <- ms.ms.model$getNodeNames(includeData = FALSE)
C.model$calculate(nonXnodes)
site <- 1
sp <- 1
nyear <- 10
max.nreps <- 7

with(C.model,
     dDynamicOccupancy(X[site, 1:nyear, 1:max.nreps, sp],
                       nrep=nrep[site, 1:nyear, sp],
                       psi1=psi[site,1,sp],
                       phi=phi[site,1:(nyear-1),sp],
                       gamma=gam[site,1:(nyear-1),sp],
                       p=p[site, 1:nyear, 1:max.nreps, sp])
     ) ## Produces NaNs, so let's debug

debugonce(dDynamicOccupancy)
with(C.model,
     dDynamicOccupancy(X[site, 1:nyear, 1:max.nreps, sp],
                       nrep=nrep[site, 1:nyear, sp],
                       psi1=psi[site,1,sp],
                       phi=phi[site,1:(nyear-1),sp],
                       gamma=gam[site,1:(nyear-1),sp],
                       p=p[site, 1:nyear, 1:max.nreps, sp])
     ) ## Produces NaNs, so let's debug

## There are phi's greater than 1 and gam's less than 0.
## I think these are supposed to be expit()ed, right?
## I am just going to do this to see if dDynamicOccupancy works
## when given reasonable values.
with(C.model,
     dDynamicOccupancy(X[site, 1:nyear, 1:max.nreps, sp],
                       nrep=nrep[site, 1:nyear, sp],
                       psi1=expit(psi[site,1,sp]),
                       phi=expit(phi[site,1:(nyear-1),sp]),
                       gamma=expit(gam[site,1:(nyear-1),sp]),
                       p=p[site, 1:nyear, 1:max.nreps, sp],
                       log = TRUE)
     ) ## Oh, log = TRUE is how it will be called from MCMC

ms.ms.model$expandNodeNames("nrep")
## configure and build mcmc
mcmc.spec <- configureMCMC(ms.ms.model,
                           print=FALSE,
                           monitors = model.input$monitors)
mcmc <- buildMCMC(mcmc.spec)
C.mcmc <- compileNimble(mcmc, project = ms.ms.model)

niter <- 1
burnin <- 0
nchains <- 1
## run model
ms.ms.nimble <- runMCMC(C.mcmc, niter=niter,
                        nchains=nchain,
                        nburnin=burnin,
                        WAIC=FALSE)


save(ms.ms.nimble, file=file.path(save.dir,
                                  sprintf('runs/nimble_bees_%s_%s.Rdata',
                                          natural.decay, include.int)))

## checking right hand side values
C.model$HRarea
C.model$natural

## ## ****************************************************************
## ## runn cppp on model
## ##
## *****************************************************************
install_github("nimble-dev/nimble",
               ref = "add_CPPP_etc",
               subdir = "packages/nimble")

source(sprintf('src/models/complete_%s_cppp.R', include.int))

load(file=file.path(save.dir,
                    sprintf('runs/nimble_bees_%s_%s.Rdata',
                            natural.decay, include.int)))

ms.ms.model <- nimbleModel(code=ms.ms.occ,
                           constants=model.input$constants,
                           data=model.input$data,
                           inits=model.input$inits,
                           check=FALSE,
                           calculate=FALSE)


likeDiscFuncGenerator <- nimbleFunction(
    setup = function(model, ...){},
    run = function(){
        output <- -calculate(model)
        returnType(double(0))
        return(output)
    },
    contains = discrepancyFunction_BASE
)

samples.4.cppp <- do.call(rbind, ms.ms.nimble)

## drop things ther were monitored that were not parameters
samples.4.cppp <- samples.4.cppp[, colnames(samples.4.cppp) %in%
                                   ms.ms.model$getNodeNames(includeData=FALSE,
                                              stochOnly=TRUE)]

model.cppp <- runCPPP(ms.ms.model,
                      dataNames= "X",
                      discrepancyFunction= likeDiscFuncGenerator,
                      nCores=1,
                      origMCMCOutput= samples.4.cppp)

## ability to deal with a multi chain input
## catch unpassed arguments eariler, origMCMCOutput, discFunction,
## dataNames
## breaks if you monitor things other than parameters, why?


