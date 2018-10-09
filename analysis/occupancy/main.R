## ************************************************************
## setwd('~/Dropbox/hedgerow_metacom')
rm(list=ls())
setwd('analysis/occupancy')
args <- commandArgs(trailingOnly=TRUE)
source('src/initialize.R')

## ************************************************************
## prep data
## ************************************************************

model.input <- prepOccModelInput(nzero=0,
                    threshold=5,
                    spec=spec,
                    sr.sched=sr.sched,
                    traits=all.traits,
                    col.name.trait1 = "r.degree",
                    col.name.trait2 = "MeanITD",
                    HRarea= hr.area.sum, #sum.dist.area, ##spstats
                    natural.mat=nat.area.sum, ## natural
                    natural.decay=natural.decay,
                    HR.decay=HR.decay,
                    veg=by.site,
                    col.name.div.type = "Div",
                    save.dir=save.dir)

plotVariablesWrapper()
model.input <- prepModel()

## ## ## ****************************************************************
## ## run model
## ## *****************************************************************
burnin <- 1e1*scale
niter <- (1e3)*scale
nthin <- 2
nchain <- 3

## build R model
ms.ms.model <- nimbleModel(code=ms.ms.occ,
                           constants=model.input$constants,
                           data=model.input$data,
                           inits=model.input$inits,
                           check=FALSE,
                           calculate=FALSE)
## compile R model
C.model <- compileNimble(ms.ms.model)

## configure and build mcmc
mcmc.spec <- configureMCMC(ms.ms.model,
                           print=FALSE,
                           monitors = model.input$monitors,
                           enableWAIC = TRUE)
mcmc <- buildMCMC(mcmc.spec,
                  enableWAIC = TRUE)
C.mcmc <- compileNimble(mcmc, project = ms.ms.model)

## run model
ms.ms.nimble <- runMCMC(C.mcmc, niter=niter,
                        nchains=nchain,
                        nburnin=burnin,
                        WAIC=TRUE)
plotPosteriors()
checkChains()
save(ms.ms.nimble, model.input, file=file.path(save.dir,
                    sprintf('runs/%s_nimble_bees_%s_%s.Rdata',
                                          data.subset,
                                          natural.decay, HR.decay)))

source('src_plotting/plotResults.R')
## ************************************************************
