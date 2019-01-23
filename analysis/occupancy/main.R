## ************************************************************
## setwd('~/Dropbox/hedgerow_metacom')
setwd('analysis/occupancy')
args <- commandArgs(trailingOnly=TRUE)
source('src/initialize.R')
print(args)

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
                                 HRarea= hr.area.sum,
                                 natural.mat=nat.area.sum,
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

## ## build R model
ms.ms.model <- nimbleModel(code=ms.ms.occ,
                           constants=model.input$constants,
                           data=model.input$data,
                           inits=model.input$inits,
                           check=FALSE,
                           calculate=FALSE)
## compile R model
C.model <- compileNimble(ms.ms.model)

## for WAIC, need to monitor all stocastic nodes but this results is
## an absolutely huge output. Monitoring top-level nodes for non-WAIC
## purposes.

## configure and build mcmc
mcmc.spec <- configureMCMC(ms.ms.model,
                           print=FALSE,
                           monitors = model.input$monitors,
                           enableWAIC = FALSE)
mcmc <- buildMCMC(mcmc.spec,
                  enableWAIC = FALSE)
C.mcmc <- compileNimble(mcmc, project = ms.ms.model)

## run model
ms.ms.nimble <- runMCMC(C.mcmc, niter=niter,
                        nchains=nchain,
                        nburnin=burnin,
                        WAIC=FALSE)
save(ms.ms.nimble, model.input, file=file.path(save.dir,
                                               sprintf('runs/%s_%s_%s.Rdata',
                                                       data.subset,
                                                       natural.decay,
                                                       HR.decay)))

source('src_plotting/plotResults.R')
checkChains()
all.samples <- plotPosteriors()
************************************************************
