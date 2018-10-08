## ************************************************************
## setwd('~/Dropbox/hedgerow_metacom')
rm(list=ls())
setwd('analysis/occupancy')
args <- commandArgs(trailingOnly=TRUE)
## args <- c("allInt","2500", "350", "filtering","all",1e2)
source('src/initialize.R')

## ************************************************************
## prep data
## ************************************************************

model.input <- prepOccModelInput(nzero=0,
                    threshold=2,
                    save.dir=save.dir,
                    spec,
                    sr.sched,
                    all.traits,
                    col.name.trait1 = "r.degree",
                    col.name.trait2 = "BodyLength",
                    HRarea= hr.area.sum, #sum.dist.area, ##spstats
                    natural.mat=nat.area.sum, ## natural
                    natural.decay=natural.decay,
                    HR.decay=HR.decay,
                    veg=by.site, #raw.flower.data,
                    load.inits=FALSE,
                    model.type=include.int,
                    col.name.div.type = "Div",## div.visits, Div
                    raw.flower=FALSE,
                    drop.li.ht=FALSE,
                    only.li.ht=FALSE)


## variables to plot
pdf.f(plotVariables,
      file=file.path(save.dir, 'figures/variables',
                     sprintf('%s%s%s.pdf', natural.decay, HR.decay, data.subset)),
      height= 6, width=3)


burnin <- 1e1*scale
niter <- (1e3)*scale
nthin <- 2
nchain <- 1

source(sprintf('src/models/complete_%s.R', include.int))

if(filtering){
    source('src/dynamicOcc.R')
    model.input$data$Z <- NULL
    model.input$inits$Z <- NULL
    ## We do not want any X element equal to NA or they will not be
    ## considered data and will be sampled.
    model.input$data$X[ is.na(model.input$data$X) ] <- -1000
    source(sprintf('src/models/complete_%s_filter.R', include.int))
}

testOccData()

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
plotPosteriors()
checkChains()
save(ms.ms.nimble, file=file.path(save.dir,
                    sprintf('runs/no_li_ht_%s_nimble_bees_%s_%s_%s.Rdata',
                                          data.subset,
                                          natural.decay, HR.decay,
                                          include.int)))

## ****************************************************************
## posterior probabilities and interaction plots
## ****************************************************************

if(!exists("ms.ms.nimble")){
    load(file=file.path(save.dir,
                        sprintf('runs/no_li_ht_%s_nimble_bees_%s_%s_%s.Rdata',
                                data.subset,
                                natural.decay, HR.decay,
                                include.int)))
}

if(is.list(ms.ms.nimble$samples)){
    samples.4.table <- do.call(rbind, ms.ms.nimble$samples)
} else{
    samples.4.table <- ms.ms.nimble$samples
}

means <- apply(samples.4.table, 2, mean)
se <- apply(samples.4.table, 2, sd)



pdf.f(f.plotInteractionsHRRemnant.k, file=file.path(save.dir,
                sprintf("figures/interactions/HRinteractions-k-%s-%s.pdf",
                                                 natural.decay, HR.decay)),
      width=6, height=11)


pdf.f(f.plotInteractionsFloralDiv, file=file.path(save.dir,
                 sprintf("figures/interactions/HRinteractions-fra-%s-%s.pdf",
                                                 natural.decay, HR.decay)),
      width=3, height=9)




pdf.f(plotInteractionsB, file=file.path(save.dir,
                 sprintf("figures/interactions/HRinteractions-B-%s-%s.pdf",
                                                 natural.decay, HR.decay)),
      width=3, height=9)


## posterior probs
if(!exists("ms.ms.model")){
    ms.ms.model <- nimbleModel(code=ms.ms.occ,
                               constants=model.input$constants,
                               data=model.input$data,
                               inits=model.input$inits,
                               check=FALSE,
                               calculate=FALSE)
}

samples.4.table <- samples.4.table[, colnames(samples.4.table) %in%
                      ms.ms.model$getNodeNames(includeData=FALSE,
                                               stochOnly=TRUE)]
## H param > 0
h1 <- apply(samples.4.table,
            2, function(x) sum(x > 0)/length(x))
## H param < 0
h2 <- apply(samples.4.table,
            2, function(x) sum(x < 0)/length(x))
posterior.probs <- cbind(h1,h2)
makeTable()


