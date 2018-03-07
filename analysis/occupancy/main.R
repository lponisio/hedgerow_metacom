## ************************************************************
## setwd('~/Dropbox/hedgerow_metacom')
rm(list=ls())
setwd('analysis/occupancy')
args <- commandArgs(trailingOnly=TRUE)
args <- c("allInt","350","filtering","control",1e2)
source('src/initialize.R')

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
                    HRarea= hr.area.sum, #sum.dist.area, ##spstats
                    natural.mat=nat.area.sum, ## natural
                    natural.decay=natural.decay,
                    veg=by.site,
                    w.ypr=w.ypr,
                    load.inits=FALSE,
                    model.type=include.int,
                    col.name.div.type = "Div") ## div.visits, Div


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
plotPosteriors()
checkChains()
save(ms.ms.nimble, file=file.path(save.dir,
                    sprintf('runs/%s_nimble_bees_%s_%s.Rdata',
                                          data.subset,
                                          natural.decay,
                                          include.int)))

## ****************************************************************
## posterior probabilities
## *****************************************************************
if(is.null(ms.ms.nimble)){
    load(file=file.path(save.dir,
                        sprintf('runs/%s_nimble_bees_%s_%s.Rdata',
                                data.subset, natural.decay, include.int)))
}

if(is.list(ms.ms.nimble$samples)){
    samples.4.table <- do.call(rbind, ms.ms.nimble$samples)
} else{
    samples.4.table <- ms.ms.nimble$samples
}

if(is.null(ms.ms.model)){
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
## H param == 0
h0 <- apply(samples.4.table,
            2, function(x) dnorm(0, mean(x), sd(x))/length(x))
## H param < 0
h2 <- apply(samples.4.table,
            2, function(x) sum(x < 0)/length(x))
posterior.probs <- cbind(h1,h0,h2)
makeTable()


## variables to plot
if(data.subset=="all"){
    by.site <- by.site[by.site$Site %in% rownames(model.input$data$fra),]
    controls <- by.site$Site[by.site$SiteStatus == "control"]
    hedgerows <- by.site$Site[by.site$SiteStatus == "mature" | by.site$SiteStatus == "maturing"]
    pdf.f(plotVariables,
          file=file.path(save.dir, 'figures/variables/all.pdf'),
          height= 9, width=3)
}

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


