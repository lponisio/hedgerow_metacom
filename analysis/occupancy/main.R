## ************************************************************
rm(list=ls())
setwd('~/Dropbox/hedgerow_metacom/analysis/occupancy')
source('src/misc.R')
library('abind')
library('nimble')
library('R2jags')
source('src/prep.R')
source('src/initialize.R')
w.rain <- FALSE
w.ypr <- FALSE


## ************************************************************
## prep data
## ************************************************************

model.input <- prep(nzero=0,
                    threshold=5,
                    phen=TRUE,
                    save.dir=save.dir,
                    spec,
                    sr.sched,
                    all.traits,
                    trait.ev1 = "r.degree",
                    trait.ev2 = "BodyLength",
                    rain=rain,
                    HRarea=sum.dist.area, ##spstats
                    natural.mat=nat.area.sum,
                    natural.decay="350",
                    veg=by.site,
                    w.rain=w.rain,
                    w.ypr=w.ypr,
                    load.inits=FALSE)

## buffer "d2000"

scale <- 1e2
burnin <- 1e1*scale
niter <- (1e3)*scale
nthin <- scale/10


source('src/complete_noRain.R')

input1 <- c(code=ms.ms.occ,
            model.input)

## *****************************************************************
## run in nimble, this works but the mid lelvel parameters do not mix
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
                                  'runs/bees_noRain.Rdata'))

## *****************************************************************
## run in jags, parent values error
## *****************************************************************
source('src/complete_jags.R')

ms.ms.jags <- ms.ms(d=model.input)

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
                    'gam.0',
                    'gam.hr.area',
                    'gam.nat.area',
                    'gam.fra',
                    'gam.traits1',
                    'gam.traits2',
                    'gam.traits1.fra')
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
## comparisons
## *****************************************************************
source('../../../occupancy/analysis/all/plotting.R')

load(file=file.path(save.dir, 'runs/crosslevel.Rdata'))
load(file=file.path(save.dir, 'runs/bees_noRain.Rdata'))


ms.ms.occ.all <- combine_MCMC_comparison_results(ms.ms.nimble[[1]],
                                                 ## ms.ms.crosslevel[[1]],
                                                 name = "ms.ms")


make_MCMC_comparison_pages(ms.ms.occ.all,
                           dir=file.path(save.dir, "figures/comparisons"))

checkChains(ms.ms.occ.all$ms.ms$samples,
            f.path = file.path(save.dir, "figures/chains/%s.pdf"))



## paramter groups
params <- dimnames(ms.ms.occ.all$ms.ms$summary)[[3]]
groups <- list()
i <- 1
while(length(params) > 0){
    id <- agrep(params[1], params, max.distance = 0.2)
    groups[[i]] <- params[id]
    params <- params[-id]
    i <- i + 1
}

f <- function(){
    layout(matrix(1:3, ncol=1))
    par(oma=c(3,1,1,1))
    cols <- rainbow(dim(ms.ms.occ.all$ms.ms$summary)[1])
    for(group in groups){
        if(length(group) > 1){
            for(i in 1:dim(ms.ms.occ.all$ms.ms$summary)[1]){
                this.samp <- ms.ms.occ.all$ms.ms$summary[i,,group]
                xs <- jitter(1:dim(this.samp)[2])
                if(i == 1){
                    plot(x=xs, y=this.samp["mean",], pch=16,
                         col=cols[i],
                         ylim=range(c(ms.ms.occ.all$ms.ms$summary[,'CI95_upp',
                                                                  group],
                                      ms.ms.occ.all$ms.ms$summary[,'CI95_low',
                                                                  group])),
                         xaxt="n",
                         ylab="Estimate",
                         xlab="")
                    abline(x=0, type="dashed")

                    axis(1, at=1:dim(this.samp)[2],
                         labels=FALSE)
                    text(x=1:dim(this.samp)[2],
                         y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
                         labels=group, srt=45, adj=1, xpd=TRUE)
                } else{
                    points(x=xs,
                           y=this.samp["mean",], pch=16, col=cols[i])
                }
                arrows(y1=this.samp['CI95_upp',],
                       y0=this.samp['CI95_low',],
                       x0=xs,
                       code=0, angle=90, length=0.02, lwd=1, col=cols[i])
            }
        }

        legend("topright", legend=dimnames(ms.ms.occ.all$ms.ms$summary)[[1]],
               pch=16, col=cols)
    }
}

pdf.f <- function(f, file, ...) {
    cat(sprintf('Writing %s\n', file))
    pdf(file, ...)
    on.exit(dev.off())
    f()
}

pdf.f(f, file=file.path(save.dir, "figures/allparams.pdf"),
      height=10, width=8.5 )

## *****************************************************************
## not using mcmc suite
## *****************************************************************


ms.ms.model <- nimbleModel(code=ms.ms.occ,
                           constants=model.input$constants,
                           data=model.input$data,
                           inits=model.input$inits,
                           check=FALSE)

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
