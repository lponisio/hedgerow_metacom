## ************************************************************
rm(list=ls())
setwd('~/Dropbox/hedgerow_metacom/analysis/occupancy')
source('src/misc.R')
library('abind')
library('nimble')
library('arm')
source('src/prep.R')
source('src/initialize.R')
load(file=file.path(save.dir, 'runs/nimble.Rdata'))

nimble.mcmc <- as.matrix(ms.ms.nimble[[1]]$samples["nimble",,])

coef.vect <- ms.ms.nimble[[1]]$summary["nimble",,][1,]
sd.vect <- ms.ms.nimble[[1]]$summary["nimble",,][3,]
short.names <- colnames(ms.ms.nimble[[1]]$summary["nimble",,])

mus <- grepl("^mu", names(coef.vect))
sigmas <- grepl("^sigma", names(coef.vect))

mu.coef <- coef.vect[mus]
sigma.coef <- coef.vect[sigmas]
mu.sd <- sd.vect[mus]
sigma.sd <- sd.vect[sigmas]
mu.names <- short.names[mus]
sigma.names <- short.names[sigmas]

coefplot(mu.coef, mu.sd,
         varnames = mu.names,
         main = "",
         xlim=c(-10, 10))

## using ggplot2 from scratch
library(ggplot2)
nimble.mat <- as.matrix(nimble.mcmc)
nimble.dat <- as.data.frame(nimble.mat)
coef.vect <- apply(nimble.dat, 1, mean)
lower.vect <- apply(nimble.dat, 1,
                    function(x) quantile(x, probs = c(0.025)))
upper.vect <- apply(nimble.dat, 1,
                    function(x) quantile(x, probs = c(0.975)))
long.names <- c("Intercept", "Diversity", "Mobility", "Deviance")
plot.dat <- data.frame(coef.vect,
                       lower.vect,
                       upper.vect,
                       long.names)[c(2,3), ]

p <- ggplot(data = plot.dat,
            aes(x = coef.vect, y = long.names)) + geom_point() + geom_segment(aes(x = lower.vect,
                   xend = upper.vect,
                   y = long.names, yend = long.names))

p <- p + geom_vline(xintercept = 0, colour = "blue", linetype = 2)+ xlab("Posterior estimates") + ylab("")
p
