
dDynamicOccupancy <- nimbleFunction(
    ## DynamicOccupancy removes the z's and muZ's from the model and computes
    ## the probability of all reps over all years for one site, species
    run = function(x = double(2),
                   nrep = double(1),
                   psi1 = double(0),
                   phi = double(1),
                   gamma = double(1),
                   p = double(2),
                   log = double(0, default = 0)) {
        ## x is a year by rep matix
        ## prob of the occupied given p
        numObs1 <- sum(x[1,1:nrep[1]])
        ## prob of observation given occupided
        ProbOccAndCount <- psi1 * exp(sum(dbinom(x[1,1:nrep[1]], size = 1, p = p[1,1:nrep[1]], log = 1)))
        ## prob of the empty site
        ProbUnoccAndCount <- (1 - psi1) * (numObs1 == 0)
        ## probably of the observed states
        ProbCount <- ProbOccAndCount + ProbUnoccAndCount
        ProbOccGivenCount <- ProbOccAndCount / ProbCount
        ## occupided and persists, or unoccupied and colonizes
        ProbOccNextTime <- ProbOccGivenCount * phi[1] +
            (1-ProbOccGivenCount) * gamma[1]
        ll <- log(ProbCount)
        nyears <- dim(x)[1]
        for(t in 2:nyears) {
            numObs <- sum(x[t,1:nrep[t]])
            ProbOccAndCount <- ProbOccNextTime *
                exp(sum(dbinom(x[t,1:nrep[t]], size = 1, p = p[t,1:nrep[t]], log = 1)))
            ProbUnoccAndCount <- (1-ProbOccNextTime) * (numObs == 0)
            ProbCount <- ProbOccAndCount + ProbUnoccAndCount
            ProbOccGivenCount <- ProbOccAndCount / ProbCount
            ll <- ll + log(ProbCount)
            if(t < nyears) ProbOccNextTime <- ProbOccGivenCount * phi[t] +
                               (1-ProbOccGivenCount) * gamma[t]
        }
        if(log) return(ll)
        else return(exp(ll))
        returnType(double(0))
    }
)

## make matching changes here (at least to the input arguments)
## and check UserManual for non-needed "r" functions
rDynamicOccupancy <- nimbleFunction(
    run = function(n = double(),
                  nrep = double(1),
                   psi1 = double(0),
                   phi = double(1),
                   gamma = double(1),
                   p = double(2),
                   log = double(0, default = 0)) {
        nyear <- dim(p)[1]
        nreps <- dim(p)[2]
        ans <- matrix(1, nrow=nyear, ncol=nreps)
        returnType(double(2))
        return(ans)
    }
)

registerDistributions(list(
    dDynamicOccupancy = list(
        BUGSdist = "dDynamicOccupancy(nrep, psi1, phi, gamma, p)",
        Rdist = "dDynamicOccupancy(nrep, psi1, phi, gamma, p)",
        types = c('value = double(2)',
                  'nrep = double(1)',
                  'psi1 = double(0)',
                  'phi = double(1)',
                  'gamma = double(1)',
                  'p = double(2)'))
))



## test.occ <- function(x,
##                      nrep,
##                      psi1,
##                      phi,
##                      gamma,
##                      p){
##     ## prob of the occupied given p
##     numObs1 <- sum(x[1,], na.rm=TRUE)
##     ## prob of observation given occupided
##     ProbOccAndCount <- psi1 * exp(sum(dbinom(x[1,], size = 1, p = p[1,], log = 1)))
##     ## prob of the empty site
##     ProbUnoccAndCount <- (1 - psi1) * (numObs1 == 0)
##     ## probably of the observed states
##     ProbCount <- ProbOccAndCount + ProbUnoccAndCount
##     ProbOccGivenCount <- ProbOccAndCount / ProbCount
##     ## occupided and persists, or unoccupied and colonizes
##     ProbOccNextTime <- ProbOccGivenCount * phi[1] +
##         (1-ProbOccGivenCount) * gamma[1]
##     ll <- log(ProbCount)
##     nyears <- dim(x)[1]
##     for(t in 2:nyears) {
##         numObs <- sum(x[t,])
##         ProbOccAndCount <- ProbOccNextTime *
##             exp(sum(dbinom(x[t,], size = 1, p = p[t,], log = 1)))
##         ProbUnoccAndCount <- (1-ProbOccNextTime) * (numObs == 0)
##         ProbCount <- ProbOccAndCount + ProbUnoccAndCount
##         ProbOccGivenCount <- ProbOccAndCount / ProbCount
##         ll <- ll + log(ProbCount)
##         if(t < nyears) ProbOccNextTime <- ProbOccGivenCount * phi[t] +
##                            (1-ProbOccGivenCount) * gamma[t]
##     }
## }

## nsite <- model.input$constants$nsite
## nsp <- model.input$constants$nsp
## nyear <- model.input$constants$nyear
## max.nreps <- model.input$constants$max.nreps
## X <- model.input$data$X
## psi <- adrop(X[,,1,,drop=FALSE],drop=3)*0.5
## phi <- psi*0.25
## gamma <- psi*0.5
## p <- X*0.05

## for(site in 1:nsite) {
##     for(sp in 1:nsp) {
##     test.occ(X[site, 1:nyear, 1:max.nreps, sp], nrep[site, 1:nyear, sp],  psi[site,1,sp],
##     phi[site,1:(nyear-1),sp], gamma[site,1:(nyear-1),sp], p[site, 1:nyear, 1:max.nreps, sp])

##     }
## }
