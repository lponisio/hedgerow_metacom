ms.ms.occ <- nimbleCode({

    ## multi-species priors
    mu.p.0     ~ dnorm(0,0.001)
    mu.p.day.1 ~ dnorm(0,0.001)
    mu.p.day.2 ~ dnorm(0,0.001)

    sigma.p.0     ~ dunif(0,10)
    sigma.p.day.1 ~ dunif(0,10)
    sigma.p.day.2 ~ dunif(0,10)

    mu.phi.0  ~ dnorm(0,0.001)
    mu.gam.0  ~ dnorm(0,0.001)

    sigma.phi.0 ~ dunif(0,10)
    sigma.gam.0 ~ dunif(0,10)

    phi.traits     ~ dnorm(0,0.001)
    gam.traits     ~ dnorm(0,0.001)
    phi.rain        ~ dnorm(0,0.001)
    gam.rain        ~ dnorm(0,0.001)
    phi.rain.traits ~ dnorm(0,0.001)
    gam.rain.traits ~ dnorm(0,0.001)
    phi.hr.area ~ dnorm(0,0.001)
    gam.hr.area ~ dnorm(0,0.001)

    ## species-specific detectability
    for(sp in 1:nsp) {
        p.0[sp]     ~ dnorm(mu.p.0,     sd=sigma.p.0)
        p.day.1[sp] ~ dnorm(mu.p.day.1, sd=sigma.p.day.1)
        p.day.2[sp] ~ dnorm(mu.p.day.2, sd=sigma.p.day.2)
        for(site in 1:nsite) {
            for(yr in 1:nyear) {
                for(rep in 1:nrep[site,yr,sp]) {
                    logit(p[site,yr,rep,sp]) <-
                        p.0[sp] +
                        p.day.1[sp]*day[site,yr,rep,sp] +
                        p.day.2[sp]*day[site,yr,rep,sp]^2
                }
            }
        }
    }

    for(sp in 1:nsp) {
        phi.0[sp] ~ dnorm(mu.phi.0,sd=sigma.phi.0)
        gam.0[sp] ~ dnorm(mu.gam.0,sd=sigma.gam.0)

        for(site in 1:nsite) {

            ## occupancy in year 1
            psi[site,1,sp] <- frac.presence[site,sp]

            ## Z[site,1,sp] ~ dbern(psi[site,1,sp])

            ## ## detectability in year 1
            ## for(rep in 1:nrep[site,1,sp]) {
            ##     mu.p[site,1,rep,sp] <- Z[site,1,sp]*p[site,1,rep,sp]
            ##     X[site,1,rep,sp] ~ dbern(mu.p[site,1,rep,sp])
            ## }

            ## occupancy and detectability in subsequent years
            for(yr in 1:(nyear-1)) {
                phi[site,yr,sp] <-
                    phi.0[sp] +
                    phi.traits*traits[sp] +
                    phi.rain*rain[site, yr] +
                    phi.rain.traits*rain[site,yr]*traits[sp] +
                    phi.hr.area*HRarea[site]

                gam[site,yr,sp] <-
                    gam.0[sp] +
                    gam.traits*traits[sp] +
                    gam.rain*rain[site,yr] +
                    gam.rain.traits*rain[site, yr]*traits[sp] +
                    gam.hr.area*HRarea[site]

                ## logit(psi[site,yr+1,sp]) <-
                ##     Z[site,yr,sp] * phi[site,yr,sp] +
                ##     (1-Z[site,yr,sp]) * gam[site,yr,sp]

                ## Z[site,yr+1,sp] ~ dbern(psi[site,yr+1,sp])

                ## for(rep in 1:nrep[site,yr+1,sp]) {
                ##     mu.p[site,yr+1,rep,sp] <- Z[site,yr+1,sp]*p[site,yr+1,rep,sp]
                ##     X[site,yr+1,rep,sp] ~ dbern(mu.p[site,yr+1,rep,sp])
                ## }
            }
        }
    }

    ## iterate over sp
    gam[1:nsite,1:nyear,sp] <-
        gam.0[sp] +
        gam.traits*traits[sp] +
        gam.rain*rain[1:nsite,1:nyear] +
        gam.rain.traits*rain[1:nsite, 1:nyear]*traits[sp] +
        gam.hr.area*HRarea[1:nsite]
        ## similar for phi


## dDynamicOccupancy <- nimbleFunction(
##   ## I've checked that this runs and compiles, but I haven't tested if
##   ## I got the logic right!
##   run = function(x = double(2),
##     nrep = double(),
##     psi1 = double(),
##     phi = double(1),
##     gamma = double(1),
##     p = double(1),
## )) ## rewrite to use x in transpose manner; to use a vector nrep (by time); psi1 could be a vector (by site); make p a matrix

    ## iterate over site, sp:
X[site, 1:nyear, 1:MAXNREP, sp] ~ dDynamicOccupancy(nrep[site, 1:nyear, sp],  psi[1:nsite,1,sp] , phi[site,1:nyear,sp], gamma[site,1:nyear,sp], p[site, 1:nyear , 1:MAXNREP , sp])


    ## oops, need Z's for state dynamics
    ## for(sp in 1:nsp) {
    ##     for(yr in 1:(nyear-1)) {
    ##         X[1:nsite, yr+1, 1:nrep, sp] ~ dBernDetectionMatrix( occProb = psi[1:nsite,yr+1,sp],
    ##                                                             detectionProb = p[1:nsite,yr+1,1:nrep,sp],
    ##                                                             numReps = nrep[1:nsite,yr+1,sp])
    ##     }
    ## }


    ## calculate some useful stuff
    for(yr in 1:nyear) {
        for(site in 1:nsite) {
            N[site,yr] <- sum(Z[site,yr, 1:nsp])
        }
    }
})

