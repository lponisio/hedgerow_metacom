ms.ms.occ <- nimbleCode({

    ## multi-species priors
    mu.p.0     ~ dnorm(0,0.01)
    mu.p.day.1 ~ dnorm(0,0.01)
    mu.p.day.2 ~ dnorm(0,0.01)

    sigma.p.0     ~ dunif(0,10)
    sigma.p.day.1 ~ dunif(0,10)
    sigma.p.day.2 ~ dunif(0,10)

    mu.phi.0  ~ dnorm(0,0.01)
    mu.gam.0  ~ dnorm(0,0.01)

    sigma.phi.0 ~ dunif(0,10)
    sigma.gam.0 ~ dunif(0,10)

    mu.phi.hr.area  ~ dnorm(0,0.01)
    mu.gam.hr.area  ~ dnorm(0,0.01)

    sigma.phi.hr.area ~ dunif(0,10)
    sigma.gam.hr.area ~ dunif(0,10)

    mu.phi.nat.area  ~ dnorm(0,0.01)
    mu.gam.nat.area  ~ dnorm(0,0.01)

    sigma.phi.nat.area ~ dunif(0,10)
    sigma.gam.nat.area ~ dunif(0,10)

    mu.phi.rain  ~ dnorm(0,0.01)
    mu.gam.rain  ~ dnorm(0,0.01)

    sigma.phi.rain ~ dunif(0,10)
    sigma.gam.rain ~ dunif(0,10)

    mu.phi.fra  ~ dnorm(0,0.01)
    mu.gam.fra  ~ dnorm(0,0.01)

    sigma.phi.fra ~ dunif(0,10)
    sigma.gam.fra ~ dunif(0,10)

    mu.phi.traits1  ~ dnorm(0,0.01)
    mu.gam.traits1  ~ dnorm(0,0.01)

    sigma.phi.traits1 ~ dunif(0,10)
    sigma.gam.traits1 ~ dunif(0,10)

    mu.phi.traits2  ~ dnorm(0,0.01)
    mu.gam.traits2  ~ dnorm(0,0.01)

    sigma.phi.traits2 ~ dunif(0,10)
    sigma.gam.traits2 ~ dunif(0,10)

    mu.phi.traits1.fra  ~ dnorm(0,0.01)
    mu.gam.traits1.fra  ~ dnorm(0,0.01)

    sigma.phi.traits1.fra ~ dunif(0,10)
    sigma.gam.traits1.fra ~ dunif(0,10)

    mu.phi.traits1.rain  ~ dnorm(0,0.01)
    mu.gam.traits1.rain  ~ dnorm(0,0.01)

    sigma.phi.traits1.rain ~ dunif(0,10)
    sigma.gam.traits1.rain ~ dunif(0,10)

    ## species-specific detectability
    for(sp in 1:nsp) {

        ## day
        p.0[sp]     ~ dnorm(mu.p.0,     sd=sigma.p.0)
        p.day.1[sp] ~ dnorm(mu.p.day.1, sd=sigma.p.day.1)
        p.day.2[sp] ~ dnorm(mu.p.day.2, sd=sigma.p.day.2)

        ## species specific
        phi.0[sp] ~ dnorm(mu.phi.0,sd=sigma.phi.0)
        gam.0[sp] ~ dnorm(mu.gam.0,sd=sigma.gam.0)

        ## hedgerow area
        phi.hr.area[sp] ~ dnorm(mu.phi.hr.area,
                                sd=sigma.phi.hr.area)
        gam.hr.area[sp] ~ dnorm(mu.gam.hr.area,
                                sd=sigma.gam.hr.area)

        ## natural habitat
        phi.nat.area[sp] ~ dnorm(mu.phi.nat.area,
                                 sd=sigma.phi.nat.area)
        gam.nat.area[sp] ~ dnorm(mu.gam.nat.area,
                                 sd=sigma.gam.nat.area)

        ## rain
        phi.rain[sp] ~ dnorm(mu.phi.rain,
                             sd=sigma.phi.rain)
        gam.rain[sp] ~ dnorm(mu.gam.rain,
                             sd=sigma.gam.rain)

        ## fra
        phi.fra[sp] ~ dnorm(mu.phi.fra,
                            sd=sigma.phi.fra)
        gam.fra[sp] ~ dnorm(mu.gam.fra,
                            sd=sigma.gam.fra)

        ## traits 1
        phi.traits1[sp] ~ dnorm(mu.phi.traits1,
                                sd=sigma.phi.traits1)
        gam.traits1[sp] ~ dnorm(mu.gam.traits1,
                                sd=sigma.gam.traits1)

        ## traits 2
        phi.traits2[sp] ~ dnorm(mu.phi.traits2,
                                sd=sigma.phi.traits2)
        gam.traits2[sp] ~ dnorm(mu.gam.traits2,
                                sd=sigma.gam.traits2)

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
        for(site in 1:nsite) {
            phi.yr.mean[site,sp] <- mean(phi[site, 1:(nyear-1), sp])
            gam.yr.mean[site,sp] <- mean(gam[site, 1:(nyear-1), sp])
            psi.1[site,sp] <- gam.yr.mean[site,sp]/
                (1 - phi.yr.mean[site,sp] + gam.yr.mean[site,sp])
        }
    }


    for(sp in 1:nsp) {
        for(site in 1:nsite) {

            ## occupancy in year 1
            psi[site,1,sp] <- psi.1[site,sp]
            Z[site,1,sp] ~ dbern(psi[site,1,sp])

            ## detectability in year 1
            for(rep in 1:nrep[site,1,sp]) {
                mu.p[site,1,rep,sp] <- Z[site,1,sp]*p[site,1,rep,sp]
                X[site,1,rep,sp] ~ dbern(mu.p[site,1,rep,sp])
            }

            ## occupancy in subsequent years
            for(yr in 1:(nyear-1)) {
                phi[site,yr,sp] <-
                    phi.0[sp] +
                    phi.traits1[sp]*traits1[sp] +
                    phi.traits2[sp]*traits2[sp] +
                    phi.rain[sp]*rain[yr] +
                    phi.hr.area[sp]*HRarea[site] +
                    phi.nat.area[sp]*natural[site, yr] +
                    phi.fra[sp]*fra[site, yr] +
                    phi.traits1.fra[sp]*fra[site, yr]*traits1[sp] +
                    phi.traits1.rain[sp]*rain[yr]*traits1[sp]

                gam[site,yr,sp] <-
                    gam.0[sp] +
                    gam.traits1[sp]*traits1[sp] +
                    gam.traits2[sp]*traits2[sp] +
                    gam.rain[sp]*rain[yr] +
                    gam.hr.area[sp]*HRarea[site] +
                    gam.nat.area[sp]*natural[site, yr] +
                    gam.fra[sp]*fra[site, yr] +
                    gam.traits1.fra[sp]*fra[site, yr]*traits1[sp] +
                    gam.traits1.rain[sp]*rain[yr]*traits1[sp]

                logit(psi[site,yr+1,sp]) <-
                    Z[site,yr,sp] * phi[site,yr,sp] +
                    (1-Z[site,yr,sp]) * gam[site,yr,sp]

                Z[site,yr+1,sp] ~ dbern(psi[site,yr+1,sp])

                ## detectability in != year 1
                for(rep in 1:nrep[site,yr+1,sp]) {
                    mu.p[site,yr+1,rep,sp] <- Z[site,yr+1,sp]*p[site,yr+1,rep,sp]
                    X[site,yr+1,rep,sp] ~ dbern(mu.p[site,yr+1,rep,sp])
                }
            }
        }
    }

    ## calculate some useful stuff
    ## for(yr in 1:nyear) {
    ##     for(site in 1:nsite) {
    ##         N[site,yr] <- sum(Z[site,yr, 1:nsp])
    ##     }
    ## }

    for(site in 1:nsite) {
        phi.yr.sp.mean[site] <- mean(phi[site, 1:(nyear-1), 1:nsp])
        gam.yr.sp.mean[site] <- mean(gam[site, 1:(nyear-1), 1:nsp])
    }

    for(sp in 1:nsp) {
        phi.sp.mean[sp] <- mean(phi[1:nsite, 1:(nyear-1), sp])
        gam.sp.mean[sp] <- mean(gam[1:nsite, 1:(nyear-1), sp])
    }


})

