ms.ms.occ <- nimbleCode({

    ## species-specific  parameters
    for(sp in 1:nsp) {
        ## day
        p.0     ~ dnorm(0, 0.01)
        p.day.1 ~ dnorm(0, 0.01)
        p.day.2 ~ dnorm(0, 0.01)

        ## species specific
        phi.0 ~ dnorm(0, 0.01)
        gam.0 ~ dnorm(0, 0.01)

        ## hedgerow area
        phi.hr.area ~ dnorm(0, 0.01)
        gam.hr.area ~ dnorm(0, 0.01)

        ## natural habitat
        phi.nat.area ~ dnorm(0, 0.01)
        gam.nat.area ~ dnorm(0, 0.01)

        ## fra
        phi.fra ~ dnorm(0, 0.01)
        gam.fra ~ dnorm(0, 0.01)

        ## traits 1
        phi.traits1 ~ dnorm(0, 0.01)
        gam.traits1 ~ dnorm(0, 0.01)

        ## traits 2
        phi.traits2 ~ dnorm(0, 0.01)
        gam.traits2 ~ dnorm(0, 0.01)

        ## traits 1 * fra interaction
        phi.traits1.fra ~ dnorm(0, 0.01)
        gam.traits1.fra ~ dnorm(0, 0.01)

        for(site in 1:nsite) {
            for(yr in 1:nyear) {
                for(rep in 1:nrep[site,yr,sp]) {
                    logit(p[site,yr,rep,sp]) <-
                        p.0  +
                        p.day.1*day[site,yr,rep,sp] +
                        p.day.2*day[site,yr,rep,sp]^2
          }
        }
      }


    }


    for(sp in 1:nsp) {
        for(site in 1:nsite) {
            ## start off at the average for each species, site across years
            phi.site.mean[site,sp] <- mean(phi[site, 1:(nyear-1), sp])
            gam.site.mean[site,sp] <- mean(gam[site, 1:(nyear-1), sp])

            psi.1[site,sp] <- gam.site.mean[site,sp]/
                (1 - phi.site.mean[site,sp] + gam.site.mean[site,sp])

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
                    phi.0 +
                    phi.traits1*traits1[sp] +
                    phi.traits2*traits2[sp] +
                    phi.hr.area*HRarea[site] +
                    phi.nat.area*natural[site] +
                    phi.fra*fra[site, yr] +
                    phi.traits1.fra*fra[site, yr]*traits1[sp]

                gam[site,yr,sp] <-
                    gam.0 +
                    gam.traits1*traits1[sp] +
                    gam.traits2*traits2[sp] +
                    gam.hr.area*HRarea[site] +
                    gam.nat.area*natural[site] +
                    gam.fra*fra[site, yr] +
                    gam.traits1.fra*fra[site, yr]*traits1[sp]

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
    for(yr in 1:nyear) {
        for(site in 1:nsite) {
            N[site,yr] <- sum(Z[site,yr, 1:nsp])
        }
    }

    for(sp in 1:nsp) {
        phi.sp.mean[sp] <- mean(phi[1:nsite, 1:(nyear-1), sp])
        gam.sp.mean[sp] <- mean(gam[1:nsite, 1:(nyear-1), sp])
    }


})

