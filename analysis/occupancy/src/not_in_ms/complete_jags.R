ms.ms <- function(d,
                  ni=1100,
                  nt=100,
                  nb=100,
                  nc=3) {

  model.jags <- function() {

    ## multi-species priors
    mu.p.0     ~ dnorm(0,0.01)
    mu.p.day.1 ~ dnorm(0,0.01)
    mu.p.day.2 ~ dnorm(0,0.01)
    sigma.p.0     ~ dunif(0,10)
    sigma.p.day.1 ~ dunif(0,10)
    sigma.p.day.2 ~ dunif(0,10)
    tau.p.0 <- 1/(sigma.p.0*sigma.p.0)
    tau.p.day.1 <- 1/(sigma.p.day.1*sigma.p.day.1)
    tau.p.day.2 <- 1/(sigma.p.day.2*sigma.p.day.2)

    mu.phi.0  ~ dnorm(0,0.01)
    mu.gam.0  ~ dnorm(0,0.01)
    sigma.phi.0 ~ dunif(0,10)
    sigma.gam.0 ~ dunif(0,10)
    tau.phi.0 <-  1/(sigma.phi.0*sigma.phi.0)
    tau.gam.0 <-  1/(sigma.gam.0*sigma.gam.0)

    mu.phi.hr.area  ~ dnorm(0,0.01)
    mu.gam.hr.area  ~ dnorm(0,0.01)
    sigma.phi.hr.area ~ dunif(0,10)
    sigma.gam.hr.area ~ dunif(0,10)
    tau.phi.hr.area <- 1/(sigma.phi.hr.area*sigma.phi.hr.area)
    tau.gam.hr.area <- 1/(sigma.gam.hr.area*sigma.gam.hr.area)

    mu.phi.nat.area  ~ dnorm(0,0.01)
    mu.gam.nat.area  ~ dnorm(0,0.01)
    sigma.phi.nat.area ~ dunif(0,10)
    sigma.gam.nat.area ~ dunif(0,10)
    tau.phi.nat.area <- 1/(sigma.phi.nat.area*sigma.phi.nat.area)
    tau.gam.nat.area <- 1/(sigma.gam.nat.area*sigma.gam.nat.area)

    mu.phi.fra  ~ dnorm(0,0.01)
    mu.gam.fra  ~ dnorm(0,0.01)
    sigma.phi.fra ~ dunif(0,10)
    sigma.gam.fra ~ dunif(0,10)
    tau.phi.fra <- 1/(sigma.phi.fra*sigma.phi.fra)
    tau.gam.fra <- 1/(sigma.gam.fra*sigma.gam.fra)

    mu.phi.traits1  ~ dnorm(0,0.01)
    mu.gam.traits1  ~ dnorm(0,0.01)
    sigma.phi.traits1 ~ dunif(0,10)
    sigma.gam.traits1 ~ dunif(0,10)
    tau.phi.traits1 <- 1/(sigma.phi.traits1*sigma.phi.traits1)
    tau.gam.traits1 <- 1/(sigma.gam.traits1*sigma.gam.traits1)

    mu.phi.traits2  ~ dnorm(0,0.01)
    mu.gam.traits2  ~ dnorm(0,0.01)
    sigma.phi.traits2 ~ dunif(0,10)
    sigma.gam.traits2 ~ dunif(0,10)
    tau.phi.traits2 <- 1/(sigma.phi.traits2*sigma.phi.traits2)
    tau.gam.traits2 <- 1/(sigma.gam.traits2*sigma.gam.traits2)

    mu.phi.traits1.fra  ~ dnorm(0,0.01)
    mu.gam.traits1.fra  ~ dnorm(0,0.01)
    sigma.phi.traits1.fra ~ dunif(0,10)
    sigma.gam.traits1.fra ~ dunif(0,10)
    tau.phi.traits1.fra <-
      1/(sigma.phi.traits1.fra*sigma.phi.traits1.fra)
    tau.gam.traits1.fra <-
      1/(sigma.gam.traits1.fra*sigma.gam.traits1.fra)


    ## species-specific  parameters
    for(sp in 1:nsp) {
      ## day
      p.0[sp]     ~ dnorm(mu.p.0,     tau.p.0)
      p.day.1[sp] ~ dnorm(mu.p.day.1, tau.p.day.1)
      p.day.2[sp] ~ dnorm(mu.p.day.2, tau.p.day.2)

      ## species specific
      phi.0[sp] ~ dnorm(mu.phi.0, tau.phi.0)
      gam.0[sp] ~ dnorm(mu.gam.0, tau.gam.0)

      ## hedgerow area
      phi.hr.area[sp] ~ dnorm(mu.phi.hr.area,
                              tau.phi.hr.area)
      gam.hr.area[sp] ~ dnorm(mu.gam.hr.area,
                              tau.gam.hr.area)

      ## natural habitat
      phi.nat.area[sp] ~ dnorm(mu.phi.nat.area,
                               tau.phi.nat.area)
      gam.nat.area[sp] ~ dnorm(mu.gam.nat.area,
                               tau.gam.nat.area)

      ## fra
      phi.fra[sp] ~ dnorm(mu.phi.fra,
                          tau.phi.fra)
      gam.fra[sp] ~ dnorm(mu.gam.fra,
                          tau.gam.fra)

      ## traits 1
      phi.traits1[sp] ~ dnorm(mu.phi.traits1,
                              tau.phi.traits1)
      gam.traits1[sp] ~ dnorm(mu.gam.traits1,
                              tau.gam.traits1)

      ## traits 2
      phi.traits2[sp] ~ dnorm(mu.phi.traits2,
                              tau.phi.traits2)
      gam.traits2[sp] ~ dnorm(mu.gam.traits2,
                              tau.gam.traits2)

      ## traits 1 * fra interaction
      phi.traits1.fra[sp] ~ dnorm(mu.phi.traits1.fra,
                                  tau.phi.traits1.fra)
      gam.traits1.fra[sp] ~ dnorm(mu.gam.traits1.fra,
                                  tau.gam.traits1.fra)
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
        ## start off at the average for each species, site across years
        logit(phi.site.mean[site,sp]) <- mean(phi[site,,sp])
        logit(gam.site.mean[site,sp]) <- mean(gam[site,,sp])
        
        psi.1[site,sp] <- gam.site.mean[site,sp]/
          (1 - phi.site.mean[site,sp] + gam.site.mean[site,sp])

        ## occupancy in year 1
        psi[site,1,sp] <- psi.1[site,sp]
        Z[site,1,sp] ~ dbern(psi.1[site,sp])
        
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
            phi.hr.area[sp]*HRarea[site] +
            phi.nat.area[sp]*natural[site] +
            phi.fra[sp]*fra[site, yr] +
            phi.traits1.fra[sp]*fra[site, yr]*traits1[sp]

          gam[site,yr,sp] <-
            gam.0[sp] +
            gam.traits1[sp]*traits1[sp] +
            gam.traits2[sp]*traits2[sp] +
            gam.hr.area[sp]*HRarea[site] +
            gam.nat.area[sp]*natural[site] +
            gam.fra[sp]*fra[site, yr] +
            gam.traits1.fra[sp]*fra[site, yr]*traits1[sp]

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
  }
  model.out <- jags(data=list(X=d$data$X,
                              day=d$data$day,
                              fra=d$data$fra,
                              HRarea=d$data$HRarea,
                              natural=d$data$natural,
                              traits1=d$data$traits1,
                              traits2=d$data$traits2,
                              nrep=d$constants$nrep,
                              nsp=d$constants$nsp,
                              nsite=d$constants$nsite,
                              nyear=d$constants$nyear),
                    inits=d$inits,
                    parameters.to.save=d$monitors,
                    model.file=model.jags,
                    n.chains=nc,
                    n.thin=nt,
                    n.iter=ni,
                    n.burnin=nb,
                    working.directory=NULL)
  return(model.out)
}

## specify the parameters to be monitored
get.params <- function(){
    c('mu.p.0',
      'mu.p.day.1',
      'mu.p.day.2',
      'sigma.p.0',
      'sigma.p.day.1',
      'sigma.p.day.2',
      'mu.phi.0',
      'sigma.phi.0',
      'mu.phi.rain',
      'sigma.phi.rain',
      'mu.phi.hr.area',
      'sigma.phi.hr.area',
      'mu.phi.nat.area',
      'sigma.phi.nat.area',
      'mu.phi.fra',
      'sigma.phi.fra',
      'mu.phi.traits1',
      'sigma.phi.traits1',
      'mu.phi.traits2',
      'sigma.phi.traits2',
      'mu.gam.0',
      'sigma.gam.0',
      'mu.gam.rain',
      'sigma.gam.rain',
      'mu.gam.hr.area',
      'sigma.gam.hr.area',
      'mu.gam.nat.area',
      'sigma.gam.nat.area',
      'mu.gam.fra',
      'sigma.gam.fra',
      'mu.gam.traits1',
      'sigma.gam.traits1',
      'mu.gam.traits2',
      'sigma.gam.traits2',
      'mu.gam.traits1.fra',
      'sigma.gam.traits1.fra',
      'phi.site.mean',
      'gam.site.mean',
      'phi.sp.mean',
      'gam.sp.mean')
}
