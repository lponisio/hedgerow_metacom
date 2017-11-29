ms.ms <- function(d,
                  ni=1100,
                  nt=100,
                  nb=100,
                  nc=3) {

  sink('model.jags')
  cat('model {

    ## multi-species priors
    mu.p.0     ~ dnorm(0,0.001)
    mu.p.day.1 ~ dnorm(0,0.001)
    mu.p.day.2 ~ dnorm(0,0.001)

    sigma.p.0     ~ dunif(0,10)
    sigma.p.day.1 ~ dunif(0,10)
    sigma.p.day.2 ~ dunif(0,10)
    tau.p.0     <- pow(sigma.p.0,    -2)
    tau.p.day.1 <- pow(sigma.p.day.1,-2)
    tau.p.day.2 <- pow(sigma.p.day.2,-2)

    mu.phi.0  ~ dnorm(0,0.001)
    mu.gam.0  ~ dnorm(0,0.001)

    sigma.phi.0 ~ dunif(0,10)
    sigma.gam.0 ~ dunif(0,10)
    tau.phi.0 <- pow(sigma.phi.0,-2)
    tau.gam.0 <- pow(sigma.gam.0,-2)

    phi.dprime     ~ dnorm(0,0.001)
    gam.dprime     ~ dnorm(0,0.001)
    phi.ypr        ~ dnorm(0,0.001)
    gam.ypr        ~ dnorm(0,0.001)
    phi.ypr.dprime ~ dnorm(0,0.001)
    gam.ypr.dprime ~ dnorm(0,0.001)

    ## species-specific detectability
    for(sp in 1:nsp) {
      p.0[sp]     ~ dnorm(mu.p.0,    tau.p.0)
      p.day.1[sp] ~ dnorm(mu.p.day.1,tau.p.day.1)
      p.day.2[sp] ~ dnorm(mu.p.day.2,tau.p.day.2)
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

      phi.0[sp] ~ dnorm(mu.phi.0,tau.phi.0)
      gam.0[sp] ~ dnorm(mu.gam.0,tau.gam.0)

      for(site in 1:nsite) {

        ## occupancy in year 1
        psi[site,1,sp] <- frac.presence[site,sp]

        Z[site,1,sp] ~ dbern(psi[site,1,sp])

        ## detectability in year 1
        for(rep in 1:nrep[site,1,sp]) {
          mu.p[site,1,rep,sp] <- Z[site,1,sp]*p[site,1,rep,sp]
          X[site,1,rep,sp] ~ dbern(mu.p[site,1,rep,sp])
        }

        ## occupancy and detectability in subsequent years
        for(yr in 1:(nyear-1)) {
          phi[site,yr,sp] <-
            phi.0[sp] +
              phi.dprime*dprime[sp] +
                phi.ypr*ypr[site,yr] +
                  phi.ypr.dprime*ypr[site,yr]*dprime[sp]

          gam[site,yr,sp] <-
            gam.0[sp] +
              gam.dprime*dprime[sp] +
                gam.ypr*ypr[site,yr] +
                  gam.ypr.dprime*ypr[site,yr]*dprime[sp]

          logit(psi[site,yr+1,sp]) <-
            Z[site,yr,sp] * phi[site,yr,sp] +
              (1-Z[site,yr,sp]) * gam[site,yr,sp]

          Z[site,yr+1,sp] ~ dbern(psi[site,yr+1,sp])

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
        N[site,yr] <- sum(Z[site,yr,])
      }
    }
  }',fill = TRUE)
  sink()

  attach(d$data)
  model.out <- jags.parallel(data=list('X', 'day', 'nrep', 'nsp',
                               'nsite', 'nyear', 'ypr', 'dprime',
                               'site.presence', 'frac.presence'),
                             inits=d$inits,
                             parameters.to.save=d$params,
                             model.file='model.jags',
                             n.chains=nc,
                             n.thin=nt,
                             n.iter=ni,
                             n.burnin=nb,
                             working.directory=NULL)
  detach(d$data)
  model.out
}

## specify the parameters to be monitored
get.params <- function()
  c('mu.p.0',
    'mu.p.day.1',
    'mu.p.day.2',
    'sigma.p.0',
    'sigma.p.day.1',
    'sigma.p.day.2',
    'p.0',
    'p.day.1',
    'p.day.2',
    'mu.phi.0',
    'sigma.phi.0',
    'phi.0',
    'phi.ypr',
    'phi.dprime',
    'phi.ypr.dprime',
    'phi',
    'mu.gam.0',
    'sigma.gam.0',
    'gam.0',
    'gam.ypr',
    'gam.dprime',
    'gam.ypr.dprime',
    'gam',
    'psi',
    'N')
