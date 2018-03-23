

model {
  
  # Modified version of Tingley et al. 2016 Pyrodiveristy paper
  
  # Hyper-priors
  # Occupancy covariates
  for (l in 1:2) {
    mu.b[l] ~ dnorm(0, 0.01)
    tau.b[l] ~ dgamma(0.1, 0.1)
  }
  # Detection covariates
  for (m in 1:3) {
    mu.a[m] ~ dnorm(0, 0.01)
    tau.a[m] ~ dgamma(0.1, 0.1)
  }
  mu.phi ~ dnorm(0, 0.01)
  tau.phi ~ dgamma(0.1, 0.1)
  
  ### Species loop
  
  for (i in 1:(n.species)) {    
    
    # Species-specific parameter priors
    for (n in 1:2) {
      b[i,n] ~ dnorm(mu.b[n], tau.b[n])
    }
    for (o in 1:3) {
      a[i,o] ~ dnorm(mu.a[o], tau.a[o])
    }
    #phi[i] ~ dnorm(mu.phi, tau.phi)
    
    # State model for time = 1
    
    for (j in 1:n.sites) {
      logit(psi[j,i]) <- b[i,1] + b[i,2]*hab[j]
      Z[j,i] ~ dbin(psi[j,i], 1)
      
      # Data model for time = 1
      
      for (k in 1:n.visits){
        logit(theta[j,k,i]) <- a[i,1] + a[i,2]*wind[j,k] + a[i,3]*effort[k]
        mu.p[j,k,i] <- theta[j,k,i]*Z[j,i]
        y[j,k,i] ~ dbin(mu.p[j,k,i], 1)
      }
      
    }
  }
}