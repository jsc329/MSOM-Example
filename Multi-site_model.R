model {
  
  # Hyper-priors
  
  # Occupancy covariates
  
  for (i in 1:n.statecovs) {
    
    mu.b[i] ~ dnorm(0, 0.01)
    
    tau.b[i] ~ dgamma(0.1, 0.1)
    
  }
  
  # Detection covariates
  
  for (i in 1:n.detcovs) {
    
    mu.a[i] ~ dnorm(0, 0.01)
    
    tau.a[i] ~ dgamma(0.1, 0.1)
    
  }
  
  
 for (i in 1:(n.sites-1)){
   mu.site[i] ~ dnorm(0, 0.01)
   tau.site[i] ~ dgamma(0.1, 0.1)
 }
  

   


### Species loop
  
  
  
  for (i in 1:(n.species)) {    
    
    
    
    # Species-specific parameter priors
    
    for (n in 1:n.statecovs) {
      
      b[i,n] ~ dnorm(mu.b[n], tau.b[n])
      
    }
    
    for (o in 1:n.detcovs) {
      
      a[i,o] ~ dnorm(mu.a[o], tau.a[o])
      
    }
      # Set the first site to 0 so we can have constrasts
      siteeff[i, 1] <- 0
    
    for (r in 1:(n.sites-1)){
      siteeff[i, r+1] ~ dnorm(mu.site[r], tau.site[r]) 
    }
    
    
    for (j in 1:n.stations) {
      
      logit(psi[j,i]) <- b[i,1] + b[i,2]*hab[j] + siteeff[i, sitenum[j]]
      
      Z[j,i] ~ dbin(psi[j,i], 1)
      
      
      
      
      for (k in 1:n.visits){
        
        logit(theta[j,k,i]) <- a[i,1] + a[i,2]*wind[j]
        
        mu.p[j,k,i] <- theta[j,k,i]*Z[j,i]
        
        y[j,k,i] ~ dbin(mu.p[j,k,i], 1)
        
      }
      
      
    }
    
  }
  
}
