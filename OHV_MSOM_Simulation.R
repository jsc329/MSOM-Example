

points <- 300
species <- 30
nspec <- species
years <- 2
visits <- 3
nrep <- visits
mean.psi <- 0.7
mean.p <- 0.5

sig.lpsi <- 0.5
mu.beta1.lpsi <- -0.5
mu.beta2.lpsi <- 0.8
sig.beta.lpsi  <- 0.2
sig.lp <- 0.2
mu.beta.lp <- -0.5
sig.beta.lp <- 0.3
mu.lsiteOHV <- -0.5
mu.lsitenonOHV <- 0.8
sig.lsite <- 0.1
mu.lyear <- -0.5
sig.lyear <- 0.2


# Make empty dataframes to hold the simulation output
y.all <- array(dim=c(points, visits, years, species), data=0)
detected.at.all <- array(dim=c(years, species))
psi <- array(dim=c(points, years, species), data=0)
z <- array(dim=c(points, years, species), data=0)
p <- array(dim=c(points, visits, years, species), data=0)
detcovs <- array(dim=c(points, years), data=rnorm(points*years))
habcovs <- array(dim=c(points, years), data=rnorm(points*years))
ohveff <- array(dim=c(species, 2), data=0)
yeareff <- array(dim=c(species, years), data=0)

habitat <- habcovs
wind <- detcovs


ohvind <- c(rep(1,150), rep(2,150))

yearnum <- NULL
for (i in 1:years){
  yearnum <- c(yearnum, rep(i, points*years))
}

pointnum <- 1:points


# Draw species-specific intercepts and slopes from their normal distributions
# Build up linear predictors for occupancy and detection
mu.lpsi <- ifelse(mean.psi == '1', 500, qlogis(mean.psi))
mu.lp <- ifelse(mean.p == '1', 500, qlogis(mean.p))

beta0 <- rnorm(species, mu.lpsi, sig.lpsi)           # occupancy intercept
beta1OHV <- rnorm(species, mu.beta1.lpsi, sig.beta.lpsi) # occupancy slope on habitat
beta1nonOHV <- rnorm(species, mu.beta2.lpsi, sig.beta.lpsi) # occupancy slope on habitat
alpha0 <- rnorm(species, mu.lp, sig.lp)              # detection intercept
alpha1 <- rnorm(species, mu.beta.lp, sig.beta.lp)    # detection slope on wind

  ohveff[,1] <- rnorm(species, mu.lsiteOHV, sig.lsite)
  ohveff[,2] <- rnorm(species, mu.lsitenonOHV, sig.lsite)


yeareff[,1] <- 0
for (t in 2:years){
  yeareff[,t] <- rnorm(species, mu.lyear, sig.lyear)
}

# Add in an extra level for years
for (t in 1:years){
  for(k in 1:nspec){
    for (j in 1:points){
      if (j <= 150){
        psi[j,t,k] <- plogis(beta0[k] + beta1OHV[k] * habitat[j,t] + ohveff[k, ohvind[j]] + yeareff[k, t])
      }
      else {
        psi[j,t,k] <- plogis(beta0[k] + beta1nonOHV[k] * habitat[j,t] + ohveff[k, ohvind[j]] + yeareff[k, t])
      }
    }
    for(r in 1:nrep){
      p[,r,t,k] <- plogis(alpha0[k] + alpha1[k] * wind[,t])
    }
  }
}
# Distribute species over sites (simulate true state)
for(k in 1:nspec){
  for (t in 1:years){
    z[,t,k] <- rbinom(points, 1, psi[,t,k])
  }
}
occurring.in.sample <- apply(z, c(2,3), max) # Presence/absence at study sites

# Measurement of presence/absence (simulate observation)
for (t in 1:years) {
for(k in 1:nspec) {
  for(i in 1:points){
    for(j in 1:nrep) {
      y.all[i,j,t,k] <- rbinom(1, z[i,t,k], p[i,j,t,k])
    }
  }
  # detected.at.all[k] <- if(any(y.all[,,k]>0)) TRUE else FALSE
  detected.at.all[t,k] <- any(y.all[, , t, k] > 0)
}
}

detected.at.all <- apply(detected.at.all, 2, all)

y.obs <- y.all[,,,detected.at.all]    # Drop species never detected
detected.at.site <- apply(y.obs>0, c(1,3,4), any)
y.sum.all <- apply(y.all, c(1,3,4), sum) # Detection frequency for all species
y.sum.obs <- y.sum.all[,,detected.at.all]# Detection frequency for obs. species
z.obs <- apply(y.all, c(1,3,4), max)   # Observed presence/absence matrix
missed.sites <- z - z.obs            # Sites where species missed
Ntotal.fs <- sum(occurring.in.sample)# Number of species in finite-sample
Ntotal.obs <- sum(detected.at.all)   # Observed species richness (all sites)
S.true <- apply(z, 1, sum)           # Vector of true local richness
S.obs <- apply(z.obs, 1, sum)        # Vector of observed local richness


