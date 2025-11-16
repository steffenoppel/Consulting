model {
  # -------------------------------------------------
  # Vital parameters:
  # s: survival probability
  # p: recapture probability
  # r: recovery probability
  # -------------------------------------------------
  # States: 
  # 1 alive 
  # 2 recently dead 
  # 3 recently dead but not recovered, or long dead (absorbing)
  # -------------------------------------------------
  # Observations:
  # 1 seen alive
  # 2 recovered dead
  # 3 neither seen nor recovered
  # -------------------------------------------------
  # Data:
  # n_occ = number of occasions; integer
  # n_ind = number of individuals; integer
  # CH = capture history matrix; dim(n_ind, n_occ)
  # AM = age matrix; dim(n_ind, n_occ)
  # TM = tag type matrix; dim(n_ind, n_occ)
  # n_tag = number of different tags in dataset (3 = ring or VHF or GPS, 2 = ring or GPS/VHF)
  # sex = vector of sex of individuals (1 = female, 2 = male, NA)
  # sex.ratio = sex ratio computed form available sexes, to be able to impute NAs in sex; vector of 2
  # first = vector of occasion of first encounter; length n_ind
  # FR = ones = vector with the value 1 repeated n_ind times
  # -------------------------------------------------
  #
  # Priors
  for(a in 1:2) { # age
    s_age[a] ~ dunif(0,1)
  }
  for(tt in 1:n_tag) { # tag type
    p_tag[tt] ~ dunif(0,1)
    r_tag[tt] ~ dunif(0,1)
  }
  #
  # Linear predictors - Constraints
  for(i in 1:n_ind){
    for(t in first[i]:(n_occ-1)){
      s[i,t] <- s_age[AM[i,t]]
      p[i,t] <- p_tag[TM[i,t]]
      r[i,t] <- r_tag[TM[i,t]]
    }
  }
  #
  # Define state-transition and observation matrices
  for (i in 1:n_ind){
    for (t in first[i]:(n_occ-1)){
      # State transition probabilities
      ps[1,i,t,1] <- s[i,t]
      ps[1,i,t,2] <- (1-s[i,t])
      ps[1,i,t,3] <- 0
      ps[2,i,t,1] <- 0
      ps[2,i,t,2] <- 0
      ps[2,i,t,3] <- 1
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 1
      # Detection probabilities
      po[1,i,t,1] <- p[i,t]
      po[1,i,t,2] <- 0
      po[1,i,t,3] <- 1-p[i,t]
      po[2,i,t,1] <- 0
      po[2,i,t,2] <- r[i,t]
      po[2,i,t,3] <- 1-r[i,t]
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- 0
      po[3,i,t,3] <- 1
    }
  }
  #
  # Likelihood - Marginalization
  for(i in 1:n_ind){
    zeta[i,first[i],1] <- 1
    zeta[i,first[i],2] <- 0
    zeta[i,first[i],3] <- 0
    for(t in first[i]:(n_occ-1)){
      for(k in 1:3) { # number of states
        zeta[i,(t+1),k] <- inprod(zeta[i,t,], ps[,i,t,k]) * po[k,i,t,CH[i,(t+1)]]
      }
    }
    lik[i] <- sum(zeta[i,n_occ,])
    ones[i] ~ dbin(lik[i], FR[i])
  }
}