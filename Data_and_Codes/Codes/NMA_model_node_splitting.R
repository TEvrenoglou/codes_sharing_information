writeLines("
model{

for (i in 1:ns) {
    w[i, 1] <- 0
    delta[i, t[i, 1]] <- 0
    j[i,1] <-0 ### NEW
    u[i] ~ dnorm(0, .0001)
    for (k in 1:na[i]) {
      y[i, t[i, k]] ~ dnorm(phi[i, t[i, k]], prec[i, t[i,k]])
      
      phi[i, t[i, k]] <-(u[i] + delta[i, t[i, k]]) * pooled.sd[i]
      
      
      index[i,k] <- split[i] * (equals(t[i,k], pair[1]) + equals(t[i,k], pair[2])) ############# we need the split it is 0 when comparison exists in study and 1 else ####
      
    }
    for (k in 2:na[i]) {
      delta[i, t[i, k]] ~ dnorm(md[i, t[i, k]], taud[i,t[i, k]])
      
      md[i, t[i, k]] <-(d[t[i, k]] - d[t[i, 1]] + sw[i,k])*(1-index[i,m[i,k]])+((1-sign)*direct*index[i,m[i,k]]-sign*direct*index[i,m[i,k]]) 
      
      j[i,k] <- k - (equals(1, split[i]) * step(k-3))
      
      taud[i, t[i, k]] <- PREC*2*(j[i,k]-1)/j[i,k]   
      
      w[i, k] <- (delta[i, t[i, k]] - d[t[i, k]] + d[t[i,1]])*(1-index[i,k])
      
      sw[i, k] <- sum(w[i, 1:(k - 1)]) / (j[i,k] - 1)
    }
}
  
  d[ref] <- 0
  
  direct~dnorm(direct.mean,direct.prec)
  
    for (k in 1:(ref-1)) {
      #prec.d[k] <- 1/pow(prior.d[k,2],2)
      prec.d[k] <- prior.d[k,2]
      m.d[k] <- prior.d[k,1]
      d[k] ~ dnorm(m.d[k],prec.d[k]) # Asssign the informative priors for the basic comparisons
      SMD.ref[k] <- d[k] - d[ref]
      predSMD.ref[k] ~ dnorm(SMD.ref[k],PREC) # Create predictive distributions
    }
    for (k in (ref+1):nt) {
      #prec.d[k] <- 1/pow(prior.d[k-1,2],2)
      prec.d[k] <- prior.d[k-1,2]
      m.d[k] <- prior.d[k-1,1]
      d[k] ~ dnorm(m.d[k],prec.d[k]) # Asssign the informative priors for the basic comparisons
      SMD.ref[k] <- d[k] - d[ref]
      predSMD.ref[k] ~ dnorm(SMD.ref[k],PREC) # Create predictive distributions
    }
  
  tau ~ dnorm(0,1)T(0,)
  
  PREC <- 1 / pow(tau, 2)
  
  for (c in 1:(nt - 1)) {
    for (k in (c + 1):nt) {
      SMD[c, k] <- d[c] - d[k]
    }
  }
  
  diff <- direct - SMD.ref[pair[1]] ##################################
  
  # # calculate p-value
  prob <- step(diff)
  
  
}",con="CA_model_consistency.txt")


modfile_consistency="CA_model_consistency.txt"

