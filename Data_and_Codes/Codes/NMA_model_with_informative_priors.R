# NMA model for the analysis of CA. The mains structure is the same as described in Dias et al. (DOI: 10.1177/0272989X12458724)

writeLines(
  "
  model{
    for(i in 1:ns) {
      u[i] ~ dnorm(0,.0001) # Noninformative prior for the baseline group
      w[i,1] <- 0
      delta[i,t[i,1]] <- 0
      for(k in 1:na[i]) {
        y[i,t[i,k]] ~ dnorm(theta[i,t[i,k]],prec[i,t[i,k]]) # Within-studies assumption
        theta[i,t[i,k]] <- (u[i]+delta[i,t[i,k]])*pooled.sd[i] # Parametrized NMA model
      }
      for (k in 2:na[i]) {
        delta[i,t[i,k]] ~ dnorm(md[i,t[i,k]],taud[i,t[i,k]]) # Across-studies assumption
        md[i,t[i,k]] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # Consistency equations with correction for multi-arms studies
        taud[i,t[i,k]] <- PREC * 2 * (k-1) / k
        w[i,k] <- (delta[i,t[i,k]] - d[t[i,k]] + d[t[i,1]])
        sw[i,k] <- sum(w[i,1:(k-1)]) / (k-1)
      }
    }
    d[ref] <- 0
    for (k in 1:(ref-1)) {
      prec.d[k] <- 1/pow(prior.d[k,2],2)
      m.d[k] <- prior.d[k,1]
      d[k] ~ dnorm(m.d[k],prec.d[k]) # Asssign the informative priors for the basic comparisons
      SMD.ref[k] <- d[k] - d[ref]
      predSMD.ref[k] ~ dnorm(SMD.ref[k],PREC) # Create predictive distributions
    }
    for (k in (ref+1):nt) {
      prec.d[k] <- 1/pow(prior.d[k-1,2],2)
      m.d[k] <- prior.d[k-1,1]
      d[k] ~ dnorm(m.d[k],prec.d[k]) # Asssign the informative priors for the basic comparisons
      SMD.ref[k] <- d[k] - d[ref]
      predSMD.ref[k] ~ dnorm(SMD.ref[k],PREC) # Create predictive distributions
    }
    tau ~ dnorm(0,1)T(0,) # weakly informative prior for heterogeneity parameter
    PREC <- 1 / pow(tau,2)
    for (c in 1:(nt-1)) {
      for (k in (c+1):nt) {
        SMD[c,k] <- d[c] - d[k] # SMD's for all pairwise comparisons in the network
        pred.SMD[c,k] ~ dnorm(SMD[c,k],PREC)
      }
    }
    # Calculate SUCRAS for ranking
    order[1:nt] <- rank(d[1:nt])
    for (k in 1:nt) {
      most.effective[k] <- equals(order[k],1)
      for (j in 1:nt) {
        effectiveness[k,j] <- equals(order[k],j)
      }
    }
    for (k in 1:nt) {
      for (j in 1:nt) {
        cum.effectiveness[k,j] <- sum(effectiveness[k,1:j])
      }
    }
    for (k in 1:nt) {
      SUCRA[k] <- sum(cum.effectiveness[k,1:(nt-1)]) / (nt-1)
    }
  }
           ",
  con = "CA.model.txt"
)

modfile_ca = 'CA.model.txt'
