# NMA model for the analysis of GP. The mains structure is the same as described in Dias et al. (DOI: 10.1177/0272989X12458724)
writeLines(
  "
  model{
    v1 ~ dbeta(3,3) # generate a value from Beta(3,3)
    v2 <- 1
    for(i in 1:ns) {
      u[i] ~ dnorm(0,.0001)
      v[i] <- ind[i]*v1+(1-ind[i])*v2 # downweight the study i if ind=1, otherwise not
      w[i,1] <- 0
      delta[i,t[i,1]] <- 0
      for(k in 1:na[i]) {
        y[i,t[i,k]] ~ dnorm(theta[i,t[i,k]],v[i]*prec[i,t[i,k]]) ### within-studies assumption with inflated variance
        theta[i,t[i,k]] <- (u[i]+delta[i,t[i,k]])*pooled.sd[i] ### parametrized NMA model for continuous outcomes
      }
      for (k in 2:na[i]) {
        delta[i,t[i,k]] ~ dnorm(md[i,t[i,k]],taud[i,t[i,k]]) ### distribution of the random effects
        md[i,t[i,k]] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # consistency equations for treatment effects with correction for multi-arm studies
        taud[i,t[i,k]] <- PREC * 2 * (k-1) / k
        w[i,k] <- delta[i,t[i,k]] - d[t[i,k]] + d[t[i,1]] # correction for multi-arm studies
        sw[i,k] <- sum(w[i,1:(k-1)]) / (k-1)
      }
    }
    d[ref] <- 0
    for (k in 1:(ref-1)) {
      d[k] ~ dnorm(0,.0001)
      SMD.ref[k] <- (d[k] - d[ref]) - (beta[k]-beta[ref]) # shift the SMD's based on the location parameters
      predSMD.ref[k] ~ dnorm(SMD.ref[k],PREC)
    }
    for (k in (ref+1):nt) {
      d[k] ~ dnorm(0,.0001)
      SMD.ref[k] <- (d[k] - d[ref]) - (beta[k]-beta[ref]) # shift the SMD's based on the location parameters
      predSMD.ref[k] ~ dnorm(SMD.ref[k],PREC)
    }
    tau ~ dnorm(0,1)T(0,) # weakly informative prior for heterogeneity parameter
    PREC <- 1 / pow(tau,2)
    for (k in 1:nt){
      prec.b[k] <- 1/pow(b[k,2],2)
      m.b[k] <- b[k,1]
      beta[k] ~ dnorm(m.b[k],prec.b[k]) # prior for the location parameter as constructed for each drug
    }
  }
           ",
  con = "GP.model.txt"
)

modfile_GP_expert_downweight = 'GP.model.txt'