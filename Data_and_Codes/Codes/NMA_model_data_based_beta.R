writeLines(
  "
  model{
    v1 ~ dbeta(3,3) # generate a value from Beta(3,3)
    v2 <- 1
for (i in 1:ns) {
  w[i, 1] <- 0
  delta[i, t[i, 1]] <- 0
  v[i] <- ind[i]*v1+(1-ind[i])*v2 # downweight the study i if ind=1, otherwise not
  u[i] ~ dnorm(0, 1e-04) # noninformative prior for the baseline group

  for (k in 1:na[i]) {
    y[i, t[i, k]] ~ dnorm(phi[i, t[i, k]],v[i]*prec[i, t[i,k]]) ### within-studies assumption

    phi[i, t[i, k]] <-(u[i] + delta[i, t[i, k]]) * pooled.sd[i] ### parameterized NMA model for continuous outcome
  }
  for (k in 2:na[i]) {
    delta[i, t[i, k]] ~ dnorm(md[i, t[i, k]]-beta[t[i, k]], taud[i,t[i, k]]) ### distribution of the random effects with the location parameter

    md[i, t[i, k]] <-d[t[i, k]] - d[t[i, 1]] + sw[i,k] # consistency equations for treatment effects with correction for multi-arm studies

    taud[i, t[i, k]] <-PREC * 2 * (k - 1) / k # correction for multi-arm studies

    w[i, k] <- (delta[i, t[i, k]] - d[t[i, k]] + d[t[i,1]])

    sw[i, k] <-sum(w[i, 1:(k - 1)]) / (k - 1)
  }
}
d[ref] <- 0
beta[ref] <- 0

for (k in 1:(ref - 1)) {
  d[k] ~ dnorm(0, 1e-04)
}
for (k in (ref + 1):nt) {
  d[k] ~ dnorm(0, 1e-04)
}
   for (k in 1:(ref-1)){
      prec.b[k] <- 1/pow(b[k,2],2)
      m.b[k] <- b[k,1]
      beta[k] ~ dnorm(m.b[k],prec.b[k])

   }

     for (k in (ref+1):nt){
      prec.b[k] <- 1/pow(b[k,2],2)
      m.b[k] <- b[k,1]
      beta[k] ~ dnorm(m.b[k],prec.b[k])

  }

tau ~ dnorm(0,1)T(0,) # weakly informative prior for heterogeneity parameter
PREC <- 1 / pow(tau, 2)

#### obtain SMD's for all treatment comparisons
for (c in 1:(nt - 1)) {
  for (k in (c + 1):nt) {
    SMD[c, k] <- d[c] - d[k]
  }
}

#### obtain SMD's only for the basic treatment comparison
for (c in 1:nt) {
  SMD.ref[c] <- d[c] - d[ref]
}

#### obtain predictive distributions for the basic treatment comparisons
for (c in 1:(ref - 1)) {
  X[c] <- d[c] - d[ref]
  predSMD.ref[c] ~ dnorm(X[c], PREC)
}
for (c in (ref + 1):nt) {
  X[c] <- d[c] - d[ref]
  predSMD.ref[c] ~ dnorm(X[c], PREC)
}

#### obtain predictive distributions for all treatment comparison
for (c in 1:(nt - 1)) {
  for (k in (c + 1):nt) {
    predSMD[c, k] ~ dnorm(SMD[c, k], PREC)
  }
}

#### obtain SUCRAS for all treatments in the network
order[1:nt] <- rank(d[1:nt])
for (k in 1:nt) {
  most.effective[k] <- equals(order[k], 1)
  for (j in 1:nt) {
    effectiveness[k, j] <- equals(order[k], j)
  }
}
for (k in 1:nt) {
  for (j in 1:nt) {
    cumeffectiveness[k, j] <- sum(effectiveness[k, 1:j])
  }
}
for (k in 1:nt) {
  SUCRA[k] <- sum(cumeffectiveness[k, 1:(nt - 1)]) / (nt -1)

}
}", con="GP.model.data.txt")

modfile_gp = 'GP.model.data.txt'

