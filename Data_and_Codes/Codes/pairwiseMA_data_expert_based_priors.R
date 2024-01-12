### pairwise meta-analysis model used for the data based beta priors
writeLines(
  "
           model {
              for (i in 1:N1) {
                 
                 w1[i] = 1/v1[i]
                 
                 y1[i] ~ dnorm(theta1[i], w1[i])
                 
                 theta1[i] ~ dnorm(mean1, prec1)

              }
                 ## model
                 for (i in 1:N2) {
                 
                 w2[i] <- 1/v2[i]
                 
                 y2[i] ~ dnorm(theta2[i], w2[i])
                 
                 theta2[i]=mean2

                 }
              beta = mean1-mean2 ## define the betas
              
              mean1 ~ dnorm(0, .0001)
              
              tau1 ~ dnorm(0,1)T(0,) # weakly informative prior for heterogeneity
              
              prec1 <- 1/pow(tau1,2)

              mean2 ~ dnorm(0, .0001)

           }
           ",
  con = "model_random"
)

modfile_meta <- 'model_random'

### parameters to save
save_ma_data <- c("beta")




### pairwise meta-analysis model used for the expert opinion based beta priors
writeLines(
  "
  model{
  ## Pool expert's opinion for each drug separately
    for(i in 1:nCA) {
    
      for(k in 1:naCA[i]) {
        
        yCA[i,tCA[i,k]] ~ dnorm(thetaCA[i,tCA[i,k]],gamma[i]*precCA[i,tCA[i,k]]) ### sd's provides by each expert are inflated with their confidence gamma
        
        thetaCA[i,tCA[i,k]] = uCA[i,tCA[i,k]]*pooled.sdCA[i] # standardize the true underlying effect
        
        uCA[i,tCA[i,k]] ~ dnorm(ksiCA[tCA[i,k]],prec.t)
      }
    }
    
    tau ~ dnorm(0,1)T(0,)
    
    prec.t <- 1/pow(tau,2)

    for(k in 1:nt){
      
      ksiCA[k] ~ dnorm(0,.0001)

      ksiGP[k] = yGP[k]/pooled.sdGP[k] # standardize the pooled result from GP

      beta[k] = ksiGP[k]-ksiCA[k] # calculate drug-speciic shift
    }
  }
  ",
  con = "betas.model.txt"
)

modfile_expert = 'betas.model.txt'

save_ma_expert <- c("beta")



