
settings.consistency <- function(data,priors,treat1,treat2,beta.setting,downweight.setting,expand=F){
  
  trts <- levels(as.factor(data$Drug))  
  
  pair1 <- which(trts==treat1)
  
  pair2 <- which(trts==treat2)
  
  if(expand==T){
    
    comparison_expand <- cbind.data.frame(trts[-pair2],trts[pair2])
    
    names(comparison_expand) <- c("treat1","treat2")
    
    comparisons_expand <- paste(comparison_expand$treat1,comparison_expand$treat2,sep = " vs ")
    
    priors_expand <- matrix(ncol = 5,nrow=length(comparisons_expand))
    
    priors_expand <- as.data.frame(priors_expand)
    
    priors_expand[,1] <- comparisons_expand
    
    priors_expand[,2] <- 0
    
    priors_expand[,3] <- 1/1000
    
    priors_expand[,4] <- "none"
    
    priors_expand[,5] <- "none"
    
    names(priors_expand) <- names(priors)
    
  }
  
  sign <- NA
  
  if(pair1<pair2){
    sign=1
  }else{
    sign=0
  }
  
  pair <- c(pair1,pair2)
  
  comp <- paste(treat1,"vs",treat2,sep=" ")  
  
  if(expand==F){
    prior <- priors %>% 
      filter(beta==beta.setting) %>% 
      filter(downweight==downweight.setting)
  }else{
    prior <- priors_expand %>% 
      filter(beta==beta.setting) %>% 
      filter(downweight==downweight.setting)  
  }
  
  
  E <- which(prior$comparison==comp)
  
  direct_mean <- prior$mean[E]
  
  direct_prec <- prior$prec[E]
  
  prior$mean[E] <- 0
  
  prior$prec[E] <- 0.0001
  
  prior <- as.data.frame(prior)
  
  res <- list("comparison" = comp,
              "priors" = prior,
              "direct_mean" = direct_mean,
              "direct_prec" = direct_prec,
              "pair" = pair,
              "sign" = sign
  )
  
  return(res)
}

run_consistency_check <- function(data,
                                  data.priors,
                                  beta,
                                  downweight,
                                  ref,
                                  treat1,
                                  treat2,
                                  posterior.samples=T,
                                  seed,expand=F){
  
  set.seed(seed)
  
  runCA.data <- list()
  
  t <- list()
  
  pair1 <- list()
  
  settings <- list()
  
  prior.d <- list()
  
  pair <- list()
  
  checkPair <- list()
  
  split <- list()
  
  na <- list()
  
  bi <- list()
  
  m <- list()
  
  si <- list()
  
  final.data <- list()
  
  node_split <- list()
  
  res <- list()
  
  direct_comparison <- list()
  
  difference <- list()
  
  P <- list()
  
  bayesian_pval <- list()
  
  basic_comparisons <- list()
  
  indirect_comparison <- list()
  
  dir <- list()
  
  L_dir <- list()
  
  U_dir <- list()
  
  median_dir <- list()
  
  ind <- list()
  
  L_ind <- list()
  
  U_ind <- list()
  
  median_ind <- list()
  
  diff <- list()
  
  L_diff <- list()
  
  U_diff <- list()
  
  median_diff <- list()
  
  pval <- list()
  
  data_results <- list()
  
  data_paper <- list()
  
  interval_direct <- list()
  
  direct_paper <- list()
  
  interval_indirect <- list()
  
  indirect_paper <- list()
  
  interval_diff <- list()
  
  difference_paper <- list()
  
  bugs_out <- list()
  
  post_direct <- list()
  
  basic_bugs_out <- list()
  
  post_indirect <- list()
  
  evidence <- list()
  
  density_data <- list()
  
  name1 <- list()
  
  name2 <- list()
  
  name <- list()
  
  namef <- list()
  
  post_name <- list()
  
  path <- list()
  
  for(i in 1:length(treat1)){
    
    runCA.data[[i]]=long2jags(Study_No,
                              Drug,
                              mean=OverallEfficM,
                              sd=OverallEfficSD,
                              n=OverallEfficN,
                              data=data)
    
    t[[i]]=runCA.data[[i]]$T
    
    settings[[i]] <- settings.consistency(data = data,
                                          priors = data.priors,
                                          treat1 = treat1[i],
                                          treat2 = treat2[i],
                                          beta.setting = beta,
                                          downweight.setting = downweight,
                                          expand=expand)
    
    prior.d[[i]] <- settings[[i]]$priors[,c(2,3)]
    
    pair[[i]] <- settings[[i]]$pair 
    
    checkPair[[i]] <- PairXY(t[[i]], pair[[i]])
    
    split[[i]]=checkPair[[i]][,"split"]
    
    na[[i]]  <- runCA.data[[i]]$n.k
    
    # Build vector bi[i] with baseline treatment: t[i, b[i]]
    bi[[i]] <- Basetreat(t[[i]], checkPair[[i]][,"b"])
    
    # Indexes to sweep non-baseline arms only
    m[[i]] <- NonbaseSweep(checkPair[[i]], na[[i]])
    
    # Build matrix si[i,k] with non-baseline treatments: t[i, m[i,k]]
    si[[i]] <- Sweeptreat(t[[i]],m[[i]])
    
    runCA.data[[i]][[10]] = split[[i]]
    
    runCA.data[[i]][[11]] = pair[[i]]
    
    runCA.data[[i]][[12]] = m[[i]]
    
    names(runCA.data[[i]])[10] = "split"
    
    names(runCA.data[[i]])[11] = "pair"
    
    names(runCA.data[[i]])[12] = "m"
    
    final.data[[i]] = list(ns = runCA.data[[i]]$k,
                           nt = runCA.data[[i]]$n,
                           na = runCA.data[[i]]$n.k,
                           t = runCA.data[[i]]$T,
                           y = runCA.data[[i]]$Y,
                           prec = runCA.data[[i]]$Pr,
                           pooled.sd = runCA.data[[i]]$pooled.sd,
                           ref = ref,
                           split = split[[i]],
                           pair = pair[[i]],
                           m = m[[i]],
                           prior.d = prior.d[[i]],
                           direct.mean = settings[[i]]$direct_mean,
                           direct.prec = settings[[i]]$direct_prec,
                           sign = settings[[i]]$sign
    )
    
    
    node_split[[i]] <- jags(data = final.data[[i]], 
                            inits = NULL,
                            parameters.to.save = c("SMD.ref","direct","diff","prob"), 
                            n.chains = 2, 
                            n.iter = 50000,
                            n.burnin = 10000,
                            DIC=F,
                            model.file = modfile_consistency)
    
    #### Gather the results
    
    res[[i]]=as.data.frame(node_split[[i]]$BUGSoutput$summary)
    
    res[[i]]$ind=as.character(rownames(res[[i]]))
    
    direct_comparison[[i]]=res[[i]] %>%
      filter(ind=="direct")
    
    difference[[i]]=res[[i]] %>%
      filter(ind=="diff")
    
    P[[i]]=res[[i]] %>%
      filter(ind=="prob")
    
    P[[i]]=P[[i]] %>%
      dplyr::select(mean)
    
    bayesian_pval[[i]]=2*min(P[[i]],1-P[[i]])
    
    #### Obtain indirect comparison#####
    basic_comparisons[[i]] = res[[i]] %>%
      filter(grepl("SMD.ref",ind))
    
    pair1[[i]] <- pair[[i]][which(settings[[i]]$pair!=ref)]
    
    indirect_comparison[[i]] = basic_comparisons[[i]] %>%
      filter(row.names(basic_comparisons[[i]]) == paste("SMD.ref[",pair1[[i]],"]",sep = ""))
    
    
    ######### Save results in a dataframe ####
    
    dir[[i]] = direct_comparison[[i]]$mean
    
    L_dir[[i]] = direct_comparison[[i]]$`2.5%`
    
    U_dir[[i]] = direct_comparison[[i]]$`97.5%`
    
    median_dir[[i]] = direct_comparison[[i]]$`50%`
    
    ind[[i]] = indirect_comparison[[i]]$mean
    
    L_ind[[i]] = indirect_comparison[[i]]$`2.5%`
    
    U_ind[[i]] = indirect_comparison[[i]]$`97.5%`
    
    median_ind[[i]] = indirect_comparison[[i]]$`50%`
    
    diff[[i]] = difference[[i]]$mean
    
    L_diff[[i]] = difference[[i]]$`2.5%`
    
    U_diff[[i]] = difference[[i]]$`97.5%`
    
    median_diff[[i]] = difference[[i]]$`50%`
    
    P[[i]] = P[[i]]$mean
    
    pval[[i]] = bayesian_pval[[i]]
    
    data_results[[i]] = cbind.data.frame(dir[[i]],
                                         L_dir[[i]],
                                         U_dir[[i]],
                                         median_dir[[i]],
                                         ind[[i]],
                                         L_ind[[i]],
                                         U_ind[[i]],
                                         median_ind[[i]],
                                         diff[[i]],
                                         L_diff[[i]],
                                         U_diff[[i]],
                                         median_diff[[i]],
                                         P[[i]],
                                         pval[[i]])
    
    names(data_results[[i]])=c("direct",
                               "lb.direct",
                               "ub.direct",
                               "median_direct",
                               "indirect",
                               "lb.indirect",
                               "ub.indirect",
                               "median_indirect",
                               "difference",
                               "lb.difference",
                               "ub.difference",
                               "median_difference",
                               "P(diff>0)",
                               "pval")
    
    ## create  properly row.names for data_results
    row.names(data_results[[i]])=paste(settings[[i]]$comparison)
    
    ##### table for paper ####
    data_paper[[i]]=format(round(data_results[[i]],digits = 2),nsmall=2)
    
    interval_direct[[i]] = paste("[",data_paper[[i]]$lb.direct,", ",data_paper[[i]]$ub.direct,"]",sep="")
    
    direct_paper[[i]] = paste(data_paper[[i]]$direct,interval_direct[[i]],sep=" ")
    
    interval_indirect[[i]] = paste("[",data_paper[[i]]$lb.indirect,", ",data_paper[[i]]$ub.indirect,"]",sep="")
    
    indirect_paper[[i]] = paste(data_paper[[i]]$indirect,interval_indirect[[i]],sep=" ")
    
    interval_diff[[i]] = paste("[",data_paper[[i]]$lb.difference,", ",data_paper[[i]]$ub.difference,"]",sep="")
    
    difference_paper[[i]] = paste(data_paper[[i]]$difference,interval_diff[[i]],sep=" ")
    
    data_paper[[i]]=cbind.data.frame(direct_paper[[i]],
                                     indirect_paper[[i]],
                                     difference_paper[[i]],
                                     data_paper[[i]]$`P(diff>0)`,
                                     data_paper[[i]]$pval)
    
    names(data_paper[[i]])=c("direct [95% CI]",
                             "indirect [95% CI]",
                             "difference [95% CI]",
                             "P(difference>0)",
                             "pval")
    
    row.names(data_paper[[i]])=row.names(data_results[[i]])
    
    
    if(posterior.samples==T){
      
      ###### Extract posterior samples from direct and indirect evidence #######################
      bugs_out[[i]] = as.data.frame(node_split[[i]]$BUGSoutput$sims.matrix)
      
      post_direct[[i]] = bugs_out[[i]]$direct
      
      if(expand==F){
        
        basic_bugs_out[[i]] = bugs_out[[i]][,1:14]
        
      }else{
        
        basic_bugs_out[[i]] = bugs_out[[i]][,1:33]
        
      }
      
      post_indirect[[i]] = basic_bugs_out[[i]][,which(names(basic_bugs_out[[i]]) == paste("SMD.ref[",pair1[[i]],"]",sep = ""))] 
      
      evidence[[i]] = c(rep("Direct evidence",length(post_direct[[i]])),rep("Indirect evidence",length(post_indirect[[i]])))
      
      density_data[[i]] = c(post_direct[[i]],post_indirect[[i]])
      
      density_data[[i]] = cbind.data.frame(density_data[[i]],evidence[[i]])
      
      names(density_data[[i]])=c("values","evidence")
      
      #### REMOVE ####
      if(beta=="data"){
        name1[[i]] <- "Data based beta"
      }
      else if(beta=="expert"){
        name1[[i]] <- "Expert opinion beta"
      }
      else if(beta=="none"){
        name1[[i]] <- "standard"
      }
      
      if(downweight=="no"){
        name2[[i]] <- "No DW"
      }
      else if(downweight=="rob"){
        name2[[i]] <- "RoB DW"
      }
      else if(downweight=="nct"){
        name2[[i]] <- "NCT DW"
      }
      else if(downweight=="none"){
        name2[[i]] <- "NMA"
      }
      
      name[[i]] <- paste(name1[[i]],name2[[i]],sep = "-")
      
      namef[[i]] <- paste(row.names(data_paper[[i]]),name[[i]],sep="_")
      
      post_name[[i]] = paste("Posterior_samples_",namef[[i]],".csv",sep="")
      
      density_data[[i]]$comparison <- row.names(data_paper[[i]])
      
      density_data[[i]]$model <- name[[i]]
      
    }
    
  }
  
  numerical_results <- bind_rows(data_paper, .id = "column_label")
  
  if(posterior.samples==T){
    final_results <- list("results"=numerical_results,
                          "posterior_samples"=density_data)
    
    final_results$results$column_label = NULL
    
  }else{
    final_results <- list("results"=numerical_results)
    
    final_results$results$column_label = NULL
  }
  
  
  return(final_results)
  
}


PairXY <- function(treat, pair)
  # Check if pair(X,Y) in row i of data 
  # and give baseline for data row i
{
  N <- nrow(treat)
  out <- cbind(split=rep(0,N), b=rep(0,N))
  for (i in 1:N) {
    # returns positions of matches to elements of pair in t[i,]
    # or zero if not present
    pos <- match(pair, treat[i,], nomatch=0)   # lenght = length(pair) = 2
    out[i,1] <- ifelse(prod(pos)>0, 1, 0)      # 1 if pair in line i, 0 o.w.
    out[i,2] <- ifelse(prod(pos)==0, 1, pos[1])
  }
  out
}


NonbaseSweep <- function(index, na)
  # gives na-1 indexes to sweep non-baseline arms only
{
  N <- NROW(na)
  C <- max(na)
  out <- matrix(nrow=N, ncol=C)
  for (i in 1:N) {
    for (k in 2:na[i]) {
      out[i,k] <- k - (index[i,"b"] >= k)
    }
  }
  out
}


Sweeptreat <- function(treat, m)
  # Builds matrix with non-baseline treatments
{
  N <- NROW(treat)
  C <- NCOL(m)
  out <- matrix(nrow=N, ncol=C)
  for (i in 1:N) {
    for (k in 2:C) {
      out[i,k] <- treat[i,m[i,k]]
    }
  }
  out
}


Basetreat <- function(treat, b)
  # Builds vector with baseline treatments
{
  N <- nrow(treat)
  out <- rep(0,N)
  for (i in 1:N) {
    out[i] <- treat[i,b[i]]
  }
  out
}


# Code obtained from Law et al. (https://doi.org/10.1186/s12874-019-0689-9)

cd <- function(ests, Cov, method="community", weighted=FALSE, quants=c(0.2, 0.4, 0.6, 0.8),
               fix.layout=if(method=="community")FALSE else TRUE, fast=FALSE, n=1000, seed=1, save=FALSE){
  
  require(mvtnorm)
  require(igraph)
  require(gplots)  
  
  #
  # ests: Vector of treatment effect estimates
  #
  # Cov: covariance matrix
  #
  # method: "community": Undertake community detection to group similar treatments WRT relative treatment effect 
  #                      estimates, SEs and p values, then plot. 
  #         "boostrap": Use bootstrap replication to obtain n sets of relative trt effect estimates and covariance matrices.
  #                     Undertake community detection on each of the n replications, and visualise the proportion of times
  #                     each pair of treatments is in the same community.
  #
  # weighted:  FALSE: Detect communities using unweighted adjacency matrices (i.e. thresholds only)
  #             TRUE: Detect communities using weighted adjacency matrices (i.e. combining thresholds and weights)
  #
  # quants: A scalar or vector of the quantiles to be used as thresholds in the creation of the adjacency matrices that form
  #         the basis of both methods. Default is quartiles (20, 40, 60, 80 percentiles).
  #
  #
  # fix.layout: For "community" method, fixes location of each node when plotting when TRUE. For "bootstrap" method, 
  #             fix.layout=TRUE returns heatmaps with treatments ordered alphabetically. When fix.layout=FALSE, the
  #             first heatmap returned has treatments ordered by R, and subsequent heatmaps with the same ordering.
  # 
  # fast: TRUE: Speeds up code by using the cluster edge betweenness algorithm to find the best set of communities. 
  #       FALSE: Searches over all possible sets of communities.
  #
  # n: Number of bootstrap replications to create ("bootstrap" method only)
  #
  # seed: Sets seed, for reproduceable plots (and results when using "bootstrap" method)
  #
  # save: If TRUE, saves plots as pdfs.
  #
  
  set.seed(seed)
  
  ############## FUNCTION TO CREATE SQUARE MATRICES OF EFFECT ESTIMATES, SEs AND P VALUES #####################
  
  fns <- function(est, covar=NULL)  {
    
    
    # Create estimates matrix:
    est <- c(0, est)
    mat <- matrix(rep(est, each=length(est)), nrow=length(est)) - matrix(rep(est, times=length(est)), nrow=length(est))
    
    if(is.null(covar)) {
      return(mat)
    }else{
      # Create SEs matrix
      trts <- nrow(covar)+1
      se_mat <- matrix(0, nrow=nrow(covar)+1, ncol=nrow(covar)+1)
      
      # First row and column is simply square root of the diagonal of the covariance matrix:
      se_mat[1, ] <- c(0, diag(covar))^0.5
      se_mat[, 1] <- c(0, diag(covar))^0.5
      
      # For remaining elements, use variance-covariance of linear combinations: var(AY)= A * Y * t(A)
      for(i in 1:nrow(covar)) {
        
        for(j in 1:ncol(covar)) {
          lin_com <- rep(0, nrow(covar))
          lin_com[i] <- lin_com[i] - 1
          lin_com[j] <- lin_com[j] + 1
          se_mat[i+1,j+1] <- (t(lin_com) %*% covar %*% lin_com)^0.5
        }
      }
      
      # Create p value matrix
      z <- mat/se_mat
      p <- 2 * pnorm(-abs(z))
      diag(p) <- 1 
      
      return(list(est=abs(mat), se=se_mat, p=1-p)) # We take absolute treatment effect estimates, and compliment of p value.
    }
  }
  
  
  ########################## END OF INITIAL FUNCTION #############################
  
  ############################### MAIN FUNCTION ##################################
  
  # Create required square matrices using above function:
  matrices <- fns(est=ests, covar=Cov)
  est <- matrices$est
  se <- matrices$se
  p <- matrices$p
  
  trts <- LETTERS[1:nrow(est)]
  
  weights <- switch(weighted + 1, NULL, TRUE) # If weighted==FALSE, weights=NULL. If weighted==TRUE, weights=TRUE.
  
  
  ################### COMMUNITY DETECTION WITHOUT SIMULATION ###############
  
  if(method=="community") 
  {
    
    rownames(est) <- trts
    colnames(est) <- trts
    rownames(se) <- trts
    colnames(se) <- trts
    rownames(p) <- trts
    colnames(p) <- trts
    
    # Create vectors for finding quantiles. Take the lower triangle of each matrix:
    est.vec <- est[c(lower.tri(est))] 
    se.vec <- se[c(lower.tri(se))]
    p.vec <- p[c(lower.tri(p))]
    
    # Find quantiles
    est.quant <- quantile(est.vec, probs=quants)
    se.quant <- quantile(se.vec, probs=quants)
    p.quant <- quantile(p.vec, probs=quants)
    
    # Create three empty lists, to be populated with the adjacency matrices:
    est.adj.mat <- vector("list", length(quants))
    se.adj.mat <- vector("list", length(quants))
    p.adj.mat <- vector("list", length(quants))
    
    names(est.adj.mat) <- paste("est", as.character(100*quants), sep="")
    names(se.adj.mat) <- paste("se", as.character(100*quants),  sep="")
    names(p.adj.mat) <- paste("p", as.character(100*quants),  sep="")
    
    #par(mfrow=c(1,3))
    
    for(i in 1:length(quants))
    {
      
      ################## UNWEIGHTED ADJACENCY MATRICES #########################     
      
      if(weighted==FALSE)
      {
        
        # Create unweighted adjacency matrices: If est/se/p is lower than threshold -> 1; otherwise -> 0
        est.adj.mat[[i]] <- 1*(est<est.quant[i]) 
        se.adj.mat[[i]] <- 1*(se<se.quant[i])
        p.adj.mat[[i]] <- 1*(p<p.quant[i])
        
      }
      
      ################## WEIGHTED ADJACENCY MATRICES #########################     
      
      if(weighted==TRUE)
      {
        # In weighted adjacency matrices, larger numbers -> more weight -> more similar, 
        # therefore take reciprocal of est, SE and (compliment of) p for non-zero adjacencies (i.e. > 0)
        
        est.adj.mat[[i]] <- est
        est.adj.mat[[i]][est.adj.mat[[i]]>=est.quant[i]] <- 0 # If estimated trt effect contrast is large, no connection
        est.adj.mat[[i]][est.adj.mat[[i]] > 0] <- 1 / est.adj.mat[[i]][est.adj.mat[[i]] > 0] 
        
        se.adj.mat[[i]] <- se
        se.adj.mat[[i]][se.adj.mat[[i]]>=se.quant[i]] <- 0 # If estimated SE is large, no connection
        se.adj.mat[[i]][se.adj.mat[[i]] > 0] <- 1 / se.adj.mat[[i]][se.adj.mat[[i]] > 0] 
        
        p.adj.mat[[i]] <- p
        p.adj.mat[[i]][p.adj.mat[[i]]>=p.quant[i]] <- 0 # If this p is large (ie original p value is small), no connection
        p.adj.mat[[i]][p.adj.mat[[i]] > 0] <- 1 / p.adj.mat[[i]][p.adj.mat[[i]] > 0]
        
      }
      
      
      # Graphs:
      # g.est <- graph_from_adjacency_matrix(est.adj.mat[[i]], mode="undirected", diag=F, weighted=weights)
      g.se <- graph_from_adjacency_matrix(se.adj.mat[[i]], mode="undirected", diag=F, weighted=weights)
      #g.p <- graph_from_adjacency_matrix(p.adj.mat[[i]], mode="undirected", diag=F, weighted=weights)
      
      if(fast)
      {
        #ceb.est <- cluster_edge_betweenness(g.est)
        ceb.se <- cluster_edge_betweenness(g.se)
        #ceb.p <- cluster_edge_betweenness(g.p)
      }
      else
      {
        #ceb.est <- cluster_optimal(g.est)
        ceb.se <- cluster_optimal(g.se)
        #ceb.p <- cluster_optimal(g.p)
      }
      
      
      ########## PLOTTING ############
      
      if(fix.layout==FALSE)
      {
        #layout.est <- layout_nicely(g.est) 
        layout.se <- layout_nicely(g.se) 
        #layout.p <- layout_nicely(g.p)
      }
      else
      {
        browser()
        #layout.est <- layout_in_circle(g.est) 
        layout.se <- layout_in_circle(g.se) 
        #layout.p <- layout_in_circle(g.p)   
      }
      
      #  plot(ceb.est, g.est, layout=layout.est)
      # title(paste("Est: ", names(est.quant[i]), " (", round(est.quant[i],2), ")", sep=""), cex.main=2)
      
      plot(ceb.se, g.se, layout=layout.se)
      title(paste("SE: ", names(se.quant[i]), " (", round(se.quant[i],2), ")", sep=""),cex.main=2)
      
      
      
      #plot(ceb.p, g.p, layout=layout.p)
      #title(paste("P: ", names(p.quant[i]), " (", round(1-p.quant[i],2), ")", sep=""), cex.main=2)
      
      
      if(i!=length(quants))
      {
        cat("Press [enter] to continue")
        line <- readline()
      }
      
      if(save==TRUE) dev.print(pdf, paste(100*quants[i], " percentile",".pdf", sep=""))
      
    } # end of for loop i in 1:length(quants)
  } # end of if statement <if method=="community">
  
  
  
  ####################### HEATMAP/BOOTSTRAPPING METHOD #######################
  
  
  if(method=="bootstrap")
  {
    
    # Create bootstrap replications (p-1 columns, n rows):
    pred <- rmvnorm(n, mean=ests, sigma=Cov) 
    
    # Create lists of length n for storing
    # - the n square matrices of relative treatment effect estimates
    # - the n vectors containing the lower triangle of those matrices
    # - the n adjacency matrices:
    
    rel.effects.list <- vector("list", n)
    est.vec.pred <- vector("list", n)
    est.adj.mat.list <-  vector("list", n)
    
    # Create a list of matrices, each matrix to hold n rows of treatment memberships:
    membership.list <- vector("list", length(quants)) 
    
    # All possible pairs of treatments:      
    pair <- combn(trts, 2)
    pair <- t(pair)
    
    for(i in 1:n)
    { 
      # Create matrix of relative effects:
      rel.effects.list[[i]] <- fns(pred[i,])
      rel.effects.list[[i]] <- abs(rel.effects.list[[i]])
      # Make vector out of lower triangle of matrix:
      est.vec.pred[[i]] <- c(rel.effects.list[[i]])[c(lower.tri(rel.effects.list[[i]]))] 
    }
    
    # Find the quantiles for each bootstrap replication:
    
    est.quant.pred <- lapply(est.vec.pred, quantile, probs=quants) 
    
    for(j in 1:length(quants))
    {
      if(weighted==TRUE) est.adj.mat.list <- rel.effects.list
      
      for(i in 1:n)
      {
        
        if(weighted==TRUE)
        {
          est.adj.mat.list[[i]][est.adj.mat.list[[i]]>=est.quant.pred[[i]][j]] <- 0 # If estimated trt effect contrast is large, no connection
          est.adj.mat.list[[i]][est.adj.mat.list[[i]] > 0] <- 1 / est.adj.mat.list[[i]][est.adj.mat.list[[i]] > 0] # Take reciprocal of non-zero elements.
        }
        else{
          est.adj.mat.list[[i]] <- 1*(rel.effects.list[[i]]<est.quant.pred[[i]][j])  # If trt contrast is below threshold, -> 1
        }
      } # end of 1:n
      
      
      # Create and store all n memberships for current quantile:
      
      ceb.list <- if(fast) lapply(est.adj.mat.list, function(x)  cluster_edge_betweenness(graph_from_adjacency_matrix(x, mode="undirected", diag=F, weighted = weights))$membership)
      else lapply(est.adj.mat.list, function(x)  cluster_optimal(graph_from_adjacency_matrix(x, mode="undirected", diag=F, weighted = weights))$membership)
      
      # Collapse n memberships into a matrix of n rows
      membership.list[[j]] <- matrix(data=unlist(ceb.list), byrow=TRUE, nrow=n,  dimnames=list(1:n, LETTERS[1:(ncol(pred)+1)]))
      
      # For the current quantile, there now exists n sets of detected communities. How often does each unique pair of trts appear in the same community?
      n.same.comm <- apply(pair, 1, function(x) sum(membership.list[[j]][,x[1]]==membership.list[[j]][,x[2]]))
      prob.same.comm <- n.same.comm/nrow(pred)
      
      # Create symmetric matrix of these proportions:
      mat <- matrix(NA, nrow=ncol(pred)+1, ncol=ncol(pred)+1, dimnames=list(trts, trts))
      mat[c(lower.tri(mat))] <- prob.same.comm
      mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
      
      
      if(mean(mat, na.rm=T)!=1) # Make sure at least some probabilities are != 1 (requirement of heatmap.2 function)
      {
        fix.rows <- ifelse(j==1 & fix.layout==FALSE, yes=TRUE, no=FALSE) # For reordered layout (used in call to heatmap.2 below).
        
        if(j>1 & fix.layout==FALSE) mat <- mat[row.order, row.order] # For reordered layout, use row order from first heatmap (ie when j==1)
        
        heat <- heatmap.2(mat, dendrogram = "none", trace="none", cellnote=round(mat,2), notecol=1,
                          density.info = "none", col=rev(heat.colors(16)), symm=TRUE,
                          breaks=seq(from=0, to=1, length=17),
                          notecex=4*(1/sqrt(ncol(mat))),
                          srtCol=0,
                          main=paste("Proportion of times\ntreatments in same community (", names(est.quant.pred[[1]][j]), " quantile)" ,sep=""),
                          Rowv=fix.rows, Colv=fix.rows, revC=fix.rows
        )
        
        if(j==1) row.order <- heat$rowInd  # Recording row order of first heatmap.
        
        
        if(j!=length(quants))
        {
          cat("Press [enter] to continue")
          line <- readline()
        }
        
        filetext <- ifelse(weighted, "weighted", "unweighted")
        if(save==TRUE) dev.print(pdf, paste("Predictive ests (", filetext, "), ", 100*quants[j], ".pdf", sep=""))
        
      }
      
      # If all probabilities = 1, don't attempt heatplot:
      else print(paste("For", names(est.quant.pred[[1]][j]), " quantile, all treatments in the same community for all simulations; no heatmap produced."),q=F)
    }
  } # end of method=="bootstrap" if statement
  
}


## modifies JAGS output and produces the league tables
get_results <- function(JAGSobject, 
                        parameter=parameter, 
                        forestplot=F, 
                        treatnames=NA, 
                        rounding=3){
  resultstable=JAGSobject$BUGSoutput$summary
  allvariablenames=rownames(resultstable)
  rowsmatching=substr(allvariablenames,1,nchar(parameter))
  rowstokeep=startsWith(rowsmatching,parameter)
  
  resultstabletokeep=resultstable[rowstokeep,c(1,3,7,2)]
  rowstokeep2=startsWith(dimnames(resultstabletokeep)[[1]], paste(parameter,"[",sep=""))
  resultstabletokeep=resultstabletokeep[rowstokeep2,]
  library(stringr)
  tosplit=unlist(strsplit(rownames(resultstabletokeep),","))
  tosplit2 <- as.numeric(str_extract(tosplit, "[0-9]+"))
  nroftreatments=max(tosplit2)
  location=matrix(tosplit2,ncol=2,byrow=T)
  meanmat=CImat=sdmat=matrix(NA,nrow=nroftreatments,ncol=nroftreatments)
  
  for(i in 1:nrow(location)){
    meanmat[location[i,1],location[i,2]]=resultstabletokeep[i,1]
    sdmat[location[i,1],location[i,2]]=resultstabletokeep[i,4]
    CImat[location[i,1],location[i,2]]=resultstabletokeep[i,3]#high CI
    CImat[location[i,2],location[i,1]]=resultstabletokeep[i,2]#low CI
    
  }
  if(forestplot){
    library(metafor)
    library(ggforestplot)
    slab1=rep(1:(nroftreatments-1),(nroftreatments-1):1)
    a=t(sapply(1:nroftreatments,rep,nroftreatments))
    slab2=a[lower.tri(a,F)]
    slab=paste(slab1,"vs",slab2,sep="")
    df=data.frame(estimate=meanmat[upper.tri(meanmat)],se=(CImat[upper.tri(CImat)]-CImat[lower.tri(CImat)])/3.92,
                  name=slab)
    forest(x=meanmat[upper.tri(meanmat)], ci.lb=CImat[lower.tri(CImat)],ci.ub=CImat[upper.tri(CImat)], slab=slab,xlab="Network meta-analysis results")
    
    
    
  }
  meanmat=round(meanmat,rounding)
  CImat=round(CImat,rounding)
  
  #create to print
  Ttreat=dim(meanmat)[1]
  toprintmat=matrix(nrow=Ttreat,ncol=Ttreat)
  
  for(i in c(1:(Ttreat-1)))
  {for (j in c((i+1):Ttreat))
  {
    toprintmat[i,j]=paste(meanmat[i,j],"(",CImat[j,i],",",CImat[i,j],")",sep="")
    toprintmat[j,i]=paste(c(-meanmat[i,j]),"(",c(-CImat[i,j]),",",c(-CImat[j,i]),")",sep="")
  }}
  if(!missing(treatnames)){
    diag(meanmat)=treatnames
    diag(CImat)=treatnames
    diag(toprintmat)=treatnames
  }
  
  list(Means=meanmat,CI=CImat, leaguetable=toprintmat)
}


## produces the density plots in appendix 1 figures 2-6 
density_plot <- function(data, treat){
  
  if(treat=="Aripiprazole"){
    data1 <- data %>% 
      filter(comparison=="Aripiprazole vs Placebo")
    
  }
  else if(treat=="Olanzapine"){
    
    data1 <- data %>% 
      filter(comparison=="Olanzapine vs Placebo")
  }
  else if(treat=="Paliperidone"){
    
    data1 <- data %>% 
      filter(comparison=="Paliperidone vs Placebo")
  }
  else if(treat=="Quetiapine"){
    
    data1 <- data %>% 
      filter(comparison=="Quetiapine vs Placebo")
  }
  else if(treat=="Risperidone"){
    
    data1 <- data %>% 
      filter(comparison=="Risperidone vs Placebo")
  }    
  
  ## ensure the order of the models in the graph  
  data1$model=factor(data1$model,
                     levels = c(
                       "Naive synthesis","Data based beta-No DW","Data based beta-RoB DW","Data based beta-NCT DW",
                       "Expert opinion beta-No DW","Expert opinion beta-RoB DW","Expert opinion beta-NCT DW",
                       "NMA with non-informative priors")
  )  
  
  g <- ggplot(data1,aes(x=values, color=evidence,fill=evidence)) +
    facet_wrap(~model,ncol = 2)+
    scale_fill_lancet()+
    scale_color_lancet()+
    theme_bw(base_size = 10)+
    geom_density(alpha=0.5,size=0.6)+
    xlab("Standardized mean difference")+
    ylab("Density")+
    ggtitle(paste(treat," vs Placebo",sep=""))+
    theme(
      legend.text=element_text(size=8),
      legend.title = element_text(size=10),
      plot.title = element_text(hjust = 0.5),
      strip.background = element_rect(colour = "lightblue",fill = "lightblue"),
      strip.text = element_text(size=8),
      
      
      
    )+
    labs(color="Type of Evidence",fill="Type of Evidence")
  
  return(g)
}


## produces the heatmap in appendix 1 figure 7
heatmap <- function(data){
  
  # ensure the order of the treatments in the final graph
  data$model_type=factor(data$model_type,
                         levels = c(
                           "Naive synthesis",
                           "Data based beta-No DW",
                           "Data based beta-RoB DW",
                           "Data based beta-NCT DW",
                           "Expert opinion beta-No DW",
                           "Expert opinion beta-RoB DW",
                           "Expert opinion beta-NCT DW",
                           "NMA with non-informative priors"
                         ))
  
  ## create order according to NMA with non-informative priors model
  naive_order <- data %>%
    filter(model_type=="NMA with non-informative priors")
  
  order <- naive_order$drug
  
  data$drug <- factor(data$drug,levels=order)
  
  # create the heatmap
  g <- ggplot(data, aes(drug, fct_rev(model_type), fill= mean)) + 
    geom_tile(color = "black") +
    geom_text(aes(label = paste(format(round(mean,digits = 2),nsmall=2) )),
              size=3,color="black") +
    
    scale_fill_gradient(low="yellow", high="red") +
    
    guides(fill = guide_colourbar(label = FALSE,
                                  ticks = FALSE,
                                  barwidth = 1.2))+           
    
    labs(x="",y="",fill = "SUCRAS")+
    scale_y_discrete(expand=c(0,0))+
    theme_bw()+
    theme(
      #bold font for legend text
      legend.text=element_text(face="bold"),
      #set thickness of axis ticks
      axis.ticks=element_line(size=1),
      #remove plot background
      plot.background=element_blank(),
      #remove plot border
      panel.border=element_blank(),
      axis.text.x=element_text(angle=0,hjust=0.5,size=8),
      
      axis.text.y=element_text(size=10),
      legend.position = "right",
      legend.direction = "vertical",
      legend.title = element_text(hjust = 0.5)
      
    )
  
  return(g)
  
}

# create an indicator parameter which is 1 when a study is of high RoB and 0 otherwise
ind_rob = function(data){
  
  u = length(unique(data$Study_No))
  
  u1 = as.data.frame(table(data$Study_No))
  
  target = unique(data$Study_No)
  
  u1 = u1[match(target, u1$Var1), ]
  
  sub = list()
  
  for (i in 1:nrow(u1)) {
    sub[[i]] = rep(i, u1$Freq[i])
  }
  
  sub = unlist(sub)
  
  data$sub = sub
  
  r = list()
  
  rob = list()
  
  ind_robf = list()
  
  for (i in 1:u) {
    r[[i]] = data %>%
      filter(sub == i)
    
    rob[[i]] = r[[i]]$overall_rob
    
    rob[[i]] = rob[[i]][1]
    
    ind_robf[[i]] = ifelse(rob[[i]] == "high", 1, 0)
    
  }
  
  res = unlist(ind_robf)
  
  return(res)
}


## creates the indicator variable used in the model settings that involve NCT downweight
ind_drugs = function(data){
  
  sub.data = subset(all.data,(all.data$GeneralPopulation==1 | all.data$ChildAdolesc==1) & !is.na(all.data$OverallEfficM) & !is.na(all.data$OverallEfficSD) & !is.na(all.data$OverallEfficN))
  
  ca.data = subset(sub.data,sub.data$ChildAdolesc==1)
  
  gp.data = subset(sub.data,sub.data$GeneralPopulation==1)
  
  unique(ca.data$Drug)
  
  unique(gp.data$Drug)
  
  #identify the common drugs for CA and GP
  drugs = intersect(unique(ca.data$Drug),unique(gp.data$Drug))
  
  u = length(unique(gp.data$Study_No))
  
  u1 = as.data.frame(table(gp.data$Study_No))
  
  target = unique(gp.data$Study_No)
  
  u1 = u1[match(target, u1$Var1),]
  
  sub = list()
  
  for(i in 1:nrow(u1)){
    sub[[i]] = rep(i,u1$Freq[i])
  }
  
  sub = unlist(sub)
  
  gp.data$sub = sub
  
  r=list()
  
  d=list()
  
  d1=list()
  
  q=list()
  
  ind=list()
  
  for(i in 1:u){
    
    r[[i]] = gp.data %>% 
      filter(sub==i)
    
    d[[i]] = r[[i]]$Drug
    
    d1[[i]] = intersect(d[[i]],drugs)
    
    q[[i]] = length(d[[i]])-length(d1[[i]])
    
    if(q[[i]]==0){
      
      ind[[i]]=0
    }else{
      ind[[i]]=1
    }
    
  }
  
  ind=unlist(ind)
  
  return(ind)
}

# produces the forest plot presented in Figure 2 of the main paper
forest_plot <- function(data){
  
  ## make sure that means and credible interval bounds are treated as numeric values  
  data$mean <- as.numeric(data$mean)
  
  data$lower <- as.numeric(data$lower)
  
  data$upper <- as.numeric(data$upper)
  
  ## assign some pseudo values for NA's (those value will never be plotted in the graph)
  data$mean=ifelse(is.na(data$mean),100,data$mean)
  
  data$lower=ifelse(is.na(data$lower),90,data$lower)
  
  data$upper=ifelse(is.na(data$upper),110,data$upper)
  
  # Create the order of the models in the final graph
  data$model_type=factor(data$model_type,levels =c("Naive synthesis","Data based beta-No DW","Data based beta-RoB DW","Data based beta-NCT DW",
                                                   "Expert opinion beta-No DW","Expert opinion beta-RoB DW","Expert opinion beta-NCT DW",
                                                   "NMA with non-informative priors","Direct comparison"))
  
  
  # Plot the estimates in a forestplot
  g=ggplot(data=data, aes(y=fct_rev(model_type), x=mean, xmin=lower, xmax=upper,colour=model_type,shape=model_type))+
    theme_bw(base_size = 13)+
    facet_wrap(~comparison,strip.position = c("top"))+
    geom_vline(xintercept = 0,linetype=2)+
    geom_point(size=3, alpha = 1) +
    geom_errorbarh(height=.3,size=0.8)+
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
          
          strip.background = element_rect(
            color="black", fill="lightblue", size=1, linetype="solid"),
          strip.text = element_text(face="bold"),
          axis.ticks.y = element_blank(),
          legend.position = "bottom",
          legend.box.background = element_rect(color="black", size=0.5),
          legend.box.margin = margin(0.1, 6, 6, 6)
    )+
    scale_colour_manual(name = "Models",
                        
                        labels = c("Naive synthesis","Data based beta-No DW","Data based beta-RoB DW","Data based beta-NCT DW",
                                   "Expert opinion beta-No DW","Expert opinion based beta-RoB DW","Expert opinion beta-NCT DW",
                                   "NMA with non-informative priors","Direct comparison"),
                        
                        values = c("black","#0073C2FF","#CD534CFF","orange","#0073C2FF","#CD534CFF","orange","darkgreen","purple"))+
    
    scale_shape_manual(name = "Models",
                       
                       labels = c("Naive synthesis","Data based beta-No DW","Data based beta-RoB DW","Data based beta-NCT DW",
                                  "Expert opinion beta-No DW","Expert opinion based beta-RoB DW","Expert opinion beta-NCT DW",
                                  "NMA with non-informative priors","Direct comparison"),
                       values = c(19,15,15,15,17,17,17,19,19)
                       
    )+
    xlab("")+
    coord_cartesian(xlim=c(-1.5, 1.5))+
    scale_x_continuous(breaks = c(-1.5, -1, -0.5, 0, 0.5, 1,1.5))
  
  return(g)
  
}


