## main sampler:
Net.fit.Murray.sampler = function(Z, X, spp_design # some trait based measurements
                    ,mcmc.iter = 10000, burn.in = 10 
                    ,vars_prop = list( beta_occu = rep(1e-5,2 * ncol(X))
                                        ,beta_graph = rep(1e-5,length(spp_design))
                                        ,eta_intra = rep(1e-5,nspp)
                                        ,eta_inter = rep(1e-5,nspp*(nspp-1)/2)
                                        ,d_intra=rep(1e-5,nspp)
                                        ,d_inter = rep(1e-5,nspp)
                                        )
                    ,para_prior = list( beta_occu = rep(1000,2 * ncol(X))
                                        ,beta_graph = rep(1e-5,length(spp_design))
                                        ,eta_intra = rep(1e-1,nspp)
                                        ,eta_inter = rep(1000,nspp*(nspp-1)/2)
                                        ,d_intra=rep(1000,nspp)
                                        ,d_inter = rep(1000,nspp)
                                        )
          					,uni_prior = F
                    ,distM,link_map
                    ,dist_mainland , link_mainland
                    ,int_range_intra="nn",int_range_inter="nn"
                    ,seed = 42,ini,thin.by = 100,report.by=100,nIter=100, Importance = F,method="CFTP"){ # ini has same formate of theta
  
  cat("Initializing...\n\n")
  require(coda)
  require(Matrix)
  require(RcppArmadillo)
  source("./R/misc.R")
  source("./R/Murray_ratio.R")
  Rcpp::sourceCpp("./src/IsingCpp_CFTP_sparse.cpp")
  set.seed(seed)
  if(uni_prior) getlogprior = getlogprior_uniform
  else getlogprior = getlogprior_normal
	Z = matrix(Z,ncol = 1)
  vars_prop = vars_prop[names(ini)]
  para_prior = para_prior[names(ini)]
  
  cat("MCMC reported every",report.by,"iterations and thinned by",thin.by,"\n\n")
  
  nsite = (nrow(distM))
  
  theta_curr = ini
  spp_neig = 1 *( spp_design[[1]]!=0 )
  
  nspp = nrow(spp_design[[1]])
  nrep = length(spp_design[[1]])
  
    #theta_tuta=ini
  ini = lapply(ini,as.matrix)
  theta.mcmc = list()
  for(i in 1:length(ini)){
    theta.mcmc[[i]] = mcmc(matrix(nrow = floor(mcmc.iter/thin.by),ncol = length(ini[[i]])),thin = thin.by)
    
  }
  
  
	
  
  names(theta.mcmc) = names(ini)
  
  
	
  
  

  n_para_group = length(ini)
  spp_mat_temp = lapply(1:length(spp_design),function(i,spp_design,beta){
    spp_design[[i]] * beta[i]  
  },spp_design,theta_curr$beta_graph)
  
  ini$spp_mat = Reduce("+",spp_mat_temp)	
	
  MRF_curr = getMRF(ini,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter)
  
  theta_curr = ini	
	
  cat("Burn in...\n")
  
  accept_theta_occu = 0


  low_acc_theta_occu = 0


  timing = proc.time()
  
  for(i in 1:burn.in){# to burn in 
    #propose theta 
    theta_prop = theta_curr

    for(j in c(1:(length(theta_curr)-1))){ 
      	theta_prop[[j]] = matrix( rnorm(length(theta_curr[[j]]),mean = 0,sd = sqrt(vars_prop[[j]])),nrow(theta_curr[[j]]),ncol(theta_curr[[j]]) )+ theta_curr[[j]]
    }
	
    spp_mat_temp = lapply(1:length(spp_design),function(i,spp_design,beta){
      spp_design[[i]] * beta[i]  
    },spp_design,theta_prop$beta_graph)
  
    theta_prop$spp_mat = Reduce("+",spp_mat_temp)	

    
    # MH ratio
    
	MRF_prop = getMRF(theta_prop,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter)
	Z_temp = rIsingOccu_multi(MRF_prop,n=1,method = method,nIter)
	
    Murray_ratio=Murray_ratio_occu_theta(MRF_curr ,MRF_prop, getlogprior(theta_prop,theta_curr,para_prior)
                        ,Z
                        ,Z_temp
                        ,para_prior
                        ,distM,link_map
                        ,dist_mainland,link_mainland
                        ,int_range_intra,int_range_inter)
    r = runif(1)
    if(is.na(Murray_ratio)){
      Murray_ratio = 0
      #cat("occuNA\n")
    }
    if(Murray_ratio<exp(-10)) low_acc_theta_occu = low_acc_theta_occu + 1
    if(r<=Murray_ratio){
      theta_curr=theta_prop
      MRF_curr = MRF_prop
      accept_theta_occu = accept_theta_occu + 1
    }
      
    
    if(i%%report.by == 0) {
      
      cat("Burn in iteration",i-report.by+1,"to",i,":\n\n")
      #cat("    # of Z proposed for imperfect detection: ",propose_Z_num,"\n")
      cat("    # of occupancy theta acceptance: " , accept_theta_occu,"\n")
      cat("    # of occupancy acceptance ratio <exp(-10): ",low_acc_theta_occu,"\n\n")
      timing = proc.time()- timing
      cat("Time used in this" ,report.by,":",timing[1],"s\n")
      cat("\n\n")
      accept_Z = 0
      accept_Z_missing_obs = 0
      accept_theta_occu = 0
      accept_theta_det = 0
      low_acc_Z = 0
      low_acc_Z_missing_obs = 0
      low_acc_theta_occu = 0
      low_acc_theta_det = 0
      propose_Z_num = 0
      propose_Z_missing_obs = 0
      timing = proc.time()
      }
  }
  cat("Start sampling...\n")
  accept_Z = 0
  accept_Z_missing_obs = 0
  accept_theta_occu = 0
  accept_theta_det = 0
  low_acc_Z = 0
  low_acc_Z_missing_obs = 0
  low_acc_theta_occu = 0
  low_acc_theta_det = 0
  propose_Z_num = 0
  propose_Z_missing_obs = 0
  timing = proc.time()  
  for(i in 1:(mcmc.iter)){ # for to save 
    #propose theta 
    theta_prop = theta_curr

    for(j in c(1:(length(theta_curr)-1))){ # no detection proposing
      	theta_prop[[j]] = matrix( rnorm(length(theta_curr[[j]]),mean = 0,sd = sqrt(vars_prop[[j]])),nrow(theta_curr[[j]]),ncol(theta_curr[[j]]) )+ theta_curr[[j]]
    }
	
    spp_mat_temp = lapply(1:length(spp_design),function(i,spp_design,beta){
      spp_design[[i]] * beta[i]  
    },spp_design,theta_prop$beta_graph)
  
    theta_prop$spp_mat = Reduce("+",spp_mat_temp)	

    
    # MH ratio
    
	MRF_prop = getMRF(theta_prop,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter)
	Z_temp = rIsingOccu_multi(MRF_prop,n=1,method = method,nIter)
	
    Murray_ratio=Murray_ratio_occu_theta(MRF_curr ,MRF_prop, getlogprior(theta_prop,theta_curr,para_prior)
                        ,Z
                        ,Z_temp
                        ,para_prior
                        ,distM,link_map
                        ,dist_mainland,link_mainland
                        ,int_range_intra,int_range_inter)
    r = runif(1)
    if(is.na(Murray_ratio)){
      Murray_ratio = 0
      #cat("occuNA\n")
    }
    if(Murray_ratio<exp(-10)) low_acc_theta_occu = low_acc_theta_occu + 1
    if(r<=Murray_ratio){
      theta_curr=theta_prop
      MRF_curr = MRF_prop
      accept_theta_occu = accept_theta_occu + 1
    }
      
    
    if(i %% thin.by==0){
      for(j in 1:(length(theta_curr)-1)){
        theta.mcmc[[j]][i/thin.by,] =as.vector( theta_curr[[j]])
      } # saving the results
    }
    
    
    
    if(i%%report.by == 0) { # reporting
      cat("Sampling iteration",i-report.by+1,"to",i,":\n\n")
      cat("    # of occupancy theta acceptance: " , accept_theta_occu,"\n")
      cat("    # of occupancy acceptance ratio <exp(-10): ",low_acc_theta_occu,"\n\n")
      timing = proc.time()-timing
      cat("Time used in this" ,report.by,":",timing[1],"s\n")
      cat("\n\n")
      accept_theta_occu = 0
      low_acc_theta_occu = 0
      timing = proc.time()
      }
  }
  
  theta.mean =lapply(theta.mcmc,function(thetaa){ apply(thetaa,2,mean,na.rm=T)})
  
  res = list(theta.mcmc = theta.mcmc
             ,means = theta.mean
             ,vars=vars_prop
             ,distM = distM
             ,dist_mainland = dist_mainland
             ,linkmap = link_map
             ,link_mainland = link_mainland
             ,interaction.range =list( inter = int_range_inter,intra = int_range_intra)
             ,envX=X, spp_design,uni_prior=uni_prior,para_prior = para_prior,method = method)
  return(res)
}

