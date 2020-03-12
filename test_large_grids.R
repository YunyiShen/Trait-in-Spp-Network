source("./R/misc.R")
source("./R/Main_Sampler.R")
source("./R/Simu_data_Sampling.R")
require(Matrix)
require(Rcpp)
Rcpp::sourceCpp("src/IsingCpp_CFTP_sparse.cpp")


## generate graph 
n_grids = 12 # 15 by 15 grid system
link_inner = adjacency.matrix(n_grids) # nearest neighborhood 
link_outer = link_inner

###### True Parameter Setting ######
nspp = 6
spp_neig = matrix(1,nspp,nspp)
diag(spp_neig)=0
spp_neig = as(spp_neig,'dsCMatrix')
envX = matrix(1,n_grids^2,1)
envX = cbind(envX ,rnorm(n_grids^2))

theta = list(beta_occu = runif(nspp*ncol(envX),-1,1),
             beta_graph = c(.3,-.5),
             eta_intra = rep(0,nspp),
             eta_inter = rep(0,nspp))

link_map = 
  list(inter = 0*link_outer,
       intra = 0*link_inner) # no spatial auto-correlation for now

nsite = n_grids^2


distM_mainland = matrix(0,nsite,1)

spp_mat_2 = matrix(runif(nspp^2,0,.5),nspp,nspp)
spp_mat_2 = spp_neig * (spp_mat_2 + t(spp_mat_2))
spp_design  = list(spp_neig,spp_mat_2)

theta$spp_mat = spp_design[[1]]*theta$beta_graph[1] + spp_design[[2]]*theta$beta_graph[2]

###### Simulate Data ######
set.seed(42)
MRF = getMRF(theta,envX,distM = 0*link_map[[1]],link_map,link_mainland = distM_mainland, dist_mainland = distM_mainland ,
			 int_range_intra="nn",int_range_inter="nn")

Z_simu = IsingSamplerCpp(1, MRF$A, MRF$thr, 1, 30, c(-1,1), F,NA+MRF$thr) ## take true occupancy

###### Run the Model! ######

vars_prop = list( beta_occu = 5e-4
                  ,beta_graph = 5e-4 # no extra det thing
                  ,eta_intra = rep(1e-10,nspp)
                  ,eta_inter = c(1e-10,1e-10))

para_prior = list( beta_occu = rep(1000,nspp * ncol(envX))
                   ,beta_graph = rep(1000,2)
                   ,eta_intra = rep(1000,nspp)
                   ,eta_inter = rep(1000,nspp)
                   )


kk = Net.fit.Murray.sampler(Z = t(Z_simu), X = envX, spp_design = spp_design
                                            , mcmc.iter = 50000, burn.in = 5000
                                            , vars_prop = vars_prop
                                            , para_prior = para_prior

                                            , uni_prior = F
                                            , distM=link_map[[1]]+1,link_map=link_map
                                            , dist_mainland = distM_mainland , link_mainland =  distM_mainland 
                                            , int_range_intra="nn",int_range_inter="nn"                                          
                                            , seed = 42
                                            , ini = theta,thin.by = 10,report.by = 100,nIter = 50,method = "CFTP")


















