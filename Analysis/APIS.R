source("./R/misc.R")
source("./R/Main_Sampler.R")
source("./R/Simu_data_Sampling.R")
require(Matrix)
require(Rcpp)

# fisher marten coyote red_fox gray_fox bobcat in kg, male
Bodysize = matrix( c(4.75,0.85,12.5,8.1,5.3,12.35))

Bodysize_mat = Bodysize %*% matrix(1,1,6) - matrix(1,6,1) %*% t(Bodysize) 
Bodysize_mat = abs(Bodysize_mat)
Bodysize_mat = Bodysize_mat/max(Bodysize_mat)

Zs = read.csv("./Data/PA_all_full.csv")
Zs = Zs[,c("Fisher","Marten","Coyote","Fox_red","Fox_gray","Bobcat")]
Z_vec = matrix( as.matrix( Zs),ncol = 1)

dist_mainland = read.csv("./Data/dist_to_mainland.csv")$X0
dist_mainland = (dist_mainland-min(dist_mainland))/(max(dist_mainland)-min(dist_mainland))

exp_mainland = exp(-2*dist_mainland)

envX = cbind(1,exp_mainland)

link_inner = Matrix(0,155,155,sparse = T)
link_outer = link_inner
link_map = 
  list(inter = 0*link_outer,
       intra = 0*link_inner) # no spatial auto-correlation for now


nspp = 6

spp_neig = Matrix(1,nspp,nspp,sparse = T)
diag(spp_neig) = 0

spp_design  = list(spp_neig,Bodysize_mat)

theta = list(beta_occu = runif(nspp*ncol(envX),-1,1),
             beta_graph = c(.3,-.5),
             eta_intra = rep(0,nspp),
             eta_inter = rep(0,nspp))


vars_prop = list( beta_occu = 5e-4
                  ,beta_graph = 5e-4 # no extra det thing
                  ,eta_intra = rep(1e-10,nspp)
                  ,eta_inter = c(1e-10,1e-10))

para_prior = list( beta_occu = rep(1000,nspp * ncol(envX))
                   ,beta_graph = rep(1000,2)
                   ,eta_intra = rep(1000,nspp)
                   ,eta_inter = rep(1000,nspp)
)


kk = Net.fit.Murray.sampler(Z = Z_vec, X = envX
                            , spp_design = spp_design
                            , mcmc.iter = 500, burn.in = 5
                            , vars_prop = vars_prop
                            , para_prior = para_prior
                            
                            , uni_prior = F
                            , distM=link_map[[1]],link_map=link_map
                            , dist_mainland = 0*dist_mainland , link_mainland = 0* dist_mainland 
                            , int_range_intra="nn",int_range_inter="nn"                                          
                            , seed = 42
                            , ini = theta,thin.by = 1,report.by = 30,nIter = 50,method = "CFTP")


