## misc helper functions
getintralayerGraph = function(distM,link_map,eta,d,int_range = "exp",spp_mat) #it can be used multiple times for interislan and intra-island
{
  # pass all graphs as sparse matrix in package Matrix
  nspp = nrow(spp_mat) # which is the interspecific neighborhood matrix
  A = list() # intralayer graphs are passed using lists
  link_map = as(as.matrix(link_map),"dgCMatrix")
  if(int_range=="arth"){
    A = lapply(1:nspp,function(i,eta,distM,d){
      eta[i]*as.matrix(1/((distM)^(2+d[i])))
    },eta,distM,d)
  }
  else{
    if(int_range=="exp"){
	  A = lapply(1:nspp,function(i,eta,d,distM,link_map){
		At = eta[i]*as.matrix(exp(-exp(d[i])*distM)) * (link_map)
	    diag(At)=0
	    return(At)
	  },eta,d,distM,link_map)
    }
    else{
      if(int_range=="nn"){
	  A = lapply(1:nspp, function(i,eta,link_map){
		    eta[i]*as.matrix((link_map))
	    },eta,link_map)
      }
      else{
        #print("int_range must be exp or arth, will assume exp")
		A = lapply(1:nspp,function(i,eta,d,distM,link_map){
		  At = eta[i]*as.matrix(exp(-exp(d[i])*distM)) * (link_map)
	      diag(At)=0
	      return(At)
	    },eta,d,distM,link_map)
      }
    }
  }
  return(A) # if link map is sparse, then A is sparse
} 
  # passed 2019/3/18

getfullGraph = function(A_ex,A_in,spp_mat){
  nspp = nrow(spp_mat)
  nsite = nrow(A_ex[[1]])
  A = Matrix(0,nspp*nsite,nspp*nsite,sparse = T)
  for(i in 2:nspp-1){
    A[1:nsite + (i-1)*nsite,1:nsite + (i-1)*nsite]=A_ex[[i]] + A_in[[i]] # diagonal part
    for(j in (i+1):nspp){
      
      diag(A[1:nsite + (i-1)*nsite,1:nsite + (j-1)*nsite])=spp_mat[i,j]
      diag(A[1:nsite + (j-1)*nsite,1:nsite + (i-1)*nsite])=spp_mat[j,i]
      
    }
  }
  i=nspp
  A[1:nsite + (i-1)*nsite,1:nsite + (i-1)*nsite]=A_ex[[i]] + A_in[[i]]
  A = as(A,'symmetricMatrix')
  return(A)
} 
  # passed 2019/3/18

mainland_thr = function(dist_mainland,link_mainland,eta,d,int_range_inter="exp"){
	A = 0*dist_mainland
	link_mainland = (as.matrix(link_mainland))
	if(int_range_inter=="arth"){
			A = eta*as.matrix(1/((dist_mainland)^(2+d)))
	}
	else{
		if(int_range_inter=="exp"){
			A = eta*as.matrix(exp(-exp(d)*dist_mainland)) * (link_mainland)
		}
	  else{
	    if(int_range_inter=="nn")
	    A = eta * (link_mainland)
	  }
	}
	return(A)
	# test for 2spp passed 3/18/2019
}


getMRF = function(theta,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra="nn",int_range_inter="exp"){
  	nsite = nrow(envX)
	beta_occu = theta$beta_occu
	eta_intra = theta$eta_intra # intra spp, intra island if apply
	d_intra = theta$d_intra
	spp_mat = theta$spp_mat
	nspp = nrow(spp_mat)
	#nrep = ncol(Z_vec)
	A_in = getintralayerGraph(distM,link_map$intra,eta_intra,d_intra,int_range = int_range_intra,spp_mat)
	eta_inter = theta$eta_inter # assume there is a 
	d_inter = theta$d_inter
	A_ex = getintralayerGraph(distM,link_map$inter,eta_inter,d_inter,int_range = int_range_inter,spp_mat) # graph among islands, if apply, distM should only contain graph 
	A=getfullGraph(A_ex,A_in,spp_mat)
    rm(A_ex,A_in)
	thr = lapply(1:nspp,
	  function(i,envX,beta_occu,dist_mainland,link_mainland,eta_inter,d_inter,int_range_inter){
	    envX %*% beta_occu[1:ncol(envX)+(i-1)*ncol(envX)] + 
		      mainland_thr(dist_mainland,link_mainland,eta_inter[i],d_inter[i],int_range_inter)
	    },envX,beta_occu,dist_mainland,link_mainland,eta_inter,d_inter,int_range_inter)
	thr = Reduce(rbind,thr)
    return(list(A = A,thr = thr))
}


Hamiltonian = function(MRF,Z_vec){
	nrep = ncol(Z_vec)
	Ham = lapply(1:nrep,function(i,Z,J,h){H(J,Z[,i],h)},Z=Z_vec,J=MRF$A,h=( MRF$thr))
  
	Ham = Reduce(rbind,Ham)
	return(Ham) # if we have repeat, just make Z_vec has two cols 
	
}


rIsingOccu_multi = function(MRF,n=1,method = "CFTP",nIter = 100){
	Z = IsingSamplerCpp(n=n,graph = MRF$A,thresholds=MRF$thr, responses = matrix( c(-1L, 1L),2,1),beta = 1,nIter=nIter,exact = (method=="CFTP"),constrain = NA + MRF$thr)
  return(t(Z))
	# test for 2spp case, passed 3/18/2019
}
  

extract_thr = function(i,thr_list){
	nspp = length(thr_list)
	thr = sapply(thr_list,function(thr1,i){t(thr1[i,])},i=i) # thr at site i for all spps, will return a matrix with ncol = nspp, nrow = nperiod
	return(thr)
}


         
write_json.IsingOccu_samples = function(x,path){
  n_sample = nrow(x$Z.mcmc)
  x$theta.mcmc = lapply(x$theta.mcmc,matrix,nrow = n_sample)
  x$Z.mcmc = matrix(x$Z.mcmc,nrow = n_sample)
  class(x) = 'list'
  jsonlite::write_json(x,path)
}
		   
		   
adjacency.matrix = function(m, n = NULL)
{
    if (missing(n))
    {
        A = Matrix(0, m^2, m^2,sparse = T)
        for (i in 1:m^2)
        {
            up = i - m
            down = i + m
            left = i - 1
            right = i + 1
            if (up > 0)
                A[i, up] = 1
            if (down <= m^2)
                A[i, down] = 1
            if (left %% m != 0)
                A[i, left] = 1
            if (i %% m != 0)
                A[i, right] = 1
        }
    }
    else
    {
        A = Matrix(0, m * n, m * n,sparse = T)
        for (i in 1:(m * n))
        {
            up = i - n
            down = i + n
            left = i - 1
            right = i + 1
            if (up > 0)
                A[i, up] = 1
            if (down <= (m * n))
                A[i, down] = 1
            if (left %% n != 0)
                A[i, left] = 1
            if (i %% n != 0)
                A[i, right] = 1
        }
    }
    A
}   
		   
		   
		   
		   