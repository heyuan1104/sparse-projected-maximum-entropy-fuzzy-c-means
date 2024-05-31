My.rkm=function(X,nclus,ndim,iterations,nstart,seed){
  X=as.matrix(X)
  #number of cluster
  k <- nclus
  
  #rank
  d <- ndim
  
  #the number of subjects
  n <- nrow(X)
  
  #the number of variables
  p <- ncol(X)
  
  #The number of iteration
  ite.num <- iterations
  
  #initial parameters of component loadings
  A <- eigen(t(X) %*% X)$vectors[,1:d]
  
  #initial indicator matrix
  U <-RKM_nstart(X,k,d,ite.num,nstart,seed)$U
  
  #initial G
  G <- ginv( t(U) %*% U ) %*% t(U) %*% X %*% A
  
  #threshold value
  epsilon <- 1e-8
  
  #matrix storing values of the objective function
  OB_mat <- as.data.frame(matrix(NA,ite.num,2))
  colnames(OB_mat) <- c("Update_U","Update_A")
  
  for(ite in 1:ite.num){
    #Update U
    for(i in 1:n){
      dummy <- rep(NA,k)
      for(l in 1:k){
        U[i,] <- 0
        U[i,l] <- 1
        dummy[l] <- RKM_OB(X,U,G,A)
      }
      U[i,] <- 0
      U[i,which(min(dummy) == dummy)[1]] <- 1
      
    }
    
    #Update G
    G <- ginv( t(U) %*% U ) %*% t(U) %*% X %*% A
    
    #Calculating the value
    OB_mat$Update_U[ite] <- RKM_OB(X,U,G,A)
    
    #Update A
    D_tmp <- t(X) %*% U %*% ginv(t(U) %*% U) %*% t(U) %*% X
    A <- eigen(D_tmp)$vectors[,1:d]
    
    #Update G
    G <- ginv( t(U) %*% U ) %*% t(U) %*% X %*% A
    
    #Calculating the value
    OB_mat$Update_A[ite] <- RKM_OB(X,U,G,A)
    if(ite > 1){
      if(OB_mat$Update_A[ite-1]- OB_mat$Update_A[ite] < epsilon) break
    }
    
  }
  
  #各クラスターのサイズ
  cluster.vec=rep(0,length=ncol(U))
  for(i in 1:nrow(U)){
    for(j in 1:ncol(U)){
      if(U[i,j]==1){
        cluster.vec[j]=cluster.vec[j]+1
      }
    }
  }
  
  #クラスター情報を得る
  cluster.vec1=rep(0,length=nrow(X))
  for(i in 1:nrow(U)){
    for(j in 1:ncol(U)){
      if(U[i,j]==1){
        cluster.vec1[i]=j
      }
    }
  }
  
  return(list("cluster size"=cluster.vec,"G1"=G, "G2"=G%*%t(A), "variable score"=X%*%A, cluster=cluster.vec1,'U'=U, 'OB'=OB_mat))
  
}

#目的関数
RKM_OB=function(X,U,G,A){
  OB <- sum((X - U %*% G %*% t(A))^2)
  return(OB)
}

#for nstart
RKM_nstart=function(X,k,d,ite.num,nstart,seed){
  
  #the number of subjects
  n <- nrow(X)
  
  #the number of variables
  p <- ncol(X)
  
  #threshold value
  epsilon <- 1e-10
  
  #matrix storing values of the objective function
  OB_mat <- as.data.frame(matrix(NA,ite.num,2))
  colnames(OB_mat) <- c("Update_U","Update_A")
  
  loss_best=Inf
  OB_value=c()
  
  for(count in 1:nstart){
    
    #initial indicator matrix
    U <- matrix(0,n,k)
    
    set.seed(seed+count)
    for(i in 1:n){
      U[i,sample(1:k,1)] <- 1
    }
    U_initial=U
    
    #initial parameters of component loadings
    A <- eigen(t(X) %*% X)$vectors[,1:d]
    
    #initial G
    G <- ginv( t(U) %*% U ) %*% t(U) %*% X %*% A
    
    for(ite in 1:ite.num){
      #Update U
      for(i in 1:n){
        dummy <- rep(NA,k)
        for(l in 1:k){
          U[i,] <- 0
          U[i,l] <- 1
          dummy[l] <- RKM_OB(X,U,G,A)
        }
        U[i,] <- 0
        U[i,which(min(dummy) == dummy)[1]] <- 1
        
      }
      
      #Update G
      G <- ginv( t(U) %*% U ) %*% t(U) %*% X %*% A
      
      #Calculating the value
      OB_mat$Update_U[ite] <- RKM_OB(X,U,G,A)
      
      #Update A
      D_tmp <- t(X) %*% U %*% ginv(t(U) %*% U) %*% t(U) %*% X
      A <- eigen(D_tmp)$vectors[,1:d]
      
      #Update G
      G <- ginv( t(U) %*% U ) %*% t(U) %*% X %*% A
      
      #Calculating the value
      OB_mat$Update_A[ite] <- RKM_OB(X,U,G,A)
      if(ite > 1){
        if(OB_mat$Update_A[ite-1]- OB_mat$Update_A[ite] < epsilon) {
          of=OB_mat$Update_A[ite]
          break
        }
      }
    }
    OB_value=append(OB_value,of)
    if(of<loss_best){
      loss_best=of
      #G_initial_best=G_initial
      U_initial_best=U_initial
    }
  }
  
  return(list('U'=U_initial_best,'OBvalue'=OB_value))
}

