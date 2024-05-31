My.spefcm.rcpp=function(X, cluster, dimension, lambda.e, lambda.l, iterations, nstart,seed){
  
  X=as.matrix(X)
  #iterations
  ite.num=iterations
  
  #define the number of variables
  n=nrow(X)
  
  #define the number of subjects
  p=ncol(X)
  
  #define the number of clusters and dimensions
  k=cluster
  d=dimension
  
  #define the threshold
  epsilon=1e-8
  
  #initial indicator matrix
  U=t(RKM_nstart(X,k,d,ite.num,nstart,seed)$U)
  
  #initial V
  Y=U %*% X
  V=matrix(0,nrow = k,ncol = p)
  for(i in 1:nrow(Y)){
    V[i,]=Y[i,]/sum(U[i,])
  }
  
  #initial A
  A <- eigen(t(X) %*% X)$vectors[,1:d]
  
  #initial B
  B <- matrix(rnorm(p*d,0,1),p,d)
  for(o2 in 1:d){
    B[,o2] <- B[,o2]/sqrt(sum(B[,o2]^2))
  }
  B <- A
  
  #matrix storing values of the objective function
  OB_mat <- as.data.frame(matrix(NA,ite.num,3))
  colnames(OB_mat) <- c("Update_U","Update_B","Update_A")
  
  for(ite in 1:ite.num){
    ##update U
    for(i in 1:k){
      for(j in 1:n){
        a=c()
        for(t in 1:k){
          a=append(a,exp(-lambda.e*sum((X[j,]-A %*% t(B) %*% V[t,])^2)))
        }
        U[i,j]=exp(-lambda.e*sum((X[j,]-A %*% t(B) %*% V[i,])^2))/sum(a)
      }
    }
    
    Y=U %*% X
    for(i in 1:nrow(Y)){
      V[i,]=Y[i,]/sum(U[i,])
    }
    #U更新後の目的関数
    OB_mat$Update_U[ite]=spefcm_OB(X,U,V,A,B,k,lambda.e,lambda.l)
    
    ##update B
    G.vec1=c()
    Z.mat1=matrix(0,(n*d*k),(p*d))
    for(i in 1:k){
      P=matrix(0,nrow = n, ncol = p)
      for(j in 1:n){
        P[j,]=V[i,]
      }
      phi=diag(U[i,])
      X1=sqrt(phi) %*% X
      X2=sqrt(phi) %*% P
      
      G <- X1 %*% A
      G.vec <- as.numeric(G)
      G.vec1 <- append(G.vec1, G.vec)
      
      Z.mat <- matrix(0,(n*d),(p*d))
      for(o in 1:d){
        Z.mat[((o-1)*n+1):(o*n),((o-1)*p+1):(o*p)] <- X2
      }
      Z.mat1[((i-1)*(n*d)+1):(i*n*d),] <- Z.mat
    }
    #coefficient vector
    B.vec <- as.numeric(B)
    
    for(j in 1:(p*d)){
      #j <- 1
      Rj <- G.vec1 - Z.mat1[,-j] %*% B.vec[-j]
      B.vec.dummy <- sum(Z.mat1[,j] * Rj)
      B.vec=Lasso(B.vec, Z.mat1, lambda.l, B.vec.dummy, j)
      B <- matrix(B.vec,p,d)
    }
    OB_mat$Update_B[ite]=spefcm_OB(X,U,V,A,B,k,lambda.e,lambda.l)
    
    ## update A
    X3=matrix(0,nrow = n*k,ncol = p)
    X4=matrix(0,nrow = n*k,ncol = p)
    for(i in 1:k){
      P=matrix(0,nrow = n, ncol = p)
      for(j in 1:n){
        P[j,]=V[i,]
      }
      phi=diag(U[i,])
      X1=sqrt(phi) %*% X
      X2=sqrt(phi) %*% P
      
      X3[((i-1)*n+1):(i*n),]=X1
      X4[((i-1)*n+1):(i*n),]=X2
    }
    ans=svd(t(X3) %*% X4 %*% B)
    A=ans$u %*% t(ans$v)
    OB_mat$Update_A[ite] <- spefcm_OB(X,U,V,A,B,k,lambda.e,lambda.l)
    
    
    if(ite > 1){
      if((OB_mat$Update_A[ite-1]- OB_mat$Update_A[ite]) < epsilon) break
    }
  }
  
  U2=matrix(0,ncol(U),nrow(U))
  for(i in 1:nrow(U2)){
    U2[i,which.max(t(U)[i,])]=1
  }
  #クラスター情報を得る
  cluster.vec1=rep(0,length=nrow(X))
  for(i in 1:nrow(U2)){
    for(j in 1:ncol(U2)){
      if(U2[i,j]==1){
        cluster.vec1[i]=j
      }
    }
  }
  
  centroid=matrix(NA,k,d)
  for(i in 1:k){
    centroid[i,]=t(B) %*% V[i,]
  }
  
  return(list('U'=t(U), 'U2'=U2, 'OB'=OB_mat, 'X'=X %*% B, 'B'=B, 'G'=V %*% B %*% t(A), 'A'=A, 'cluster'=cluster.vec1, 'centroid'=centroid))
}

spefcm_OB=function(X, U, V, A, B, k,lambda.e, lambda.l){
  n=nrow(X)
  S=matrix(0,nrow = n, ncol = k)
  for(i in 1:k){
    for(l in 1:n){
      S[l,i]=sum((X[l,]-A %*% t(B) %*% V[i,])^2)
    }
  }
  L1=sum(U * t(S))
  L2=lambda.e^(-1)*sum(U * log(U))
  L3=lambda.l*sum(abs(B))
  L=L1+L2+L3
  return(L)
}