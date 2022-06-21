AmatrixEst <- function(X,theta,w0,w1, doublezero){
  
  n <- nrow(X)
  Qvec <- expit(c(cbind(1,X)%*%theta))
  eta <- etafunct(w0, w1, Qvec, doublezero)
  ans <- (1/n)*t(cbind(1,X))%*%diag(c(Qvec*(1-Qvec)*eta))%*%cbind(1,X)
  return(ans)
  
}

AmatrixEst.MLE <- function(X,theta){
  
  n <- nrow(X)
  Qvec <- expit(c(cbind(1,X)%*%theta))
  ans <- (1/n)*t(cbind(1,X))%*%diag(c(Qvec*(1-Qvec)))%*%cbind(1,X)
   
  return(ans)
  
}

BmatrixEst <- function(y, X,theta,w0,w1, doublezero){
  
  n <- nrow(X)
  w <- y*w1 + (1-y)*w0
  Qvec <- expit(c(cbind(1,X)%*%theta))
  resid <- y - Qvec - Cfunct(t = w0, u = w1, v = Qvec,doublezero=doublezero)
  ans <- (1/n)*t(cbind(1,X)) %*% diag(c( (w^2)*(resid^2) )) %*% cbind(1,X)
  return(ans)
  
}

BmatrixEst.MLE <- function(X,y,theta){
  
  n <- nrow(X)
  Qvec <- expit(c(cbind(1,X)%*%theta))
  resid <- y - Qvec
  ans <- (1/n)*t(cbind(1,X))%*%diag(c(resid^2))%*%cbind(1,X)
  return(ans)
  
}
# BmatrixEst2 <- function(y, X,theta,w0,w1, doublezero){
#   
#   n <- nrow(X)
#   w <- y*w1 + (1-y)*w0
#   Qvec <- expit(c(cbind(1,X)%*%theta))
#   resid <- y - Qvec - Cfunct(t = w0, u = w1, v = Qvec,doublezero=doublezero)
#   #ans <- (1/n)*t(cbind(1,X)) %*% diag(c( (w^2)*(resid^2) )) %*% cbind(1,X)
#   Multiplier.Vec <- (w^2)*(resid^2)
#   ans <- matrix(0,nrow=length(theta),ncol=length(theta))
#   for(i in 1:n){
#     
#     zrow <- matrix(c(1,X[i,]),ncol=1)
#     ans <- ans + Multiplier.Vec[i]*(zrow %*% t(zrow))
#     
#   }
#   ans <- ans/n
#   
#   return(ans)
#   
# }

etafunct <- function(t,u,v,doublezero){
  
  ans <- (1-v)*u + v*t - (u-t)*Cfunct(t,u,v,doublezero)
  return(ans)
  
}


Cfunct <- function(t,u,v, doublezero){
  
  ret.vec <- v*(1-v)*(u-t)/((1-v)*t+v*u)
  if(length(doublezero)>0){
    
    ret.vec[doublezero] <- 0
    
  }
  
  return(ret.vec)
  
}

dCdv <- function(t,u,v, doublezero){
  
  retvec <- (u-t)*( t*(1-v)^2 - u*v^2 )/(( (1-v)*t + v*u )^2)
  if(length(doublezero)>0){
    
    retvec[doublezero] <- 0
    
  }
  return(retvec)
}

expit <- function(z){
  
  return(1/(1+exp(-z)))
  
}

Carroll.w <- function(d,h,p){
  
  v <- sqrt(d/(p-1))
  psi <- v*(1-(v/h)^2)^3*(abs(v)<=h)
  weights.out <- psi/v
  return(weights.out)
  
}

MainFunction <- function(w0, w1, y, X, theta, doublezero){
  
  w <- y*w1 + (1-y)*w0
  Qvec <- expit(c(cbind(1,X)%*%theta))
  resid <- y - Qvec - Cfunct(t = w0, u = w1, v = Qvec, doublezero=doublezero)
  ModMallowsX <- apply(cbind(1,X),2,function(s){s*c(w*resid)})
  ans <- colSums(ModMallowsX)
  return(ans)
  
}

DerivativeFunction <- function(w0, w1, y, X, theta, doublezero){
  
  w <- y*w1 + (1-y)*w0
  Qvec <- expit(c(cbind(1,X)%*%theta))
  resid <- y - Qvec - Cfunct(t = w0, u = w1, v = Qvec,doublezero=doublezero)
  DerivF.wt <- w*(1+dCdv(t=w0,u=w1, v = Qvec,doublezero=doublezero))*(-1)*Qvec*(1-Qvec)
  ans <- t(cbind(1,X)) %*% diag(DerivF.wt) %*% cbind(1,X)
  return(ans)
}


CalcModMallows <- function(y,X,theta.init,h=8, force.weights1 = NA){
  
   n <- length(y)
   theta.curr <- theta.init
   err.obj <- 10
   err.tol <- 10^(-6)
   mu1 <- covRob(data=X[y==1,], estim = "mcd", quan = 0.95, ntrial = 1000)$center
   mu0 <- covRob(data=X[y==0,], estim = "mcd", quan = 0.9, ntrial = 1000)$center
   Sigma1 <- covRob(data=X[y==1,], estim = "mcd", quan = 0.95, ntrial = 1000)$cov
   #Sigma1 <- var(X[y==1,])
   #Sigma0 <- var(X[y==0,])
   Sigma0 <- covRob(data=X[y==0,], estim = "mcd", quan = 0.9, ntrial = 1000)$cov
   d1 <- mahalanobis(x=X, center=mu1, cov = Sigma1)
   d0 <- mahalanobis(x=X, center=mu0, cov = Sigma0)
   w0 <- Carroll.w(d = d0, h = h, p = 3)
   w1 <- Carroll.w(d = d1, h = h, p = 3)
   #w1 <- rep(1,n)
   #w0 <- rep(1,n)
   Double.zero.id <- which( (w0==0)&(w1==0) )
   #browser()
   #w1 <- Pregibon.w(v=mahalanobis(x=X, center=mu1, cov = Sigma1), h=h)
   #w0 <- Pregibon.w(v=mahalanobis(x=X, center=mu0, cov = Sigma0), h=h)
   if(length(force.weights1)>1){
     
     w1 <- force.weights1
     
   }
    #browser()
   w <- y*w1 + (1-y)*w0
   #browser()
   STOP.iter <- FALSE
   while(!STOP.iter){
     
     #browser()
     V <- DerivativeFunction(w0 = w0, w1 = w1, y = y, X = X, theta = theta.curr, doublezero=Double.zero.id)
     theta.new <- theta.curr - solve(V) %*% MainFunction(w0=w0, w1=w1, y=y, X=X, theta = theta.curr, doublezero=Double.zero.id)
     err.obj <- sum((theta.new-theta.curr)^2)
     
     if(err.obj < err.tol){
       
       STOP.iter <- TRUE
       
     }
     if(norm(V)<= 10^(-5)){
       
       STOP.iter <- TRUE
       
     }
     theta.curr <- theta.new
     
   }
   
   A <- AmatrixEst(X = X, theta = theta.curr, w0 = w0, w1 = w1, doublezero = Double.zero.id)
   B <- BmatrixEst(y=y, X=X, theta = theta.curr, w0 = w0, w1 = w1, doublezero = Double.zero.id)
   Ainv <- solve(A)
   SigmaEst <- (1/n)*Ainv%*%B%*%t(Ainv)
  #browser()
   return(list(theta.curr=theta.curr, w=w, SigmaEst = SigmaEst))
  
}