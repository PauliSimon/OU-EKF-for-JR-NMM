sigm <- function(x)nuMax/(1+exp(r*(nu0-x)))
sigmDeriv <- function(x){
  if(abs(x-nu0)>500)return(0)
  exp(r*(nu0+x))*nuMax*r/((exp(nu0*r)+exp(r*x))^2)
}
M0 <- function(m=jr_dim/2)matrix(0,m,m)
Idm <- function(m=jr_dim/2)diag(m)
theta <- function(dampingMat,t)diag(exp(-diag(dampingMat)*t))%*%(Idm()+dampingMat*t)
thetaStrich <- function(dampingMat,t)-(dampingMat%^%2)%*%diag(exp(-diag(dampingMat)*t))*t
kappa <- function(dampingMat,t)diag(exp(-diag(dampingMat)*t))*t
kappaStrich <- function(dampingMat,t)diag(exp(-diag(dampingMat)*t))%*%(Idm()-dampingMat*t)


variableSigmaCovarMat <- function(t,Sigma,dampingMat){
  rbind(cbind(1/4*(solve(dampingMat)%^%(3))%*%(Sigma%^%2)%*%(Idm()+kappa(dampingMat,t)%*%thetaStrich(dampingMat,t)-(theta(dampingMat,t)%^%2)),
              1/2*(Sigma%^%2)%*%(kappa(dampingMat,t)%^%2)),
        cbind(1/2*(Sigma%^%2)%*%(kappa(dampingMat,t)%^%2),
              1/4*(solve(dampingMat))%*%(Sigma%^%2)%*%(Idm()+kappa(dampingMat,t)%*%thetaStrich(dampingMat,t)-(kappaStrich(dampingMat,t)%^%2)))
  )
}




OUcovMat <- function(dampingMat,t){
  res <- rbind(cbind(1/4*(solve(dampingMat)%^%(3))%*%(Sigma.%^%2)%*%(Idm()+kappa(dampingMat,t)%*%thetaStrich(dampingMat,t)-(theta(dampingMat,t)%^%2)),
                                                1/2*(Sigma.%^%2)%*%(kappa(dampingMat,t)%^%2)),
                                          cbind(1/2*(Sigma.%^%2)%*%(kappa(dampingMat,t)%^%2),
                                                1/4*(solve(dampingMat))%*%(Sigma.%^%2)%*%(Idm()+kappa(dampingMat,t)%*%thetaStrich(dampingMat,t)-(kappaStrich(dampingMat,t)%^%2)))
  )
  return(res)
  return(Matrix(res,sparse))
  
}


WIEcovMat <- function(sigma){
  stepMatHalf.%*%diag(c(rep(0,3),tau.,sigma,tau.)^2)%*%t(stepMatHalf.)*Delta.
}

expMt <- function(dampingMat,t){
  return(rbind(cbind(theta(dampingMat,t),kappa(dampingMat,t)),
               cbind(thetaStrich(dampingMat,t),kappaStrich(dampingMat,t))))
}



JrParams <- function(A=A.,a=a.,B=B.,b=b.,G=G.,g=g.,C=C1,obsErr=obsStd,Delta=Delta.,
                     sigma=sigma.,mu=mu.,tau=tau.,maxT=maxT.){
  sigmas <- c(tau,sigma,tau,tau,tau)[jr_dim_seq_half]
  Sigma <- diag(sigmas)   ##ONLY IF CONSTANT
  I6 <- diag(jr_dim)
  dampingMat <- diag(c(a,a,b,g,b)[jr_dim_seq_half])
  Mmat <- rbind(cbind(M0(),Idm()),cbind(-(dampingMat%^%2),-2*dampingMat))
  
  covMat_WIE <- Delta*(Sigma%^%2) 
  
  baseInits <- rep(0,jr_dim)
  
  
  # covMatOU <- OUcovMat(dampingMat,Delta) 
  covMatOU <- variableSigmaCovarMat(Delta,Sigma,dampingMat)
  
  
  
  x1StepMat <- expMt(dampingMat,Delta)
  
  # muFixed <- c(0,mu,0,0,0)[jr_dim_seq_half]#(0,220,0)
  # G1fixed <- c(A*a,A*a,B*b,G*g,B*b)[jr_dim_seq_half]*muFixed
  # 
  # G1intDeltaHalf <- G1fixed*Delta/2
  
  g1int <- A*a*mu*Delta/2
  
  if(jr_dim==10){
    
    return(list(stepMat=x1StepMat,covMat = covMatOU,A=A,a=a,B=B,b=b,C=C,g1int=g1int,G=G,g=g,baseInits=rep(0,jr_dim),
                Delta=Delta,N=ceiling(maxT/Delta)))
  }
  
  return(list(stepMat=x1StepMat,covMat = covMatOU,A=A,a=a,B=B,b=b,C=C,g1int=g1int,baseInits=rep(0,jr_dim),
              Delta=Delta,N=ceiling(maxT/Delta)))
}

NetworkJrParams <- function(...){
  allPars <- list(...)
  nPop <- length(allPars)
  stepMats <- array(0,c(jr_dim,jr_dim,nPop))
  covMats <- array(0,c(jr_dim,jr_dim,nPop))
  As <- NULL;as=NULL;Bs=NULL;bs=NULL;Cs=NULL;g1ints=NULL;baseInits=matrix(0,nPop,jr_dim);Gs=NULL;gs=NULL
  compDelta <- allPars[[1]]$Delta;compN <- allPars[[1]]$N
  for(curPop in 1:nPop){
    jrpar = allPars[[curPop]]
    stepMats[,,curPop] <- jrpar$stepMat
    covMats[,,curPop] <- jrpar$covMat
    As <- c(As,jrpar$A);as <- c(as,jrpar$a);Bs <- c(Bs,jrpar$B);bs <- c(bs,jrpar$b)
    Cs <- c(Cs,jrpar$C);g1ints=c(g1ints,jrpar$g1int)
    baseInits[curPop,] <- jrpar$baseInits
    if(jr_dim==10){
      Gs <- c(Gs,jrpar$G);gs <- c(gs,jrpar$g)
    }
    if(jrpar$Delta!=compDelta|jrpar$N!=compN) stop("Please input only objects with same Delta and N")
  }
  
  if(jr_dim==10){
    return(list(stepMats=stepMats,covMats=covMats,As=As,as=as,Bs=Bs,bs=bs,Cs=Cs,Gs=Gs,gs=gs,
                g1Integrals=g1ints,initialStates=baseInits,Delta=compDelta,N=compN))
  }
  
  return(list(stepMats=stepMats,covMats=covMats,As=As,as=as,Bs=Bs,bs=bs,Cs=Cs,
              g1Integrals=g1ints,initialStates=baseInits,Delta=compDelta,N=compN))
}

Simulation_Prep <- function(netPars,adjMat,coupling_strength){
  nPop <- length(netPars$As)
  randoms <- array(0,dim=c(netPars$N-1,jr_dim,nPop))
  for(pop in 1:nPop)randoms[,,pop] <- rmvn(netPars$N-1,rep(0,jr_dim),netPars$covMats[,,pop])
  netPars[["covMats"]] <- NULL;netPars$nPop <- nPop
  netPars$coupling_strength <- coupling_strength
  netPars$adjacencyMat <- adjMat;netPars$strangRandomness <- randoms
  return(netPars)
}
