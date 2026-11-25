#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(expm)
library(data.table)
library(mvnfast)
library(rlist)
library(purrr)
library(mvnfast)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("JansenAndRit.cpp")
sourceCpp("WienerIntegrator.cpp")
sourceCpp("EulerMaruyama.cpp")
# sourceCpp("mm_filter_patch.cpp")
# sourceCpp("mm_wie_em_patch.cpp")
source("paramPrepper.R")
source("ekf_likelihood.R")
source("ciTests.R")
#source("profile_Fix.R")

# source("ekf_likelihood_legacy.R")
# source("ciTests_legacy.R")
# source("mm_patch.R")
# source("Jr estimation.R")
# source("Substep patch.R")
# source("mm_wie_em_patch.R")

# TODO simplify g1integral in non-coupling
# TODO OOP based network setup

jr_dim <- 6 # 6 or 10
jr_dim_seq <- 1:jr_dim
jr_dim_seq_half <- 0:(jr_dim/2)
network_simulator <-simulate_network_jr
if(jr_dim==10)network_simulator <- simulate_network_jr_dim10
A <- 3.25#3.25
a <- 100#100
B <- 22#22
b <- 50#20#50
G <- 10#10
g <- 500#500
A. <- A;a. <- a;B. <- B;b. <- b;G. <- G;g. <- g
C. <- 135#70#135
C1 <- C.
C2 <- .8*C.
C3 <- .25*C.
C4 <- .25*C.
C5 <- .3*C.;C6 <- .1*C.;C7 <- .8*C.
obsErr <- 0#sqrt(.5)#sqrt(1.6)#0#sqrt(1.65)
obsStd <- obsErr
Delta. <- .002#1/173.61#.002#1/256#.002#0.002 is (more than) enough for strang #173.61
sigma. <- 2000#2000
tau. <- 1
sigmas. <- c(tau.,sigma.,tau.,tau.,tau.)[jr_dim_seq_half]
Sigma. <- diag(sigmas.)   ##ONLY IF CONSTANT
muFixed. <- c(0,220,0,0,0)[jr_dim_seq_half]#(0,220,0)
mu. <- muFixed.[2]
maxT. <- 4097/173.61#10
N. <- ceiling(maxT./Delta.)
span_val<-5*maxT.
nPop <- 4
couplingStrength <- 500
adjMat <- matrix(1,nPop,nPop)-diag(nPop)
adjMat <- matrix(0,nPop,nPop);adjMat[row(adjMat)==(col(adjMat)-1)] <- 1



# xseq <- seq(0,maxT,length.out=N)
nu0 <- 6#6
nuMax <- 5#5
r <- .56#.56
sigm <- function(x)nuMax/(1+exp(r*(nu0-x)))
# sigm <- function(x)5/(1+exp(.56*(6-x)))
I6 <- diag(jr_dim)


dampingMat <- diag(c(a,a,b,g,b)[jr_dim_seq_half])
Mmat <- rbind(cbind(M0(),Idm()),cbind(-(dampingMat%^%2),-2*dampingMat))

m0 <- diag(rep(0,3))

baseInits <- rep(0,jr_dim)


covMatOU. <- OUcovMat(dampingMat,Delta.)


x1StepMat. <- Matrix(expMt(dampingMat,Delta.),sparse=TRUE)
stepMatHalf. <- expMt(dampingMat,Delta./2)


covMat_WIE <- WIEcovMat(sigma.)



G1fixed. <- c(A*a,A*a,B*b,G*g,B*b)[jr_dim_seq_half]*muFixed.

G1intDelta <- G1fixed.*Delta.
G1intDeltaHalf <- G1fixed.*Delta./2



firstMeanCrossing <- function(x,comp_mean,comp_sd)which(abs(x-comp_mean)<=comp_sd)[1]

findBurnin <- function(x,comp_mean=NULL,comp_sd=NULL){
  if(is.null(comp_mean)) comp_mean <- mean(x)
  if(is.null(comp_sd)) comp_sd <- sd(x)
  minIdx <- firstMeanCrossing(x,comp_mean,comp_sd)
  if(is.na(minIdx)|minIdx==1)return(1:length(x))
  return(-c(1:(minIdx-1)))
}



simulate_jr_double_normals <- function(Delta=.001,maxT=500,stepMat=expMt(Delta/2),C1=135,baseInits=rep(0,6),mu=220,sigmas=c(10,1000,10),removeBurnin=findBurnin){
  N <- as.integer(maxT/Delta)
  if(is.null(removeBurnin))findBurnin <- function(x)1:N
  G1integral <- c(A*a,A*a,B*b)*c(0,mu,0)*Delta
  covarMat <- variableSigmaCovarMat(Delta/2,diag(sigmas))
  
  randoms <- rmvn(2*N-2,rep(0,6),covarMat)
  simul <- strang_integrator_double_normals(N,Delta,stepMat,randoms,A,a,B,b,C1,.8*C1,.25*C1,.25*C1,G1integral,baseInits)
  rm(randoms);gc()
  return(simul[findBurnin(simul)])
}

simulate_jr <- function(Delta=Delta.,maxT=maxT.,stepMat=expMt(dampingMat,Delta),C1=C.,sigma=sigma.,
                        baseInits=rep(0,jr_dim),
                        covarMat=variableSigmaCovarMat(Delta,diag(c(tau.,sigma,tau.)),dampingMat),
                        g1dh=G1fixed./2*Delta,obsSd=obsStd,removeBurnin=NULL,fullPath=FALSE){
  N <- as.integer(maxT/Delta)+1
  
  if(is.null(removeBurnin))findBurnin <- function(x)1:N
  
  if(sum(covarMat)==0){randoms <- matrix(0,N-1,6)
  }else{
    rands <- diag(covarMat)!=0
    if(!all(rands)){
      randT <- rmvn(N-1,rep(0,sum(rands)),covarMat[rands,rands])
      randoms <- matrix(0,N-1,6)
      randoms[,rands] <- randT;rm(randT);gc()
    }else{
      randoms <- rmvn(N-1,rep(0,jr_dim),covarMat)
    }
    }
  if(!fullPath){
    if(jr_dim==6)simul <- strang_integrator(N,Delta,stepMat,randoms,A,a,B,b,C1,g1dh,baseInits)
    else simul <- strang_integrator_dim10(N,Delta,stepMat,randoms,A,a,B,b,G,g,C1,g1dh,baseInits)
    rm(randoms);gc()
    curObs <- simul[findBurnin(simul)]
    return(curObs+rnorm(length(curObs),0,obsSd))
  }
  simul <- strang_integrator_full_path(N,Delta,stepMat,randoms,A,a,B,b,C1,g1dh,baseInits)
  rm(randoms);gc()
  obs <- simul[,2]-simul[,3]
  indices <- findBurnin(obs)[-1]
  return(simul[indices,]+rnorm(length(indices),0,obsSd))
}


simulate_em_baseline <- function(subFac = 100,Delta=Delta.,maxT=maxT.,stepMat=expMt(dampingMat,Delta/2),sigma=sigma.,
                                 baseInits=rep(0,jr_dim),
                                 covarMat=Delta*diag(c(tau.,sigma,tau.)^2),obsSd=obsStd,type="em"){
  Delta <- Delta/subFac
  N <- as.integer(maxT/Delta)+1
  

  if(type=="em"){
    randoms <- rmvn(N-1,rep(0,3),covarMat)
    curObs <- jr_em(N,baseInits,Delta,A,a,B,b,C1,C2,C3,C4,muFixed.,randoms)
  }else{
    if(type=="wi"){
      randoms <- rmvn(N-1,rep(0,6),diag(6))%*%expm::sqrtm(stepMat%*%diag(c(rep(0,3),tau.,sigma,tau.)^2)%*%t(stepMat)*Delta)
      curObs <- wiener_simulate(randoms,stepMat,N,Delta,mu.,A,a,B,b,C1,nuMax,nu0,r)
  } else{
      randoms <- rmvn(N-1,rep(0,6),variableSigmaCovarMat(Delta,diag(c(tau.,sigma,tau.)),dampingMat))
      curObs <- strang_integrator(N,Delta,expMt(dampingMat,Delta),randoms,A,a,B,b,C1,G1fixed./2*Delta,baseInits)
    }
  }
  rm(randoms);gc()
  curObs <- curObs[seq(1,N,subFac)]
  return(curObs[-1]+rnorm(length(curObs)-1,0,obsSd))
}



simulate_comp <- function(){
  
  randoms <- rmvn(N.-1,rep(0,jr_dim),diag(6))
  
  # rand_ou <- t(sapply(1:(N.-1),function(i)expm::sqrtm(covMatOU.)%*%randoms[i,]))
  rand_ou <- randoms%*%expm::sqrtm(covMatOU.)
  
  rand_wie <- (randoms[,4:6])%*%(sqrt(Delta.)*Sigma.)
  x <- wiener_simulate(rand_wie,stepMatHalf.,N.,Delta.,mu.,A,a,B,b,C1,nuMax,nu0,r)
  y <- strang_integrator(N.,Delta.,as.matrix(x1StepMat.),rand_ou,A,a,B,b,C1,G1fixed./2*Delta.,baseInits)
  plot(x,type="l")
  lines(y,col="red")
  
  plot(x[1:1000],type="l")
  lines(y[1:1000],col="red")
  
  plot(x[500:510],type="l")
  lines(y[500:510],col="red")
  
  
  jr_inference_with_ci(c("C1","mu","sigma"),x)
  
  
  xf <- wiener_simulate_fullpath(rand_wie,stepMatHalf.,N.,Delta.,mu.,A,a,B,b,C1,nuMax,nu0,r)
  yf <- strang_integrator_full_path(N.,Delta.,as.matrix(x1StepMat.),rand_ou,A,a,B,b,C1,G1fixed./2*Delta.,baseInits)
  plot(xf[1:1000,5],type="l")
  lines(yf[1:1000,5],col="red")
}

simulate_variable <- function(...){
  kwargs <- list(...)
  vals <- modifyList(.defaults, kwargs)
  
  # actually how to make stuff hierarchical?
  if(TRUE)stop("Nah")
  N <- as.integer(vals$maxT/vals$Delta)
  randoms <- rmvn(N-1,rep(0,6),vals$covarMat)
  # if(jr_dim==6)
  simul <- strang_integrator(N,vals$Delta,vals$stepMat,randoms,vals$A,vals$a,vals$B,vals$b,vals$C1,vals$g1dh,vals$baseInits)
  # else simul <- strang_integrator_dim10(N,Delta,stepMat,randoms,A,a,B,b,G,g,C1,g1dh,baseInits)
  rm(randoms);gc()
  curObs <- simul[findBurnin(simul)]
  return(curObs+rnorm(length(curObs),0,obsSd))
}




simulate_coupling <- function(Delta=Delta.,maxT=maxT.,stepMat=expMt(dampingMat,Delta),C1=C.,baseInits=rep(0,jr_dim),covarMat=covMatOU.,couplingStrength=couplingStrength,obVar=obsStd^2,removeBurnin=NULL,fullPath=FALSE){
  # Delta=Delta.;maxT=maxT.;stepMat=expMt(dampingMat,Delta);C1=C.;baseInits=rep(0,jr_dim);covarMat=covMatOU.;maxT <- 2
  # N <- as.integer(maxT/Delta)
  if(!is.null(covarMat))stop("needs fixing in general?covMats ftw? but covar mat in params")
  randoms <- array(0,dim=c(N-1,jr_dim,nPop))
  for(pop in 1:nPop)randoms[,,pop] <- rmvn(N-1,rep(0,6),covMats[,,pop])
  simul <- simulate_network_jr(N,Delta,stepMats,randoms,As,as,Bs,bs,Cs,
                               g1integrals,nPop,500,adjMat,matrix(0,4,jr_dim))
  
  indices <- findBurnin(simul[,1])
  return(simul[indices,]+rmvn(1,rep(0,nPop),diag(rep(obVar,nPop))))
}



simulate_lie_trotter <-  function(Delta=Delta.,maxT=maxT.,stepMat=expMt(dampingMat,Delta),C1=C.,baseInits=rep(0,6),covarMat=covMatOU.,removeBurnin=NULL){
  N <- as.integer(maxT/Delta)
  if(is.null(removeBurnin))findBurnin <- function(x)1:N
  
  if(sum(covarMat)==0){randoms <- matrix(0,N-1,6)
  }else{ randoms <- rmvn(N-1,rep(0,6),covarMat)}
  simul <- lie_trotter_full_path(N,Delta,stepMat,randoms,A,a,B,b,C1,.8*C1,.25*C1,.25*C1,G1intDelta,baseInits)
  rm(randoms);gc()
  obs <- simul[,2]-simul[,3]
  indices <- findBurnin(obs)
  return(simul[indices,])
}




H <- c(0,1,-1,0,0,0)
t1 <- function(x,q=0) x1StepMat.%*%x+q
t2 <- function(x,stepSize=Delta./2) x+stepSize*c(rep(0,3),A*a*c(sigm(x[2]-x[3]),C2*sigm(C1*x[1])+mu.),B*b*C4*sigm(C3*x[1]))
t2Deriv <- function(x){
  fullMat <- diag(6)
  stuff <- A*a*sigmDeriv(x[2]-x[3])*Delta./2
  fullMat[4,2] <- stuff;fullMat[4,3] <- -stuff
  fullMat[5,1] <- A*a*C1*C2*sigmDeriv(C1*x[1])*Delta./2
  fullMat[6,1] <- B*b*C3*C4*sigmDeriv(C3*x[1])*Delta./2
  fullMat
}

fFull <- function(x,q=0)t2(t1(t2(x),q))
F_jac <- function(x)t2Deriv(x1StepMat.%*%t2(x))%*%x1StepMat.%*%t2Deriv(x)

L_jac <- function(m){ 
  fullMat <- diag(6)
  temp <- Delta./2*A*a*sigmDeriv((x1StepMat.[2,]%*%t2(m))[1]-(x1StepMat.[3,]%*%t2(m))[1])
  fullMat[4,2] <- temp;fullMat[4,3] <- -temp
  fullMat[5,1] <-  Delta./8*A*a*C1*C2*nuMax*r*(pracma::sech(r/2*(nu0-C1*x1StepMat.[1,]%*%c(m[1],m[2],m[3],A*a*Delta.*nuMax/(2*(1+exp(r*(nu0-m[2]+m[3]))))+m[4],
                                                                                     Delta./2*A*a*(mu.+(C2*nuMax)/(1+exp(r*(nu0-C1*m[1]))))+m[5],
                                                                                     B*b*C4*Delta.*nuMax/(2*(1+exp(r*(nu0-C3*m[1]))))+m[6])))^2)
  fullMat[6,1] <- Delta./8*B*b*C3*C4*nuMax*r*(pracma::sech(r/2*(nu0-C3*x1StepMat.[1,]%*%c(m[1],m[2],m[3],A*a*Delta.*nuMax/(2*(1+exp(r*(nu0-m[2]+m[3]))))+m[4],
                                                                                    Delta./2*A*a*(mu.+(C2*nuMax)/(1+exp(r*(nu0-C1*m[1]))))+m[5],
                                                                                    Delta.*B*b*C4*nuMax/(2*(1+exp(r*(nu0-C3*m[1]))))+m[6])))^2)
  return(fullMat)
}

network_deriv_x <- function(x){
  
}


# ###############CHECKING JACOBIANS
# xs <- rep(1,6);fFull_q <- function(q)fFull(xs,q)
# any(abs(numDeriv::jacobian(fFull,xs)-F_jac(xs))>1e-8)
# # any(abs(numDeriv::jacobian(fFull,xs)-full_step_derivative(xs,Delta./2,x1StepMat.,mu.,C1,C2,C3,C4,A,a,B,b))>1e-8)
# any(abs(numDeriv::jacobian(fFull_q,rep(0,6))-L_jac(xs))>1e-8)
# # any(abs(numDeriv::jacobian(fFull_q,rep(0,6))-random_derivatives(xs,Delta./2,x1StepMat.,mu.,C1,C2,C3,C4,A,a,B,b))>1e-8)



param_strings <- c("observation","stepMat","covMat","obsStd","nObs","Delta","mu","A","a","B","b","C1","sigma","nuMax","nu0","r")
usual_param_value <- list(NA,as.matrix(x1StepMat.),covMatOU.,obsErr,NA,Delta.,mu.,A,a,B,b,C1,sigma.,nuMax,nu0,r)
# initialGuesses <- c(obsStd=2,mu=200,A=5,a=80,B=30,b=45,C1=140,sigma=sigma.*.9,nuMax=4,nu0=5,r=1/2)
initialGuesses <- c(obsStd=.5,mu=200,A=5,a=80,B=30,b=45,C1=140,sigma=sigma.*1.1,nuMax=4,nu0=5,r=1/2)
if(jr_dim==10){
  param_strings <- c(param_strings,"G","g")
  usual_param_value <- c(usual_param_value,G,g)
  initialGuesses <- c(initialGuesses,14,600)
}
usuals <- usual_param_value
names(usuals) <- param_strings

param_strings     <- c(param_strings, "burnin")
usual_param_value <- c(usual_param_value, list(500L))
