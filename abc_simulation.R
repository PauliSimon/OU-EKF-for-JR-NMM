true_theta <- c(C1,sigma.,mu.)
data_obs   <- simulate_jr()
n_iter <- num.cores*1

simulate <- function(the) {
  simulate_jr(C1=the[1],covarMat=variableSigmaCovarMat(Delta.,diag(c(tau.,the[2],tau.)),dampingMat),
              g1dh= c(0,the[3]*A*a*Delta./2,0))
}


summary_extraction <- function(obs,bw0){

  den <- density(obs,
                 bw=bw0,
                 from=-40,
                 to=40,
                 n=1024)$y

  sp <- spectrum(obs,
                 plot=FALSE,
                 log="no",
                 span=span_val)$spec
  list(den=den,sp=sp)
}


perturbate <- function(prev_thetas,covarMat){
  mvnfast::rmvn(1,prev_thetas,covarMat)[1,]
}

dperturbate <- function(new_theta,prev_theta,covarMat){
  mvnfast::dmvn(new_theta,prev_theta,covarMat)
}

rprior <- function(){
  #C mu sigma

  c(runif(1,100,200),
    runif(1,1500,3000),
    runif(1,130,280))
}

min_lim <- c(100,1500,130);max_lim <- c(200,3000,280)

prior_numerator <- prod(max_lim-min_lim)

simulate_summaries <- function(the,bw0){
  summary_extraction(simulate(the),bw0)
}
get_dist_profile <- function(num_tries,bw0,true_summaries){
  dists_sp <- numeric(num_tries)
  dists_den <- numeric(num_tries)
  for(i in 1:num_tries){
    thetas <- rprior()
    sim_obs <- simulate(thetas)
    summs <- summary_extraction(sim_obs,bw0)
    dists_sp[i] <- sum((summs$sp-true_summaries$sp)^2)
    dists_den[i] <- sum((summs$den-true_summaries$den)^2)
  }
  return(list(dists_sp=dists_sp,dists_den=dists_den))
}

get_dist_from_theta <- function(the,bw0,true_summaries,weight_den,weight_sp){
  sm <- summary_extraction(simulate(the),bw0)
  sum(weight_den*(true_summaries$den-sm$den)^2)+sum(weight_sp*(true_summaries$sp-sm$sp)^2)
}

M <- 500
set.seed(123)
max_gen <- 20
nsim_pilot <- 1e+3+1

res <- mclapply(1:n_iter,function(i){

  start <- Sys.time()
  true_obs <- simulate_jr()
  bw0 <- bw.nrd0(true_obs)
  true_summaries <- summary_extraction(true_obs,bw0)
  rm(true_obs)

  dists_res <- get_dist_profile(nsim_pilot,bw0,true_summaries)

  weight_sp <- 1/median(dists_res$dists_sp)
  weight_den <- 1/median(dists_res$dists_den)
  dists <- quantile(rowSums(do.call(cbind,dists_res)),probs=seq(0,1,.05))

  sbdf <- do.call(cbind,dists_res)
  sbdf <- sbdf[order(rowSums(sbdf)),]


  delta1 <- sbdf[nsim_pilot/2,1]*weight_sp+ sbdf[nsim_pilot/2,2]*weight_den

  rm(dists_res,sbdf)


  ra <- 1;nsim <- 0;
  theta1 <- matrix(0,M,3)
  prev_dists <- numeric(M)
  internal_nsim <- 0
  for(h in 1:M){
    D <- Inf
    while(D>delta1){
      internal_nsim <- internal_nsim+1
      thetas <- rprior()

      D <- get_dist_from_theta(thetas,bw0,true_summaries,weight_den,weight_sp)
    }
    prev_dists[h] <- D
    theta1[h,] <- thetas
  }
  prev_theta <- theta1
  theta_new <- matrix(0,M,3)
  weights <- rep(1/M,M)

  all_dists <- matrix(0,M,max_gen)
  all_theta <- array(0,c(M,max_gen,3))
  all_delta <- numeric(max_gen)
  all_weights <- matrix(0,M,max_gen)
  nsims <- integer(max_gen)

  all_dists[,1] <- prev_dists
  all_theta[,1,] <- prev_theta
  nsims[1] <- internal_nsim
  all_weights[,1] <- weights
  all_delta[1] <- delta1


  while(ra<max_gen&(M/internal_nsim>.01)&difftime(Sys.time(),start,units="mins")<(230*60)){
    ra <- ra+1
    delta_cur <- median(prev_dists)
    prev_dists <- c()
    new_weights <- weights
    covMat <- 2*cov.wt(prev_theta, wt = weights)$cov
    internal_nsim <- 0
    for(h in 1:M){
      D=Inf
      while(D>delta_cur){
        internal_nsim <- internal_nsim+1
        theta_draw <- prev_theta[sample(1:M,1,prob=weights),]

        theta_cur <- perturbate(theta_draw,covMat)
        while(any((theta_cur-min_lim<0)|(max_lim-theta_cur<0))){
          theta_cur <- perturbate(theta_draw,covMat)
        }

        D <- get_dist_from_theta(theta_cur,bw0,true_summaries,weight_den,weight_sp)
      }
      prev_dists[h] <- D
      theta_new[h,] <- theta_cur
      new_weights[h] <- prior_numerator/sum(weights*apply(prev_theta,1,function(prt)dperturbate(theta_cur,prt,covMat)))
    }
    weights <- new_weights/sum(new_weights)
    prev_theta <- theta_new

    all_dists[,ra] <- prev_dists
    all_theta[,ra,] <- prev_theta
    nsims[ra] <- internal_nsim
    all_weights[,ra] <- weights
    all_delta[ra] <- delta_cur

  }

  return(list(all_dists=all_dists,all_theta=all_theta,nsims=nsims,
              all_weights=all_weights,all_delta=all_delta,time_passed=Sys.time()-start))

},mc.cores=num.cores)


save(res,file="abcLongTime.RData")
