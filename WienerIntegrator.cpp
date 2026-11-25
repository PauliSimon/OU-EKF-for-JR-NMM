// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cstring>
#include <cmath>
using namespace Rcpp;
using namespace arma;



static const vec::fixed<6> H = {0,1,-1,0,0,0};
static inline double path_to_observation(vec x){
  return x(1)-x(2);
}

// inline sigmoid; TODO replace if we want to alter these parameters (do we though?)

static inline double sigm(double x,double nuMax = 5.0, double nu0 = 6.0, double r = 0.56) {
  return nuMax / (1.0 + std::exp(r * (nu0 - x)));
}


static inline double sigm_deriv(double x,double nuMax = 5.0, double nu0 = 6.0, double r = 0.56){
  // if(abs(x-nuMax)>500)return(0)
  //   exp(r*(nu0+x))*nuMax*r/((exp(nu0*r)+exp(r*x))^2)
  
  if(((x-nu0<-500)||(x-nu0>500))) return 0;
  return std::exp(r*(nu0+x))*nuMax*r/pow((std::exp(nu0*r)+std::exp(r*x)),2);
}


static inline double sech(double x){
  return 1.0/std::cosh(x);
}

// [[Rcpp::export]]
vec wi_propagate(const arma::vec& half,
                 const arma::mat& stepMat,
                 const arma::mat& covar_mat,
                 const double obs_var,
                 const int N, const double Delta, const double mu, const double A, const int a,
                 const int B, const int b, const double C1,
                 const double nuMax,const  double nu0, const double r){
  
  const double C2 = 0.8*C1;
  const double C3 = 0.25*C1;
  const double C4 = 0.25*C1;
  
  vec::fixed<6> wi_nonlinear;
  wi_nonlinear.zeros();
  wi_nonlinear[3] = A*a*sigm(half[1]-half[2],nuMax,nu0,r);
  wi_nonlinear[4] = A*a*(mu + C2*sigm(C1*half[0],nuMax,nu0,r));
  wi_nonlinear[5] = B*b*C4*sigm(C3*half[0],nuMax,nu0,r);

  return stepMat*half+stepMat*(wi_nonlinear*Delta);
 }



// [[Rcpp::export]]
mat wi_jacobi_function(const arma::vec& half,
                 const arma::mat& stepMat,
                 const arma::mat& covar_mat,
                 const double obs_var,
                 const int N, const double Delta, const double mu, const double A, const int a,
                 const int B, const int b, const double C1,
                 const double nuMax,const  double nu0, const double r){
  
  const double C2 = 0.8*C1;
  const double C3 = 0.25*C1;
  const double C4 = 0.25*C1;
  mat::fixed<6,6> res; res.eye();
  
  double inter = Delta*A*a*sigm_deriv(half[1]-half[2],nuMax,nu0,r);
  res(3,1) = inter;
  res(3,2) = -inter;
  res(4,0) = Delta*A*a*C1*C2*sigm_deriv(C1*half[0],nuMax,nu0,r);
  res(5,0) = Delta*B*b*C3*C4*sigm_deriv(C3*half[0],nuMax,nu0,r);
  
  // res = res*stepMat;
  // for(int i=0;i<6;++i) res(i,i) += 1;
  
  
  return stepMat*res*stepMat;
  
}

// [[Rcpp::export]]
double wi_ekf_nll(const arma::vec& observations,
                  const arma::mat& stepMat,
                  const arma::mat& covar_mat,
                  const double obs_var,
                  const int N, const double Delta, const double mu, const double A, const int a,
                  const int B, const int b, const double C1,
                  const double nuMax,const  double nu0, const double r,const uword burnin = 2){
  
  
  using vec6 = arma::vec::fixed<6>;
  using mat66 = arma::mat::fixed<6,6>;
  
  mat66 doubleStep = stepMat*stepMat;
  
  vec6 mkCur;
  mkCur.zeros();
  mat66 PkCur;
  PkCur.eye();
  double neg_llik = 0.0;
  
  
  mat66 PkBar,function_jacobi,joseph_form;
  vec6 Kk,mkBar,half;
  
  // mat66 random_jacobi;random_jacobi.eye();
  
  double vk,sk;
  
  for(int i=0;i<N;++i){
  //   
    half= stepMat*mkCur;
    mkBar = wi_propagate(half,stepMat,covar_mat,obs_var, N, Delta,  mu,  A,  a,
                         B, b,C1,nuMax, nu0, r);
    function_jacobi = wi_jacobi_function(half,stepMat,covar_mat,obs_var,
                                         N, Delta, mu, A, a, B, b, C1, nuMax,nu0,r);
    PkBar = function_jacobi*PkCur*function_jacobi.t()+covar_mat;

    vk = observations(i)-path_to_observation(mkBar);
    sk = as_scalar(H.t()*PkBar*H)+obs_var;

    Kk = PkBar*H/sk;

    mkCur = mkBar+ Kk*vk;

    joseph_form.eye();
    joseph_form -= Kk*H.t();
    PkCur = joseph_form*PkBar*joseph_form.t()+Kk*obs_var*Kk.t();


    if(i>=burnin) neg_llik += 0.5*(log(2*M_PI*sk)+vk*vk/(sk));
  }
  return neg_llik;
}

vec wiener_nonlinear(const vec& cur, const double Delta, const double mu, const double A, const int a,
                     const int B, const int b, const double C1,
                     const double nuMax,const  double nu0, const double r){
  
  const double C2 = 0.8*C1;
  const double C3 = 0.25*C1;
  const double C4 = 0.25*C1;
  
  vec::fixed<6> res;
  res.zeros();
  
  res[3] = A*a*sigm(cur[1]-cur[2],nuMax,nu0,r);
  res[4] = A*a*(mu + C2*sigm(C1*cur[0],nuMax,nu0,r));
  res[5] = B*b*C4*sigm(C3*cur[0],nuMax,nu0,r); 
  
  return res;
  
}




// IMPORTANT!!!! This is with randoms that are Delta*Sigma^2 var, no stepmat yet
// [[Rcpp::export]]
vec wiener_simulate(const mat& randoms,
                    const arma::mat& stepMat,
                    const int N, const double Delta, const double mu, const double A, const int a,
                    const int B, const int b, const double C1,
                    const double nuMax,const  double nu0, const double r){
  
  using vec6 = arma::vec::fixed<6>;
  using mat66 = arma::mat::fixed<6,6>;
  
  //TODO check if randoms 3 col
  
  
  
  vec6 cur;
  
  vec obs(N);
  
  cur.zeros();
  obs(0)=0;
  for(int cur_iter=0;cur_iter<(N-1);++cur_iter){
    cur = stepMat*cur;
    cur += wiener_nonlinear(cur,Delta, mu,A, a,
                                       B,b,C1,nuMax, nu0, r)*Delta;
    for(int i = 0;i<3;++i)cur[i+3] += randoms(cur_iter,i); 
    cur = stepMat*cur;
    obs(cur_iter+1) = cur(1)-cur(2);
  }
  return obs;
}


// [[Rcpp::export]]
mat wiener_simulate_fullpath(const mat& randoms,
                    const arma::mat& stepMat,
                    const int N, const double Delta, const double mu, const double A, const int a,
                    const int B, const int b, const double C1,
                    const double nuMax,const  double nu0, const double r){
  
  using vec6 = arma::vec::fixed<6>;
  using mat66 = arma::mat::fixed<6,6>;
  
  vec6 cur;
  
  mat obs(N,6);
  
  cur.zeros();
  obs.row(0).zeros();
  for(int cur_iter=0;cur_iter<(N-1);++cur_iter){
    cur = stepMat*cur;
    cur += wiener_nonlinear(cur,Delta, mu,A, a,
                            B,b,C1,nuMax, nu0, r)*Delta;
    for(int i = 0;i<3;++i)cur[i+3] += randoms(cur_iter,i);
    cur = stepMat*cur;
    obs.row(cur_iter+1) = cur.t();
  }
  return obs;
}






// [[Rcpp::export]]
List wi_ekf_reconstruct(const arma::vec& observations,
                  const arma::mat& stepMat,
                  const arma::mat& covar_mat,
                  const double obs_var,
                  const int N, const double Delta, const double mu, const double A, const int a,
                  const int B, const int b, const double C1,
                  const double nuMax,const  double nu0, const double r){
  
  
  using vec6 = arma::vec::fixed<6>;
  using mat66 = arma::mat::fixed<6,6>;
  
  arma::mat mk(N, 6, arma::fill::zeros);
  arma::cube Pk(6,6,N, arma::fill::zeros);
  
  mat66 doubleStep = stepMat*stepMat;
  
  vec6 mkCur;
  mkCur.zeros();
  mat66 PkCur;
  PkCur.eye();
  double neg_llik = 0.0;
  
  
  mat66 PkBar,function_jacobi,joseph_form;
  vec6 Kk,mkBar,half;
  
  // mat66 random_jacobi;random_jacobi.eye();
  
  double vk,sk;
  
  for(int i=0;i<N;++i){
    //   
    half= stepMat*mkCur;
    mkBar = wi_propagate(half,stepMat,covar_mat,obs_var, N, Delta,  mu,  A,  a,
                         B, b,C1,nuMax, nu0, r);
    function_jacobi = wi_jacobi_function(half,stepMat,covar_mat,obs_var,
                                         N, Delta, mu, A, a, B, b, C1, nuMax,nu0,r);
    PkBar = function_jacobi*PkCur*function_jacobi.t()+covar_mat;
    
    vk = observations(i)-path_to_observation(mkBar);
    sk = as_scalar(H.t()*PkBar*H)+obs_var;
    
    Kk = PkBar*H/sk;
    
    mkCur = mkBar+ Kk*vk;
    
    joseph_form.eye();
    joseph_form -= Kk*H.t();
    PkCur = joseph_form*PkBar*joseph_form.t()+Kk*obs_var*Kk.t();
    
    mk.row(i) = mkCur.t();
    Pk.slice(i) = PkCur;
  }
  
  List ret;
  ret["mk"]=mk;
  ret["Pk"]=Pk;
  
  return ret;
}

// [[Rcpp::export]]
arma::mat wi_ekf_terms(const arma::vec& observations,
                       const arma::mat& stepMat,
                       const arma::mat& covar_mat,
                       const double obs_var,
                       const int N, const double Delta, const double mu, const double A, const int a,
                       const int B, const int b, const double C1,
                       const double nuMax, const double nu0, const double r,
                       const int burnin = 2){
  
  using vec6 = arma::vec::fixed<6>;
  using mat66 = arma::mat::fixed<6,6>;
  
  vec6 mkCur;
  mkCur.zeros();
  mat66 PkCur;
  PkCur.eye();
  
  mat66 PkBar, function_jacobi, joseph_form;
  vec6 Kk, mkBar, half;
  double vk, sk;
  
  arma::mat out(N, 2, arma::fill::zeros);   // col 0: nll term, col 1: standardized innovation
  
  for(int i = 0; i < N; ++i){
    half = stepMat * mkCur;
    mkBar = wi_propagate(half, stepMat, covar_mat, obs_var, N, Delta, mu, A, a,
                         B, b, C1, nuMax, nu0, r);
    function_jacobi = wi_jacobi_function(half, stepMat, covar_mat, obs_var,
                                         N, Delta, mu, A, a, B, b, C1, nuMax, nu0, r);
    PkBar = function_jacobi * PkCur * function_jacobi.t() + covar_mat;
    
    vk = observations(i) - path_to_observation(mkBar);
    sk = as_scalar(H.t() * PkBar * H) + obs_var;
    
    Kk = PkBar * H / sk;
    mkCur = mkBar + Kk * vk;
    
    joseph_form.eye();
    joseph_form -= Kk * H.t();
    PkCur = joseph_form * PkBar * joseph_form.t() + Kk * obs_var * Kk.t();
    
    if(i >= burnin){
      out(i, 0) = 0.5 * (std::log(2.0 * M_PI * sk) + vk * vk / sk);
    }
    out(i, 1) = vk / std::sqrt(sk);
  }
  return out;
}