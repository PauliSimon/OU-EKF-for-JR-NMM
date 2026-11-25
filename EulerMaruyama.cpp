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
arma::vec jr_em(int N,
                const arma::vec &x0,
                double Delta,
                double A, double a,
                double B, double b,
                double C1, double C2, double C3, double C4,
                const arma::vec &muFixed,
                const arma::mat &randomness,
                double nuMax = 5.0, double nu0 = 6.0, double r = 0.56) {
  
  if (N <= 0) Rcpp::stop("N must be positive.");
  if (x0.n_elem != 6) Rcpp::stop("x0 must be length 6.");
  if (muFixed.n_elem != 3) Rcpp::stop("muFixed must be length 3.");
  if (N > 1 && (randomness.n_rows != (std::size_t)(N-1) || randomness.n_cols != 3))
    Rcpp::stop("randomness must be (N-1) x 3 when N > 1.");
  
  arma::mat Xfull((std::size_t)N, 6);
  Xfull.zeros();
  
  Xfull.row(0) = x0.t();
  
  
  for (int i = 1; i < N; ++i) {
    const arma::rowvec prev = Xfull.row(i - 1); // length 6
    
    double x1 = prev(0) + Delta * prev(3);
    double x2 = prev(1) + Delta * prev(4);
    double x3 = prev(2) + Delta * prev(5);
    
    double r1 = (N > 1) ? randomness(i - 1, 0) : 0.0;
    double r2 = (N > 1) ? randomness(i - 1, 1) : 0.0;
    double r3 = (N > 1) ? randomness(i - 1, 2) : 0.0;
    
    double v1 = prev(3) + Delta * ( A * a * ( muFixed(0) + sigm(prev(1) - prev(2),nuMax,nu0,r) )
                                    - 2.0 * a * prev(3)
                                    - a * a * prev(0) )
    + r1;
    
    double v2 = prev(4) + Delta * ( A * a * ( muFixed(1) + C2 * sigm(C1 * prev(0),nuMax,nu0,r) )
                                    - 2.0 * a * prev(4)
                                    - a * a * prev(1) )
    + r2;
    
    double v3 = prev(5) + Delta * ( B * b * ( muFixed(2) + C4 * sigm(C3 * prev(0),nuMax,nu0,r) )
                                    - 2.0 * b * prev(5)
                                    - b * b * prev(2) )
    + r3;
    
    Xfull(i, 0) = x1;
    Xfull(i, 1) = x2;
    Xfull(i, 2) = x3;
    Xfull(i, 3) = v1;
    Xfull(i, 4) = v2;
    Xfull(i, 5) = v3;
  }
  
  return conv_to<vec>::from(Xfull.col(1)-Xfull.col(2));
}

vec em_propagate(const arma::vec& prev,
                 const arma::mat& stepMat,
                 const arma::mat& covar_mat,
                 const double obs_var,
                 const int N, const double Delta, const double mu, const double A, const int a,
                 const int B, const int b, const double C1,
                 const double nuMax,const  double nu0, const double r){
  
  vec::fixed<6> res;
  const double C2 = 0.8*C1;
  const double C3 = 0.25*C1;
  const double C4 = 0.25*C1;
 
  double x1 = prev(0) + Delta * prev(3);
  double x2 = prev(1) + Delta * prev(4);
  double x3 = prev(2) + Delta * prev(5);
  
  double v1 = prev(3) + Delta * ( A * a * ( sigm(prev(1) - prev(2),nuMax,nu0,r) )
                                  - 2.0 * a * prev(3)
                                  - a * a * prev(0) );
  
  double v2 = prev(4) + Delta * ( A * a * ( mu + C2 * sigm(C1 * prev(0),nuMax,nu0,r) )
                                  - 2.0 * a * prev(4)
                                  - a * a * prev(1) );
  
  double v3 = prev(5) + Delta * ( B * b * (  C4 * sigm(C3 * prev(0),nuMax,nu0,r) )
                                  - 2.0 * b * prev(5)
                                  - b * b * prev(2) );
  
  res(0) = x1;
  res(1) = x2;
  res(2) = x3;
  res(3) = v1;
  res(4) = v2;
  res(5) = v3;
  
  return res;                                                                                                
}

mat em_function_jacobi(const arma::vec& prev,
                       const arma::mat& stepMat,
                       const arma::mat& covar_mat,
                       const double obs_var,
                       const int N, const double Delta, const double mu, const double A, const int a,
                       const int B, const int b, const double C1,
                       const double nuMax,const  double nu0, const double r){
  mat::fixed<6,6> res;
  res.eye();
  const double C2 = 0.8*C1;
  const double C3 = 0.25*C1;
  const double C4 = 0.25*C1;
  for(int i = 0;i<3;++i) res(i,i+3)=Delta;
  res(3,0) = -Delta*a*a;
  res(3,1) = Delta*A*a*sigm_deriv(prev(1)-prev(2),nuMax,nu0,r);
  res(3,2) = -Delta*A*a*sigm_deriv(prev(1)-prev(2),nuMax,nu0,r);
  res(3,3) -= 2.0*a*Delta;
  
  res(4,0) = Delta*A*a*C1*C2*sigm_deriv(C1*prev(0),nuMax,nu0,r);
  res(4,1) = -Delta*a*a;
  res(4,4) -= 2.0*a*Delta;
  
  res(5,0) = Delta*B*b*C3*C4*sigm_deriv(C3*prev(0),nuMax,nu0,r);
  res(5,2) = -Delta*b*b;
  res(5,5) -= 2*b*Delta;
  
  return res;
}



// [[Rcpp::export]]
double em_ekf_nll(const arma::vec& observations,
                  const arma::mat& stepMat,
                  const arma::mat& covar_mat,
                  const double obs_var,
                  const int N, const double Delta, const double mu, const double A, const int a,
                  const int B, const int b, const double C1,
                  const double nuMax,const  double nu0, const double r,const uword burnin=2){
  
  
  using vec6 = arma::vec::fixed<6>;
  using mat66 = arma::mat::fixed<6,6>;
  
  vec6 mkCur;
  mkCur.zeros();
  mat66 PkCur;
  PkCur.eye();
  
  const double C2 = 0.8*C1;
  const double C3 = 0.25*C1;
  const double C4 = 0.25*C1;
  double neg_llik = 0.0;
  
  
  mat66 PkBar,function_jacobi,joseph_form;
  vec6 Kk,mkBar;
  
  // mat66 random_jacobi;random_jacobi.eye();
  
  double vk,sk;
  
  for(int i=0;i<N;++i){
    mkBar = em_propagate(mkCur,stepMat,covar_mat,obs_var, N, Delta,  mu,  A,  a,
                         B, b,C1,nuMax, nu0, r);
    function_jacobi = em_function_jacobi(mkCur,stepMat,covar_mat,obs_var,
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

// [[Rcpp::export]]
List em_ekf_reconstruct(const arma::vec& observations,
                  const arma::mat& stepMat,
                  const arma::mat& covar_mat,
                  const double obs_var,
                  const int N, const double Delta, const double mu, const double A, const int a,
                  const int B, const int b, const double C1,
                  const double nuMax,const  double nu0, const double r){
  
  
  using vec6 = arma::vec::fixed<6>;
  using mat66 = arma::mat::fixed<6,6>;
  
  vec6 mkCur;
  mkCur.zeros();
  mat66 PkCur;
  PkCur.eye();
  
  arma::mat mk(N, 6, arma::fill::zeros);
  arma::cube Pk(6,6,N, arma::fill::zeros);
  
  const double C2 = 0.8*C1;
  const double C3 = 0.25*C1;
  const double C4 = 0.25*C1;
  double neg_llik = 0.0;
  
  
  mat66 PkBar,function_jacobi,joseph_form;
  vec6 Kk,mkBar;
  
  // mat66 random_jacobi;random_jacobi.eye();
  
  double vk,sk;
  
  for(int i=0;i<N;++i){
    mkBar = em_propagate(mkCur,stepMat,covar_mat,obs_var, N, Delta,  mu,  A,  a,
                         B, b,C1,nuMax, nu0, r);
    function_jacobi = em_function_jacobi(mkCur,stepMat,covar_mat,obs_var,
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
arma::mat em_ekf_terms(const arma::vec& observations,
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
  vec6 Kk, mkBar;
  double vk, sk;
  
  arma::mat out(N, 2, arma::fill::zeros);   // col 0: nll term, col 1: standardized innovation
  
  for(int i = 0; i < N; ++i){
    mkBar = em_propagate(mkCur, stepMat, covar_mat, obs_var, N, Delta, mu, A, a,
                         B, b, C1, nuMax, nu0, r);
    function_jacobi = em_function_jacobi(mkCur, stepMat, covar_mat, obs_var,
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