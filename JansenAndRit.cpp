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
static inline double sigm_second_deriv(
    double x,
    double nuMax = 5.0,
    double nu0 = 6.0,
    double r = 0.56
) {
    const double s = sigm(x, nuMax, nu0, r);
    const double u = s / nuMax;

    return r*r*s*(1.0-u)*(1.0-2.0*u);
}

static inline double sech(double x){
  return 1.0/std::cosh(x);
}


static inline vec nonlinear_step(vec x,double stepSize, double mu, double C1, double C2, double C3, double C4,
                                 double A, double a, double B, double b,double nuMax = 5.0, double nu0 = 6.0, double r = 0.56){
  //x+stepSize*c(rep(0,3),A*a*c(sigm(x[2]-x[3]),C2*sigm(C1*x[1])+220),B*b*C4*sigm(C3*x[1])
  arma::vec::fixed<6> rhs=  {0,0,0,A*a*sigm(x(1)-x(2),nuMax,nu0,r),A*a*(C2*sigm(C1*x(0),nuMax,nu0,r)+mu),B*b*C4*sigm(C3*x(0),nuMax,nu0,r)};
  return x+stepSize*rhs;
}


static inline vec nonlinear_step_derivative(const vec& x,double stepSize, double mu, double C1, double C2, double C3, double C4,
                                            double A, double a, double B, double b,double nuMax = 5.0, double nu0 = 6.0, double r = 0.56){
  return {A*a*sigm_deriv(x(1)-x(2),nuMax,nu0,r)*stepSize,
    A*a*C1*C2*sigm_deriv(C1*x(0),nuMax,nu0,r)*stepSize,B*b*C3*C4*sigm_deriv(C3*x(0),nuMax,nu0,r)*stepSize};
}

static inline mat full_step_derivative(vec x,double stepSize,mat stepMat, double mu, double C1, double C2, double C3, double C4,
                                       double A, double a, double B, double b,double nuMax = 5.0, double nu0 = 6.0, double r = 0.56){
  vec::fixed<3> half_step_deriv = nonlinear_step_derivative(stepMat*nonlinear_step(x,stepSize,mu,C1,C2,C3,C4,A,a,B,b,nuMax,nu0,r),
                            stepSize,mu,C1,C2,C3,C4,A,a,B,b,nuMax,nu0,r);
  
  vec::fixed<3> deriv = nonlinear_step_derivative(x,stepSize,mu,C1,C2,C3,C4,A,a,B,b,nuMax,nu0,r);
  
  
  double temp;
  
  // for(int j=0;j<6;++j){
  //   Q(3,j) = Q(3,j) + (Q(1,j)-Q(2,j))*half_step_deriv(0);
  //   Q(4,j)
  // }
  // IMPORTANT!! THIS DEPENDS ON THE SPARSE STRUCTURE OF exp(M delta) i hope that this is not wrong lmao
  
  
  for(int i=4;i<6;++i){
    for(int j=0;j<4;j+=3)stepMat(i,j) += stepMat(0,j)*half_step_deriv(i-3);
  }
  for(int j=1;j<5;j+=3)stepMat(3,j) += stepMat(1,j)*half_step_deriv(0);
  for(int j=2;j<6;j+=3)stepMat(3,j) -=  stepMat(2,j)*half_step_deriv(0);
  
  stepMat(0,1) = stepMat(0,3)*deriv(0);
  stepMat(0,2) = -stepMat(0,3)*deriv(0);
  
  stepMat(1,0) = stepMat(1,4)*deriv(1);
  
  stepMat(2,0) = stepMat(2,5)*deriv(2);
  
  for(int i=3;i<6;++i){
    temp = stepMat(i,3)*deriv(0);
    stepMat(i,1)+=temp;
    stepMat(i,2)-=temp;
  }
  
  stepMat(3,0) += stepMat(3,4)*deriv(1)+stepMat(3,5)*deriv(2);
  
  stepMat(4,0) += stepMat(4,4)*deriv(1);
  
  stepMat(5,0) += stepMat(5,5)*deriv(2);
  
  return stepMat;
  
}



static inline mat random_derivatives(vec x,double stepSize,const mat& stepMat, double mu, double C1, double C2, double C3, double C4,
                                double A, double a, double B, double b,double nuMax = 5.0, double nu0 = 6.0, double r = 0.56){
  
  mat::fixed<6,6> result_mat;
  result_mat.eye();
  
  const vec innerVec = {x(0),x(1),x(2),A*a*stepSize*nuMax/(1+std::exp(r*(nu0-x(1)+x(2))))+x(3),
                    stepSize*A*a*(mu+(C2*nuMax)/(1+std::exp(r*(nu0-C1*x(0)))))+x(4),
                    stepSize*B*b*C4*nuMax/(1.0+std::exp(r*(nu0-C3*x(0))))+x(5)};
  
  const vec nStep = nonlinear_step(x,stepSize,mu,C1,C2,C3,C4,A,a,B,b,nuMax,nu0,r);
  
  const double temp = stepSize*A*a*sigm_deriv(as_scalar(stepMat.row(1)*nStep)-as_scalar(stepMat.row(2)*nStep),nuMax,nu0,r);
  
  result_mat(3,1)=temp;
  result_mat(3,2)=-temp;
  
  result_mat(4,0) = stepSize/4.0*A*a*C1*C2*nuMax*r*pow(sech(r/2.0*(nu0-C1*as_scalar(stepMat.row(0)*innerVec))),2);
  result_mat(5,0) = stepSize/4.0*B*b*C3*C4*nuMax*r*pow(sech(r/2.0*(nu0-C3*as_scalar(stepMat.row(0)*innerVec))),2);
  
  return result_mat;
}


static inline vec G2cur(const vec& Qcur,
                 double A, double a, double C2, double C1,
                 double B, double b, double C4, double C3,
                 double nuMax = 5.0, double nu0 = 6.0, double r = 0.56) {
  vec out(3);
  out(0) = A * a * sigm(Qcur(1) - Qcur(2),nuMax,nu0,r);          
  out(1) = A * a * C2 * sigm(C1 * Qcur(0),nuMax,nu0,r);          
  out(2) = B * b * C4 * sigm(C3 * Qcur(0),nuMax,nu0,r);         
  return out;
}

// [[Rcpp::export]]
arma::vec strang_integrator(
    const int N,
    const double Delta, 
    const arma::mat& stepMatStrang,       
    const arma::mat& strangRandomness,    
    const double A, const double a, const double B, const double b,
    double C1,
    const arma::vec& G1integral,          
    const arma::rowvec& initialState,      
    const double nuMax = 5.0,const  double nu0 = 6.0, const double r = 0.56 // ,const bool dim10 = false
) {
  const double C2 = 0.8*C1;
  const double C3 = 0.25*C1;
  // if(dim10) return strang_integrator_dim10(N,Delta,stepMatStrang,strangRandomness,A,a,B,b,C1,C2,C3,
  //    .3*C1,0.1*C1,0.8*C1,G1integral,initialState,nuMax,nu0,r);
  if (N <= 0) Rcpp::stop("N must be positive.");
  const int N_anc = N-1;
  if ((int)initialState.n_elem != 6) Rcpp::stop("Initial rows must be length 6.");
  if ((int)stepMatStrang.n_rows != 6 || (int)stepMatStrang.n_cols != 6) Rcpp::stop("stepMatStrang must be 6x6.");
  if ((int)strangRandomness.n_cols != 6) Rcpp::stop("strangRandomness must have 6 columns.");
  if ((int)G1integral.n_elem != 3) Rcpp::stop("G1integral must be length 3.");
  if ((int)strangRandomness.n_rows != N_anc) {
    Rcpp::stop("strangRandomness must have %d rows (N-1).", N_anc);
  }
  
  arma::vec observables(N);
  observables(0) = initialState[1] - initialState[2];
  
  
  double step[36];
  std::memcpy(step, stepMatStrang.memptr(), 36 * sizeof(double));
  
  const double* randPtr = strangRandomness.memptr();
  const arma::uword rand_nrows = strangRandomness.n_rows;
  const double* colptr[6];
  for (int c = 0; c < 6; ++c) colptr[c] = randPtr + (arma::uword)c * rand_nrows;
  
  double prev[3];
  double x2half[3];
  double x1step[6];
  double q2cur0, q2cur1, q2cur2;
  double G20, G21, G22;
  double s,s1,s2,s3;

  const double g1_0 = G1integral(0), g1_1 = G1integral(1), g1_2 = G1integral(2);
  const double DeltaHalf = Delta/2;
  
  q2cur0 = initialState[0];q2cur1 = initialState[1];q2cur2=initialState[2];
  for (int i = 3; i < 6; ++i) prev[i-3] = initialState[i];
  
  for (int curIdx = 0; curIdx < N_anc; ++curIdx) {


    s1 = sigm(q2cur1 - q2cur2,nuMax,nu0,r);
    s2 = sigm(C1 * q2cur0,nuMax,nu0,r);
    s3 = sigm(C3 * q2cur0,nuMax,nu0,r);
    G20 = A * a * s1;
    G21 = A * a * C2 * s2;
    G22 = B * b * C3 * s3;
    
    
    x2half[0] = prev[0] + G20 * DeltaHalf + g1_0;
    x2half[1] = prev[1] + G21 * DeltaHalf + g1_1;
    x2half[2] = prev[2] + G22 * DeltaHalf + g1_2;
    
    for (int row = 0; row < 6; ++row) {
      s = 0.0;
      s += step[row] * q2cur0;
      s += step[6 + row] * q2cur1;
      s += step[12 + row] * q2cur2;
      s += step[18 + row] * x2half[0];
      s += step[24 + row] * x2half[1];
      s += step[30 + row] * x2half[2];
      x1step[row] = s;
    }
    
    
    x1step[0] += colptr[0][curIdx];
    x1step[1] += colptr[1][curIdx];
    x1step[2] += colptr[2][curIdx];
    x1step[3] += colptr[3][curIdx];
    x1step[4] += colptr[4][curIdx];
    x1step[5] += colptr[5][curIdx];
    
    q2cur0 = x1step[0];
    q2cur1 = x1step[1];
    q2cur2 = x1step[2];
    
    
    s1 = sigm(q2cur1 - q2cur2,nuMax,nu0,r);
    s2 = sigm(C1 * q2cur0,nuMax,nu0,r);
    s3 = sigm(C3 * q2cur0,nuMax,nu0,r);
    G20 = A * a * s1;
    G21 = A * a * C2 * s2;
    G22 = B * b * C3 * s3;
    
    prev[0] = x1step[3] + G20 * DeltaHalf + g1_0;
    prev[1] = x1step[4] + G21 * DeltaHalf + g1_1;
    prev[2] = x1step[5] + G22 * DeltaHalf + g1_2;
    
    
    observables(curIdx+1) = q2cur1 - q2cur2;
  }
  
  return observables;
}



// [[Rcpp::export]]
arma::mat strang_integrator_full_path(
    int N,
    double Delta, // TODO maybe instead specify 1/Delta to avoid numerical problems?
    const arma::mat& stepMatStrang,       
    const arma::mat& strangRandomness, 
    double A, double a, double B, double b,
    double C1, const arma::vec& G1integral,         
    const arma::rowvec& initialState,     
    const double nuMax = 5.0,const  double nu0 = 6.0, const double r = 0.56
) {
  if (N <= 0) Rcpp::stop("N must be positive.");
  const int N_anc = N-1;
  if ((int)initialState.n_elem != 6) Rcpp::stop("Initial rows must be length 6.");
  if ((int)stepMatStrang.n_rows != 6 || (int)stepMatStrang.n_cols != 6) Rcpp::stop("stepMatStrang must be 6x6.");
  if ((int)strangRandomness.n_cols != 6) Rcpp::stop("strangRandomness must have 6 columns.");
  if ((int)G1integral.n_elem != 3) Rcpp::stop("G1integral must be length 3.");
  if ((int)strangRandomness.n_rows != N_anc) {
    Rcpp::stop("strangRandomness must have %d rows (N-1).", N_anc);
  }
  
  
  const double C2 = 0.8*C1;
  const double C3 = 0.25*C1;
  const double C4 = 0.25*C1;
  
  arma::mat trajectory(N,6);
  double step[36];
  std::memcpy(step, stepMatStrang.memptr(), 36 * sizeof(double));
  const double* randPtr = strangRandomness.memptr();
  const arma::uword rand_nrows = strangRandomness.n_rows;
  const double* colptr[6];
  for (int c = 0; c < 6; ++c) colptr[c] = randPtr + (arma::uword)c * rand_nrows;
  
  
  double prev[3];
  double x2half[3];
  double x1step[6];
  double q2cur0, q2cur1, q2cur2;
  double G20, G21, G22;
  double s,s1,s2,s3;
  int curIdx_succ;
  
  const double g1_0 = G1integral(0), g1_1 = G1integral(1), g1_2 = G1integral(2);
  const double DeltaHalf = Delta/2;
  
  q2cur0 = initialState[0];q2cur1 = initialState[1];q2cur2=initialState[2];
  for (int i = 3; i < 6; ++i) prev[i-3] = initialState[i];
  
  for (int curIdx = 0; curIdx < N_anc; ++curIdx) {
    curIdx_succ = curIdx+1;
    
    s1 = sigm(q2cur1 - q2cur2,nuMax,nu0,r);
    s2 = sigm(C1 * q2cur0,nuMax,nu0,r);
    s3 = sigm(C3 * q2cur0,nuMax,nu0,r);
    G20 = A * a * s1;
    G21 = A * a * C2 * s2;
    G22 = B * b * C4 * s3;
    
    
    x2half[0] = prev[0] + G20 * DeltaHalf + g1_0;
    x2half[1] = prev[1] + G21 * DeltaHalf + g1_1;
    x2half[2] = prev[2] + G22 * DeltaHalf + g1_2;
    
    for (int row = 0; row < 6; ++row) {
      s = 0.0;
      s += step[row] * q2cur0;
      s += step[6 + row] * q2cur1;
      s += step[12 + row] * q2cur2;
      s += step[18 + row] * x2half[0];
      s += step[24 + row] * x2half[1];
      s += step[30 + row] * x2half[2];
      x1step[row] = s;
    }
    
    
    x1step[0] += colptr[0][curIdx];
    x1step[1] += colptr[1][curIdx];
    x1step[2] += colptr[2][curIdx];
    x1step[3] += colptr[3][curIdx];
    x1step[4] += colptr[4][curIdx];
    x1step[5] += colptr[5][curIdx];
    
    // q2cur = x1half[0..2]
    q2cur0 = x1step[0];
    q2cur1 = x1step[1];
    q2cur2 = x1step[2];
    
    
    s1 = sigm(q2cur1 - q2cur2,nuMax,nu0,r);
    s2 = sigm(C1 * q2cur0,nuMax,nu0,r);
    s3 = sigm(C3 * q2cur0,nuMax,nu0,r);
    G20 = A * a * s1;
    G21 = A * a * C2 * s2;
    G22 = B * b * C4 * s3;
    
    // build rhs: [q2cur; x1half[3..5] + G2*Delta + G1integral]
    prev[0] = x1step[3] + G20 * DeltaHalf + g1_0;
    prev[1] = x1step[4] + G21 * DeltaHalf + g1_1;
    prev[2] = x1step[5] + G22 * DeltaHalf + g1_2;
    
    trajectory(curIdx_succ,0)=q2cur0;
    trajectory(curIdx_succ,1)=q2cur1;
    trajectory(curIdx_succ,2)=q2cur2;
    trajectory(curIdx_succ,3)=prev[0];
    trajectory(curIdx_succ,4)=prev[1];
    trajectory(curIdx_succ,5)=prev[2];
    
  }
  
  return trajectory;
}



// [[Rcpp::export]]
List ekf_ou(const arma::vec& observations,
                 const arma::mat& stepMat,
                 const arma::mat& covar_mat,
                 const double obs_var,
                 const int N, const double Delta, const double mu, const double A, const int a,
                 const int B, const int b, const double C1,
                 const double nuMax,const  double nu0, const double r){

  using vec6 = arma::vec::fixed<6>;
  using mat66 = arma::mat::fixed<6,6>;
  double halfStep = Delta/2.0;

  arma::mat mk(N, 6, arma::fill::zeros);
  arma::cube Pk(6,6,N, arma::fill::zeros);

  vec6 mkCur;
  mkCur.zeros();
  mat66 PkCur;
  PkCur.eye();
  
  const double C2 = 0.8*C1;
  const double C3 = 0.25*C1;
  const double C4 = 0.25*C1;


  mat66 PkBar,function_jacobi,random_jacobi,joseph_form;
  vec6 Kk,mkBar;
  double vk,sk;

  for(int i=0;i<N;++i){
    mkBar = nonlinear_step(stepMat*nonlinear_step(mkCur,halfStep,mu,C1,C2,C3,C4,A,a,B,b,nuMax,nu0,r),halfStep,mu,C1,C2,C3,C4,A,a,B,b,nuMax,nu0,r);
    function_jacobi = full_step_derivative(mkCur,halfStep,stepMat,mu,C1,C2,C3,C4,A,a,B,b,nuMax,nu0,r);
    random_jacobi = random_derivatives(mkCur,halfStep,stepMat,mu,C1,C2,C3,C4,A,a,B,b,nuMax,nu0,r);
    PkBar = function_jacobi*PkCur*function_jacobi.t()+random_jacobi*covar_mat*random_jacobi.t();
    
    vk = observations(i)-path_to_observation(mkBar);
    sk = as_scalar(H.t()*PkBar*H)+obs_var;
    Kk = PkBar*H/sk;
    
    mkCur = mkBar+ Kk*vk;
    // PkCur = PkBar - Kk*sk*Kk.t();
    // PkCur = PkBar-Kk*H.t()*PkBar-PkBar*H*Kk.t()+Kk*sk*Kk.t(); //Joseph Form, todo test whether this is quicker than the compressed I-HK thing
    joseph_form.eye();
    joseph_form -= Kk*H.t();
    PkCur = joseph_form*PkBar*joseph_form.t()+Kk*obs_var*Kk.t();
    
    
    mk.row(i)=mkCur.t();
    Pk.slice(i)=PkCur;
  }

  List ret;
  ret["mk"]=mk;
  ret["Pk"]=Pk;
  return ret;
}




// [[Rcpp::export]]
List ekf_innovations(const arma::vec& observations,
            const arma::mat& stepMat,
            const arma::mat& covar_mat,
            const double obs_var,
            const int N, const double Delta, const double mu, const double A, const int a,
            const int B, const int b, const double C1,
            const double nuMax,const  double nu0, const double r){
  
  using vec6 = arma::vec::fixed<6>;
  using mat66 = arma::mat::fixed<6,6>;
  double halfStep = Delta/2.0;
  
  vec6 mk;
  mk.zeros();
  mat66 Pk;
  Pk.eye();
  
  const double C2 = 0.8*C1;
  const double C3 = 0.25*C1;
  const double C4 = 0.25*C1;
  
  vec vk(N),sk(N);
  mat66 PkBar,function_jacobi,random_jacobi,joseph_form;
  vec6 Kk,mkBar;
  double vkCur,skCur;
  
  for(int i=0;i<N;++i){
    mkBar = nonlinear_step(stepMat*nonlinear_step(mk,halfStep,mu,C1,C2,C3,C4,A,a,B,b,nuMax,nu0,r),halfStep,mu,C1,C2,C3,C4,A,a,B,b,nuMax,nu0,r);
    function_jacobi = full_step_derivative(mk,halfStep,stepMat,mu,C1,C2,C3,C4,A,a,B,b,nuMax,nu0,r);
    random_jacobi = random_derivatives(mk,halfStep,stepMat,mu,C1,C2,C3,C4,A,a,B,b,nuMax,nu0,r);
    PkBar = function_jacobi*Pk*function_jacobi.t()+random_jacobi*covar_mat*random_jacobi.t();
    
    vkCur = observations(i)-path_to_observation(mkBar);
    skCur = as_scalar(H.t()*PkBar*H)+obs_var;
    Kk = PkBar*H/skCur;
    
    mk = mkBar+ Kk*vkCur;
    // PkCur = PkBar - Kk*sk*Kk.t();
    // PkCur = PkBar-Kk*H.t()*PkBar-PkBar*H*Kk.t()+Kk*sk*Kk.t(); //Joseph Form, todo test whether this is quicker than the compressed I-HK thing
    joseph_form.eye();
    joseph_form -= Kk*H.t();
    Pk = joseph_form*PkBar*joseph_form.t()+Kk*obs_var*Kk.t();
    
    
    vk(i)=vkCur;
    sk(i)=skCur;
  }
  
  List ret;
  ret["vk"]=vk;
  ret["sk"]=sk;
  return ret;
}



// [[Rcpp::export]]
double construct_ekf_likelihood(const arma::vec& observations,
            const arma::mat& stepMat,
            const arma::mat& covar_mat,
            const double obs_var,
            const int N, const double Delta, const double mu, const double A, const int a,
            const int B, const int b, const double C1,
            const double nuMax,const  double nu0, const double r,const int burnin=250){
  
  using vec6 = arma::vec::fixed<6>;
  using mat66 = arma::mat::fixed<6,6>;
  double halfStep = Delta/2.0;
  
  vec6 mkCur;
  mkCur.zeros();
  mat66 PkCur;
  PkCur.eye();
  
  const double C2 = 0.8*C1; //TODO NO
  const double C3 = 0.25*C1;
  const double C4 = 0.25*C1;
  double neg_llik = 0.0;
  
  
  mat66 PkBar,function_jacobi,random_jacobi,joseph_form;
  vec6 Kk,mkBar;
  double vk,sk;
  
  for(int i=0;i<N;++i){
    mkBar = nonlinear_step(stepMat*nonlinear_step(mkCur,halfStep,mu,C1,C2,C3,C4,A,a,B,b,nuMax,nu0,r),halfStep,mu,C1,C2,C3,C4,A,a,B,b,nuMax,nu0,r);
    function_jacobi = full_step_derivative(mkCur,halfStep,stepMat,mu,C1,C2,C3,C4,A,a,B,b,nuMax,nu0,r);
    random_jacobi = random_derivatives(mkCur,halfStep,stepMat,mu,C1,C2,C3,C4,A,a,B,b,nuMax,nu0,r);
    PkBar = function_jacobi*PkCur*function_jacobi.t()+random_jacobi*covar_mat*random_jacobi.t();
    
    vk = observations(i)-path_to_observation(mkBar);
    sk = as_scalar(H.t()*PkBar*H)+obs_var;
    
    Kk = PkBar*H/sk;
    
    mkCur = mkBar+ Kk*vk;
    // PkCur = PkBar - Kk*sk*Kk.t();
    // PkCur = PkBar-Kk*H.t()*PkBar-PkBar*H*Kk.t()+Kk*sk*Kk.t(); //Joseph Form, todo test whether this is quicker than the compressed I-HK thing
    joseph_form.eye();
    joseph_form -= Kk*H.t();
    PkCur = joseph_form*PkBar*joseph_form.t()+Kk*obs_var*Kk.t();
    // curVar = PkCur(1,1)+PkCur(2,2)-2*PkCur(1,2);
    // if(curVar>0)neg_llik += 0.5*(log(2*M_PI*curVar)+pow(observations(i)-path_to_observation(mkCur),2)/(curVar));
    // else neg_llik += 10000;
    if(i>=burnin)neg_llik += 0.5*(log(2*M_PI*sk)+vk*vk/(sk));
  }
    return neg_llik;
}



// [[Rcpp::export]]
double construct_ekf_likelihood_initials_submitted(const arma::vec& observations,
                                                   vec& mkCur,
            const arma::mat& stepMat,
            const arma::mat& covar_mat,
            const double obs_var,
            const int N, const double Delta, const double mu, const double A, const int a,
            const int B, const int b, const double C1,
            const double nuMax,const  double nu0, const double r,const int burnin=250){
  
  using vec6 = arma::vec::fixed<6>;
  using mat66 = arma::mat::fixed<6,6>;
  double halfStep = Delta/2.0;
  
  // vec6 mkCur;
  // mkCur.zeros();
  mat66 PkCur;
  PkCur.eye();
  
  const double C2 = 0.8*C1; //TODO NO
  const double C3 = 0.25*C1;
  const double C4 = 0.25*C1;
  double neg_llik = 0.0;
  
  
  mat66 PkBar,function_jacobi,random_jacobi,joseph_form;
  vec6 Kk,mkBar;
  double vk,sk;
  
  for(int i=0;i<N;++i){
    mkBar = nonlinear_step(stepMat*nonlinear_step(mkCur,halfStep,mu,C1,C2,C3,C4,A,a,B,b,nuMax,nu0,r),halfStep,mu,C1,C2,C3,C4,A,a,B,b,nuMax,nu0,r);
    function_jacobi = full_step_derivative(mkCur,halfStep,stepMat,mu,C1,C2,C3,C4,A,a,B,b,nuMax,nu0,r);
    random_jacobi = random_derivatives(mkCur,halfStep,stepMat,mu,C1,C2,C3,C4,A,a,B,b,nuMax,nu0,r);
    PkBar = function_jacobi*PkCur*function_jacobi.t()+random_jacobi*covar_mat*random_jacobi.t();
    
    vk = observations(i)-path_to_observation(mkBar);
    sk = as_scalar(H.t()*PkBar*H)+obs_var;

    Kk = PkBar*H/sk;
    
    mkCur = mkBar+ Kk*vk;
    joseph_form.eye();
    joseph_form -= Kk*H.t();
    PkCur = joseph_form*PkBar*joseph_form.t()+Kk*obs_var*Kk.t();
    if(i>=burnin)neg_llik += 0.5*(log(2*M_PI*sk)+vk*vk/(sk));
  }
    return neg_llik;
}


// [[Rcpp::export]]
double construct_ekf_combined_likelihood(
    const arma::vec& observations,
    const arma::mat& stepMatFine,
    const arma::mat& covar_matFine,
    const double obs_var,
    const int N,
    const double Delta,
    const double mu,
    const double A,
    const int a,
    const int B,
    const int b,
    const double C1,
    const double nuMax,
    const double nu0,
    const double r,
    const int nSubsteps = 1,
    const int kStep = 1,
    const double kWeight = 1.0,
    const int burnin = 250
){
  
  using vec6  = arma::vec::fixed<6>;
  using mat66 = arma::mat::fixed<6,6>;
  
  // Internal integration step
  const double delta = Delta / static_cast<double>(nSubsteps);
  const double halfStep = delta / 2.0;
  
  vec6 mkCur;
  mkCur.zeros();
  
  mat66 PkCur;
  PkCur.eye();
  
  const double C2 = 0.8 * C1;
  const double C3 = 0.25 * C1;
  const double C4 = 0.25 * C1;
  
  double neg_llik = 0.0;
  
  auto predict_interval =
    [&](const vec6& mStart,
        const mat66& PStart,
        vec6& mOut,
        mat66& POut){
      
      mOut = mStart;
      POut = PStart;
      
      for(int j = 0; j < nSubsteps*kStep; ++j){
        
        vec6 mBefore = mOut;
        
        mOut =
          nonlinear_step(
            stepMatFine *
              nonlinear_step(
                mBefore,
                halfStep,
                mu, C1, C2, C3, C4,
                A, a, B, b,
                nuMax, nu0, r
              ),
              halfStep,
              mu, C1, C2, C3, C4,
              A, a, B, b,
              nuMax, nu0, r
          );
        
        mat66 Fx =
          full_step_derivative(
            mBefore,
            halfStep,
            stepMatFine,
            mu, C1, C2, C3, C4,
            A, a, B, b,
            nuMax, nu0, r
          );
        
        arma::mat Fq =
          random_derivatives(
            mBefore,
            halfStep,
            stepMatFine,
            mu, C1, C2, C3, C4,
            A, a, B, b,
            nuMax, nu0, r
          );
        
        POut =
          Fx * POut * Fx.t()
          + Fq * covar_matFine * Fq.t();
        
        POut = 0.5 * (POut + POut.t());
      }
    };
    
    for(int i = 0; i < N; i+=kStep){
      
      vec6 mOne;
      mat66 POne;
      
      predict_interval(
        mkCur,
        PkCur,
        mOne,
        POne
      );
      
      double vOne =
        observations(i)
        - path_to_observation(mOne);
      
      double sOne =
      as_scalar(H.t() * POne * H)
        + obs_var;
      
      
      vec6 KOne = POne * H / sOne;
      
      mkCur = mOne + KOne * vOne;
      
      mat66 joseph_form;
      joseph_form.eye();
      joseph_form -= KOne * H.t();
      
      PkCur =
        joseph_form * POne * joseph_form.t()
        + KOne * obs_var * KOne.t();
      
      PkCur = 0.5 * (PkCur + PkCur.t());
      
      
      if(i >= burnin){
        
        neg_llik +=
          0.5 *
          (
              std::log(2.0 * M_PI * sOne)
        + vOne * vOne / sOne
          );
      }
    }
    
    return neg_llik;
}



// [[Rcpp::export]]
double get_ekf_likelihood(mat m, cube P,vec y,int nObs,double obs_var,bool log_likelihood=true){
  double llik = 0.0;
  double var;
  mat Pcur;
  for(int i=0;i<nObs;++i){
    Pcur = P.slice(i);
    var = Pcur(1,1)+Pcur(2,2)-2*Pcur(1,2)+obs_var;
    if(var<=0)llik -= 10000.0;
    else{
      llik -= 0.5*log(2*M_PI*var)+pow(y(i)-m(i,1)+m(i,2),2)/(2.0*var);
    }
  }
  if(log_likelihood) return llik;
  return exp(llik);
}

























// [[Rcpp::export]]
arma::mat lie_trotter_full_path(
    int N,
    double Delta, // TODO maybe instead specify 1/Delta to avoid problems?
    const arma::mat& stepMatLie,      
    const arma::mat& lieRandomness,    
    double A, double a, double B, double b,
    double C1, double C2, double C3, double C4,
    const arma::vec& G1integral,         
    const arma::rowvec& initialState,    
    const double nuMax = 5.0,const  double nu0 = 6.0, const double r = 0.56
) {
  if (N <= 0) Rcpp::stop("N must be positive.");
  const int N_anc = N-1;
  if ((int)initialState.n_elem != 6) Rcpp::stop("Initial rows must be length 6.");
  if ((int)stepMatLie.n_rows != 6 || (int)stepMatLie.n_cols != 6) Rcpp::stop("stepMatLie must be 6x6.");
  if ((int)lieRandomness.n_cols != 6) Rcpp::stop("lieRandomness must have 6 columns.");
  if ((int)G1integral.n_elem != 3) Rcpp::stop("G1integral must be length 3.");
  if ((int)lieRandomness.n_rows != N_anc) {
    Rcpp::stop("lieRandomness must have %d rows (N-1).", N_anc);
  }
  
  arma::mat trajectory(N,6);
 
  double step[36];
  std::memcpy(step, stepMatLie.memptr(), 36 * sizeof(double));
 
  const double* randPtr = lieRandomness.memptr();
  const arma::uword rand_nrows = lieRandomness.n_rows;
  const double* colptr[6];
  for (int c = 0; c < 6; ++c) colptr[c] = randPtr + (arma::uword)c * rand_nrows;
  
  double prev[3];
  double x2half[3];
  double x1step[6];
  double q2cur0, q2cur1, q2cur2;
  double G20, G21, G22;
  double s,s1,s2,s3;
  int curIdx_succ;
  
  const double g1_0 = G1integral(0), g1_1 = G1integral(1), g1_2 = G1integral(2);
  
  q2cur0 = initialState[0];q2cur1 = initialState[1];q2cur2=initialState[2];
  for (int i = 3; i < 6; ++i) prev[i-3] = initialState[i];
  
  for (int curIdx = 0; curIdx < N_anc; ++curIdx) {
    curIdx_succ = curIdx+1;
    
    s1 = sigm(x1step[1] - x1step[2],nuMax,nu0,r);
    s2 = sigm(C1 * x1step[0],nuMax,nu0,r);
    s3 = sigm(C3 * x1step[0],nuMax,nu0,r);
    G20 = A * a * s1;
    G21 = A * a * C2 * s2;
    G22 = B * b * C4 * s3;
    
    
    x2half[0] = x1step[3] + G20 * Delta + g1_0;
    x2half[1] = x1step[4] + G21 * Delta + g1_1;
    x2half[2] = x1step[5] + G22 * Delta + g1_2;
    
    // x1half = step * prev
    for (int row = 0; row < 6; ++row) {
      s = 0.0;
      s += step[row] * x1step[0];
      s += step[6 + row] * x1step[1];
      s += step[12 + row] * x1step[2];
      s += step[18 + row] * x2half[0];
      s += step[24 + row] * x2half[1];
      s += step[30 + row] * x2half[2];
      x1step[row] = s;
    }
    
// + randoms
    x1step[0] += colptr[0][curIdx];
    x1step[1] += colptr[1][curIdx];
    x1step[2] += colptr[2][curIdx];
    x1step[3] += colptr[3][curIdx];
    x1step[4] += colptr[4][curIdx];
    x1step[5] += colptr[5][curIdx];
    
    trajectory(curIdx_succ,0)=x1step[0];
    trajectory(curIdx_succ,1)=x1step[1];
    trajectory(curIdx_succ,2)=x1step[2];
    trajectory(curIdx_succ,3)=x1step[3];
    trajectory(curIdx_succ,4)=x1step[4];
    trajectory(curIdx_succ,5)=x1step[5];
    
   
    
  }
  
  return trajectory;
}


// [[Rcpp::export]]
mat simulate_network_jr(
    const uword N,
    const double Delta, // TODO maybe instead specify 1/Delta to avoid problems?
    const arma::cube& stepMats,       
    const arma::cube& strangRandomness,    
    const vec& As, const vec& as, const vec& Bs,const  vec& bs,
    const vec& Cs,
    const vec& g1Integrals,
    const uword nPop,
    const int coupling_strength,
    const mat& adjacencyMat,         
    const arma::mat& initialStates,     
    const double nuMax = 5.0,const  double nu0 = 6.0, const double r = 0.56
){
  arma::mat observables(N,nPop);
  
  
  // double step[36];
  // std::memcpy(step, stepMatStrang.memptr(), 36 * sizeof(double));
  
  // const double* randPtr = strangRandomness.memptr();
  // const arma::uword rand_nrows = strangRandomness.n_rows;
  // const double* colptr[6];
  // for (int c = 0; c < 6; ++c) colptr[c] = randPtr + (arma::uword)c * rand_nrows;
  
  double s;
  vec prevq0(nPop);
  mat prev(3,nPop);
  mat x2half(3,nPop);
  mat x1step(6,nPop);
  mat q2cur(3,nPop);
  mat G2(3,nPop);
  mat sMat(3,nPop);
  const int N_anc = N-1;
  const double DeltaHalf = Delta/2;
  
  for(int i=0;i<nPop;++i){observables(0,i) = initialStates(i,1) - initialStates(i,2);
  q2cur(0,i) = initialStates(i,0);q2cur(1,i) = initialStates(i,1);q2cur(2,i)=initialStates(i,2);
  for (int j = 3; j < 6; ++j) prev(j-3,i) = initialStates(i,j);}
  
  for (int curIdx = 0; curIdx < N_anc; ++curIdx) {
    for(int curPop=0;curPop<nPop;++curPop){
      prevq0[curPop] = q2cur(0,curPop); 
      sMat(0,curPop) = sigm(q2cur(1,curPop) - q2cur(2,curPop),nuMax,nu0,r);
      sMat(1,curPop) = sigm(Cs[curPop] * q2cur(0,curPop),nuMax,nu0,r);
      sMat(2,curPop) = sigm(Cs[curPop]*0.25 * q2cur(0,curPop),nuMax,nu0,r);
      G2(0,curPop) = As[curPop] * as[curPop] * sMat(0,curPop);
      G2(1,curPop) = As[curPop] * as[curPop] * Cs[curPop]*0.8* sMat(1,curPop);
      G2(2,curPop) = Bs[curPop] * bs[curPop] * Cs[curPop]*0.25 * sMat(2,curPop);
      
      
      x2half(0,curPop) = prev(0,curPop) + G2(0,curPop) * DeltaHalf;
      x2half(1,curPop) = prev(1,curPop) + G2(1,curPop) * DeltaHalf + g1Integrals[curPop];
      x2half(2,curPop) = prev(2,curPop) + G2(2,curPop) * DeltaHalf;
      
      s=0.0;
      for(int otherPop=0;otherPop<nPop;++otherPop){
        if(otherPop!=curPop && adjacencyMat(otherPop,curPop)==1){
          // s += sigm(observables(curIdx,otherPop)); // ESCUAIN-POOLE COUPLING
          if(otherPop>curPop)s += q2cur(0,otherPop); // IRENE COUPLING
          else s += prevq0(otherPop);
        }
      }
      x2half(1,curPop) += DeltaHalf*As[curPop]*as[curPop]*coupling_strength*s;
      
      
      
      
      // x1half = step * prev + randoms
      for (int row = 0; row < 6; ++row) {
        // s = 0.0;
        // 
        // for(int halfIdx = 0;halfIdx<3;++halfIdx){
        //   s += stepMats(row,halfIdx,curPop)*q2cur(halfIdx,curPop)+ stepMats(row,halfIdx+3,curPop)*x2half(halfIdx,curPop);
        // }
        // x1step(row,curPop) = s;
        
        
        x1step(row,curPop) = stepMats(row,row%3,curPop)*q2cur(row%3,curPop) + stepMats(row,row%3+3,curPop)*x2half(row%3,curPop);
      }
      
      
      x1step(0,curPop) += strangRandomness(curIdx,0,curPop);
      x1step(1,curPop) += strangRandomness(curIdx,1,curPop);
      x1step(2,curPop) += strangRandomness(curIdx,2,curPop);
      x1step(3,curPop) += strangRandomness(curIdx,3,curPop);
      x1step(4,curPop) += strangRandomness(curIdx,4,curPop);
      x1step(5,curPop) += strangRandomness(curIdx,5,curPop);
      
      // q2cur = x1half[0..2]
      q2cur(0,curPop) = x1step(0,curPop);
      q2cur(1,curPop) = x1step(1,curPop);
      q2cur(2,curPop) = x1step(2,curPop);
      
      
      sMat(0,curPop) = sigm(q2cur(1,curPop) - q2cur(2,curPop),nuMax,nu0,r);
      sMat(1,curPop) = sigm(Cs[curPop] * q2cur(0,curPop),nuMax,nu0,r);
      sMat(2,curPop) = sigm(Cs[curPop]*0.25 * q2cur(0,curPop),nuMax,nu0,r);
      G2(0,curPop) = As[curPop] * as[curPop] * sMat(0,curPop);
      G2(1,curPop) = As[curPop] * as[curPop] * Cs[curPop]*0.8 * sMat(1,curPop);
      G2(2,curPop) = Bs[curPop] * bs[curPop] * Cs[curPop]*0.25 * sMat(2,curPop);
      
      // build rhs: [q2cur; x1half[3..5] + G2*Delta + G1integral]
      prev(0,curPop) = x1step(3,curPop) + G2(0,curPop) * DeltaHalf;
      prev(1,curPop) = x1step(4,curPop) + G2(1,curPop) * DeltaHalf + g1Integrals[curPop];
      prev(2,curPop) = x1step(5,curPop) + G2(2,curPop) * DeltaHalf;
      
      s=0.0;
      for(int otherPop=0;otherPop<nPop;++otherPop){
        if(otherPop!=curPop && adjacencyMat(otherPop,curPop)==1){
          // s += sigm(observables(curIdx,otherPop)); // ESCUAIN-POOLE COUPLING
          if(otherPop>curPop)s += q2cur(0,otherPop); // IRENE COUPLING
          else s += prevq0(otherPop);
        }
      }
      prev(1,curPop) += DeltaHalf*As[curPop]*as[curPop]*coupling_strength*s;
      
      
      
      observables(curIdx+1,curPop) = q2cur(1,curPop) - q2cur(2,curPop);
    }
  }
  
  return observables; 
}

























// [[Rcpp::export]]
void jacobian_by_x(const sp_mat& jac_rt2x,
                   const sp_mat& stepmat,
                   sp_mat& jac_x,const uword nPop) {
  // const arma::uword nPop = static_cast<arma::uword>(nPop);

  
  uword ind2;
  for(uword ind1=0;ind1<jac_x.n_rows;ind1+=6){
    ind2=ind1+5;
    jac_x.rows(ind1,ind2) = stepmat*jac_x.rows(ind1,ind2);
  }
  
  jac_x = jac_rt2x*jac_x;
}


// [[Rcpp::export]]
arma::sp_mat random_network_jacobian_VECTOR(const uword nPop,
                                     const double Delta,
                                     const arma::vec& A,
                                     const arma::vec& a,
                                     const arma::vec& B,
                                     const arma::vec& b,
                                     const arma::vec& C,
                                     const arma::mat& couplingStrength,
                                     const vec& Xs) {
  
  
  const uword p = 3 * nPop;
  const uword m = 6 * nPop;
  const double scale = 0.5 * Delta;
  
  std::vector<uword> rows;
  std::vector<uword> cols;
  std::vector<double> vals;
  
  rows.reserve(nPop * nPop + 9 * nPop);
  cols.reserve(nPop * nPop + 9 * nPop);
  vals.reserve(nPop * nPop + 9 * nPop);
  
  auto add_entry = [&](uword r, uword col, double v) {
    if (v != 0.0) {
      rows.push_back(r);
      cols.push_back(col);
      vals.push_back(v);
    }
  };
  
  for (uword k = 0; k < nPop; ++k) {
    const uword base = 6 * k;
    
    const double Ak  = A[k];
    const double ak  = a[k];
    const double Bk  = B[k];
    const double bk  = b[k];
    const double C1k = C[k];
    const double C2k = 0.8  * C1k;
    const double C3k = 0.25 * C1k;
    const double C4k = 0.25 * C1k;
    
    const double x1 = Xs(base + 0);
    const double x2 = Xs(base + 1);
    const double x3 = Xs(base + 2);
    
    // Identity block
    add_entry(base + 0, base + 0, 1.0);
    add_entry(base + 1, base + 1, 1.0);
    add_entry(base + 2, base + 2, 1.0);
    add_entry(base + 3, base + 3, 1.0);
    add_entry(base + 4, base + 4, 1.0);
    add_entry(base + 5, base + 5, 1.0);
    
    // Diagonal nonlinear block
    const double s23 = sigm_deriv(x2 - x3);
    const double s1  = sigm_deriv(C1k * x1);
    const double s3  = sigm_deriv(C3k * x1);
    
    // row 3 = p1 update
    add_entry(base + 3, base + 1, scale * ( Ak * ak * s23));
    add_entry(base + 3, base + 2, scale * (-Ak * ak * s23));
    
    // row 4 = p2 update
    add_entry(base + 4, base + 0, scale * (Ak * ak * C1k * C2k * s1));
    
    // row 5 = p3 update
    add_entry(base + 5, base + 0, scale * (Bk * bk * C3k * C4k * s3));
    
    // Off-diagonal coupling: only into p2, exactly as in t2_network_step()
    for (arma::uword j = 0; j < nPop; ++j) {
      if (j == k) continue;
      
      const double v = scale * (Ak * ak * couplingStrength(j, k));
      add_entry(base + 4, 6 * j + 0, v);
    }
  }
  
  umat loc(2, vals.size());
  vec  x(vals.size());
  
  for (uword i = 0; i < vals.size(); ++i) {
    loc(0, i) = rows[i];
    loc(1, i) = cols[i];
    x(i) = vals[i];
  }
  sp_mat J(loc, x);
  return J;
}

// [[Rcpp::export]]
arma::sp_mat random_network_jacobian(const uword nPop,
                             const double Delta,
                             const arma::vec& A,
                             const arma::vec& a,
                             const arma::vec& B,
                             const arma::vec& b,
                             const arma::vec& C,
                             const arma::mat& couplingStrength,
                             const mat& Xs) {
  if (A.n_elem != nPop || a.n_elem != nPop || B.n_elem != nPop || b.n_elem != nPop ||
      C.n_elem != nPop) {
    Rcpp::stop("A, a, B, b, and C must all have length N.");
  }
  if (Xs.n_rows < 3 || Xs.n_cols != nPop) {
    if(Xs.n_rows==6*nPop&&Xs.n_cols==1){
      return random_network_jacobian_VECTOR(nPop,Delta,A,a,B,b,C,couplingStrength,Xs);
    } else {
    Rcpp::stop("Xs must be at least 3 x N (but it should probably be 6xN).");}
  }

  
  const uword p = 3 * nPop;
  const uword m = 6 * nPop;
  const double scale = 0.5 * Delta;
  
  std::vector<uword> rows;
  std::vector<uword> cols;
  std::vector<double> vals;
  
  rows.reserve(nPop * nPop + 9 * nPop);
  cols.reserve(nPop * nPop + 9 * nPop);
  vals.reserve(nPop * nPop + 9 * nPop);
  
  auto add_entry = [&](uword r, uword col, double v) {
    if (v != 0.0) {
      rows.push_back(r);
      cols.push_back(col);
      vals.push_back(v);
    }
  };
  
  for (uword k = 0; k < nPop; ++k) {
    const uword q = 3 * k;
    
    const double Ak  = A[k];
    const double ak  = a[k];
    const double Bk  = B[k];
    const double bk  = b[k];
    const double C1k = C[k];
    const double C2k = 0.8  * C1k;
    const double C3k = 0.25 * C1k;
    const double C4k = 0.25 * C1k;
    
    const double x1 = Xs(0, k);
    const double x2 = Xs(1, k);
    const double x3 = Xs(2, k);
    
    //Id mats
    add_entry(q + 0, q + 0, 1.0);
    add_entry(q + 1, q + 1, 1.0);
    add_entry(q + 2, q + 2, 1.0);
    add_entry(p + q + 0, p + q + 0, 1.0);
    add_entry(p + q + 1, p + q + 1, 1.0);
    add_entry(p + q + 2, p + q + 2, 1.0);
    
    // Diagonal block
    const double s23 = sigm_deriv(x2 - x3);
    const double s1  = sigm_deriv(C1k * x1);
    const double s3  = sigm_deriv(C3k * x1);
    
    add_entry(p + q + 0, q + 1, scale * ( Ak * ak * s23));
    add_entry(p + q + 0, q + 2, scale * (-Ak * ak * s23));
    
    add_entry(p + q + 1, q + 0, scale * (Ak * ak * C1k * C2k * s1));
    add_entry(p + q + 2, q + 0, scale * (Bk * bk * C3k * C4k * s3));
    
    // Off-diagonal blocks
    for (arma::uword j = 0; j < nPop; ++j) {
      if (j == k) continue;
      
      const double v = scale * (Ak * ak * couplingStrength(j, k));
      add_entry(p + q + 1, 3 * j + 0, v);
    }
  }
  
  umat loc(2, vals.size());
  vec  x(vals.size());
  
  for (uword i = 0; i < vals.size(); ++i) {
    loc(0, i) = rows[i];
    loc(1, i) = cols[i];
    x(i) = vals[i];
  }
  
  sp_mat J(loc, x);
  return J;
}


// [[Rcpp::export]]
mat t2_network_step(const vec& X_slice,const double DeltaHalf,const uword nPop,const double mu,const vec& As,
                const vec& as,const vec& Bs,const vec& bs, const vec& Cs,const mat& coupling_strength,
                const double nuMax, const double nu0, const double r_for_sigm){
  vec::fixed<6> popX;
  mat resultMat(6,nPop); double coup;
  for(int i=0;i<nPop;++i){
    popX = X_slice.subvec(6*i,6*i+5); 
    for(int j=0;j<6;++j) resultMat(j,i) = popX(j);
    resultMat(3,i) += DeltaHalf*As(i)*as(i)*sigm(popX(1)-popX(2), nuMax, nu0, r_for_sigm);
    resultMat(4,i) += DeltaHalf*As(i)*as(i)*(Cs(i)*0.8*sigm(Cs(i)*popX(0), nuMax, nu0, r_for_sigm)+mu);
    resultMat(5,i) += DeltaHalf*Bs(i)*bs(i)*Cs(i)*0.25*sigm(Cs(i)*0.25*popX(0), nuMax, nu0, r_for_sigm);
    
    coup = 0.0;
    for(int otherPop=0;otherPop<nPop;++otherPop){
      if(otherPop==i)continue;
      coup += coupling_strength(otherPop,i)*X_slice(6*otherPop);
    }
    
    resultMat(4,i) += DeltaHalf*As(i)*as(i)*coup;
  }
  return resultMat;
}

// THIS CALCULATES blockdiag(A,A,A,A)*B
// [[Rcpp::export]]
static inline sp_mat blockdiag_times_mat(const sp_mat& A,const sp_mat& B){
  const uword dim_A = A.n_rows;const uword dim_B=B.n_rows;const uword multi = dim_B/dim_A;
  uword ind2;
  sp_mat C(dim_B,dim_B);
  for (uword ind1=0;ind1<dim_B;ind1+=dim_A) {
    ind2 = ind1+dim_A-1;
    C.rows(ind1,ind2) = A*B.rows(ind1, ind2);
  }
  return C;
}





// [[Rcpp::export]]
double ekf_neg_llik_network(const mat& observations, 
                const uword N,
                const double Delta, // TODO maybe instead specify 1/Delta to avoid problems?
                const sp_mat& stepMat,       
                const sp_mat& covar_mat,   
                const vec& As, const vec& as, const vec& Bs,const  vec& bs,
                const vec& Cs,const double mu,
                double obs_var,
                const uword nPop,
                const mat& coupling_strength,        
                const double nuMax = 5.0,const  double nu0 = 6.0, const double r = 0.56){

  const double halfStep = Delta/2.0;
  const uword dim = 6*nPop;
  
  vec curEst(nPop),mkCur(dim),mkBar(dim),vk(nPop),xProcVec(dim);
  mkCur.zeros();


  double neg_llik = 0.0;
  
  
  sp_mat PkBar(dim,dim),joseph_form(dim,dim),H_by_x(nPop,dim),
    X_slice, PkCur(dim,dim), random_jacobi(dim,dim),function_jacobi(dim,dim),
    Kk(dim,nPop);
  mat Sk(nPop,nPop),SkInv(nPop,nPop);
  // if(obs_var<1)obs_var=1;
  sp_mat obs_vars = obs_var * speye<sp_mat>(nPop, nPop);

  

  H_by_x.zeros();
  for(int i=0;i<nPop;++i){
    H_by_x(i,i*6+1)=1;
    H_by_x(i,i*6+2)=-1;
  }
  
  mat xProcess(6,nPop);
  
  
  PkCur = eye(dim,dim);
  for(int i=0;i<N;++i){
    
    xProcess = t2_network_step(mkCur, halfStep,nPop,mu, As,
                             as, Bs, bs,  Cs,coupling_strength,
                             nuMax, nu0, r);
    xProcVec = vectorise(stepMat*xProcess);
    
    random_jacobi = random_network_jacobian_VECTOR(nPop,
      Delta, As, as,Bs,bs,Cs,coupling_strength,xProcVec);
    function_jacobi = random_network_jacobian_VECTOR(nPop,
      Delta, As, as,Bs,bs,Cs,coupling_strength,mkCur);
    jacobian_by_x(random_jacobi,stepMat,function_jacobi,nPop);
    mkBar = vectorise(t2_network_step(xProcVec,halfStep,nPop,mu,As,
                               as,Bs,bs,Cs,coupling_strength,
                               nuMax,nu0,r));
    
    
    PkBar = function_jacobi*PkCur*function_jacobi.t() + random_jacobi*blockdiag_times_mat(covar_mat,random_jacobi.t());
    // for(int ah=0;ah<nPop;++ah)curEst(ah) = mkBar(ah*6+1)-mkBar(ah*6+2);
    curEst = H_by_x*mkBar;
    vk = observations.row(i).t()-curEst;
    Sk = H_by_x*PkBar*H_by_x.t() + obs_vars;
    
    
    //Try to ensure symmetric Sk
    Sk = (Sk+Sk.t())/2.0;
    
    
    // arma::mat L;
    // bool is_spd = arma::chol(L, Sk);
    // if(!is_spd)Rcout<<covar_mat;
    
    SkInv = inv_sympd(Sk);
    Kk = PkBar*H_by_x.t()*SkInv;
    // Kk = solve(Sk, (H_by_x * PkBar).t()).t();
    mkCur = mkBar + Kk*vk;
    joseph_form.eye();
    joseph_form -= Kk*H_by_x;
    // PkCur = PkBar-Kk*Sk*Kk.t();
    // if(i==5)return neg_llik;
    PkCur = joseph_form*PkBar*joseph_form.t()+Kk*obs_var*Kk.t();
    neg_llik += 0.5*(nPop*log(2*M_PI)+log_det_sympd(Sk)+as_scalar(vk.t()*SkInv*vk)); //TODO check correctness
  }
  return neg_llik;
  
}







///////////////////////////////////////////////////////////////////////////////
//////////////////////////////// CISTUFF ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////



// [[Rcpp::export]]
arma::mat construct_ekf_terms(
    const arma::vec& observations,
    const arma::mat& stepMat,        // expMt(dampingMat, Delta/nSubsteps)
    const arma::mat& covar_mat,      // covariance for ONE substep of size Delta/nSubsteps
    const double obs_var,
    const int N,
    const double Delta,
    const double mu,
    const double A,
    const int a,
    const int B,
    const int b,
    const double C1,
    const double nuMax,
    const double nu0,
    const double r,
    const int nSubsteps = 1,
    const int burnin = 250
){
  using vec6  = arma::vec::fixed<6>;
  using mat66 = arma::mat::fixed<6,6>;
  
  const double delta    = Delta / static_cast<double>(nSubsteps);
  const double halfStep = delta / 2.0;
  
  vec6 mkCur; mkCur.zeros();
  mat66 PkCur; PkCur.eye();
  
  const double C2 = 0.8  * C1;
  const double C3 = 0.25 * C1;
  const double C4 = 0.25 * C1;
  
  arma::mat out(N, 2, arma::fill::zeros);   // col 0: nll term, col 1: standardized innovation
  
  for(int i = 0; i < N; ++i){
    
    vec6  mOut = mkCur;
    mat66 POut = PkCur;
    
    for(int j = 0; j < nSubsteps; ++j){
      vec6 mBefore = mOut;
      
      mOut = nonlinear_step(
        stepMat * nonlinear_step(mBefore, halfStep, mu, C1, C2, C3, C4, A, a, B, b, nuMax, nu0, r),
        halfStep, mu, C1, C2, C3, C4, A, a, B, b, nuMax, nu0, r);
      
      mat66 Fx = full_step_derivative(mBefore, halfStep, stepMat, mu, C1, C2, C3, C4, A, a, B, b, nuMax, nu0, r);
      arma::mat Fq = random_derivatives(mBefore, halfStep, stepMat, mu, C1, C2, C3, C4, A, a, B, b, nuMax, nu0, r);
      
      POut = Fx * POut * Fx.t() + Fq * covar_mat * Fq.t();
      POut = 0.5 * (POut + POut.t());
    }
    
    double vk = observations(i) - path_to_observation(mOut);
    double sk = as_scalar(H.t() * POut * H) + obs_var;
    
    vec6 Kk = POut * H / sk;
    mkCur = mOut + Kk * vk;
    
    mat66 joseph_form; joseph_form.eye();
    joseph_form -= Kk * H.t();
    PkCur = joseph_form * POut * joseph_form.t() + Kk * obs_var * Kk.t();
    PkCur = 0.5 * (PkCur + PkCur.t());
    
    if(i >= burnin){
      out(i, 0) = 0.5 * (std::log(2.0 * M_PI * sk) + vk * vk / sk);
    }
    out(i, 1) = vk / std::sqrt(sk);
  }
  return out;
}