#include "fem.h"
#define GAMMA 1.4

int main(int argc, char **argv){

  int N = 4;
  int K = 1e6;
  int Np = 44; // 44 for N = 4
  int Nq = Np;
  int Nfields = 4;
  
  // ========================== set up OCCA application

  App *app = new App;
  //  app->device.setup("mode: 'Serial'");
  app->device.setup("mode: 'CUDA', device_id: 0");

  // define macros for OCCA kernels
  app->props["defines/p_Nq"] = Nq;
  
  // switch dfloat type (double/float) in types.h
  if (sizeof(dfloat)==4){
    app->props["defines/USE_DOUBLE"] = 0;
  }else{
    app->props["defines/USE_DOUBLE"] = 1;
  }

  // Euler-specific values
  app->props["defines/p_gamma"] = GAMMA;
  app->props["defines/p_Nfields"] = Nfields;

    // build occa kernels  
  string path = "okl/EulerTest.okl";

  //testing
  occa::kernel volume;
  volume = app->device.buildKernel(path.c_str(),"volume",app->props);

  // =================== set initial condition

  MatrixXd rho = (MatrixXd::Random(Np,K)).array().abs();
  MatrixXd u = MatrixXd::Random(Np,K);
  MatrixXd v = MatrixXd::Random(Np,K);
  MatrixXd w = MatrixXd::Random(Np,K);
  MatrixXd p = (MatrixXd::Random(Np,K)).array().abs();

  MatrixXd rhou = rho.array()*u.array();
  MatrixXd rhov = rho.array()*v.array();
  MatrixXd rhow = rho.array()*w.array();    
  MatrixXd E = p.array()/(GAMMA-1.0) + .5*rho.array()*(u.array().square() + v.array().square() + w.array().square());

  MatrixXd Q(Nfields*Np,K);
  Q << rho,rhou,rhov,E;

  // ================= construct discretization matrices

  // fake SBP ops
  MatrixXd Qr = MatrixXd::Random(Nq,Nq);  
  MatrixXd Qs = MatrixXd::Random(Nq,Nq);
  MatrixXd Qt = MatrixXd::Random(Nq,Nq);  
  MatrixXd QNr = Qr - Qr.transpose();
  MatrixXd QNs = Qs - Qs.transpose();
  MatrixXd QNt = Qt - Qt.transpose();  

  occa::memory o_Q; // solution    
  setOccaArray(app, Q, o_Q);

  occa::memory o_rhs;  // timestep stuff
  setOccaArray(app,MatrixXd::Zero(Np*Nfields,K),o_rhs);

  occa::memory o_QNr, o_QNs, o_QNt;

  // set operators
  setOccaArray(app, QNr, o_QNr);
  setOccaArray(app, QNs, o_QNs);
  setOccaArray(app, QNt, o_QNt);  
  
  // ============== run RK solver ==================

  int Nsteps = 100;
  for (int i = 0; i < Nsteps; ++i){      
    volume(K, o_QNr, o_QNs, o_QNt, o_Q, o_rhs);
  }
  
}
