#include "fem.h"
#include "types.h"

// may want to generalize later to multiple devices...

// sets up mesh-based occa parameters
void setupOccaMesh2d(Mesh *mesh, App *app){

  app->device.setup("mode: 'Serial'");
  // app->device.setup("mode: 'CUDA', deviceID: 0");

  app->props = occa::getKernelProperties();

  app->props["defines/p_Np"] = mesh->Np; // number of dims
  app->props["defines/p_Np1"] = mesh->N + 1;  // (N+1)
  
  app->props["defines/p_Nfaces"] = mesh->Nfaces;
  app->props["defines/p_Nfp"] = mesh->Nfp;
  app->props["defines/p_NfpNfaces"] = mesh->Nfp * mesh->Nfaces;  

  int Nvgeo = 5; 
  int Nfgeo = 3; 
  app->props["defines/p_Nvgeo"] = Nvgeo;
  app->props["defines/p_Nfgeo"] = Nfgeo;

  // switch dfloat type (double/float) in types.h
  if (sizeof(dfloat)==4){
    app->props["defines/USE_DOUBLE"] = 0;
  }else{
    app->props["defines/USE_DOUBLE"] = 1;
  }

  // allocate space for geofaces
  int K = mesh->K;
  int Np = mesh->Np;
  int NfpNfaces = (mesh->Nfp)*(mesh->Nfaces);  
  MatrixXd vgeo(Nvgeo*(mesh->Vq.rows()),K);

  // interpolate to quad pts and store
  MatrixXd Vq   = mesh->Vq;
  MatrixXd rxJq = Vq*(mesh->rxJ);
  MatrixXd sxJq = Vq*(mesh->sxJ);
  MatrixXd ryJq = Vq*(mesh->ryJ);
  MatrixXd syJq = Vq*(mesh->syJ);
  MatrixXd Jq   = Vq*(mesh->J);
  
  vgeo << rxJq,
    ryJq,
    sxJq,
    syJq,
    Jq;	  

  // normals computed at face quad points already
  MatrixXd fgeo(Nfgeo*(mesh->Vf.rows()),K);
  fgeo << mesh->nxJ,
    mesh->nyJ,
    mesh->sJ;

  MatrixXd rhs(Np,K), rhsf(NfpNfaces,K), res(Np,K);
  rhs.fill(0.0);
  rhsf.fill(0.0);  
  res.fill(0.0);  

  setOccaArray(app, rhs, app->o_rhs);
  setOccaArray(app, rhsf,  app->o_rhsf);
  setOccaArray(app, res, app->o_res);  
  
  setOccaArray(app, vgeo, app->o_vgeo);
  setOccaArray(app, fgeo, app->o_fgeo);
  setOccaIntArray(app, mesh->mapPq, app->o_mapPq);

  setOccaArray(app, mesh->D1D, app->o_D1D);
  setOccaArray(app, mesh->Vf1D, app->o_Vf1D);  
  setOccaArray(app, mesh->wq1D, app->o_wq1D);
}
  
// set occa array:: cast to dfloat
void setOccaArray(App *app, MatrixXd A, occa::memory &c_A){
  int r = A.rows();
  int c = A.cols();
  dfloat *f_A = (dfloat*)malloc(r*c*sizeof(dfloat));
  Map<MatrixXdf >(f_A,r,c) = A.cast<dfloat>();
  c_A = app->device.malloc(r*c*sizeof(dfloat),f_A);
  free(f_A);
}

// set occa array for int matrices
void setOccaIntArray(App *app, MatrixXi A, occa::memory &c_A){
  int r = A.rows();
  int c = A.cols();
  int *f_A = (int*)malloc(r*c*sizeof(int));
  Map<MatrixXi >(f_A,r,c) = A;
  c_A = app->device.malloc(r*c*sizeof(int),f_A);
  free(f_A);
}

// get occa array:: cast to double
void getOccaArray(App *app, occa::memory c_A, MatrixXd &A){
  int r = A.rows();
  int c = A.cols();
  dfloat *f_A = (dfloat*)malloc(r*c*sizeof(dfloat));
  c_A.copyTo(f_A);
  for (int i = 0; i < r*c; ++i){
    A(i) = f_A[i]; // somewhat inefficient
  }

  free(f_A);
}
		 

