#include "fem.h"

#define GAMMA 1.4f

int main(int argc, char **argv){

  //occa::printModeInfo();

  int N = 2;
  int K1D = 2;
  
  Mesh *mesh = (Mesh*) calloc(1, sizeof(Mesh));  
  HexMesh3d(mesh,K1D,K1D,K1D); // make Cartesian mesh

#if 0
  // [0,20] x [-5,5] for vortex
  double Lx = 5;
  double Ly = 10;
  double Lz = 5;  
  mesh->VX = (mesh->VX.array()+1.0)*Lx;
  mesh->VY = (mesh->VY.array()+1.0)*Ly;
  mesh->VZ = (mesh->VZ.array()+1.0)*Lz;  
#endif
  
  InitRefData3d(mesh, N);

  int dim = 3;
  ConnectElems(mesh,dim);  

  MapNodes3d(mesh); // low order mapping
  GeometricFactors3d(mesh);
  Normals3d(mesh);
  
  MatrixXd xf = (mesh->Vf)*(mesh->x);
  MatrixXd yf = (mesh->Vf)*(mesh->y);
  MatrixXd zf = (mesh->Vf)*(mesh->z);  
  MatrixXi mapPq;
  BuildFaceNodeMaps(mesh,xf,yf,zf,mapPq);
  
  double DX = mesh->VX.maxCoeff()-mesh->VX.minCoeff();
  double DY = mesh->VY.maxCoeff()-mesh->VY.minCoeff();
  double DZ = mesh->VZ.maxCoeff()-mesh->VZ.minCoeff();
  MakeNodeMapsPeriodic3d(mesh,xf,yf,zf,DX,DY,DZ,mapPq);
  mesh->mapPq = mapPq;  

  // ============ problem dependent stuff ============

  int Np = mesh->Np;
  int NfpNfaces = mesh->Nfp * mesh->Nfaces;  
  int K = mesh->K;

  int Nfields = 5; 
  mesh->Nfields = Nfields;
  
  App *app = (App*) calloc(1, sizeof(App));
  setupOccaMesh3d(mesh,app); // build mesh geofacs

  app->props["defines/p_gamma"] = GAMMA;
  app->props["defines/p_Nfields"] = Nfields;
  app->props["defines/p_tau"] = 1.0;

  MatrixXd Q(Nfields*Np,K);
  MatrixXd rho(Np,K),rhou(Np,K),rhov(Np,K),rhow(Np,K),E(Np,K);
  rho.fill(2.0);
  rhou.fill(.10);
  rhov.fill(.20);
  rhow.fill(.30);
  E.fill(2.0);
  Q << rho,
    rhou,
    rhov,
    rhow,
    E;
  setOccaArray(app, Q, app->o_Q); 

  // build occa kernels
  string path = "okl/Euler3D.okl";
  app->volume = app->device.buildKernel(path.c_str(),"volume",app->props);
  app->surface = app->device.buildKernel(path.c_str(),"surface",app->props);
  app->update = app->device.buildKernel(path.c_str(),"update",app->props);
  app->eval_surface = app->device.buildKernel(path.c_str(),"eval_surface",app->props);  
  
#if 1
  // testing
  app->eval_surface(K,app->o_Vf1D,app->o_Q,app->o_Qf);

  app->volume(K,app->o_vgeo, app->o_vfgeo,
	      app->o_D1D, app->o_Vf1D, app->o_Lf1D,
	      app->o_Q,app->o_Qf,app->o_rhs,app->o_rhsf);
  
  app->surface(K,app->o_vgeo,app->o_fgeo,app->o_mapPq,app->o_Lf1D,
	       app->o_Qf,app->o_rhsf,app->o_rhs);
  
  MatrixXd Qf(Nfields*NfpNfaces,K);
  getOccaArray(app,app->o_Q,Q);
  getOccaArray(app,app->o_Qf,Qf);
  MatrixXd rhs(Nfields*Np,K);    
  getOccaArray(app,app->o_rhs,rhs);  
  //  cout << "Q = " << endl << Q << endl;
  //cout << "Qf = " << endl << Qf << endl;
  cout << "rhs = " << endl << rhs << endl;   
  return 0;
#endif
  
  // set initial condition  
  MatrixXd xq = mesh->Vq*mesh->x;
  MatrixXd yq = mesh->Vq*mesh->y;
  MatrixXd zq = mesh->Vq*mesh->z;
  /*cout << "xq = [" << xq << "];" << endl;
  cout << "yq = [" << yq << "];" << endl;
  cout << "zq = [" << zq << "];" << endl; */

  
  return 0;
  
}
