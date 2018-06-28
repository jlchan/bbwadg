#include "fem.h"

#define GAMMA 1.4f

static void VortexSolution3d(MatrixXd x, MatrixXd y, MatrixXd z, double t,
			     MatrixXd &rho, MatrixXd &rhou, MatrixXd &rhov, MatrixXd &rhow, MatrixXd &E){
  // 3D vortex on [0,10] x [0,20] x [0,10]

  rho.resize(x.rows(),x.cols());
  rhou.resize(x.rows(),x.cols());
  rhov.resize(x.rows(),x.cols());
  rhow.resize(x.rows(),x.cols());  
  E.resize(x.rows(),x.cols()); 
  
  double x0 = 5.0;
  double y0 = 5.0;
  MatrixXd xt = x.array() - x0;
  MatrixXd yt = y.array() - y0 - t;

  // cross(X,[0,0,1]) = [-y,x,0]
  MatrixXd rx = -yt;
  MatrixXd ry = xt;
  MatrixXd rz = MatrixXd::Zero(x.rows(),x.cols());
  MatrixXd r2 = rx.array().square()
    + ry.array().square() + rz.array().square();

  double rho0 = 1.0;
  double p0 = 1.0 / GAMMA;
  double Lmax = .4;
  double gm1 = GAMMA-1.0;
  
  MatrixXd L = Lmax * (.5 * (1.0 - r2.array())).exp();
  MatrixXd tmp = (1.0 - .5 * gm1 * L.array() * L.array()).pow(1.0 / gm1);
  rho = rho0 * tmp;
  rhou = rho.array()*(0.0 + rx.array() * L.array());
  rhov = rho.array()*(1.0 + ry.array() * L.array());
  rhow = rho.array()*(0.0 + rz.array() * L.array());
  E = p0 / gm1 * (1.0 + tmp.array().pow(GAMMA)) +
    .5 * (rhou.array().square() + rhov.array().square() + rhow.array().square())/rho.array();  
}


int main(int argc, char **argv){

  //occa::printModeInfo();

  int N = 5;
  int K1D = 4;
  
  Mesh *mesh = (Mesh*) calloc(1, sizeof(Mesh));  
  HexMesh3d(mesh,K1D,K1D,K1D); // make Cartesian mesh

#if 1
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

  MatrixXd xq = mesh->Vq*mesh->x;
  MatrixXd yq = mesh->Vq*mesh->y;
  MatrixXd zq = mesh->Vq*mesh->z;

  MatrixXd rho, rhou, rhov, rhow, E;
  double time = 0.0;
  VortexSolution3d(xq,yq,zq,time,rho,rhou,rhov,rhow,E);

  MatrixXd Q(Nfields*Np,K);  
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
  

  // ============== solver ===================

  double h = mesh->J.maxCoeff() / mesh->sJ.maxCoeff(); // J = O(h^d), Jf = O(h^{d-1}) in d dims
  double CN = dim * (double)((N+1)*(N+2))/2.0; // trace constant for GQ hexes
  double CFL = .25;
  double dt = CFL * h / CN;
  
  double FinalTime = 5.0;
  int Nsteps = (int) ceil(FinalTime/dt);
  dt = FinalTime/(double) Nsteps;

  printf("dt = %f, FinalTime = %f, Nsteps = %d\n",dt,FinalTime,Nsteps);  

  // interp to surface to start
  app->eval_surface(K, app->o_Vf1D,app->o_Q, app->o_Qf);
    
  int interval = ceil(Nsteps/10);
  int NINT = mesh->rk4a.size();
  for (int i = 0; i < Nsteps; ++i){
    for (int INTRK = 0; INTRK < NINT; ++INTRK){

      const dfloat fdt = (dfloat) dt;
      const dfloat fa  = (dfloat) mesh->rk4a[INTRK];
      const dfloat fb  = (dfloat) mesh->rk4b[INTRK];

      app->volume(K,app->o_vgeo, app->o_vfgeo,
		  app->o_D1D, app->o_Vf1D, app->o_Lf1D,
		  app->o_Q,app->o_Qf,app->o_rhs,app->o_rhsf);
      
      app->surface(K,app->o_vgeo,app->o_fgeo,app->o_mapPq,app->o_Lf1D,
		   app->o_Qf,app->o_rhsf,app->o_rhs);
      
      app->update(K, fa, fb, fdt,
		  app->o_Q, app->o_rhs, app->o_res);
      
      app->eval_surface(K, app->o_Vf1D,app->o_Q, app->o_Qf);
    }
    
    if (i % interval == 0){
      printf("on timestep %d out of %d\n",i,Nsteps);
    }
  }
  getOccaArray(app,app->o_Q,Q);
  rho = Q.middleRows(0,Np);
  rhou = Q.middleRows(Np,Np);
  rhov = Q.middleRows(2*Np,Np);
  E = Q.middleRows(3*Np,Np);    

  MatrixXd rhoex,rhouex,rhovex,rhowex,Eex;
  VortexSolution3d(xq,yq,zq,FinalTime,rhoex,rhouex,rhovex,rhowex,Eex);

#if 0
  cout << "xq = [" << xq << "];" << endl;
  cout << "yq = [" << yq << "];" << endl;
  cout << "rho = [" << rho << "];" << endl;
  cout << "rhoex = [" << rhoex << "];" << endl;
  return 0;
#endif

  MatrixXd wJq = mesh->wq.asDiagonal() * (mesh->Vq*mesh->J);  
  //  MatrixXd rhouex = rhoex.array()*u.array();
  //  MatrixXd rhovex = rhoex.array()*v.array();
  //  MatrixXd Eex = p.array()/(GAMMA-1.0) + .5*rhoex.array()*(u.array().square() + v.array().square());
  MatrixXd werr = wJq.array()*(rhoex - rho).array().square();
  printf("L2 error for rho = %g\n",sqrt(werr.sum()));
 
  return 0;
  
}
