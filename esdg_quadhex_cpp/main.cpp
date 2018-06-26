#include "fem.h"

#define GAMMA 1.4f

static void VortexSolution2d(MatrixXd x,MatrixXd y,double t,
		      MatrixXd &rho, MatrixXd &u, MatrixXd &v, MatrixXd &p){
  double x0 = 5.0;
  double y0 = 0.0;
  double beta = 5.0;
  MatrixXd r2 = (x.array() - x0 - t).square() + (y.array()-y0).square();

  rho.resize(x.rows(),x.cols());
  u.resize(x.rows(),x.cols());
  v.resize(x.rows(),x.cols());
  p.resize(x.rows(),x.cols());

  MatrixXd exp_r2 = (1.0-r2.array()).exp();
  u   = 1.0 - beta*exp_r2.array()*(y.array()-y0)/(2.0*PI);
  v   = beta*exp_r2.array()*(x.array()-x0-t)/(2.0*PI);
  rho = 1.0 - (1.0/(8.0*GAMMA*PI*PI))*(GAMMA-1.0)/2.0*(beta*exp_r2.array()).square();
  rho = rho.array().pow(1.0/(GAMMA-1.0));
  p   = rho.array().pow(GAMMA);

}

int main(int argc, char **argv){

  //occa::printModeInfo();

  int N = 4;
  int K1D = 8;
  
  Mesh *mesh = (Mesh*) calloc(1, sizeof(Mesh));  
  QuadMesh2d(mesh,2*K1D,K1D); // make Cartesian mesh

  // [0,20] x [-5,5] for vortex
  double Lx = 10;
  double Ly = 5;  
  mesh->VX = (mesh->VX.array()+1.0)*Lx;
  mesh->VY = (mesh->VY.array())*Ly;

  // ============ physics independent stuff ===========

  InitRefData2d(mesh, N);
  GeometricFactors2d(mesh);
  Normals2d(mesh); 

  int dim = 2;
  ConnectElems(mesh,dim);

  MatrixXd xf = (mesh->Vf)*(mesh->x);
  MatrixXd yf = (mesh->Vf)*(mesh->y);
  MatrixXd zf(xf.rows(),xf.cols()); zf.fill(0.0);
  MatrixXi mapPq;
  BuildFaceNodeMaps(mesh,xf,yf,zf,mapPq);

  double DX = mesh->VX.maxCoeff()-mesh->VX.minCoeff();
  double DY = mesh->VY.maxCoeff()-mesh->VY.minCoeff();
  MakeNodeMapsPeriodic2D(mesh,xf,yf,DX,DY,mapPq);
  mesh->mapPq = mapPq;

  //  cout << "xf = [" << xf << "];" << endl;
  //  cout << "yf = [" << yf << "];" << endl;
  //  cout << "mapPq = [" << mapPq << "];" << endl;
  //  return 0;
  

  // ============ problem dependent stuff ============

  int Np = mesh->Np;
  int NfpNfaces = mesh->Nfp * mesh->Nfaces;  
  int K = mesh->K;

  int Nfields = 4; 
  mesh->Nfields = Nfields;
  
  App *app = (App*) calloc(1, sizeof(App));
  setupOccaMesh2d(mesh,app); // build mesh geofacs

  app->props["defines/p_gamma"] = GAMMA;
  app->props["defines/p_Nfields"] = Nfields;
  app->props["defines/p_tau"] = 1.0;  
  
  // build occa kernels
  string path = "okl/Euler2D.okl";
  app->volume = app->device.buildKernel(path.c_str(),"volume",app->props);
  app->surface = app->device.buildKernel(path.c_str(),"surface",app->props);
  app->update = app->device.buildKernel(path.c_str(),"update",app->props);
  app->eval_surface = app->device.buildKernel(path.c_str(),"eval_surface",app->props);  

  // set initial condition
  MatrixXd xq = mesh->Vq*mesh->x;
  MatrixXd yq = mesh->Vq*mesh->y;  

  // vortex
  MatrixXd rho,u,v,p;
  double time = 0.0;
  VortexSolution2d(xq,yq,time,rho,u,v,p);

  MatrixXd rhou = rho.array()*u.array();
  MatrixXd rhov = rho.array()*v.array();
  MatrixXd E = p.array()/(GAMMA-1.0) + .5*rho.array()*(u.array().square() + v.array().square());

  /*
  // const sol
  rho.fill(1.0);
  rhou.fill(1.0);
  rhov.fill(1.0);
  E.fill(2.0);  

  //testing
  rho.col(0).fill(1.0);
  rho.col(1).fill(2.0);  
  */

  MatrixXd Q(Nfields*Np,K);
  Q << rho,
    rhou,
    rhov,
    E;

  setOccaArray(app, Q, app->o_Q);


  // ============== run RK solver ==================

  double h = 2.0 / (double) K1D; //min(DX,DY) / (double) K1D;
  double CN = dim * (double)((N+1)*(N+2))/2.0; // trace constant for GQ hexes
  double CFL = .5;
  double dt = CFL * h / CN;
  
  double FinalTime = 1.0;
  int Nsteps = (int) ceil(FinalTime/dt);
  dt = FinalTime/(double) Nsteps;

  printf("dt = %f, FinalTime = %f, Nsteps = %d\n",dt,FinalTime,Nsteps);  

  // interp to surface to start
  app->eval_surface(K, app->o_Vf1D,app->o_Q, app->o_Qf);

#if 0
  MatrixXd Qf(Nfields*NfpNfaces,K);
  getOccaArray(app,app->o_Q,Q);
  getOccaArray(app,app->o_Qf,Qf);
  cout << "Q = " << endl << Q << endl;
  cout << "Qf = " << endl << Qf << endl; 
  return 0;
#endif
  
#if 0

  // test rhs eval
  app->volume(K, app->o_vgeo, app->o_vfgeo,
	      app->o_D1D, app->o_Vf1D, app->o_Lf1D,
	      app->o_Q, app->o_Qf,
	      app->o_rhs, app->o_rhsf);

  getOccaArray(app,app->o_rhs,Q);
  //cout << "vol only rhs = " << endl << Q << endl;
  
  app->surface(K, app->o_fgeo, app->o_mapPq,
	       app->o_Lf1D, app->o_Qf, app->o_rhsf,
	       app->o_rhs); 
  
  getOccaArray(app,app->o_rhs,Q);
  //  cout << "rhs = " << endl << Q << endl;
  
  return 0;
#endif 
  
  int interval = ceil(Nsteps/10);
  int NINT = mesh->rk4a.size();
  for (int i = 0; i < Nsteps; ++i){
    for (int INTRK = 0; INTRK < NINT; ++INTRK){

      const dfloat fdt = (dfloat) dt;
      const dfloat fa  = (dfloat) mesh->rk4a[INTRK];
      const dfloat fb  = (dfloat) mesh->rk4b[INTRK];

      app->volume(K, app->o_vgeo, app->o_vfgeo,
		  app->o_D1D, app->o_Vf1D, app->o_Lf1D,
		  app->o_Q, app->o_Qf,
		  app->o_rhs, app->o_rhsf);

      app->surface(K, app->o_fgeo, app->o_mapPq,
		   app->o_Lf1D, app->o_Qf, app->o_rhsf,
		   app->o_rhs); 
      
      //      getOccaArray(app,app->o_rhs,Q);
      //      cout << "rhs = " << endl << Q << endl;
      
      app->update(K, fa, fb, fdt, app->o_vgeo, app->o_Vf1D, 
		  app->o_Q, app->o_Qf, app->o_rhs, app->o_res);      
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

  MatrixXd rhoex;
  VortexSolution2d(xq,yq,FinalTime,rhoex,u,v,p);

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
