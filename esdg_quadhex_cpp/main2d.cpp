#include "fem.h"

#define GAMMA 1.4

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

  /*  
  rho.fill(2.0);
  u.fill(.1);
  v.fill(.2);
  p.fill(1.0);  
  */
}

int main(int argc, char **argv){

  //occa::printModeInfo();   return 0;
  int N = 3;
  int K1D = 16;

  //printf("argc = %d\n",argc);
  if (argc == 3){
    N = atoi(argv[1]);
    K1D = atoi(argv[2]);
    printf("setting N = %d, K1D = %d\n",N,K1D);
  }

  double FinalTime = 5.0;
  double CFL = .5;  
  double a = .25; // curved warping

  if (argc > 5){
    CFL = atof(argv[3]);
    FinalTime = atof(argv[4]);
    a = atof(argv[5]);
    printf("setting CFL = %f, T = %f, curved warping a =%f\n",CFL,FinalTime,a);
  }
  
  //Mesh *mesh = (Mesh*) calloc(1, sizeof(Mesh));
  Mesh *mesh = new Mesh;
  QuadMesh2d(mesh,2*K1D,K1D); // make Cartesian mesh

  // [0,20] x [-5,5] for vortex
  double Lx = 10;
  double Ly = 5;  
  mesh->VX = (mesh->VX.array()+1.0)*Lx;
  mesh->VY = (mesh->VY.array())*Ly;

  // ============ physics independent stuff ===========

  InitRefData2d(mesh, N);
  MapNodes2d(mesh); // low order mapping

  // curvilinear meshing
  MatrixXd x = mesh->x;
  MatrixXd y = mesh->y;
  MatrixXd dx = (.5*PI*(x.array()-Lx)/Lx).cos()*(1.5*PI*y.array()/Ly).cos();
  x = x + a*Lx*dx;
  MatrixXd dy = (2.0*PI*(x.array()-Lx)/Lx).sin()*(.5*PI*y.array()/Ly).cos();
  y = y + a*Ly*dy;

  mesh->x = x;
  mesh->y = y;  

  /*
  cout << "x = [" << mesh->x << "];" << endl;
  cout << "y = [" << mesh->y << "];" << endl;
  return 0;
  */
  
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
  MakeNodeMapsPeriodic2d(mesh,xf,yf,DX,DY,mapPq);
  mesh->mapPq = mapPq;
  //return 0;
  /*
  cout << "xf = [" << xf << "];" << endl;
  cout << "yf = [" << yf << "];" << endl;
  cout << "mapPq = [" << mapPq << "];" << endl;
  return 0;
  */

  // ============ problem dependent stuff ============

  int Np = mesh->Np;
  int NfpNfaces = mesh->Nfp * mesh->Nfaces;  
  int K = mesh->K;

  int Nfields = 4; 
  mesh->Nfields = Nfields;
  
  //App *app = (App*) calloc(1, sizeof(App));
  App *app = new App;
  //app->device.setup("mode: 'Serial'");
  app->device.setup("mode: 'CUDA', device_id: 0");  
  setupOccaMesh2d(mesh,app); // build mesh geofacs

  app->props["defines/p_gamma"] = GAMMA;
  app->props["defines/p_Nfields"] = Nfields;
  if (sizeof(dfloat)==4){
    app->props["defines/p_tau"] = 1.f;
  }else{
    app->props["defines/p_tau"] = 1.0;
  }
  
  // build occa kernels  
  string path = "okl/Euler2D.okl";

  //testing
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

  //double h = 2.0 / (double) K1D; //min(DX,DY) / (double) K1D;
  double h = mesh->J.maxCoeff() / mesh->sJ.maxCoeff(); // J = O(h^d), Jf = O(h^{d-1}) in d dims
  double CN = dim * (double)((N+1)*(N+2))/2.0; // trace constant for GQ hexes
  double dt = CFL * h / CN;
  
  int Nsteps = (int) ceil(FinalTime/dt);
  dt = FinalTime/(double) Nsteps;

  printf("dt = %f, FinalTime = %f, Nsteps = %d\n",dt,FinalTime,Nsteps);  

  // interp to surface to start
  app->eval_surface(K, app->o_Vf1D,app->o_Q, app->o_Qf);

#if 0
  MatrixXd Qf(Nfields*NfpNfaces,K);
  getOccaArray(app,app->o_Q,Q);
  getOccaArray(app,app->o_Qf,Qf);
  /*cout << "xq = [" << endl << xq <<  "];" << endl;
  cout << "yq = [" << endl << yq << "];" << endl;
  cout << "xf = [" << endl << xf <<  "];" << endl;
  cout << "yf = [" << endl << yf << "];" << endl;*/
  cout << "rho = [" << endl << Q.middleRows(0,Np) <<  "];" << endl;
  cout << "rhof = [" << endl << Qf.middleRows(0,NfpNfaces) << "];" << endl; 
  return 0;
#endif
  
#if 0

  // test rhs eval
  app->volume(K, app->o_vgeo, app->o_vfgeo,
	      app->o_D1D, app->o_Vf1D, app->o_Lf1D,
	      app->o_Q, app->o_Qf,
	      app->o_rhs, app->o_rhsf);

  getOccaArray(app,app->o_rhs,Q);
  cout << "vol only rhs = " << endl << Q.middleRows(0,Np) << endl;
  
  app->surface(K, app->o_vgeo, app->o_fgeo, app->o_mapPq,
	       app->o_Lf1D, app->o_Qf, app->o_rhsf,
	       app->o_rhs); 
  
  getOccaArray(app,app->o_rhs,Q);
    cout << "vol+surf rhs = " << endl << Q.middleRows(0,Np) << endl;
  
  return 0;
#endif 
  
  int interval = max((int) ceil(Nsteps/10),1);  
  printf("Interval = %d\n",interval);

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

      app->surface(K, app->o_vgeo, app->o_fgeo, app->o_mapPq,
		   app->o_Lf1D, app->o_Qf, app->o_rhsf,
		   app->o_rhs); 
      
      //      getOccaArray(app,app->o_rhs,Q);
      //      cout << "rhs = " << endl << Q << endl;
      
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

  // finer quadrature for error eval
  VectorXd rq1D2, wq1D2;
  JacobiGQ(N+1, 0, 0, rq1D2, wq1D2);
  int Np1 = rq1D2.size();
  int Np2 = Np1*Np1;
  VectorXd rq2(Np2),sq2(Np2),wq2(Np2);
  int sk = 0;
  for (int i = 0; i < Np1; ++i){
    for (int j = 0; j < Np1; ++j){
      rq2(sk) = rq1D2(i);
      sq2(sk) = rq1D2(j);      
      wq2(sk) = wq1D2(i)*wq1D2(j);
      ++sk;
    }
  }
  MatrixXd Vqtmp = Vandermonde2DQuad(N,rq2,sq2);
  MatrixXd Vq = Vandermonde2DQuad(N,mesh->rq,mesh->sq);  
  MatrixXd Vq2 = mrdivide(Vqtmp,Vq);
  MatrixXd xq2 = Vq2 * (mesh->Vq*mesh->x);
  MatrixXd yq2 = Vq2 * (mesh->Vq*mesh->y);

  //cout << "Vq2 = "  << Vq2 << endl;
  
  MatrixXd rhoex;
  VortexSolution2d(xq2,yq2,FinalTime,rhoex,u,v,p);

#if 0
  cout << "xq = [" << xq << "];" << endl;
  cout << "yq = [" << yq << "];" << endl;
  cout << "rho = [" << rho << "];" << endl;
  cout << "rhoex = [" << rhoex << "];" << endl;
  return 0;
#endif

  MatrixXd wJq = wq2.asDiagonal() * (Vq2*(mesh->Vq*mesh->J));  
  MatrixXd rhouex = rhoex.array()*u.array();
  MatrixXd rhovex = rhoex.array()*v.array();
  MatrixXd Eex = p.array()/(GAMMA-1.0) + .5*rhoex.array()*(u.array().square() + v.array().square());

  MatrixXd werr = wJq.array()*((rhoex - Vq2*rho).array().square() +
			       (rhouex - Vq2*rhou).array().square() +
			       (rhovex - Vq2*rhov).array().square() +
			       (Eex - Vq2*E).array().square());
  printf("L2 error for rho = %g\n",sqrt(werr.sum()));

  return 0;
  
}
