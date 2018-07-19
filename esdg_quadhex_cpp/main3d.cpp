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

/*
    rho.fill(4.0);
    rhou.fill(.1);
    rhov.fill(.25);
    rhow.fill(.5);
    E.fill(2.0);  
    */

  /*
    rho = 2.0 + .5*(PI*x.array()).sin()*(PI*y.array()).sin()*(PI*z.array()).sin();
    rhou = 1.0 + .5*(PI*x.array()).sin()*(PI*y.array()).sin()*(PI*z.array()).sin();
    rhov = .0 + .5*(PI*x.array()).sin()*(PI*y.array()).sin()*(PI*z.array()).sin();
    rhow = -1.0 + .5*(PI*x.array()).sin()*(PI*y.array()).sin()*(PI*z.array()).sin();
    E = 2.0 + .5*(PI*x.array()).sin()*(PI*y.array()).sin()*(PI*z.array()).sin();
  
    int K = x.cols();
    int Np = x.rows();  
    for (int e = 0; e < K; ++e){
    for (int i = 0; i < Np; ++i){
    if (y(i,e) > 10.0){
    rho(i,e) += 2.0;
    rhou(i,e) += 1.0;
    rhov(i,e) += -1.0;
    rhow(i,e) += 1.0;		
    E(i,e) += 1.0;
    }
    }
    }
  */
}

static void TaylorGreen3d(MatrixXd x, MatrixXd y, MatrixXd z, double t,
			  MatrixXd &rho, MatrixXd &rhou, MatrixXd &rhov, MatrixXd &rhow, MatrixXd &E){

  // 3D vortex on [0,10] x [0,20] x [0,10]
  rho.resize(x.rows(),x.cols());
  rhou.resize(x.rows(),x.cols());
  rhov.resize(x.rows(),x.cols());
  rhow.resize(x.rows(),x.cols());  
  E.resize(x.rows(),x.cols());

  double gm1 = GAMMA - 1.0;
  MatrixXd u =  x.array().sin() * y.array().cos() * z.array().cos();
  MatrixXd v = -x.array().cos() * y.array().sin() * z.array().cos();
  MatrixXd w = MatrixXd::Zero(x.rows(),x.cols());
  MatrixXd p = 100.0/GAMMA + (1.0/16.0) * ((2.0*x).array().cos() + (2.0*y).array().cos()) * (2.0 + (2.0*z.array()).cos());

  rho.fill(1.0);
  rhou = rho.array()*u.array();
  rhov = rho.array()*v.array();
  rhow = rho.array()*w.array();  
  E = p.array() / (GAMMA-1.0) + .5*rho.array()*(u.array().square() + v.array().square() + w.array().square());
}

int main(int argc, char **argv){

  //occa::printModeInfo();

  int N = 3;
  int K1D = 8;
  double CFL = 1.0; 
  double FinalTime = 1.0;
  double a = .50;  
  if (argc > 2){
    N = atoi(argv[1]);
    K1D = atoi(argv[2]);
    printf("setting N = %d, K1D = %d\n",N,K1D);
  }
  if (argc > 5){
    CFL = atof(argv[3]);
    FinalTime = atof(argv[4]);
    a = atof(argv[5]);
    printf("setting CFL = %f, T = %f, curved warping a =%f\n",CFL,FinalTime,a);
  }

  
  Mesh *mesh = (Mesh*) calloc(1, sizeof(Mesh));  

#define TEST 0
#define VORTEX 1
#define TAYLOR_GREEN 0

#if TEST
  HexMesh3d(mesh,K1D,K1D,K1D); // make Cartesian mesh  
#endif
  
#if VORTEX   // isentropic vortex
  // [0,10] x [0,20] x [0,10] for vortex  
  HexMesh3d(mesh,K1D,2*K1D,K1D); // make Cartesian mesh
  double Lx = 5;
  double Ly = 10;
  double Lz = 5;  
  mesh->VX = (mesh->VX.array()+1.0)*Lx;
  mesh->VY = (mesh->VY.array()+1.0)*Ly;
  mesh->VZ = (mesh->VZ.array()+1.0)*Lz;  
#endif  

#if TAYLOR_GREEN   // Taylor Green on [-pi,pi]^3
  HexMesh3d(mesh,K1D,K1D,K1D); // uniform Cartesian mesh
  mesh->VX = (mesh->VX.array())*PI;
  mesh->VY = (mesh->VY.array())*PI;
  mesh->VZ = (mesh->VZ.array())*PI;  
#endif
  
  InitRefData3d(mesh, N);
  int dim = 3;
  ConnectElems(mesh,dim);  
  MapNodes3d(mesh); // low order mapping of nodes

  MatrixXd x = mesh->x;
  MatrixXd y = mesh->y;
  MatrixXd z = mesh->z;

#if TEST
  MatrixXd d = (.125*((x+y+z).array())).exp();
  x = x + a*d;
  y = y + a*d;
  z = z + a*d;
  mesh->x = x;
  mesh->y = y;
  mesh->z = z;
#endif
  
#if VORTEX

  MatrixXd xx = PI*(x.array()-5.0)/5.0;
  MatrixXd yy = 2.0*PI*(y.array()-10.0)/10.0;
  MatrixXd zz = PI*(z.array()-5.0)/5.0;
  MatrixXd d = xx.array().sin()*yy.array().sin()*zz.array().sin();
  MatrixXd dx = d;
  MatrixXd dy = d;
  MatrixXd dz = d;  
  /*
  MatrixXd xx1 = .5*PI*(x.array()-5.0)/5.0;
  MatrixXd yy1 = 1.5*PI*(y.array()-10.0)/10.0;
  MatrixXd zz1 = .5*PI*(z.array()-5.0)/5.0;
  MatrixXd xx2 = 1.5*PI*(x.array()-5.0)/5.0;
  MatrixXd yy2 = .5*PI*(y.array()-10.0)/10.0;
  MatrixXd zz2 = 1.5*PI*(z.array()-5.0)/5.0;
  MatrixXd dx = xx1.array().cos()*yy1.array().cos()*zz1.array().cos();
  MatrixXd dy = xx2.array().cos()*yy2.array().cos()*zz2.array().cos();
  MatrixXd dz = xx1.array().cos()*yy1.array().cos()*zz1.array().cos();
  */
  x = x + a*dx;
  y = y + a*dy;
  z = z + a*dz;
  mesh->x = x;
  mesh->y = y;
  mesh->z = z;
#endif

  /*
  MatrixXd xyz(x.rows(),x.cols()*3);
  xyz << x,y,z;
  cout << "xyz = [" << endl << xyz << "];" << endl;
  return 0;
  */  
  
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
  
  //App *app = (App*) calloc(1, sizeof(App));
  App *app = new App;
  //app->device.setup("mode: 'Serial'");
  app->device.setup("mode: 'CUDA', device_id: 0");
  //app->device.setup("mode: 'OpenCL', platformID : 0, deviceID: 0");

  setupOccaMesh3d(mesh,app); // build mesh geofacs
 
  app->props["defines/p_gamma"] = GAMMA;
  app->props["defines/p_Nfields"] = Nfields;
  app->props["defines/p_tau"] = 1.0;

  MatrixXd Vq = mesh->Vq;  
  MatrixXd xq = Vq*mesh->x;
  MatrixXd yq = Vq*mesh->y;
  MatrixXd zq = Vq*mesh->z;

#if TEST
  // test 3d geofacs

  MatrixXd Dr = Vq * mrdivide(mesh->Dr,Vq);
  MatrixXd Ds = Vq * mrdivide(mesh->Ds,Vq);
  MatrixXd Dt = Vq * mrdivide(mesh->Dt,Vq);
  MatrixXd rxJ = Vq * mesh->rxJ;
  MatrixXd sxJ = Vq * mesh->sxJ;
  MatrixXd txJ = Vq * mesh->txJ;
  MatrixXd J = Vq * mesh->J;  

  MatrixXd u = xq; //xq.array().exp();
  MatrixXd duex = 1.0 + 0*xq.array();
  
  MatrixXd d1 = rxJ.array()*(Dr*u).array() +
    sxJ.array()*(Ds*u).array() +
    txJ.array()*(Dt*u).array();
  MatrixXd d2 = Dr*(u.array()*rxJ.array()).matrix() +
    Ds*(u.array()*sxJ.array()).matrix() +
    Dt*(u.array()*txJ.array()).matrix();

  MatrixXd dudx = .5*(d1 + d2).array() / J.array();

  MatrixXd geo(rxJ.rows(),rxJ.cols()*3);
  geo << rxJ, sxJ, txJ;

  //cout << "Dt = " << Dt << endl;
  cout << setprecision(11);
  //cout << "geo = [" << endl << geo << "];" << endl;
  // printf("Vq size = %d, %d\n",Vq.rows(),Vq.cols());
  cout << "uc = [" << u << "];" << endl;
  cout << "d1c = [" << d1 << "];" << endl;
  cout << "d2c = [" << d2 << "];" << endl;  
  cout << "dudxc = [" << dudx << "];" << endl;    

  printf("diff = %g\n",(dudx.array()-duex.array()).abs().sum());

  return 0;
#endif

  
  double time = 0.0;
  MatrixXd rho, rhou, rhov, rhow, E;
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

  //  occa::kernel test = app->device.buildKernel(path.c_str(),"test_exclusive",app->props);
  //  test(1);  return 0;
  
  MatrixXd wJq = mesh->wq.asDiagonal() * (mesh->Vq*mesh->J);

#if 0
  app->eval_surface(K, app->o_Vf1D,app->o_Q, app->o_Qf);
  MatrixXd Qf(Nfields*NfpNfaces,K);

  //  getOccaArray(app,app->o_Qf,Qf);
  //  MatrixXd rhof = Qf.middleRows(0,NfpNfaces);
  //  cout << "rhof = " << rhof << endl;
  
  // test rhs eval
  app->volume(K, app->o_vgeo, app->o_vfgeo,
	      app->o_D1D, app->o_Vf1D, app->o_Lf1D,
	      app->o_Q, app->o_Qf,
	      app->o_rhs, app->o_rhsf);

  getOccaArray(app,app->o_rhs,Q);
  MatrixXd rhs1  = Q.middleRows(0,Np);  
  //  cout << "vol only rhs = " << endl << rhs1 << endl;
  
  app->surface(K, app->o_vgeo, app->o_fgeo, app->o_mapPq,
	       app->o_Lf1D, app->o_Qf, app->o_rhsf,
	       app->o_rhs); 
  
  getOccaArray(app,app->o_rhs,Q);
  rhs1 = Q.middleRows(0,Np);
  //cout << "vol+surf rhs = " << endl << rhs1 << endl;

  return 0;
  
#endif 

  // ============== solver ===================

  double h = mesh->J.maxCoeff() / mesh->sJ.maxCoeff(); // J = O(h^d), Jf = O(h^{d-1}) in d dims
  double CN = dim * (double)((N+1)*(N+2))/2.0; // trace constant for GQ hexes
  double dt = CFL * h / CN;
  int Nsteps = (int) ceil(FinalTime/dt);
  dt = FinalTime/(double) Nsteps;

  printf("dt = %f, FinalTime = %f, Nsteps = %d\n",dt,FinalTime,Nsteps);  

  // interp to surface to start
  app->eval_surface(K, app->o_Vf1D,app->o_Q, app->o_Qf);
    
  int interval = ceil(Nsteps/10);
  if (Nsteps < 10){
    interval = 1;
  }
  
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
  rhow = Q.middleRows(3*Np,Np);  
  E = Q.middleRows(4*Np,Np);

  // finer quadrature for error eval
  VectorXd rq1D2, wq1D2;
  JacobiGQ(N+1, 0, 0, rq1D2, wq1D2);
  int Np1 = rq1D2.size();
  int Nq3 = Np1*Np1*Np1;
  VectorXd rq2(Nq3),sq2(Nq3),tq2(Nq3),wq2(Nq3);
  int sk = 0;
  for (int i = 0; i < Np1; ++i){
    for (int j = 0; j < Np1; ++j){
      for (int k = 0; k < Np1; ++k){
	rq2(sk) = rq1D2(i);
	sq2(sk) = rq1D2(j);
	tq2(sk) = rq1D2(k);
	wq2(sk) = wq1D2(i)*wq1D2(j)*wq1D2(k);
	++sk;
      }
    }
  } 
  MatrixXd Vqtmp = Vandermonde3DHex(N,rq2,sq2,tq2);
  MatrixXd VqN = Vandermonde3DHex(N,mesh->rq,mesh->sq,mesh->tq);  
  MatrixXd Vq2 = mrdivide(Vqtmp,VqN);
  MatrixXd xq2 = Vq2 * (mesh->Vq*mesh->x);
  MatrixXd yq2 = Vq2 * (mesh->Vq*mesh->y);
  MatrixXd zq2 = Vq2 * (mesh->Vq*mesh->z);  
  
  MatrixXd rhoex,rhouex,rhovex,rhowex,Eex;
  VortexSolution3d(xq2,yq2,zq2,FinalTime,rhoex,rhouex,rhovex,rhowex,Eex);

  MatrixXd wJq2 = wq2.asDiagonal() * (Vq2*(mesh->Vq*mesh->J));  
  MatrixXd werr = wJq2.array()*((rhoex - Vq2*rho).array().square() +
				(rhouex - Vq2*rhou).array().square() +
				(rhovex - Vq2*rhov).array().square() +
				(rhowex - Vq2*rhow).array().square() +
				(Eex - Vq2*E).array().square());
  printf("L2 error for rho = %g\n",sqrt(werr.sum()));
 
  return 0;
  
}
