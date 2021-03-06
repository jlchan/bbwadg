#include "fem.h"

#define GAMMA 1.4

// 3D vortex on [0,10] x [0,20] x [0,10]
static void VortexSolution3d(MatrixXd x, MatrixXd y, MatrixXd z, double t,
			     MatrixXd &rho, MatrixXd &rhou, MatrixXd &rhov, MatrixXd &rhow, MatrixXd &E){
  rho.resize(x.rows(),x.cols());
  rhou.resize(x.rows(),x.cols());
  rhov.resize(x.rows(),x.cols());
  rhow.resize(x.rows(),x.cols());  
  E.resize(x.rows(),x.cols()); 
  
  double x0 = 7.5;
  double y0 = 7.5;
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

  // DEBUGGING
#if 0
  rho.fill(4.0);    rhou.fill(.1);    rhov.fill(.25);    rhow.fill(.5);    E.fill(2.0);  
#endif
  
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
  E = p.array() / (GAMMA-1.0) + .5*rho.array()*(u.array()*u.array() + v.array()*v.array() + w.array()*w.array());

  /*
  rho.fill(2.0);
  rhou.fill(.50);
  rhov.fill(.20);
  rhow.fill(.10);
  E.fill(1.0);
  */
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

#define VORTEX 1
#define TAYLOR_GREEN 0

#if VORTEX   // isentropic vortex
  printf("Running isentropic vortex\n");
  // [0,15] x [0,20] x [0,5] for vortex  
  HexMesh3d(mesh,3*K1D,4*K1D,K1D); // make Cartesian mesh
  double Lx = 15;
  double Ly = 20;
  double Lz = 5;  
  mesh->VX = .5*(mesh->VX.array()+1.0)*Lx;
  mesh->VY = .5*(mesh->VY.array()+1.0)*Ly;
  mesh->VZ = .5*(mesh->VZ.array()+1.0)*Lz;  
#endif  

#if TAYLOR_GREEN   // Taylor Green on [-pi,pi]^3
  printf("Running Taylor-Green\n");  
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

  // apply curved warping
#if VORTEX
  MatrixXd xx = (x.array()-.5*Lx)/Lx;
  MatrixXd yy = (y.array()-.5*Ly)/Ly;
  MatrixXd zz = (z.array()-.5*Lz)/Lz;
  MatrixXd dy = (3.*PI*xx).array().cos() * (PI*yy).array().cos() * (PI*zz).array().cos();
  y = y + a*Ly*dy;
  
  yy = (y.array()-.5*Ly)/Ly;    
  MatrixXd dx = (PI*xx).array().cos() * (4.*PI*yy).array().sin() * (PI*zz).array().cos();
  x = x + a*Lx*dx;

  xx = (x.array()-.5*Lx)/Lx;
  MatrixXd dz = (3.*PI*xx).array().cos() * (4.*PI*yy).array().sin() * (PI*zz).array().cos();
  //MatrixXd dz = (PI*xx).array().cos() * (2.*PI*yy).array().sin() * (PI*zz).array().cos();
  z = z + a*Lz*dz;  
#endif

#if TAYLOR_GREEN
  MatrixXd d = x.array().sin()*y.array().sin()*z.array().sin();
  MatrixXd dx = d;
  MatrixXd dy = d;
  MatrixXd dz = d;
  x = x + a*dx;
  y = y + a*dy;
  z = z + a*dz;
  
#endif
  
  mesh->x = x;
  mesh->y = y;
  mesh->z = z;

  /*
  MatrixXd xyz(x.rows(),x.cols()*3);
  xyz << x,y,z;
  cout << "xyz = [" << endl << xyz << "];" << endl;
  return 0;
  */  
  
  // 2Ngeo-2 <= N, so Ngeo <= floor(N/2) + 1
  //int Ngeo = (N/2) + 1; // int divide auto floors
  //printf("Using geometry order %d\n",Ngeo);
  //GeometricFactors3d_Ngeo(mesh,Ngeo);
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
  //app->device.setup("mode      : 'CUDA', "
  //		    "device_id : 0");
  
  //app->device.setup("mode: 'OpenCL', platformID : 0, deviceID: 0");

  setupOccaMesh3d(mesh,app); // build mesh geofacs
 
  if (sizeof(dfloat)==4){
    app->props["defines/p_gamma"] = 1.4f;
  }else{
    app->props["defines/p_gamma"] = 1.4;
  }
  app->props["defines/p_Nfields"] = Nfields;
  app->props["defines/p_tau"] = 1.0;

  MatrixXd Vq = mesh->Vq;  
  MatrixXd xq = Vq*mesh->x;
  MatrixXd yq = Vq*mesh->y;
  MatrixXd zq = Vq*mesh->z;

  double time = 0.0;
  MatrixXd rho, rhou, rhov, rhow, E;
#if VORTEX
  VortexSolution3d(xq,yq,zq,time,rho,rhou,rhov,rhow,E);
#endif
#if TAYLOR_GREEN
  TaylorGreen3d(xq,yq,zq,time,rho,rhou,rhov,rhow,E);
#endif

  MatrixXd Q(Nfields*Np,K);  
  Q << rho,
    rhou,
    rhov,
    rhow,
    E;
  setOccaArray(app, Q, app->o_Q); 
#if VORTEX
  app->props["defines/p_ceil2Nq"] = 1; // dummy var
#endif
  
#if TAYLOR_GREEN
  // for KE computation
  int log2Nq = (int)ceil(log2(Np));
  int ceil2Nq = (int) pow(2,log2Nq);
  app->props["defines/p_ceil2Nq"] = ceil2Nq;

  // interp to gauss points: may be diff from mesh->rq,sq,tq points
  VectorXd rq1D, wq1D;
  JacobiGQ(N, 0, 0, rq1D, wq1D); 
  int Np1 = (N+1);
  int Np3 = Np1*Np1*Np1;
  VectorXd rq(Np3),sq(Np3),tq(Np3),wq(Np3);  
  int sk = 0;
  for (int k = 0; k < Np1; ++k){
    for (int i = 0; i < Np1; ++i){
      for (int j = 0; j < Np1; ++j){
	rq(sk) = rq1D(i);
	sq(sk) = rq1D(j);
	tq(sk) = rq1D(k);

	wq(sk) = wq1D(i)*wq1D(j)*wq1D(k);	
	++sk;
      }
    }
  }
  MatrixXd Vqtmp = Vandermonde3DHex(N, mesh->rq, mesh->sq, mesh->tq);
  MatrixXd VqGauss = mrdivide(Vqtmp,mesh->V); // V = GLL VDM
  MatrixXd wJq = wq.asDiagonal() * (VqGauss * mesh->J);

  // 1D interp matrix
  MatrixXd Vq1D_new = Vandermonde1D(N,rq1D);
  MatrixXd Vq1D_old = Vandermonde1D(N,mesh->rq1D);
  MatrixXd V1D = mrdivide(Vq1D_new,Vq1D_old);

  //  MatrixXd V1D = MatrixXd::Identity(N+1,N+1);  
  //cout << "V1D = " << V1D << endl;

  MatrixXd KE = MatrixXd::Zero(K,1);  
  occa::memory o_V1D, o_wJq, o_KE;
  setOccaArray(app, V1D, o_V1D);
  setOccaArray(app, wJq, o_wJq);
  setOccaArray(app, KE, o_KE);
#endif
 
  //app->eval_surface = app->device.buildKernel("okl/test3D.okl","eval_surface",app->props);return 0;
  
  // build occa kernels
  string path = "okl/Euler3D.okl";
  app->volume = app->device.buildKernel(path.c_str(),"volume",app->props);
  app->surface = app->device.buildKernel(path.c_str(),"surface",app->props); 
  app->update = app->device.buildKernel(path.c_str(),"update",app->props);
  app->eval_surface = app->device.buildKernel(path.c_str(),"eval_surface",app->props);
  occa::kernel compute_aux = app->device.buildKernel(path.c_str(),"compute_aux",app->props);

  // DEBUGGING
#if 0
  app->eval_surface(K, app->o_Vf1D,app->o_Q, app->o_Qf);
  /* 
     MatrixXd Qf(Nfields*NfpNfaces,K);
     getOccaArray(app,app->o_Qf,Qf);
     MatrixXd rhof = Qf.middleRows(0,NfpNfaces);
     cout << "rhof = " << rhof << endl;
     return 0;
  */
  
  // test rhs eval  
  app->volume(K, app->o_vgeo, app->o_vfgeo,app->o_fgeo,
	      app->o_D1D, app->o_Vf1D, app->o_Lf1D,
	      app->o_Q, app->o_Qf,
	      app->o_rhs, app->o_rhsf);  
  getOccaArray(app,app->o_rhs,Q);
  MatrixXd rhs1  = Q.middleRows(0,Np);  
  //  cout << "vol only rhs = " << endl << rhs1 << endl;

  return 0;
  
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

  printf("h = %f, CN = %f, CFL = %f, dt = %f, FinalTime = %f, Nsteps = %d\n",h,CN,CFL,dt,FinalTime,Nsteps);  

  // interp to surface to start
  app->eval_surface(K, app->o_Vf1D,app->o_Q, app->o_Qf);
    
  int interval = ceil(Nsteps/10);
  if (Nsteps < 10){
    interval = 1;
  }

  MatrixXd KE_time = MatrixXd::Zero(Nsteps,1);
  
  int NINT = mesh->rk4a.size();
  for (int i = 0; i < Nsteps; ++i){
    for (int INTRK = 0; INTRK < NINT; ++INTRK){

      const dfloat fdt = (dfloat) dt;
      const dfloat fa  = (dfloat) mesh->rk4a[INTRK];
      const dfloat fb  = (dfloat) mesh->rk4b[INTRK];

      app->volume(K,app->o_vgeo, app->o_vfgeo, app->o_fgeo,
		  app->o_D1D, app->o_Vf1D, app->o_Lf1D,
		  app->o_Q,app->o_Qf,app->o_rhs,app->o_rhsf);
      
      app->surface(K,app->o_vgeo,app->o_fgeo,app->o_mapPq,app->o_Lf1D,
		   app->o_Qf,app->o_rhsf,app->o_rhs);
      
      app->update(K, fa, fb, fdt,
		  app->o_Q, app->o_rhs, app->o_res);
      
      app->eval_surface(K, app->o_Vf1D,app->o_Q, app->o_Qf);
    }

#if TAYLOR_GREEN    
    compute_aux(K, o_V1D, o_wJq, app->o_Q, o_KE);
    getOccaArray(app,o_KE,KE);
    KE_time(i) = KE.sum();
    if (i % interval == 0){
      printf("on timestep %d out of %d: KE = %g\n",i,Nsteps,KE_time(i));
    }    
#else
    if (i % interval == 0){
      printf("on timestep %d out of %d\n",i,Nsteps);
    }
#endif
    
  }
  getOccaArray(app,app->o_Q,Q);
  rho = Q.middleRows(0,Np);
  rhou = Q.middleRows(Np,Np);
  rhov = Q.middleRows(2*Np,Np);
  rhow = Q.middleRows(3*Np,Np);  
  E = Q.middleRows(4*Np,Np);

#if TAYLOR_GREEN

  VectorXd dkedt(KE_time.rows()-1);
  int t_int = max(1,dkedt.rows()/1000);
  double scale = 1.0/pow(2*PI,3);
  cout << "scale = " << scale << endl;
  printf("dkedt = [");
  for (int i = 0; i < KE_time.rows()-1; ++i){
    dkedt(i) = scale*(KE_time(i+1)-KE_time(i))/dt;
    if (i % t_int==0){
      printf(" %f", dkedt(i));
    }
  }
  printf("];\n");
  //Map<VectorXd,0,InnerStride<10> > dkedt2(dkedt.data(), dkedt.size()/10);
  //cout << "dkedt = [" << dkedt2 << "];" << endl;
  //cout << "KE_time = [" << KE_time << "];" << endl;
#endif

#if VORTEX  
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
  printf("L2 error = %g\n",sqrt(werr.sum()));
#endif 

  return 0;
  
}
