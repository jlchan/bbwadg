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


#if 0  // testing
  rho.fill(2.0);
  u.fill(.1);
  v.fill(.2);
  p.fill(1.0);  
#endif
  
}

int main(int argc, char **argv){

  int N = 4;
  int K1D = 8;

  //occa::printModeInfo();   return 0;
  if (argc > 2){
    N = atoi(argv[1]);
    K1D = atoi(argv[2]);
    printf("setting N = %d, K1D = %d\n",N,K1D);
  }

  double FinalTime = .5;
  double CFL = .25;
  double a = 0*.25; // curved warping

  if (argc > 4){
    CFL = atof(argv[3]);
    FinalTime = atof(argv[4]);
    a = atof(argv[5]);
    printf("setting CFL = %f, T = %f, curved warping a =%f\n",CFL,FinalTime,a);
  }
  
  Mesh *mesh = new Mesh;
 
  TriMesh2d(mesh,2*K1D,K1D); // make Cartesian mesh
  double Lx = 10;
  double Ly = 5;
  mesh->VX = (mesh->VX.array() + 1.0).array()*Lx;
  mesh->VY = (mesh->VY)*Ly;
  /*
  cout << "VX = " << mesh->VX << endl;
  cout << "VY = " << mesh->VY << endl;
  */

  // ============ physics independent stuff ===========

  printf("N = %d, K = %d\n",N,mesh->K);
  
  InitRefTri(mesh, N);

  // make physical nodes + geofacs
  MapTriNodes(mesh); // low order mapping
  GeometricFactors2d(mesh); 
  Normals2d(mesh);

  // makes EToE, EToF
  ConnectElems(mesh);

  // data
  int Nfields = 4; 
  mesh->Nfields = Nfields;
  int K = mesh->K;
  int Np = mesh->Np;
  int Nq = mesh->Nq;    
  int NfqNfaces = mesh->Nfq * mesh->Nfaces;
  
  // build node maps
  MatrixXd xf = (mesh->Vf)*(mesh->x);
  MatrixXd yf = (mesh->Vf)*(mesh->y);
  MatrixXd zf(xf.rows(),xf.cols()); zf.fill(0.0);
  MatrixXi mapPq;
  BuildFaceNodeMaps(mesh,xf,yf,zf,mapPq);

  // builds periodicity into the node maps
  double DX = mesh->VX.maxCoeff()-mesh->VX.minCoeff();
  double DY = mesh->VY.maxCoeff()-mesh->VY.minCoeff();
  MakeNodeMapsPeriodic2d(mesh,xf,yf,DX,DY,mapPq);

  // correct mapPq for Nfields > 1  
  MatrixXi mapPqNfields(NfqNfaces,K);
  for (int i = 0; i < NfqNfaces*K; ++i){
    int idP = mapPq(i);
    int e = idP / NfqNfaces;
    int fidP = idP % NfqNfaces; 
    mapPqNfields(i) = fidP + e*NfqNfaces*Nfields;
  }  
  mesh->mapPq = mapPqNfields;

  // ========================== set up OCCA application

  App *app = new App;
  //  app->device.setup("mode: 'Serial'");
  app->device.setup("mode: 'CUDA', device_id: 0");

  // define macros for OCCA kernels
  app->props["defines/p_Np"] = mesh->Np; // number of dims
  app->props["defines/p_Nq"] = Nq;
  app->props["defines/p_NfqNfaces"] = NfqNfaces;  
  app->props["defines/p_NqT"] = Nq + NfqNfaces; // total quadrature point 
  app->props["defines/p_T"] = max(mesh->Nfq * mesh->Nfaces,mesh->Np);
  
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

  // Euler-specific values
  app->props["defines/p_gamma"] = GAMMA;
  app->props["defines/p_Nfields"] = Nfields;
  if (sizeof(dfloat)==4){
    app->props["defines/p_tau"] = 1.f;
  }else{
    app->props["defines/p_tau"] = 1.0;
  }

  // interpolate to quad pts and store  
  //  // remove eventually...  
  //  setupOccaMesh2d(mesh,app); // store and pack mesh geofacs in OCCA mem

    // build occa kernels  
  string path = "okl/Euler2DTri.okl";

  //testing
  occa::kernel volume, surface, update, project;
  volume = app->device.buildKernel(path.c_str(),"volume",app->props);
  surface = app->device.buildKernel(path.c_str(),"surface",app->props);
  update = app->device.buildKernel(path.c_str(),"update",app->props);
  project = app->device.buildKernel(path.c_str(),"project",app->props);

  // =================== set initial condition

  MatrixXd x = mesh->x;
  MatrixXd y = mesh->y;  
 
  // smooth isentropic vortex
  MatrixXd rho,u,v,p;  
  double time = 0.0;  
  VortexSolution2d(x,y,time,rho,u,v,p);

  MatrixXd rhou = rho.array()*u.array();
  MatrixXd rhov = rho.array()*v.array();
  MatrixXd E = p.array()/(GAMMA-1.0) + .5*rho.array()*(u.array().square() + v.array().square());

  MatrixXd Q(Nfields*Np,K);
  Q << rho,rhou,rhov,E;
  
  // ================= set node maps + boundary condition
  
  MatrixXi bcFlag(xf.rows(),xf.cols());  
  bcFlag.fill(0); // no BCs

  occa::memory o_bcFlag, o_mapPq;
  setOccaIntArray(app, bcFlag, o_bcFlag);
  setOccaIntArray(app, mesh->mapPq, o_mapPq);  

  // ================= construct discretization matrices

  MatrixXd Vq = mesh->Vq;
  MatrixXd Vf = mesh->Vf;
  MatrixXd M = Vq.transpose()*mesh->wq.asDiagonal()*Vq;
  MatrixXd VqTW = Vq.transpose()*mesh->wq.asDiagonal();
  MatrixXd VfTWf = Vf.transpose()*mesh->wf.asDiagonal();
  MatrixXd Pq = mldivide(M,VqTW);
  MatrixXd Lf = mldivide(M,VfTWf);

  // hybridized SBP ops
  int NqT = Nq+NfqNfaces;
  MatrixXd Ef = Vf*Pq;
  MatrixXd Wf = mesh->wf.asDiagonal();
  MatrixXd Br = Wf * mesh->nrJ.asDiagonal();
  MatrixXd Bs = Wf * mesh->nsJ.asDiagonal();
  MatrixXd Qr = Pq.transpose()*(M*mesh->Dr)*Pq;
  MatrixXd Qs = Pq.transpose()*(M*mesh->Ds)*Pq;  
  MatrixXd Zf = MatrixXd::Zero(NfqNfaces,NfqNfaces);

  // skew form of the operators
  MatrixXd QNr(NqT,NqT);
  QNr << Qr-Qr.transpose(), Ef.transpose()*Br, -Br*Ef, Zf;
  //QNr << Qr-.5*Ef.transpose()*Br*Ef, .5*Ef.transpose()*Br, -.5*Br*Ef, .5*Br;
  MatrixXd QNs(NqT,NqT);
  QNs << Qs-Qs.transpose(), Ef.transpose()*Bs, -Bs*Ef, Zf;
  //QNs << Qs-.5*Ef.transpose()*Bs*Ef, .5*Ef.transpose()*Bs, -.5*Bs*Ef, .5*Bs;  

  MatrixXd VN(NqT,Np);
  VN << Vq,Vf;
  MatrixXd VNT = VN.transpose();
  MatrixXd VNP = VN*Pq;
  MatrixXd PN = mldivide(M,VNT); 
  
  occa::memory o_Q, o_Qv, o_Qf; // solution   
  occa::memory o_vgeo, o_fgeo; // geofacs

  MatrixXd vgeo(Nvgeo*Nq,K);  
  vgeo << Vq*(mesh->rxJ),Vq*(mesh->ryJ),
    Vq*(mesh->sxJ),Vq*(mesh->syJ),Vq*mesh->J;

  setOccaArray(app,vgeo,o_vgeo);

  MatrixXd fgeo(Nfgeo*NfqNfaces,K);  
  fgeo << mesh->nxJ,mesh->nyJ,mesh->sJ; // already at quad pts
  setOccaArray(app,fgeo,o_fgeo);
  
  occa::memory o_rhs, o_res;  // timestep stuff
  setOccaArray(app,MatrixXd::Zero(Np*Nfields,K),o_rhs);
  setOccaArray(app,MatrixXd::Zero(Np*Nfields,K),o_res);  

  occa::memory o_VNP, o_QNr, o_QNs, o_PN, o_Lf, o_Vq; // operators

  // set solution arrays
  setOccaArray(app, Q, o_Q);
  MatrixXd Qv(Nfields*Nq,K);
  Qv << Vq*rho,Vq*rhou,Vq*rhov,Vq*E;
  
  setOccaArray(app, Qv, o_Qv);
  setOccaArray(app, MatrixXd::Zero(Nfields*NfqNfaces,K), o_Qf);    

  // set operators
  setOccaArray(app, Vq, o_Vq);
  setOccaArray(app, Lf, o_Lf);
  setOccaArray(app, QNr, o_QNr);
  setOccaArray(app, QNs, o_QNs);
  setOccaArray(app, PN, o_PN);
  setOccaArray(app, VNP, o_VNP);      


  // ============== testing RHS

  /*
  //cout << "vgeo(:,1) = " << endl << vgeo.col(0) << endl;
 
  project(K, o_VNP, o_Qv, o_Qf);

  volume(K, o_vgeo,
	 o_QNr, o_QNs, o_PN,
	 o_Qv, o_Qf,
	 o_rhs);

  surface(K, o_fgeo,
	  o_mapPq, o_bcFlag,
	  o_Lf, o_Qf, 
	  o_rhs); 

  MatrixXd rhs(Np*Nfields,K);
  getOccaArray(app,o_rhs,rhs);
  cout << "rhs = [" << endl << rhs << endl <<  "];" << endl;   
  return 0;
  */

  
  // ============== run RK solver ==================

  //double h = 2.0 / (double) K1D; //min(DX,DY) / (double) K1D;
  double h = mesh->J.maxCoeff() / mesh->sJ.maxCoeff(); // J = O(h^d), Jf = O(h^{d-1}) in d dims
  double CN = (double)((N+1)*(N+2))/2.0; // trace constant for GQ hexes
  double dt = CFL * h / CN;
  
  int Nsteps = (int) ceil(FinalTime/dt);
  dt = FinalTime/(double) Nsteps;

  printf("dt = %f, FinalTime = %f, Nsteps = %d\n",dt,FinalTime,Nsteps);  

  int interval = max((int) ceil(Nsteps/10),1);  
  printf("Interval = %d\n",interval);

  int NINT = mesh->rk4a.size(); // num RK steps
  for (int i = 0; i < Nsteps; ++i){
    for (int INTRK = 0; INTRK < NINT; ++INTRK){

      const dfloat fdt = (dfloat) dt;
      const dfloat fa  = (dfloat) mesh->rk4a[INTRK];
      const dfloat fb  = (dfloat) mesh->rk4b[INTRK];

      // entropy projection
      project(K, o_VNP, o_Qv, o_Qf);
      
      volume(K, o_vgeo,
	     o_QNr, o_QNs, o_PN,
	     o_Qv, o_Qf,
	     o_rhs);
      
      surface(K, o_fgeo,
	      o_mapPq, o_bcFlag,
	      o_Lf, o_Qf, 
	      o_rhs); 
      
      update(K, fa, fb, fdt,
	     o_vgeo,o_Vq,
	     o_Q, o_Qv,
	     o_rhs, o_res);     

    }

    if (i % interval == 0){
      printf("on timestep %d out of %d\n",i,Nsteps);
    }
  }
  getOccaArray(app,o_Q,Q);
  rho = Q.middleRows(0,Np);
  rhou = Q.middleRows(Np,Np);
  rhov = Q.middleRows(2*Np,Np);
  E = Q.middleRows(3*Np,Np);

  // should really use finer quadrature for error eval 
  MatrixXd rhoex;
  MatrixXd xq = Vq*x;
  MatrixXd yq = Vq*y;
  VortexSolution2d(xq,yq,FinalTime,rhoex,u,v,p);
  
  MatrixXd wJq = mesh->wq.asDiagonal() * (Vq*mesh->J);  
  MatrixXd rhouex = rhoex.array()*u.array();
  MatrixXd rhovex = rhoex.array()*v.array();
  MatrixXd Eex = p.array()/(GAMMA-1.0) + .5*rhoex.array()*(u.array().square() + v.array().square());

  MatrixXd werr = wJq.array()*((rhoex - Vq*rho).array().square() +
			       (rhouex - Vq*rhou).array().square() +
			       (rhovex - Vq*rhov).array().square() +
			       (Eex - Vq*E).array().square());
  printf("L2 error for rho = %g\n",sqrt(werr.sum()));


  return 0;
  
}
