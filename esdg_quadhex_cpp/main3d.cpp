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

#define pfun(rho, u, v, w, E)                                   \
  ((GAMMA - 1.) * (E - .5f * rho * (u * u + v * v + w * w)))
#define beta(rho, u, v, w, E)                           \
  (rho / (2. * pfun(rho, u, v, w, E))) // inverse temp                                                                                                                                                                                     
// map conservation to entropy vars
#define pfun(rho, u, v, w, E)                                   \
  ((GAMMA - 1.) * (E - .5f * rho * (u * u + v * v + w * w)))
#define rhoeU(rho, rhou, rhov, rhow, E)                         \
  (E - .5f * (rhou * rhou + rhov * rhov + rhow * rhow) / rho)
#define sU(rho, rhou, rhov, rhow, E)                            \
  (log((GAMMA - 1.) * rhoeU(rho, rhou, rhov, rhow, E) /    \
         pow(rho, GAMMA)))

// map entropy to conservation vars                                                                                                                                                                                                        
#define sV(V1, V2, V3, V4, V5)                                  \
  (GAMMA - V1 + (V2 * V2 + V3 * V3 + V4 * V4) / (2. * V5))
#define rhoeV(V1, V2, V3, V4, V5)                                       \
  (pow((GAMMA - 1.) / pow(-V5, GAMMA), 1. / (GAMMA - 1.)) * \
   exp(-sV(V1, V2, V3, V4, V5) / (GAMMA - 1.)))

static void VU(dfloat rho, dfloat rhou, dfloat rhov, dfloat rhow, dfloat E,
        dfloat &V1, dfloat &V2, dfloat &V3, dfloat &V4, dfloat &V5)
{
  const dfloat rhoe = rhoeU(rho, rhou, rhov, rhow, E);
  const dfloat invrhoe = 1./rhoe;
  V1 = (-E + rhoe * (GAMMA + 1. - sU(rho, rhou, rhov, rhow, E))) * invrhoe;
  V2 = rhou * invrhoe;
  V3 = rhov * invrhoe;
  V4 = rhow * invrhoe;
  V5 = (-rho) * invrhoe;
}

static void UV(dfloat V1, dfloat V2, dfloat V3, dfloat V4, dfloat V5,
        dfloat &rho, dfloat &rhou, dfloat &rhov, dfloat &rhow, dfloat &E)
{
  const dfloat rhoe = rhoeV(V1, V2, V3, V4, V5);
  rho = rhoe * (-V5);
  rhou = rhoe * (V2);
  rhov = rhoe * (V3);
  rhow = rhoe * (V4);
  E = rhoe * (1. - (V2 * V2 + V3 * V3 + V4 * V4) / (2. * V5));
}


int main(int argc, char **argv){

  //occa::printModeInfo();

  int N = 3;
  int K1D = 8;
  double CFL = .5;  
  double FinalTime = .5;
  
  Mesh *mesh = (Mesh*) calloc(1, sizeof(Mesh));  
  HexMesh3d(mesh,K1D,K1D,K1D); // make Cartesian mesh

#if 1
  // [0,10] x [0,20] x [0,10] for vortex
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

  MatrixXd x = mesh->x;
  MatrixXd y = mesh->y;
  MatrixXd z = mesh->z;
  MatrixXd dx = .5*(PI*x.array()/10).sin()*(2.0*PI*y.array()/20).sin()*(PI*z.array()/10).sin();
  MatrixXd dy = -(2.0*PI*x.array()/10).sin()*(PI*y.array()/20).sin()*(2.0*PI*z.array()/10).sin();
  MatrixXd dz = dx;

  /*  MatrixXd d = (PI*(1+x.array())/2).sin()*(PI*(1+y.array())/2).sin()*(PI*(1+z.array())/2).sin();
  dx.fill(0.0);
  dy.fill(0.0); 
  dz.fill(0.0);

  dx = d;
  // dy = d;
  //  dz = d;  
  */
  
  double a = .10;
  x = x + a*dx;
  y = y + a*dy;
  z = z + a*dz;

#if 0
  cout << "x = [" << x << "];" << endl;
  cout << "y = [" << y << "];" << endl;
  cout << "z = [" << z << "];" << endl;
  return 0;
#endif
  
  mesh->x = x;
  mesh->y = y;
  mesh->z = z;  
  
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
  app->device.setup("mode: 'CUDA', deviceID: 0");
  //app->device.setup("mode: 'OpenCL', platformID : 0, deviceID: 0");

  setupOccaMesh3d(mesh,app); // build mesh geofacs
 
  app->props["defines/p_gamma"] = GAMMA;
  app->props["defines/p_Nfields"] = Nfields;
  app->props["defines/p_tau"] = 0.0;
  
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


  MatrixXd wJq = mesh->wq.asDiagonal() * (mesh->Vq*mesh->J);

#if 0
  app->eval_surface(K, app->o_Vf1D,app->o_Q, app->o_Qf);
  
  // test rhs eval
  app->volume(K, app->o_vgeo, app->o_vfgeo,
	      app->o_D1D, app->o_Vf1D, app->o_Lf1D,
	      app->o_Q, app->o_Qf,
	      app->o_rhs, app->o_rhsf);

  getOccaArray(app,app->o_rhs,Q);
  MatrixXd rhs1  = Q.middleRows(0,Np);
  
  //cout << "vol only rhs = " << endl << rhs1 << endl;
  
  app->surface(K, app->o_vgeo, app->o_fgeo, app->o_mapPq,
	       app->o_Lf1D, app->o_Qf, app->o_rhsf,
	       app->o_rhs); 
  
  getOccaArray(app,app->o_rhs,Q);
  rhs1 = Q.middleRows(0,Np);
  cout << "vol+surf rhs = " << endl << rhs1 << endl;

  /*
MatrixXd V1(Np,K);
  MatrixXd V2(Np,K);
  MatrixXd V3(Np,K);
  MatrixXd V4(Np,K);
  MatrixXd V5(Np,K);
  double UVerr = 0.0;  
  for (int e = 0; e < K; ++e){
    for (int i = 0; i < Np; ++i){
      double V1i,V2i,V3i,V4i,V5i;
      VU(rho(i,e),rhou(i,e),rhov(i,e),rhow(i,e),E(i,e),V1i,V2i,V3i,V4i,V5i);
      V1(i,e) = V1i;
      V2(i,e) = V2i;
      V3(i,e) = V3i;
      V4(i,e) = V4i;
      V5(i,e) = V5i;
      double rhoi,rhoui,rhovi,rhowi,Ei;
      UV(V1i,V2i,V3i,V4i,V5i,rhoi,rhoui,rhovi,rhowi,Ei);

      UVerr += fabs(rho(i,e)-rhoi) + fabs(rhou(i,e)-rhoui) + fabs(rhov(i,e)-rhovi) + fabs(rhow(i,e)-rhowi) + fabs(E(i,e)-Ei);
    }
  }

  MatrixXd rhsrho  = Q.middleRows(0,Np);
  MatrixXd rhsrhou = Q.middleRows(Np,Np);
  MatrixXd rhsrhov = Q.middleRows(2*Np,Np);
  MatrixXd rhsrhow = Q.middleRows(3*Np,Np);  
  MatrixXd rhsE    = Q.middleRows(4*Np,Np);
  double rhstest = (wJq.array()*((V1.array()*rhsrho.array() + V2.array()*rhsrhou.array() + V3.array()*rhsrhov.array() + V4.array()*rhsrhow.array() + V5.array()*rhsE.array()))).matrix().sum();
  printf("uverr = %g, rhstest = %g\n",UVerr,rhstest);
  */  
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
  MatrixXd Vq = Vandermonde3DHex(N,mesh->rq,mesh->sq,mesh->tq);  
  MatrixXd Vq2 = mrdivide(Vqtmp,Vq);
  MatrixXd xq2 = Vq2 * (mesh->Vq*mesh->x);
  MatrixXd yq2 = Vq2 * (mesh->Vq*mesh->y);
  MatrixXd zq2 = Vq2 * (mesh->Vq*mesh->z);  
  
  MatrixXd rhoex,rhouex,rhovex,rhowex,Eex;
  VortexSolution3d(xq2,yq2,zq2,FinalTime,rhoex,rhouex,rhovex,rhowex,Eex);

  MatrixXd wJq2 = wq2.asDiagonal() * (Vq2*(mesh->Vq*mesh->J));  
  //  MatrixXd wJq = mesh->wq.asDiagonal() * (mesh->Vq*mesh->J);  
  //  MatrixXd rhouex = rhoex.array()*u.array();
  //  MatrixXd rhovex = rhoex.array()*v.array();
  //  MatrixXd Eex = p.array()/(GAMMA-1.0) + .5*rhoex.array()*(u.array().square() + v.array().square());
  MatrixXd werr = wJq2.array()*(rhoex - Vq2*rho).array().square();
  printf("L2 error for rho = %g\n",sqrt(werr.sum()));
 
  return 0;
  
}
