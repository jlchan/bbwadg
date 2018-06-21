#include "fem.h"

int main(int argc, char **argv){

  occa::printModeInfo();

  int N = 1;
  int K1D = 1;
  
  Mesh *mesh = (Mesh*) calloc(1, sizeof(Mesh));  
  QuadMesh2d(mesh,K1D,K1D); // make Cartesian mesh

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
  mesh->mapPq = mapPq;

  double DX = mesh->VX.maxCoeff()-mesh->VX.minCoeff();
  double DY = mesh->VY.maxCoeff()-mesh->VY.minCoeff();
  MakeNodeMapsPeriodic2D(mesh,xf,yf,DX,DY,mapPq);

  //  cout << "xf = [" << xf << "];" << endl;
  //  cout << "yf = [" << yf << "];" << endl;    
  //  cout << "mapPq = " << mapPq << endl;

  // ============ problem dependent stuff ============

  int Np = mesh->Np;
  int NfpNfaces = mesh->Nfp * mesh->Nfaces;  
  int K = mesh->K;

  int Nfields = 1; 
  mesh->Nfields = Nfields;
  
  App *app = (App*) calloc(1, sizeof(App));
  setupOccaMesh2d(mesh,app); // build mesh geofacs

  app->props["defines/p_gamma"] = 1.4f;
  app->props["defines/p_Nfields"] = Nfields;
  
  // build occa kernels
  string path = "okl/Euler2D.okl";
  app->volume = app->device.buildKernel(path.c_str(),"volume",app->props);
  app->surface = app->device.buildKernel(path.c_str(),"surface",app->props);
  app->update = app->device.buildKernel(path.c_str(),"update",app->props);

  // set initial condition
  MatrixXd xq = mesh->Vq*mesh->x;
  MatrixXd yq = mesh->Vq*mesh->y;  
  VectorXd rq = mesh->rq;
  VectorXd sq = mesh->sq;

  MatrixXd Q(Nfields*Np,K);
  Q.col(0) = rq; // Nfields = 1
  Q.col(0) = sq; // Nfields = 1  
  cout << "rq = " << rq << endl;
  cout << "sq = " << sq << endl;
  cout << "D1D = " << mesh->D1D << endl;

  
  MatrixXd Qf(Nfields*NfpNfaces,K);
  Qf.fill(0.0);
  setOccaArray(app, Q, app->o_Q);
  setOccaArray(app, Qf, app->o_Qf);  
  
   
  // ============== run RK solver ==================
  
  double h = 2.0 / (double)K1D;
  double CN = dim*(double)((N+1)*(N+2))/2.0; // trace constant
  double CFL = .25;
  double dt = CFL*h / CN;
  
  double FinalTime = .01;
  int Nsteps = std::ceil(FinalTime/dt);  
  dt = FinalTime/(double) Nsteps;
  
  app->volume(K, app->o_vgeo, app->o_fgeo,
	      app->o_D1D, app->o_Vf1D, app->o_wq1D,
	      app->o_Q, app->o_Qf, app->o_rhs, app->o_rhsf);
  
#if 0  
  Nsteps = 1;
  int NINT = 1;  
  for (int i = 0; i < Nsteps; ++i){
    for (int INTRK=0;INTRK < NINT; ++INTRK){
      const dfloat fdt = (dfloat) dt;
      const dfloat fa  = (dfloat) mesh->rk4a[INTRK];
      const dfloat fb  = (dfloat) mesh->rk4b[INTRK];

      app->volume(K, app->o_vgeo, app->o_fgeo,
		  app->o_D1D, app->o_Vf1D, app->o_wq1D,
		  app->o_Q, app->o_Qf, app->o_rhs, app->o_rhsf);
      app->surface(K, app->o_vgeo, app->o_fgeo,
		   app->o_D1D, app->o_Vf1D, app->o_wq1D,
		   app->o_Q, app->o_Qf, app->o_rhs, app->o_rhsf);
      app->update(K, fa, fb, fdt, app->o_vgeo,
		  app->o_Vf1D, app->o_wq1D,
		  app->o_Q, app->o_Qf, app->o_rhs, app->o_rhsf);      
    }    
  }
  getOccaArray(app,app->o_Q,Q);
  
#endif

  return 0;
  
}
