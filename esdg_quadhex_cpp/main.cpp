#include "fem.h"

#define TEST 0

int main(int argc, char **argv){

#if TEST // test eigen
  test_solve();
  test_basis();
  return 0;
#endif

  occa::printModeInfo();

  Mesh *mesh = (Mesh*) calloc(1, sizeof(Mesh));
  
  int N = 3;
  int Nfields = 4; // euler in 2D

  int K1D = 2;
  QuadMesh2d(mesh,K1D,K1D); // make Cartesian mesh

  InitRefData2d(mesh, N, Nfields); 

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

  App *app = (App*) calloc(1, sizeof(App));
  setupOccaMesh2d(mesh,app); // build mesh geofacs

  app->props["defines/p_gamma"] = 1.4f;
  
  // build occa kernels
  string path = "okl/Euler2D.okl";
  app->volume = app->device.buildKernel(path.c_str(),"volume",app->props);
  app->surface = app->device.buildKernel(path.c_str(),"surface",app->props);
  app->update = app->device.buildKernel(path.c_str(),"update",app->props);  
   
  // ============== run RK solver ==================
  
  double h = 2.0 / (double)K1D;
  double CN = dim*(double)((N+1)*(N+2))/2.0; // trace constant
  double CFL = .25;
  double dt = CFL*h / CN;
  
  double FinalTime = 1.0;
  int Nsteps = std::ceil(FinalTime/dt);
  dt = FinalTime/(double) Nsteps;

  for (int i = 0; i < Nsteps; ++i){
    for (int INTRK=0;INTRK < 5; ++INTRK){
      const dfloat fdt = (dfloat) dt;
      const dfloat fa  = (dfloat) mesh->rk4a[INTRK];
      const dfloat fb  = (dfloat) mesh->rk4b[INTRK];

      //app->volume(mesh->K,o_vgeo,o_fgeo,o_D,o_Vf
    }
  }
  

  return 0;

  /* end game */
  exit(0);
}
