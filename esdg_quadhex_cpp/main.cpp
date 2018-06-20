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
  
  int N = 1;
  int Nfields = 4; // euler in 2D

  QuadMesh2d(mesh,2,1); // make Cartesian mesh

  InitRefData2d(mesh, N, Nfields); 

  GeometricFactors2d(mesh); 
  
  int dim = 2;
  ConnectElems(mesh,dim);
 
  cout << "Vf = " << endl << mesh->Vf << endl;
  MatrixXd xf = (mesh->Vf)*(mesh->x);
  MatrixXd yf = (mesh->Vf)*(mesh->y);
  MatrixXd zf(xf.rows(),xf.cols()); zf.fill(0.0);
  cout << "xf = [" << xf << "];" << endl;
  cout << "yf = [" << yf << "];" << endl;  
  MatrixXi mapP;
  BuildFaceNodeMaps(mesh,xf,yf,zf,mapP);
  
  return 0;


  /* end game */
  exit(0);
}
