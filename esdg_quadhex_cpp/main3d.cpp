#include "fem.h"

#define GAMMA 1.4f

int main(int argc, char **argv){

  //occa::printModeInfo();

  int N = 1;
  int K1D = 2;
  
  Mesh *mesh = (Mesh*) calloc(1, sizeof(Mesh));  
  HexMesh3d(mesh,K1D,K1D,K1D); // make Cartesian mesh

  InitRefData3d(mesh, N);
  return 0;
  
  int dim = 3;
  ConnectElems(mesh,dim);

  //cout << "EToF = [" << endl << mesh->EToF << "];" << endl;
  cout << "EToE = [" << endl << mesh->EToE << "];" << endl;
  
  return 0;
  
}
