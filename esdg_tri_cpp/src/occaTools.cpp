#include "dfloat.h"
#include "fem.h"
//#include "types.h"

// set occa array:: cast to dfloat
void setOccaArray(App *app, MatrixXd A, occa::memory &c_A){
  int r = A.rows();
  int c = A.cols();
  dfloat *f_A = (dfloat*)malloc(r*c*sizeof(dfloat));
  Map<MatrixXdf >(f_A,r,c) = A.cast<dfloat>();
  c_A = app->device.malloc(r*c*sizeof(dfloat),f_A);
  free(f_A);  
}

// set occa array for int matrices
void setOccaIntArray(App *app, MatrixXi A, occa::memory &c_A){
  int r = A.rows();
  int c = A.cols();
  int *f_A = (int*)malloc(r*c*sizeof(int));
  Map<MatrixXi >(f_A,r,c) = A;
  c_A = app->device.malloc(r*c*sizeof(int),f_A);
  free(f_A);
}

// get occa array:: cast to double
void getOccaArray(App *app, occa::memory c_A, MatrixXd &A){
  int r = A.rows();
  int c = A.cols();
  dfloat *f_A = (dfloat*)malloc(r*c*sizeof(dfloat));
  c_A.copyTo(f_A);
  for (int i = 0; i < r*c; ++i){
    A(i) = f_A[i]; // somewhat inefficient
  }

  free(f_A);
}
		 

