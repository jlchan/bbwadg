#ifndef _FEM_HEADER
#define _FEM_HEADER

#include <strings.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "dfloat.h"

#include "Mesh.h"
#include "Basis.h"

/* prototypes for storage functions (Utils.c) */
dfloat **BuildMatrix(int Nrows, int Ncols);
dfloat  *BuildVector(int Nrows);
int    **BuildIntMatrix(int Nrows, int Ncols);
int     *BuildIntVector(int Nrows);

dfloat **DestroyMatrix(dfloat **);
dfloat  *DestroyVector(dfloat *);
int    **DestroyIntMatrix(int **);
int     *DestroyIntVector(int *);

void PrintMatrix(char *message, dfloat **A, int Nrows, int Ncols);
void SaveMatrix(char *filename, dfloat **A, int Nrows, int Ncols);

// for more general meshes
Mesh *ReadGmshHybrid(char *filename);
void *FacePairHybrid(Mesh *mesh);

// compare entries of length 4 array (for general face pairing)
static int compare4(const void *a, const void *b){
  int *aI = (int*) a;   int *bI = (int*) b;
  int a1 = aI[0], a2 = aI[1], a3 = aI[2], a4 = aI[3];
  int b1 = bI[0], b2 = bI[1], b3 = bI[2], b4 = bI[3];

  if(b4>a4) return -1;
  if(a4>b4) return 1;

  if(b3>a3) return -1;
  if(a3>b3) return 1;

  if(b2>a2) return -1;
  if(a2>b2) return 1;

  if(b1>a1) return -1;
  if(a1>b1) return 1;

  return 0;
}


/* geometric/mesh functions */
Mesh *ReadGmsh3d(char *filename);

void PrintMesh ( Mesh *mesh );

void Normals3d(Mesh *mesh, int k,
	       double *nx, double *ny, double *nz, double *sJ);


void GeometricFactors3d(Mesh *mesh, int k,
			double *drdx, double *dsdx, double *dtdx,
			double *drdy, double *dsdy, double *dtdy,
			double *drdz, double *dsdz, double *dtdz,
			double *J);

// start up
void StartUp3d(Mesh *mesh);
void BuildMaps3d(Mesh *mesh);
void FacePair3d(Mesh *mesh);

// set initial condition
void WaveSetU0(Mesh *mesh, dfloat *Q, dfloat time, int field,
	       double(*uexptr)(double,double,double,double));
void WaveSetData3d(dfloat *Q);
void WaveGetData3d(Mesh *mesh, dfloat *Q);

// solver
void test_kernels(Mesh *mesh);
void time_kernels(Mesh *mesh);
void time_curved_kernels(Mesh *mesh,int nsteps);
void Wave_RK(Mesh *mesh, dfloat FinalTime, dfloat dt, int useWADG);

// for BB vs NDG stability (paper result only)
void Wave_RK_sample_error(Mesh *mesh, dfloat FinalTime, dfloat dt,
                          double(*uexptr)(double,double,double,double));

void RK_step(Mesh *mesh, dfloat rka, dfloat rkb, dfloat fdt);
void compute_error(Mesh *mesh, double time, dfloat *Q,
		   double(*uexptr)(double,double,double,double),
		   double &L2err, double &relL2err);

// planar elements + heterogeneous media
void RK_step_WADG_subelem(Mesh *mesh, dfloat rka, dfloat rkb, dfloat fdt, dfloat time);

// curvilinear and WADG-based
void InitQuadratureArrays(Mesh *mesh);
void WaveProjectU0(Mesh *mesh, dfloat *Q, dfloat time,int field,
		   double(*uexptr)(double,double,double,double));

void BuildFaceNodeMaps(Mesh *mesh, MatrixXd xf, MatrixXd yf, MatrixXd zf,
		       MatrixXi &mapP);

void GordonHallSphere(Mesh *mesh);
void InitWADG_curved(Mesh *mesh); // for curved elements (+ optional subelem)
void InitWADG_subelem(Mesh *mesh,double(*c2_ptr)(double,double,double)); // for planar elements + subelem
void checkCurvedGeo(Mesh *mesh);
void RK_step_WADG(Mesh *mesh, dfloat rka, dfloat rkb, dfloat fdt);

void test_RK(Mesh *mesh, int KblkU);
void setOccaArray(MatrixXd A, occa::memory &B); // assumes matrix is double
void setOccaIntArray(MatrixXi A, occa::memory &B); // assumes matrix is int

// cruft - can remove, but may be useful in future
void setupCG(Mesh *mesh);

void writeVisToGMSH(string fileName,Mesh *mesh, dfloat *Q, int iField, int Nfields);

#endif
