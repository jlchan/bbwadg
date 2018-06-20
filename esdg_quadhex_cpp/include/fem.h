#ifndef _FEM_HEADER
#define _FEM_HEADER

#include <strings.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "dfloat.h"

#include "Mesh.h"
#include "Basis.h"


/* geometric/mesh functions */
Mesh *ReadGmsh3d(char *filename);

void InitRefData2d(Mesh *mesh, int N, int Nfields);
void QuadMesh2d(Mesh *mesh, int Nx, int Ny);
void ConnectElems(Mesh *mesh, int dim); // computes EToE, EToF
void GeometricFactors2d(Mesh *mesh);
void BuildFaceNodeMaps(Mesh *mesh, MatrixXd xf, MatrixXd yf, MatrixXd zf, MatrixXi &mapP);

void ReadGmsh2d(Mesh *mesh); // Gmsh




/*
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
*/
void test_RK(Mesh *mesh, int KblkU);
void setOccaArray(MatrixXd A, occa::memory &B); // assumes matrix is double
void setOccaIntArray(MatrixXi A, occa::memory &B); // assumes matrix is int

void writeVisToGMSH(string fileName,Mesh *mesh, dfloat *Q, int iField, int Nfields);

#endif
