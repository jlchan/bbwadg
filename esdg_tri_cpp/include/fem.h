#ifndef _FEM_HEADER
#define _FEM_HEADER

#include <strings.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "dfloat.h"

#include "Mesh.h"


/* geometric/mesh functions */
Mesh *ReadGmsh3d(char *filename);

// 2d tri routines
void TriMesh2d(Mesh *mesh, int Nx, int Ny);
void InitRefTri(Mesh *mesh, int N);
void MapTriNodes(Mesh *mesh);
void ConnectTriElems(Mesh *mesh); // computes EToE, EToF

// 2d quad routines
void InitRefQuad(Mesh *mesh, int N);
void QuadMesh2d(Mesh *mesh, int Nx, int Ny);
void MapQuadNodes(Mesh *mesh);
void MakeNodeMapsPeriodic2d(Mesh *mesh, MatrixXd xf, MatrixXd yf, double DX, double DY, MatrixXi &mapP);

// general 2D routines
void GeometricFactors2d(Mesh *mesh);
void Normals2d(Mesh *mesh);
void ReadGmsh2d(Mesh *mesh); // Gmsh



// 3d routines
void HexMesh3d(Mesh *mesh, int Nx, int Ny, int Nz);
void InitRefHex(Mesh *mesh, int N);
void MapHexNodes(Mesh *mesh);
void GeometricFactors3d(Mesh *mesh);
void GeometricFactors3d_Ngeo(Mesh *mesh,int Ngeo);
void Normals3d(Mesh *mesh);
void MakeNodeMapsPeriodic3d(Mesh *mesh, MatrixXd xf, MatrixXd yf, MatrixXd zf,
			    double DX, double DY, double DZ, MatrixXi &mapP);

// dimension independent routines
void BuildFaceNodeMaps(Mesh *mesh, MatrixXd xf, MatrixXd yf, MatrixXd zf, MatrixXi &mapP);
void ConnectElems(Mesh *mesh); // computes EToE, EToF

// occa setup 
void setupOccaMesh2d(Mesh *mesh, App *app);
void setupOccaMesh3d(Mesh *mesh, App *app);
void setOccaArray(App *app, MatrixXd A, occa::memory &c_A);
void setOccaIntArray(App *app, MatrixXi A, occa::memory &c_A);
void getOccaArray(App *app, occa::memory c_A, MatrixXd &A);

void test_RK(Mesh *mesh, int KblkU);
void setOccaArray(MatrixXd A, occa::memory &B); // assumes matrix is double
void setOccaIntArray(MatrixXi A, occa::memory &B); // assumes matrix is int

void writeVisToGMSH(string fileName,Mesh *mesh, dfloat *Q, int iField, int Nfields);

#endif
