#ifndef _MESH_HEADER
#define _MESH_HEADER

#include <occa.hpp>
#include "Basis.h"

/*
===========================================================
Contains information about the physical mesh 
and reference elements.
===========================================================
 */



// default order
#ifndef p_N
#define p_N 6
#endif

#define NODETOL   1e-11
/*
#define p_Nfp     ((p_N+1)*(p_N+2)/2)
#define p_Np      ((p_N+1)*(p_N+2)*(p_N+3)/6)
#define p_Nfaces  6
#define p_Nfields 5 // Euler equation
#define p_Nfields 9 // elastic wave equation
*/

#define max(a,b)  ( (a>b)?a:b )
#define min(a,b)  ( (a<b)?a:b )

/* list of vertices on each edge */

typedef struct foo {

  // =============== mesh stuff ===================

  int Nv;     /* number of mesh vertices */
  int K;      /* number of mesh elements */
  VectorXd VX,VY,VZ;  
  MatrixXi EToV;  
  MatrixXi EToE;
  MatrixXi EToF;
  
  // =============== reference elem stuff ===================

  int Nverts; /* number of vertices per element */
  int Nfaces; /* number of faces per element */
  int Nedges; /* number of edges per element (3d only) */

  int N, Np, Nfp;
  int Nq, Nfq;  
  int Nfields;
  
  /* high order node info */
  MatrixXi Fmask;

  // 1D operators (for tensor product)
  MatrixXd D1D, Vf1D, Lf1D;
  VectorXd wq1D;
  
  // nodal points (GLL for quads/hexes)
  VectorXd r,s,t;
  MatrixXd V, Dr, Ds, Dt, LIFT; // Eigen-based matrices: cubature, nodal deriv/lift  
  
  // quad points (GQ for quads/hexes)
  VectorXd wq, rq, sq, tq;
  MatrixXd Vq, Vf;

  VectorXd rf,sf,tf,wf;

  // =============== global DG stuff ===================

  MatrixXd x,y,z; // xyz-coordinates of element nodes (3d)
  MatrixXd rxJ,ryJ,rzJ,sxJ,syJ,szJ,txJ,tyJ,tzJ,J; // geofacs
  MatrixXd nxJ,nyJ,nzJ,sJ; // surface geofacs

  MatrixXi mapPq; // node map
  
  // time stepping constants
  VectorXd rk4a, rk4b, rk4c;  

  // mesh size for dt - computed  
  double hMax; 
  double hMin;
 
}Mesh;


// struct for occa arrays, kernels, etc
typedef struct foobar {

  occa::device device;
  occa::properties props; 
  occa::kernel volume, surface, update, eval_surface; 

  occa::memory o_D1D, o_Vf1D, o_Lf1D; // operators

  occa::memory o_Q, o_Qf; // solution and flux vals
  occa::memory o_vgeo, o_vfgeo, o_fgeo; // geometric terms
  occa::memory o_rhs, o_rhsf, o_res; // rhs and RK residual
  occa::memory o_mapPq; // node map
  
}App;

#endif
