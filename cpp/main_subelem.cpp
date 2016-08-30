#include "fem.h"

double CavitySolution(double x, double y, double z, double time){
  return cos(M_PI*x)*cos(M_PI*y)*cos(M_PI*z)*cos(sqrt(3)*M_PI*time);
};

double PulseInitialCondition(double x, double y, double z, double time){
  double rad = sqrt(x*x + y*y + z*z);
  return exp(-16.0*rad);
};

double WaveField(double x, double y, double z){
  //return 1.0 + .25*cos(2.0*M_PI*x);
  //return 1.0 + x;
  return 1.0; // const for verification
};


int main(int argc, char **argv){

  printf("running heterogeneous subelem main\n");

  Mesh *mesh;
  int k,n, sk=0;

  // read GMSH file
  mesh = ReadGmsh3d(argv[1]);

  int KblkV = 1, KblkS = 1, KblkU = 1, KblkQ = 1, KblkQf = 1;
  if(argc >= 8){
    KblkV = atoi(argv[3]);
    KblkS = atoi(argv[4]);
    KblkU = atoi(argv[5]);
  }
  printf("for N = %d, Kblk (V,S,U,Q,Qf) = %d,%d,%d,%d,%d\n",
	 p_N,KblkV,KblkS,KblkU,KblkQ,KblkQf);

  // find element-element connectivity
  FacePair3d(mesh);

  // initialize routines and stuff
  StartUp3d(mesh);

  printf("%d elements in mesh\n", mesh->K);

  // initialize OCCA info
  dfloat WaveInitOCCA3d(Mesh *mesh,int KblkV, int KblkS, int KblkU, int KblkQ, int KblkQf);
  dfloat dt = WaveInitOCCA3d(mesh,KblkV,KblkS,KblkU,KblkQ,KblkQf);
  printf("dt = %17.15f\n", dt);

  InitQuadratureArrays(mesh); // initialize quadrature-based geofacs on host

  double (*c2_ptr)(double,double,double) = &WaveField;
  InitWADG_subelem(mesh,c2_ptr);

  //  =========== field storage (dfloat) ===========

  dfloat *Q = (dfloat*) calloc(p_Nfields*mesh->K*p_Np, sizeof(dfloat));   // 4 fields

  double (*uexptr)(double,double,double,double) = NULL;
  uexptr = &CavitySolution; // if domain = [-1,1]^3
  //uexptr = &PulseInitialCondition;

  WaveSetU0(mesh,Q,0.0,uexptr); // interpolate
  //WaveProjectU0(mesh,Q,0.0,uexptr); // L2 projection

  //writeVisToGMSH("p0.msh",mesh,Q,0,p_Nfields); return 0;

  double L2err, relL2err;
  compute_error(mesh, 0.0, Q, uexptr, L2err, relL2err);
  printf("Initial L2 error = %6.6e\n",L2err);

  // === test ===

  test_RK(mesh); return 0;


  // =============  run solver  ================

  // load data onto GPU
  WaveSetData3d(Q);

  dfloat FinalTime = .1;
  if (argc >= 2){
    FinalTime = atof(argv[2]);
  }

  printf("Running until time = %f\n",FinalTime);

  Wave_RK(mesh,FinalTime,dt,1); // run w/planar elems and WADG update

  WaveGetData3d(mesh, Q);  // unload data from GPU

  compute_error(mesh, FinalTime, Q,
                uexptr, L2err, relL2err);
  printf("N = %d, Mesh size = %f, ndofs = %d, L2 error at time %f = %6.6e\n",
	 p_N,mesh->hMax,mesh->K*p_Np,FinalTime,L2err);

  return 0;

  /* end game */
  exit(0);
}
