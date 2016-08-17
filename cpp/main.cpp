#include "fem.h"

#define TIMING 0
#define TEST 0
#define USE_WADG 1 // init WADG data structures - expensive if N >> 1
#define GORDON_HALL 1 // use sphere + curve elements

// exact solution (or initial condition if t = 0)
double BesselSolution(double x, double y, double z, double time){
  double rad = sqrt(x*x + y*y + z*z);
  double uex = (rad < 1e-10) ? 1.0 : sin(M_PI*rad)/(M_PI*rad)*cos(M_PI*time); // spherical Bessel solution
  //uex = exp(.1*(x + y + z));
  return uex;
};

// exact solution (or initial condition if t = 0)
double GaussianPulse(double x, double y, double z, double time){
  double zo = 0.5;
  double rad = sqrt(x*x + y*y + (z-zo)*(z-zo));
  double uex = exp(-100.0*rad);
  //uex = exp(.1*(x + y + z));
  return uex;
};


double CavitySolution(double x, double y, double z, double time){
  return cos(M_PI*x)*cos(M_PI*y)*cos(M_PI*z)*cos(sqrt(3)*M_PI*time);
};


int main(int argc, char **argv){

#if TEST // test eigen
  //test_solve();
  test_basis();
  return 0;
#endif

  Mesh *mesh;
  int k,n, sk=0;

  // read GMSH file
  mesh = ReadGmsh3d(argv[1]);

  int KblkV = 1, KblkS = 1, KblkU = 1, KblkQ = 1, KblkQf = 1;
  if(argc >= 8){
    KblkV = atoi(argv[3]);
    KblkS = atoi(argv[4]);
    KblkU = atoi(argv[5]);
    KblkQ = atoi(argv[6]);
    KblkQf = atoi(argv[7]);
  }
  printf("Kblk (V,S,U,Q,Qf) = %d,%d,%d,%d,%d\n",KblkV,KblkS,KblkU,KblkQ,KblkQf);

  // find element-element connectivity
  FacePair3d(mesh);

  // perform start up
  StartUp3d(mesh);

  printf("%d elements in mesh\n", mesh->K);

  dfloat dt, gdt;
  /* initialize OCCA info */
  dfloat WaveInitOCCA3d(Mesh *mesh,int KblkV, int KblkS, int KblkU, int KblkQ, int KblkQf);
  dt = WaveInitOCCA3d(mesh,KblkV,KblkS,KblkU,KblkQ,KblkQf);

  double (*uexptr)(double,double,double,double) = NULL;

  // default to cavity solution
  uexptr = &CavitySolution; // if domain = [-1,1]^3

  //setupCG(mesh);

#if GORDON_HALL
  printf("initializing WADG - Gordon Hall + WADG quadrature\n");
  GordonHallSphere(mesh); printf("Using Gordon-Hall for sphere\n");
  uexptr = &BesselSolution;  // if domain = sphere at (0,0) with radius 1
#endif

  //uexptr = &GaussianPulse;  // for fancy pic

  InitQuadratureArrays(mesh);
  InitWADG_curved(mesh);  checkCurvedGeo(mesh);
  //return 0; // check properties of Gordon-Hall

  printf("dt = %17.15f\n", dt);

  //  =========== field storage (dfloat) ===========

  dfloat *Q = (dfloat*) calloc(p_Nfields*mesh->K*p_Np, sizeof(dfloat));   // 4 fields

  // set u0
  WaveSetU0(mesh,Q,0.0,uexptr);
  //WaveProjectU0(mesh,Q,0.0,uexptr);

  double L2err, relL2err;
  compute_error(mesh, 0.0, Q, uexptr, L2err, relL2err);
  printf("Initial L2 error = %6.6e\n",L2err);

  //writeVisToGMSH("u0.msh",mesh,Q,0,p_Nfields);
  //return 0;

  // =============  run solver  ================

  // load data onto GPU
  WaveSetData3d(Q);

#if TIMING
  printf("Running in TIMING mode\n");
  time_curved_kernels(mesh,10);
  //time_curved_kernels(mesh,5); // for profiling
  return 0;
#endif

  dfloat FinalTime = .1;
  if (argc >= 2){
    FinalTime = atof(argv[2]);
  }

  printf("Running until time = %f\n",FinalTime);

  //test_RK(mesh);return 0;
#if USE_WADG
  Wave_RK(mesh,FinalTime,dt,2); // run wadg curvi
#else
  Wave_RK(mesh,FinalTime,dt,0); // run standard kernels
#endif

  // unload data from GPU
  WaveGetData3d(mesh, Q);
  //writeVisToGMSH("p.msh",mesh,Q,0,p_Nfields); return 0;

  compute_error(mesh, FinalTime, Q,
                uexptr, L2err, relL2err);
  printf("N = %d, Mesh size = %f, ndofs = %d, L2 error at time %f = %6.6e\n",
	 p_N,mesh->hMax,mesh->K*p_Np,FinalTime,L2err);

  return 0;


  /* end game */
  exit(0);
}
