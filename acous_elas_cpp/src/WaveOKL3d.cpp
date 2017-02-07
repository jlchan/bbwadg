#include <stdio.h>
#include "fem.h"

#include <occa.hpp>

// switches b/w nodal and bernstein bases
#define USE_SKEW 1
#define USE_SLICE_LIFT 0 // switch on for faster behavior if N > 6

int ngeo, nvgeo, nfgeo; // number of geometric factors

// OCCA device
occa::device device;
occa::kernelInfo dgInfo;

// OCCA array for geometric factors
occa::memory c_vgeo;
occa::memory c_fgeo;

// OCCA array for index of positive trace variables wrt face nodes
occa::memory c_vmapP;
occa::memory c_Fmask;

// rhs for RK kernels
occa::memory c_rhsQ;
occa::memory c_resQ;
occa::memory c_Q;

// OCCA arrays for nodal derivative/lift matrices
occa::memory c_Dr;
occa::memory c_Ds;
occa::memory c_Dt;
occa::memory c_LIFT;

// Bernstein derivative operators
occa::memory c_Dvals4;
occa::memory c_D_ids1;
occa::memory c_D_ids2;
occa::memory c_D_ids3;
occa::memory c_D_ids4;

// Bernstein LIFT decomposition
occa::memory c_EEL_vals;
occa::memory c_EEL_ids;
occa::memory c_L0_vals;
occa::memory c_L0_ids;
occa::memory c_cEL;
occa::memory c_slice_ids; // permutation of rows

// nodal kernels
occa::kernel rk_volume;
occa::kernel rk_surface;
// bernstein kernels
occa::kernel rk_volume_bern;
occa::kernel rk_surface_bern;
// used for both nodal and Bernstein
occa::kernel rk_update;

// WADG for heterogeneous media w/planar elems
occa::kernel rk_update_WADG;
occa::memory c_Vq_reduced;
occa::memory c_Pq_reduced;

// transversely isotropic media
occa::memory c_rhoq;
occa::memory c_lambdaq;
occa::memory c_muq;
occa::memory c_c11;
occa::memory c_c12;

// ricker ptsrc
occa::memory c_fsrc;
occa::kernel rk_volume_elas;
occa::kernel rk_surface_elas;
occa::kernel rk_update_elas;

// curvilinear WADG arrays
occa::memory c_Pq, c_Prq, c_Psq, c_Ptq, c_Pfq; // V*q'*diag(wq)
occa::memory c_VfqFace, c_Vfq, c_Vq; // interp: face nodes to quad, nodes to surf + vol quad
occa::memory c_Vrq, c_Vsq, c_Vtq; // deriv VDMs
occa::memory c_vgeoq; // 9*Nq x K - array of rstxyz
occa::memory c_fgeoq; // 4*Nfq*Nfaces x K - nxyz and sJq
occa::memory c_Jq;  // J at quadrature points
occa::memory c_Jq_reduced;  // reduced quadrature for update step

// float4 packed
occa::memory c_VPstrong, c_Vskew, c_Pskew;

occa::memory c_mapPq; // map for flux rhs
occa::memory c_Qtmp; // stores rhs at quad points
occa::memory c_Qf; // stores solution values at face quad points

// ================ curvilinear specific kernels ====================

occa::kernel rk_volume_WADG;
occa::kernel rk_surface_WADG;
occa::kernel rk_volume_planar;
occa::kernel rk_surface_planar;

// index arrays for curved/planar elements
occa::memory c_KlistAll;
occa::memory c_KlistCurved;
occa::memory c_KlistCurvedPlusNbrs;
occa::memory c_KlistPlanar;
occa::memory c_Jplanar;

// new + skew form kernels - todo
occa::kernel rk_volume_WADG_Qtmp;
occa::kernel rk_surface_WADG_Qf;
occa::kernel kernel_write_quad_pts;
occa::kernel kernel_write_vol_quad4;
occa::kernel kernel_write_vol_quad6;
occa::kernel kernel_write_surf_quad;
occa::kernel rk_volume_WADG_skew;
occa::kernel rk_surface_WADG_skew;

// blocked quadrature kernels
occa::kernel rk_volume_WADG_skew_combine;
occa::kernel rk_surface_WADG_skew_combine;
occa::kernel rk_surface_WADG_face;
occa::kernel rk_surface_WADG_skew_face;

// block sizes for optimization of diff kernels
int KblkV, KblkS, KblkU, KblkQ, KblkQf;

// runtime counters
double timeV = 0.0, timeS = 0.0, timeU=0.0, timeQ = 0.0, timeQf = 0.0;

template <typename T>
void diagnose_array(Mesh *mesh,const char *message, occa::memory &c_a, int N){

  device.finish();

  T *h_a = (T*) calloc(N, sizeof(T));

  c_a.copyTo(h_a, N*sizeof(T));

  dfloat suma = 0;
  for(int n=0;n<N;++n){
    suma += h_a[n];
  }

  printf("%s: sum = %17.15g\n", message, suma);

  free(h_a);
}

// init WADG for subelem + curvilinear vars - converts arrays to quadrature arrays
// assumes Gordon Hall already has been performed and that xyz = curvilinear
void InitWADG_curved(Mesh *mesh){

  int Nq = mesh->Nq;
  int Nfq = mesh->Nfq;
  int Nfaces = mesh->Nfaces;

  // interpolate to face cubature points
  MatrixXd Vfq = mesh->Vfq;
  MatrixXd xfq,yfq,zfq;
  xfq = Vfq*mesh->x;
  yfq = Vfq*mesh->y;
  zfq = Vfq*mesh->z;

  //// store solution at quad points to reduce shared mem use
  MatrixXi mapPq;
  BuildFaceNodeMaps(mesh,xfq,yfq,zfq,mapPq);
  for (int e = 0; e < mesh->K; ++e){
    for (int i = 0; i < mesh->Nfq*p_Nfaces; ++i){
      int idP = mapPq(i,e);

      // correct mapPq for Nfields > 1
      int f = i/mesh->Nfq;
      int eNbr = mesh->EToE(e,f);
      idP -= mesh->Nfq*p_Nfaces* eNbr; // decrement
      idP += mesh->Nfq*p_Nfaces* p_Nfields*eNbr; // re-increment

      mapPq(i,e) = idP;
    }
  }

  MatrixXd Vq = mesh->Vq;
  MatrixXd Vrq = mesh->Vrq;
  MatrixXd Vsq = mesh->Vsq;
  MatrixXd Vtq = mesh->Vtq;

  MatrixXd rxq,sxq,txq,ryq,syq,tyq,rzq,szq,tzq,Jq;
  Jq = mesh->Jq;
  rxq = mesh->rxq;  ryq = mesh->ryq;  rzq = mesh->rzq;
  sxq = mesh->sxq;  syq = mesh->syq;  szq = mesh->szq;
  txq = mesh->txq;  tyq = mesh->tyq;  tzq = mesh->tzq;

  // interp to face matrix
  MatrixXd VfqFace = mesh->VfqFace;
  MatrixXd nxq,nyq,nzq,sJq;
  nxq = mesh->nxq;  nyq = mesh->nyq;  nzq = mesh->nzq;
  sJq = mesh->sJq;

  // packing vgeoq with geofacs
  //printf("nvgeo = %d\n",nvgeo);
  MatrixXd vgeoq(Nq * nvgeo,mesh->K);
  MatrixXd fgeoq(Nfq * Nfaces * nfgeo, mesh->K);
  for(int e = 0; e < mesh->K; ++e){
    int off = 0;
    vgeoq.col(e).middleRows(off,Nq) = rxq.col(e);      off += Nq;
    vgeoq.col(e).middleRows(off,Nq) = sxq.col(e);      off += Nq;
    vgeoq.col(e).middleRows(off,Nq) = txq.col(e);      off += Nq;
    vgeoq.col(e).middleRows(off,Nq) = ryq.col(e);      off += Nq;
    vgeoq.col(e).middleRows(off,Nq) = syq.col(e);      off += Nq;
    vgeoq.col(e).middleRows(off,Nq) = tyq.col(e);      off += Nq;
    vgeoq.col(e).middleRows(off,Nq) = rzq.col(e);      off += Nq;
    vgeoq.col(e).middleRows(off,Nq) = szq.col(e);      off += Nq;
    vgeoq.col(e).middleRows(off,Nq) = tzq.col(e);

    int foffgeo = 0; // packed global storage increment
    int foff = 0; // nxyz and sJ increment
    for (int f = 0; f < Nfaces; ++f){
      fgeoq.col(e).middleRows(foffgeo,Nfq) = nxq.col(e).middleRows(foff,Nfq);      foffgeo += Nfq;
      fgeoq.col(e).middleRows(foffgeo,Nfq) = nyq.col(e).middleRows(foff,Nfq);      foffgeo += Nfq;
      fgeoq.col(e).middleRows(foffgeo,Nfq) = nzq.col(e).middleRows(foff,Nfq);      foffgeo += Nfq;
      fgeoq.col(e).middleRows(foffgeo,Nfq) = sJq.col(e).middleRows(foff,Nfq);      foffgeo += Nfq;
      foff += Nfq; // increment face offset
    }
  }

  MatrixXd invM = mesh->V * mesh->V.transpose();
  MatrixXd Pq   = invM*Vq.transpose()*mesh->wq.asDiagonal();
  MatrixXd Prq  = invM*Vrq.transpose()*mesh->wq.asDiagonal();
  MatrixXd Psq  = invM*Vsq.transpose()*mesh->wq.asDiagonal();
  MatrixXd Ptq  = invM*Vtq.transpose()*mesh->wq.asDiagonal();
  MatrixXd Pfq  = invM*Vfq.transpose()*mesh->wfq.asDiagonal();

  dgInfo.addDefine("p_Nq",mesh->Nq); // vol quadrature
  dgInfo.addDefine("p_Nfq",mesh->Nfq); // surf quadrature for one face
  dgInfo.addDefine("p_NfqNfaces",mesh->Nfq * Nfaces); // surf quadrature

  if(!strcmp(device.mode().c_str(), "CUDA")){
    cout << " Adding CUDA optimization " << endl;
    dgInfo.addCompilerFlag("--ftz=true");
    dgInfo.addCompilerFlag("--prec-div=false");
    dgInfo.addCompilerFlag("--prec-sqrt=false");
    dgInfo.addCompilerFlag("--use_fast_math");
    dgInfo.addCompilerFlag("--fmad=true"); // compiler option for cuda
    dgInfo.addCompilerFlag("-Xptxas --def-store-cache=cs");
    dgInfo.addCompilerFlag("-Xptxas --force-store-cache=cs");
  }

  if(!strcmp(device.mode().c_str(), "OpenCL")){
    cout << " Adding OpenCL optimization " << endl;
    dgInfo.addCompilerFlag("-cl-strict-aliasing");
    dgInfo.addCompilerFlag("-cl-mad-enable");
    dgInfo.addCompilerFlag("-cl-no-signed-zeros");
    dgInfo.addCompilerFlag("-cl-unsafe-math-optimizations");
    dgInfo.addCompilerFlag("-cl-finite-math-only");
    dgInfo.addCompilerFlag("-cl-fast-relaxed-math");
  }


  // reference operators
  setOccaArray(mesh,Vq,c_Vq);
  setOccaArray(mesh,Vrq,c_Vrq);
  setOccaArray(mesh,Vsq,c_Vsq);
  setOccaArray(mesh,Vtq,c_Vtq);


  setOccaArray(mesh,Pq,c_Pq);
  setOccaArray(mesh,Prq,c_Prq);
  setOccaArray(mesh,Psq,c_Psq);
  setOccaArray(mesh,Ptq,c_Ptq);

  // float4 packing
  VectorXd Vqvec(Map<VectorXd>(Vq.data(), Vq.cols()*Vq.rows()));
  VectorXd Vrqvec(Map<VectorXd>(Vrq.data(), Vq.cols()*Vq.rows()));
  VectorXd Vsqvec(Map<VectorXd>(Vsq.data(), Vq.cols()*Vq.rows()));
  VectorXd Vtqvec(Map<VectorXd>(Vtq.data(), Vq.cols()*Vq.rows()));
  VectorXd Pqvec(Map<VectorXd>(Pq.data(), Pq.cols()*Pq.rows()));
  VectorXd Prqvec(Map<VectorXd>(Prq.data(), Pq.cols()*Pq.rows()));
  VectorXd Psqvec(Map<VectorXd>(Psq.data(), Pq.cols()*Pq.rows()));
  VectorXd Ptqvec(Map<VectorXd>(Ptq.data(), Pq.cols()*Pq.rows()));
  VectorXd Vskew(4*Vq.cols()*Vq.rows());    Vskew    << Vrqvec,Vsqvec,Vtqvec,Vqvec;
  VectorXd Pskew(4*Pq.cols()*Pq.rows());    Pskew    << Prqvec,Psqvec,Ptqvec,Pqvec;
  setOccaArray(mesh,Vskew,c_Vskew);
  setOccaArray(mesh,Pskew,c_Pskew);

  setOccaArray(mesh,Pfq,c_Pfq);

  // for face flux points
  setOccaIntArray(mesh,mapPq,c_mapPq);
  MatrixXd Qf(p_Nfields*Nfq*Nfaces,mesh->K);
  setOccaArray(mesh,Qf,c_Qf);

  // save grad p, u,v,w @ quad pts
  MatrixXd Qtmp(6*mesh->Nq,mesh->K);
  setOccaArray(mesh,Qtmp,c_Qtmp);

  // geometric factors
  setOccaArray(mesh,vgeoq,c_vgeoq);
  setOccaArray(mesh,fgeoq,c_fgeoq);
  setOccaArray(mesh,Jq,c_Jq);

  // specific vol/surface routines for planar elements
  VectorXd Jplanar(mesh->K);
  for (int e = 0; e < mesh->K; ++e){
    Jplanar(e) = Jq(0,e);
  }
  VectorXi KlistAll(mesh->K);
  for (int e = 0; e < mesh->K; ++e){
    KlistAll(e) = e;
  }
  setOccaArray(mesh,Jplanar,c_Jplanar);
  setOccaIntArray(mesh,KlistAll,c_KlistAll);
  setOccaIntArray(mesh,mesh->KlistCurved,c_KlistCurved);
  setOccaIntArray(mesh,mesh->KlistPlanar,c_KlistPlanar);
  setOccaIntArray(mesh,mesh->KlistCurvedPlusNbrs,c_KlistCurvedPlusNbrs);

  // ======= use reduced strength quadrature for WADG update step

  VectorXd rq,sq,tq,wq;
  tet_cubature(min(21,2*p_N+1),rq,sq,tq,wq); // 21 = max quadrature degree
  MatrixXd Vqtmp = Vandermonde3D(p_N,rq,sq,tq);
  MatrixXd Vq_reduced = mrdivide(Vqtmp,mesh->V);
  dgInfo.addDefine("p_Nq_reduced",Vq_reduced.rows()); // vol quadrature

  MatrixXd Vrqtmp,Vsqtmp,Vtqtmp;
  GradVandermonde3D(p_N,rq,sq,tq,Vrqtmp,Vsqtmp,Vtqtmp);
  MatrixXd Vrq_reduced = mrdivide(Vrqtmp,mesh->V);
  MatrixXd Vsq_reduced = mrdivide(Vsqtmp,mesh->V);
  MatrixXd Vtq_reduced = mrdivide(Vtqtmp,mesh->V);

  MatrixXd Pq_reduced = invM*Vq_reduced.transpose()*wq.asDiagonal();
  MatrixXd Prq_reduced  = invM*Vrq_reduced.transpose()*wq.asDiagonal();
  MatrixXd Psq_reduced  = invM*Vsq_reduced.transpose()*wq.asDiagonal();
  MatrixXd Ptq_reduced  = invM*Vtq_reduced.transpose()*wq.asDiagonal();

  //MatrixXd Jq_reduced = Vq_reduced * (Pq*Jq); // project J onto PN + interp to reduced quad. note: this is not exact!!

  MatrixXd Jq_reduced(Vq_reduced.rows(),mesh->K);
  for (int e = 0; e < mesh->K; ++e){
    MatrixXd vgeo = vgeofacs3d(mesh->x.col(e),
			       mesh->y.col(e),
			       mesh->z.col(e),
			       Vrq_reduced,Vsq_reduced,Vtq_reduced);
    Jq_reduced.col(e)  = vgeo.col(9);
  }

  // scale Jq by varying coeff c^2
  MatrixXd c2q(Vq_reduced.rows(),mesh->K);
  c2q.fill(1.0);
#if 0
  for (int e = 0; e < mesh->K; ++e){
    VectorXd xq = Vq_reduced*mesh->x.col(e);
    VectorXd yq = Vq_reduced*mesh->y.col(e);
    VectorXd zq = Vq_reduced*mesh->z.col(e);
    for (int i = 0; i < mesh->Nq; ++i){
      double x = xq(i);
      double y = yq(i);
      double z = zq(i);
      double rad = sqrt(x*x + y*y + z*z);
      double c2 = 1 + .5*cos(4*M_PI*rad); // for pretty pictures...
      c2q(i,e) *= c2;
    }
  }
#endif
  //setOccaArray(mesh,c2q,c_c2q);
  setOccaArray(mesh,Jq_reduced,c_Jq_reduced);
  setOccaArray(mesh,Pq_reduced,c_Pq_reduced);
  setOccaArray(mesh,Vq_reduced,c_Vq_reduced);
  setOccaArray(mesh,Vfq,c_Vfq);
  setOccaArray(mesh,VfqFace,c_VfqFace);

  // ============ build kernels

  std::string src = "okl/WaveKernelsWADG.okl";
  printf("Building WADG kernels from %s\n",src.c_str());

  // strong form curved kernels
  rk_volume_WADG  = device.buildKernelFromSource(src.c_str(), "rk_volume_WADG", dgInfo);
  rk_surface_WADG = device.buildKernelFromSource(src.c_str(), "rk_surface_WADG", dgInfo);
  rk_update_WADG  = device.buildKernelFromSource(src.c_str(), "rk_update_WADG", dgInfo);
  //rk_update_WADG  = device.buildKernelFromSource(src.c_str(), "rk_update_WADG_block", dgInfo);
  //rk_update_WADG  = device.buildKernelFromSource(src.c_str(), "rk_update_WADG_old", dgInfo);

  // variants for strong form
  kernel_write_vol_quad4 = device.buildKernelFromSource(src.c_str(), "kernel_write_vol_quad4", dgInfo);
  rk_volume_WADG_Qtmp = device.buildKernelFromSource(src.c_str(), "rk_volume_WADG_Qtmp", dgInfo);
  rk_surface_WADG_Qf = device.buildKernelFromSource(src.c_str(), "rk_surface_WADG_Qf", dgInfo);
  kernel_write_surf_quad = device.buildKernelFromSource(src.c_str(), "kernel_write_surf_quad", dgInfo);
  rk_surface_WADG_face = device.buildKernelFromSource(src.c_str(), "rk_surface_WADG_face", dgInfo);

  // build strong-weak form kernels: split form
  kernel_write_quad_pts = device.buildKernelFromSource(src.c_str(), "kernel_write_quad_pts", dgInfo);
  kernel_write_vol_quad6 = device.buildKernelFromSource(src.c_str(), "kernel_write_vol_quad6", dgInfo);
  rk_volume_WADG_skew = device.buildKernelFromSource(src.c_str(), "rk_volume_WADG_skew", dgInfo);
  rk_surface_WADG_skew = device.buildKernelFromSource(src.c_str(), "rk_surface_WADG_skew", dgInfo);

  // monolithic strong-weak kernels
  rk_volume_WADG_skew_combine = device.buildKernelFromSource(src.c_str(), "rk_volume_WADG_skew_combine", dgInfo);
    rk_surface_WADG_skew_combine = device.buildKernelFromSource(src.c_str(), "rk_surface_WADG_skew_combine", dgInfo);

  // ================================= planar kernels
  rk_volume_planar  = device.buildKernelFromSource(src.c_str(), "rk_volume_planar", dgInfo);
  rk_surface_planar = device.buildKernelFromSource(src.c_str(), "rk_surface_planar", dgInfo);

  printf("initialized WADG for curvilinear + heterogeneous media.\n");
}

void InitWADG_subelem(Mesh *mesh,double(*c2_ptr)(double,double,double)){

  int Nq = mesh->Nq;

  // ======= use reduced strength quadrature for WADG update

  int Nq1D = p_N+1; // N+1 for GQ(1,0) in vert, N+1 for GQ(0,0)
  int Nq2 = Nq1D*Nq1D;
  int Nq3 = Nq2*Nq1D;

  VectorXd rq,sq,tq,wq;
  //  if (0){
  if (2*p_N+1 <= 15){
    tet_cubature(min(21,2*p_N+1),rq,sq,tq,wq); // 21 = max quadrature degree
  }else{

    // default to TP quadrature for now
    VectorXd a1D,wa,c1D,wc;
    JacobiGQ(p_N,0,0,a1D,wa);
    JacobiGQ(p_N,2,0,c1D,wc);
    rq.resize(Nq3);    sq.resize(Nq3);
    tq.resize(Nq3);    wq.resize(Nq3);
    int sk = 0;
    for (int k = 0; k < Nq1D; ++k){
      for (int i = 0; i < Nq1D; ++i){
        for (int j = 0; j < Nq1D; ++j){
          rq(sk) = a1D(i); sq(sk) = a1D(j);
	  tq(sk) = c1D(k); wq(sk) = wa(i)*wa(j)*wc(k);
          ++sk;
        }
      }
    }
  }
  MatrixXd Vqtmp = Vandermonde3D(p_N,rq,sq,tq);
  MatrixXd Vq_reduced = mrdivide(Vqtmp,mesh->V);

  // array of wavefield at quadrature points
  MatrixXd rhoq(Vq_reduced.rows(),mesh->K);
  MatrixXd lambdaq(Vq_reduced.rows(),mesh->K);
  MatrixXd muq(Vq_reduced.rows(),mesh->K);
  MatrixXd c11(Vq_reduced.rows(),mesh->K);
  MatrixXd c12(Vq_reduced.rows(),mesh->K);

  MatrixXd fsrcq(Vq_reduced.rows(),mesh->K);

  for (int e = 0; e < mesh->K; ++e){

    // locally interpolate to cubature points
    VectorXd xq(Vq_reduced.rows()), yq(Vq_reduced.rows()), zq(Vq_reduced.rows());
    xq = Vq_reduced*(mesh->x.col(e));
    yq = Vq_reduced*(mesh->y.col(e));
    zq = Vq_reduced*(mesh->z.col(e));

    for (int i = 0; i < Vq_reduced.rows(); ++i){
      double weight = (*c2_ptr)(xq(i),yq(i),zq(i));
      double mu = 1; double lambda = 1;
      rhoq(i,e) = 1.0;
      lambdaq(i,e) = lambda;
      muq(i,e) = mu + weight; // constant mu

#if 0
      double A = 2*mu+lambda;
      double B = lambda;

      // transverse isotropy
      c11(i,e) = A + weight;
      c12(i,e) = B + weight;
      if (zq(i) > 0){ // vertical part
	c11(i,e) *= 1.0/3.0; // 2*mu + lambda
	c12(i,e) *= 0.5;
      }
#else
      // isotropic but discontinuous
      double BB = 2.0;
      if (zq(i) < 0){
        muq(i,e) = BB + weight;
        lambdaq(i,e) = BB;
      }

#endif
      // smoothed ricker pulse
      double a = 100.0;
      double x0 = 0.0;
      double y0 = 0.0;
      double z0 = .1;
      double dx = xq(i) - x0;
      double dy = yq(i) - y0;
      double dz = zq(i) - z0;
      double r2 = dx*dx + dy*dy + dz*dz;
      fsrcq(i,e) = exp(-a*a*r2);
    }
  }

#if 0 // take local average of coeffs

  for (int e = 0; e < mesh->K; ++e){
    VectorXd wJq = wq.array()*Jq_reduced.col(e).array();
    double muavg = wJq.dot(muq.col(e));
  }
#endif 0

  //  cout << "c11" << endl  << c11.col(0) << endl;
  //  cout << "c12" << endl  << c12.col(0) << endl;
  //  cout << "mu" << endl  << muq.col(0) << endl;
  //  cout << "lambda" << endl  << lambdaq.col(0) << endl;

  dgInfo.addDefine("p_tau_v",1.0); // velocity penalty
  dgInfo.addDefine("p_tau_s",1.0);

  MatrixXd invM = mesh->V*mesh->V.transpose();
  MatrixXd Pq_reduced = invM*Vq_reduced.transpose()*wq.asDiagonal();

  dgInfo.addDefine("p_Nq_reduced",Vq_reduced.rows()); // for update step quadrature

  // not used in WADG subelem - just to compile other kernels
  dgInfo.addDefine("p_NfqNfaces",mesh->Nfq * p_Nfaces); // surf quadrature
  dgInfo.addDefine("p_Nq",mesh->Nq); // vol quadrature
  dgInfo.addDefine("p_Nfq",mesh->Nfq); // surf quadrature for one face

  setOccaArray(mesh,rhoq,c_rhoq);
  setOccaArray(mesh,lambdaq,c_lambdaq);
  setOccaArray(mesh,muq,c_muq);
  setOccaArray(mesh,c11,c_c11);
  setOccaArray(mesh,c12,c_c12);

  // smoothed ricker src
  fsrcq *= 1.0/fsrcq.array().abs().maxCoeff(); // normalize to 1
  MatrixXd fsrc = Pq_reduced * fsrcq;
  setOccaArray(mesh,fsrc,c_fsrc);

  setOccaArray(mesh,Pq_reduced,c_Pq_reduced);
  setOccaArray(mesh,Vq_reduced,c_Vq_reduced);

  // ======== build kernels

  std::string src = "okl/ElasKernelsWADG.okl";
  printf("Building heterogeneous wave propagation WADG kernel from %s\n",src.c_str());
  rk_update_elas  = device.buildKernelFromSource(src.c_str(), "rk_update_elas", dgInfo);
  rk_volume_elas  = device.buildKernelFromSource(src.c_str(), "rk_volume_elas", dgInfo);
  rk_surface_elas = device.buildKernelFromSource(src.c_str(), "rk_surface_elas", dgInfo);

}

// applies Gordon Hall to a sphere - just changes xyz coordinates.
//                                   run InitWADG to recompute geofacs
void GordonHallSphere(Mesh *mesh){

  //cout << "Fmask = " << endl << Fmask << endl;
  std::ofstream myfile;
#if 0
  myfile.open("matlabPreCurvedMesh.m");
  myfile << "EToE = zeros(" << mesh->K << "," << p_Nfaces << ");" << std::endl;
  myfile << "EToF = zeros(" << mesh->K << "," << p_Nfaces << ");" << std::endl;
  for(int e=0; e<mesh->K; ++e){
    myfile << "EToE(" << e+1 << ",:) = [";
    for(int f = 0; f< p_Nfaces; ++f){
      myfile << mesh->EToE(e,f)+1 << " ";
    }
    myfile << "];" << std::endl;

    myfile << "EToF(" << e+1 << ",:) = [";
    for(int f = 0; f< p_Nfaces; ++f){
      myfile << mesh->EToF(e,f)+1 << " ";
    }
    myfile << "];" << std::endl;
  }

  myfile << "x = [" << mesh->x << "];" << endl;
  myfile << "y = [" << mesh->y << "];" << endl;
  myfile << "z = [" << mesh->z << "];" << endl;
  myfile << "K = " << mesh->K << ";"<< std::endl;

  printf("wrote matlab pre-curved mesh file\n");
  myfile.close();
#endif

  double tol = 1e-6;

  MatrixXi eids(p_N+1,6);
  VectorXi ske(6); ske.fill(0);
  for(int i = 0; i < p_Np; ++i){
    double r = mesh->r(i);
    double s = mesh->s(i);
    double t = mesh->t(i);
    int edge = -1;
    if (fabs(s+1.0) + fabs(t+1.0) < tol){
      edge = 0;
      eids(ske(edge),edge) = i;      ++ske(edge);
    }
    if (fabs(r+s) < tol){
      edge = 1;
      eids(ske(edge),edge) = i;      ++ske(edge);
    }
    if (fabs(r+1.0) + fabs(t+1.0) < tol){
      edge = 2;
      eids(ske(edge),edge) = i;      ++ske(edge);
    }
    if (fabs(r+1.0) + fabs(s+1.0) < tol){
      edge = 3;
      eids(ske(edge),edge) = i;      ++ske(edge);
    }
    if (fabs(r+t) < tol){
      edge = 4;
      eids(ske(edge),edge) = i;      ++ske(edge);
    }
    if (fabs(s+t) < tol){
      edge = 5;
      eids(ske(edge),edge) = i;      ++ske(edge);
    }
  }

  int NpFace = max(0,(p_N-1)*(p_N-2)/2); // interior face dofs
  MatrixXi fids(NpFace,p_Nfaces); fids.fill(0);
  int fsk1 = 0,fsk2 = 0, fsk3 = 0, fsk4 = 0;
  int sk = 0;
  for (int k = 0; k <= p_N; ++k){
    for (int j = 0; j <= p_N-k; ++j){
      for (int i = 0; i <= p_N-j-k; ++i){

	// save interior face ids
	if (k==0 && i > 0 && j > 0 && i+j<p_N){
	  fids(fsk1,0) = sk; ++fsk1;
	}
	if (j==0 && i > 0 && k > 0 && i+k<p_N){
	  fids(fsk2,1) = sk; ++fsk2;
	}
	if (i+j+k==p_N && i > 0 && j > 0 && i+j<p_N){
	  fids(fsk3,2) = sk; ++fsk3;
	}
	if (i==0 && k > 0 && j > 0 && j+k < p_N){
	  fids(fsk4,3) = sk; ++fsk4;
	}
	++sk;

      }
    }
  }
  //cout << "eids = " << endl << eids << endl;
  //cout << "fids = " << endl << fids << endl;

  // hier-modal basis
  MatrixXd vertVfull, edgeVfull, faceVfull;
  VandermondeHier(p_N,mesh->r,mesh->s,mesh->t,
		  vertVfull,edgeVfull,faceVfull);
  MatrixXd vertV, edgeV, faceV;
  edgeV.resize(edgeVfull.cols(),edgeVfull.cols());
  faceV.resize(faceVfull.cols(),faceVfull.cols());

  // extract submatrix rows corresponding
  // to interior edge or interior face nodes
  sk = 0;
  for (int edge = 0; edge < 6; ++edge){
    for (int i = 0; i < p_N-1; ++i){
      int eid = eids(i+1,edge); // ignore 0, N+1 = vertex nodes
      //printf("eid(%d,%d) = %d\n",edge,i,eid);
      edgeV.row(sk) = edgeVfull.row(eid);
      ++sk;
    }
  }

  sk = 0;
  for (int face = 0; face < p_Nfaces; ++face){
    for (int i = 0; i < (p_N-2)*(p_N-1)/2; ++i){
      int fid = fids(i,face);
      faceV.row(sk) = faceVfull.row(fid);
      ++sk;
    }
  }

  // compute blending matrices for edges/faces
  MatrixXd edgeBlend = mrdivide(edgeVfull,edgeV);
  MatrixXd faceBlend = mrdivide(faceVfull,faceV);

  //  cout << "edgeV = " << endl << edgeV << endl;
  //  cout << "faceV = " << endl << faceV << endl;
  //  cout << "edgeB = " << endl << edgeBlend << endl;
  //  cout << "faceB = " << endl << faceBlend << endl;

  VectorXi edgeCurvedElems(mesh->K);
  edgeCurvedElems.fill(0);

  MatrixXd xc = mesh->x;
  MatrixXd yc = mesh->y;
  MatrixXd zc = mesh->z;
  // snap edge nodes to sphere and blend
  for (int e = 0; e < mesh->K; ++e){

    for (int edge = 0; edge < 6; ++edge){
      VectorXi ev(2);
      ev(0) = eids(0,edge);      ev(1) = eids(p_N,edge);
      VectorXd xb(2),yb(2),zb(2);
      xb(0) = mesh->x(ev(0),e);  xb(1) = mesh->x(ev(1),e);
      yb(0) = mesh->y(ev(0),e);  yb(1) = mesh->y(ev(1),e);
      zb(0) = mesh->z(ev(0),e);  zb(1) = mesh->z(ev(1),e);
      VectorXd r2 = xb.array()*xb.array() +  yb.array()*yb.array() +  zb.array()*zb.array();
      VectorXd radiu = r2.array().sqrt();
      VectorXd rerr = (radiu.array() - 1.0).array().abs();
      if (rerr.maxCoeff() < tol){ // if edge is on sphere
	edgeCurvedElems(e) = 1;

	for (int i = 0; i < p_N + 1; ++i){
	  int eid = eids(i,edge);
	  double xb = mesh->x(eid,e);
	  double yb = mesh->y(eid,e);
	  double zb = mesh->z(eid,e);
	  double rb = sqrt(xb*xb + yb*yb + zb*zb);

	  // snap nodes to unit radius sphere
	  xc(eid,e) = xb/rb;
	  yc(eid,e) = yb/rb;
	  zc(eid,e) = zb/rb;
	}

      }//if on sphere

    } //edge
  }// e

  // blend edge displacements
  int NpIntEdge = max(p_N-1,0);
  for (int e = 0; e < mesh->K; ++e){
    if (edgeCurvedElems(e)==1){
      // blend displacements into interior of elements
      VectorXd dx(NpIntEdge*6);
      VectorXd dy(NpIntEdge*6);
      VectorXd dz(NpIntEdge*6);
      int sk = 0;
      for (int edge = 0; edge < 6; ++edge){
	for (int i = 0; i < NpIntEdge; ++i){
	  int eid = eids(i+1,edge); // skip first/last edge id
	  dx(sk) = xc(eid,e) - mesh->x(eid,e);
	  dy(sk) = yc(eid,e) - mesh->y(eid,e);
	  dz(sk) = zc(eid,e) - mesh->z(eid,e);
	  ++sk;
	}
      }

      //printf("edgeBlend size = %d, %d, dx size = %d\n", edgeBlend.rows(),edgeBlend.cols(),dx.rows());
      VectorXd newx = mesh->x.col(e) + edgeBlend * dx;
      VectorXd newy = mesh->y.col(e) + edgeBlend * dy;
      VectorXd newz = mesh->z.col(e) + edgeBlend * dz;
      mesh->x.col(e) = newx;
      mesh->y.col(e) = newy;
      mesh->z.col(e) = newz;
    }
  }

  // if there are face interior nodes
  VectorXi faceCurvedElems(mesh->K);
  faceCurvedElems.fill(0);

  if (p_N > 2){
    faceCurvedElems.fill(0);
    xc = mesh->x;
    yc = mesh->y;
    zc = mesh->z;
    // snap face nodes to sphere
    for (int e = 0; e < mesh->K; ++e){
      for (int face = 0; face < p_Nfaces; ++face){
	VectorXi ev(3);
	ev(0) = mesh->Fmask(0,face); ev(1) = mesh->Fmask(p_N,face); ev(2) = mesh->Fmask(p_Nfp-1,face);
	VectorXd xb(3),yb(3),zb(3);
	xb(0) = mesh->x(ev(0),e);  xb(1) = mesh->x(ev(1),e);  xb(2) = mesh->x(ev(2),e);
	yb(0) = mesh->y(ev(0),e);  yb(1) = mesh->y(ev(1),e);  yb(2) = mesh->y(ev(2),e);
	zb(0) = mesh->z(ev(0),e);  zb(1) = mesh->z(ev(1),e);  zb(2) = mesh->z(ev(2),e);
	VectorXd r2 = xb.array()*xb.array() +  yb.array()*yb.array() +  zb.array()*zb.array();
	VectorXd radiu = r2.array().sqrt();
	VectorXd rerr = (radiu.array() - 1.0).array().abs();
	if (rerr.maxCoeff() < tol){ // if face lies on sphere
	  faceCurvedElems(e) = 1;
	  for (int i = 0; i < p_Nfp; ++i){
	    int fid = mesh->Fmask(i,face);
	    double xb = mesh->x(fid,e);
	    double yb = mesh->y(fid,e);
	    double zb = mesh->z(fid,e);
	    double rb = sqrt(xb*xb + yb*yb + zb*zb);

	    // snap nodes to unit sphere
	    xc(fid,e) = xb/rb;
	    yc(fid,e) = yb/rb;
	    zc(fid,e) = zb/rb;
	  }
	}

      }
    }

    // blend face displacements
    int NpIntFace = fids.rows();
    for (int e = 0; e < mesh->K; ++e){
      if (faceCurvedElems(e)==1){
	// blend displacements into interior of elements
	VectorXd dx(NpIntFace * p_Nfaces);
	VectorXd dy(NpIntFace * p_Nfaces);
	VectorXd dz(NpIntFace * p_Nfaces);

	int sk = 0;
	for (int f = 0; f < p_Nfaces; ++f){
	  for (int i = 0; i < fids.rows(); ++i){
	    int fid = fids(i,f);
	    dx(sk) = xc(fid,e) - mesh->x(fid,e);
	    dy(sk) = yc(fid,e) - mesh->y(fid,e);
	    dz(sk) = zc(fid,e) - mesh->z(fid,e);
	    ++sk;
	  }
	}

	// apply zero BCs to non-boundary face displacements
	sk = 0;
	for (int f = 0; f < p_Nfaces; ++f){
	  for (int i = 0; i < fids.rows(); ++i){
	    // if not boundary, don't displace nodes after edge blend
	    if (mesh->EToF(e,f)!=f){
	      dx(sk) = 0.0;
	      dy(sk) = 0.0;
	      dz(sk) = 0.0;
	    }
	    ++sk;
	  }
	}

	VectorXd newx = mesh->x.col(e) + faceBlend * dx;
	VectorXd newy = mesh->y.col(e) + faceBlend * dy;
	VectorXd newz = mesh->z.col(e) + faceBlend * dz;
	mesh->x.col(e) = newx;
	mesh->y.col(e) = newy;
	mesh->z.col(e) = newz;
      }
    }
  }

  // list curved/noncurved elements
  VectorXi isCurved(mesh->K); isCurved.fill(0);
  VectorXi isCurvedPlusNbrs(mesh->K); isCurvedPlusNbrs.fill(0);

  for (int e = 0; e < mesh->K; ++e){
    if (edgeCurvedElems(e)==1 || faceCurvedElems(e)==1){
      isCurved(e)=1;
    }
    for (int f = 0; f < p_Nfaces; ++f){
      int enbr = mesh->EToE(e,f);
      isCurvedPlusNbrs(enbr) = 1;
    }
  }
  unsigned int KCurved = isCurved.sum();
  unsigned int KPlanar = mesh->K - KCurved;
  unsigned int KCurvedPlusNbrs = isCurvedPlusNbrs.sum();

  mesh->KCurved = KCurved;
  mesh->KPlanar = KPlanar;
  mesh->KCurvedPlusNbrs = KCurvedPlusNbrs;

  VectorXi KlistCurved(max(1,KCurved)); KlistCurved.fill(0);
  VectorXi KlistPlanar(max(1,KPlanar)); KlistPlanar.fill(0);
  VectorXi KlistCurvedPlusNbrs(max(1,KCurvedPlusNbrs)); KlistCurvedPlusNbrs.fill(0);
  int sknc = 0, skc = 0, skcnbrs = 0;
  for (int e = 0; e < mesh->K; ++e){
    if (isCurved(e)){
      KlistCurved(skc) = e;
      ++skc;
    }else{
      KlistPlanar(sknc) = e;
      ++sknc;
    }
    if (isCurvedPlusNbrs(e)){
      KlistCurvedPlusNbrs(skcnbrs) = e;
      ++skcnbrs;
    }
  }
  printf("%d curved elems, %d curved + nbrs, %d planar elems, %d total elems\n",
	 KCurved,KCurvedPlusNbrs,KPlanar,mesh->K);

  mesh->KlistCurved = KlistCurved;
  mesh->KlistPlanar = KlistPlanar;
  mesh->KlistCurvedPlusNbrs = KlistCurvedPlusNbrs;

#if 0

  myfile.open("matlabCurvedMesh.m");
  myfile << "EToE = zeros(" << mesh->K << "," << p_Nfaces << ");" << std::endl;
  myfile << "EToF = zeros(" << mesh->K << "," << p_Nfaces << ");" << std::endl;
  for(int e=0; e<mesh->K; ++e){
    myfile << "EToE(" << e+1 << ",:) = [";
    for(int f = 0; f< p_Nfaces; ++f){
      myfile << mesh->EToE(e,f)+1 << " ";
    }
    myfile << "];" << std::endl;

    myfile << "EToF(" << e+1 << ",:) = [";
    for(int f = 0; f< p_Nfaces; ++f){
      myfile << mesh->EToF(e,f)+1 << " ";
    }
    myfile << "];" << std::endl;
  }

  myfile << "x = [" << mesh->x << "];" << endl;
  myfile << "y = [" << mesh->y << "];" << endl;
  myfile << "z = [" << mesh->z << "];" << endl;
  myfile << "K = " << mesh->K << ";"<< std::endl;

  printf("wrote matlab curved mesh file\n");
  myfile.close();
#endif
}

void checkCurvedGeo(Mesh *mesh){

  MatrixXd x = mesh->x;
  MatrixXd y = mesh->y;
  MatrixXd z = mesh->z;

  MatrixXd Vq = mesh->Vq;
  MatrixXd Vrq = mesh->Vrq;
  MatrixXd Vsq = mesh->Vsq;
  MatrixXd Vtq = mesh->Vtq;

  MatrixXd rxq = mesh->rxq;  MatrixXd sxq = mesh->sxq;  MatrixXd txq = mesh->txq;
  MatrixXd ryq = mesh->ryq;  MatrixXd syq = mesh->syq;  MatrixXd tyq = mesh->tyq;
  MatrixXd rzq = mesh->rzq;  MatrixXd szq = mesh->szq;  MatrixXd tzq = mesh->tzq;

  // check volume of geometry ~ 4/3 * pi
  double vol = 4.0/3.0 * 3.141592653589793;
  MatrixXd wq = mesh->wq;
  double errK = 0.0;
  for (int e = 0; e < mesh->K; ++e){
    MatrixXd Jq = mesh->Jq.col(e);
    errK += (wq.array()*Jq.array()).sum();
  }
  printf("volume of geom = %f, exact = %f, err = %g\n",errK,vol,fabs(errK-vol));

  // check to make sure linears contained in space
  MatrixXd u = 2.0*x + 3.0*y + 4.0*z;
  MatrixXd dudx =
    rxq.array()*(Vrq*u).array() +
    sxq.array()*(Vsq*u).array() +
    txq.array()*(Vtq*u).array();
  MatrixXd dudy =
    ryq.array()*(Vrq*u).array() +
    syq.array()*(Vsq*u).array() +
    tyq.array()*(Vtq*u).array();
  MatrixXd dudz =
    rzq.array()*(Vrq*u).array() +
    szq.array()*(Vsq*u).array() +
    tzq.array()*(Vtq*u).array();

  u.fill(2.0);
  double err = (dudx - u).norm();
  u.fill(3.0);
  err += (dudy - u).norm();
  u.fill(4.0);
  err += (dudz - u).norm();
  printf("err in derivative of 2x + 3y + 4z = %g\n", err);

  double maxMinDist = 0.0;
  int emax = -1, fmax = -1;
  for (int e = 0; e < mesh->K; ++e){
    for (int f = 0; f < p_Nfaces; ++f){
      int enbr = mesh->EToE(e,f);
      int fnbr = mesh->EToF(e,f);

      for (int i = 0; i < p_Nfp; ++i){
	double mindist = 100.0;
	int fid1 = mesh->Fmask(i,f);
	for (int j = 0; j < p_Nfp; ++j){
	  int fid2 = mesh->Fmask(j,fnbr);
	  double dx = x(fid1,e) - x(fid2,enbr);
	  double dy = y(fid1,e) - y(fid2,enbr);
	  double dz = z(fid1,e) - z(fid2,enbr);

	  double dist = sqrt(dx*dx + dy*dy + dz*dz);
	  mindist = min(mindist,dist);
	}
	if (max(maxMinDist,mindist) > maxMinDist){
	  emax = e; fmax = f;
	}
	maxMinDist = max(maxMinDist,mindist);

      }
    }
  }
  printf("max min nodal distance = %g on elem %d, face %d\n",maxMinDist,emax,fmax);

}


// set occa array:: cast to dfloat.
void setOccaArray(Mesh *mesh, MatrixXd A, occa::memory &c_A){
  int r = A.rows();
  int c = A.cols();
  dfloat *f_A = (dfloat*)malloc(r*c*sizeof(dfloat));
  Map<MatrixXdf >(f_A,r,c) = A.cast<dfloat>();
  c_A = device.malloc(r*c*sizeof(dfloat),f_A);
  free(f_A);
}

// set occa array:: cast to dfloat
void setOccaIntArray(Mesh *mesh,MatrixXi A, occa::memory &c_A){
  int r = A.rows();
  int c = A.cols();
  int *f_A = (int*)malloc(r*c*sizeof(int));
  Map<MatrixXi >(f_A,r,c) = A;
  c_A = device.malloc(r*c*sizeof(int),f_A);
  free(f_A);
}


dfloat WaveInitOCCA3d(Mesh *mesh, int KblkVin, int KblkSin,
		      int KblkUin, int KblkQin, int KblkQfin){

  KblkV = KblkVin;
  KblkS = KblkSin;
  KblkU = KblkUin;
  KblkQ = KblkQin;
  KblkQf = KblkQfin;

  occa::printAvailableDevices();

  //device.setup("mode = OpenCL, platformID = 0, deviceID = 0");
  //device.setup("mode = OpenMP, platformID = 0, deviceID = 0");
  device.setup("mode = Serial");
  //device.setup("mode = CUDA, platformID = 0, deviceID = 2");

  //device.setCompiler("nvcc"); device.setCompilerFlags("--use_fast_math"); device.setCompilerFlags("--fmad=true");

  printf("KblkV = %d, KblkS = %d, KblkU = %d, KblkQ (skew only) = %d\n",
	 KblkV,KblkS,KblkU,KblkQ);

  int K = mesh->K;
  int NpK = K*p_Np;

  int sz;

  // nodal operators
  sz = p_Np*p_Np*sizeof(dfloat);
  dfloat *f_Dr = (dfloat*)malloc(p_Np*p_Np*sizeof(dfloat));
  dfloat *f_Ds = (dfloat*)malloc(p_Np*p_Np*sizeof(dfloat));
  dfloat *f_Dt = (dfloat*)malloc(p_Np*p_Np*sizeof(dfloat));
  dfloat *f_LIFT = (dfloat*)malloc(p_Np*p_Nfp*p_Nfaces*sizeof(dfloat));

  // alternative way to cast matrices to dfloat
  Map<MatrixXdf >(f_Dr,p_Np,p_Np) = mesh->Dr.cast<dfloat>();
  Map<MatrixXdf >(f_Ds,p_Np,p_Np) = mesh->Ds.cast<dfloat>();
  Map<MatrixXdf >(f_Dt,p_Np,p_Np) = mesh->Dt.cast<dfloat>();
  Map<MatrixXdf >(f_LIFT,p_Np,p_Nfp*p_Nfaces) = mesh->LIFT.cast<dfloat>();
  c_Dr = device.malloc(sz,f_Dr);
  c_Ds = device.malloc(sz,f_Ds);
  c_Dt = device.malloc(sz,f_Dt);
  c_LIFT = device.malloc(p_Np*p_Nfp*p_Nfaces*sizeof(dfloat),f_LIFT);

  c_Fmask = device.malloc(p_Nfp*p_Nfaces*sizeof(int),mesh->FmaskC[0]);

  // ==================== BERNSTEIN STUFF ==========================

  // bernstein Dmatrices (4 entries per row)
  sz = 4*p_Np*sizeof(int);
  c_Dvals4 = device.malloc(4*p_Np*sizeof(dfloat), mesh->D_vals[0]);

  // barycentric deriv indices organized for ILP
  int *D_ids1 = (int*) malloc(p_Np*4*sizeof(int));
  int *D_ids2 = (int*) malloc(p_Np*4*sizeof(int));
  int *D_ids3 = (int*) malloc(p_Np*4*sizeof(int));
  int *D_ids4 = (int*) malloc(p_Np*4*sizeof(int));
  for(int i = 0; i < p_Np; ++i){
    D_ids1[4*i + 0] = mesh->D1_ids[i][0];
    D_ids1[4*i + 1] = mesh->D2_ids[i][0];
    D_ids1[4*i + 2] = mesh->D3_ids[i][0];
    D_ids1[4*i + 3] = mesh->D4_ids[i][0];

    D_ids2[4*i + 0] = mesh->D1_ids[i][1];
    D_ids2[4*i + 1] = mesh->D2_ids[i][1];
    D_ids2[4*i + 2] = mesh->D3_ids[i][1];
    D_ids2[4*i + 3] = mesh->D4_ids[i][1];

    D_ids3[4*i + 0] = mesh->D1_ids[i][2];
    D_ids3[4*i + 1] = mesh->D2_ids[i][2];
    D_ids3[4*i + 2] = mesh->D3_ids[i][2];
    D_ids3[4*i + 3] = mesh->D4_ids[i][2];

    D_ids4[4*i + 0] = mesh->D1_ids[i][3];
    D_ids4[4*i + 1] = mesh->D2_ids[i][3];
    D_ids4[4*i + 2] = mesh->D3_ids[i][3];
    D_ids4[4*i + 3] = mesh->D4_ids[i][3];
  }
  c_D_ids1 = device.malloc(sz,D_ids1);
  c_D_ids2 = device.malloc(sz,D_ids2);
  c_D_ids3 = device.malloc(sz,D_ids3);
  c_D_ids4 = device.malloc(sz,D_ids4);


  dfloat *h_EEL_vals = (dfloat*) malloc(p_Np*mesh->EEL_nnz*sizeof(dfloat));
  int *h_EEL_ids = (int*) malloc(p_Np*mesh->EEL_nnz*sizeof(int));
  for (int i = 0; i < p_Np; ++i){
    for (int j = 0; j < mesh->EEL_nnz; ++j){
      h_EEL_vals[i + j*p_Np] = mesh->EEL_vals[i][j];
      h_EEL_ids[i + j*p_Np] = mesh->EEL_ids[i][j];
    }
  }

  int L0_nnz = min(p_Nfp,7);
  mesh->L0_nnz = L0_nnz;
  dfloat *h_L0_vals = (dfloat*) malloc(p_Nfp*L0_nnz*sizeof(dfloat));
  int *h_L0_ids = (int*) malloc(p_Nfp*L0_nnz*sizeof(int));
  for (int i = 0; i < p_Nfp; ++i){
    for (int j = 0; j < L0_nnz; ++j){
      h_L0_vals[i + j*p_Nfp] = mesh->L0_vals[i][j];
      h_L0_ids[i + j*p_Nfp] = mesh->L0_ids[i][j];
    }
  }

  c_L0_vals = device.malloc(p_Nfp*L0_nnz*sizeof(dfloat),h_L0_vals);
  c_L0_ids = device.malloc(p_Nfp*L0_nnz*sizeof(int),h_L0_ids);

#if (USE_SLICE_LIFT)  // should use for N > 5 (faster)
  // store reduction matrices for orders 0,...,N
  setOccaIntArray(mesh,mesh->EEL_id_vec,c_EEL_ids);
  setOccaArray(mesh,mesh->EEL_val_vec,c_EEL_vals);
#else
  c_EEL_vals = device.malloc(p_Np*mesh->EEL_nnz*sizeof(dfloat),h_EEL_vals);
  c_EEL_ids = device.malloc(p_Np*mesh->EEL_nnz*sizeof(int),h_EEL_ids);
#endif

  //cout << "cEL = " << endl << mesh->cEL << endl;
  setOccaArray(mesh,mesh->cEL,c_cEL); // for slice-by-slice kernel

  int *h_slice_ids = (int*) malloc(p_Np*4*sizeof(int));
  for (int i = 0; i < p_Np; ++i){
    h_slice_ids[4*i+0] = mesh->slice_ids(i,0);
    h_slice_ids[4*i+1] = mesh->slice_ids(i,1);
    h_slice_ids[4*i+2] = mesh->slice_ids(i,2);
    h_slice_ids[4*i+3] = mesh->slice_ids(i,3);
  }
  //  setOccaIntArray(mesh,vol_ids,c_vol_ids); // arranged for int4 storage
  c_slice_ids = device.malloc(p_Np*4*sizeof(int),h_slice_ids);



  // =====================  geofacs ==================================

  double drdx, dsdx, dtdx;
  double drdy, dsdy, dtdy;
  double drdz, dsdz, dtdz, J;

  int sk = 0;
  int skP = -1;
  double *nxk = (double*) calloc(mesh->Nfaces,sizeof(double));//BuildVector(mesh->Nfaces);
  double *nyk = (double*) calloc(mesh->Nfaces,sizeof(double));//BuildVector(mesh->Nfaces);
  double *nzk = (double*) calloc(mesh->Nfaces,sizeof(double));//BuildVector(mesh->Nfaces);
  double *sJk = (double*) calloc(mesh->Nfaces,sizeof(double));//BuildVector(mesh->Nfaces);

  // [JC] packed geo + surface
  nvgeo = 9; // rst/xyz
  nfgeo = 4; // Fscale, (3)nxyz,
  ngeo = nfgeo*p_Nfaces + nvgeo; // nxyz + tau + Fscale (faces), G'*G (volume)
  dfloat *geo = (dfloat*) malloc(mesh->K*ngeo*sizeof(dfloat));
  dfloat *vgeo = (dfloat*) malloc(K*nvgeo*sizeof(dfloat));
  dfloat *fgeo = (dfloat*) malloc(K*nfgeo*p_Nfaces*sizeof(dfloat));

  dfloat FscaleMax = 0.f;
  for(int k=0;k<mesh->K;++k){

    GeometricFactors3d(mesh, k,
		       &drdx, &dsdx, &dtdx,
		       &drdy, &dsdy, &dtdy,
		       &drdz, &dsdz, &dtdz, &J);

    Normals3d(mesh, k, nxk, nyk, nzk, sJk);

    vgeo[k*nvgeo + 0] = drdx;    vgeo[k*nvgeo + 1] = drdy;    vgeo[k*nvgeo + 2] = drdz;
    vgeo[k*nvgeo + 3] = dsdx;    vgeo[k*nvgeo + 4] = dsdy;    vgeo[k*nvgeo + 5] = dsdz;
    vgeo[k*nvgeo + 6] = dtdx;    vgeo[k*nvgeo + 7] = dtdy;    vgeo[k*nvgeo + 8] = dtdz;

#if 0
    if (k==0){
      printf("rxyz = %f, %f, %f\n",drdx,drdy,drdz);
      printf("sxyz = %f, %f, %f\n",dsdx,dsdy,dsdz);
      printf("txyz = %f, %f, %f\n",dtdx,dtdy,dtdz);
      printf("J = %f\n\n",J);
    }
#endif

    for(int f=0;f<mesh->Nfaces;++f){

      dfloat Fscale = sJk[f]/J; //sJk[f]/(2.*J);
      dfloat nx = nxk[f];
      dfloat ny = nyk[f];
      dfloat nz = nzk[f];

      // for dt
      FscaleMax = max(FscaleMax,Fscale);

      fgeo[k*nfgeo*p_Nfaces + f*nfgeo + 0] = Fscale; // Fscale
      fgeo[k*nfgeo*p_Nfaces + f*nfgeo + 1] = nx;
      fgeo[k*nfgeo*p_Nfaces + f*nfgeo + 2] = ny;
      fgeo[k*nfgeo*p_Nfaces + f*nfgeo + 3] = nz;
    }
  }
  int *h_vmapP = (int*) malloc(mesh->K*p_Nfp*p_Nfaces*sizeof(int));
  for (int e = 0; e < mesh->K; ++e){
    for (int i = 0; i < p_Nfp*p_Nfaces; ++i){
      int f = i/p_Nfp;
      int idP = mesh->vmapP[i + p_Nfp*p_Nfaces*e];

      // correct vmapP for Nfields > 1
      int eNbr = mesh->EToE(e,f);
      idP -= p_Np*eNbr; // decrement
      idP += p_Np*p_Nfields*eNbr; // re-increment

      h_vmapP[i+p_Nfp*p_Nfaces*e] = idP;
    }
  }

  // storage for solution variables
  dfloat *f_Q    = (dfloat*) calloc(mesh->K*p_Nfields*p_Np, sizeof(dfloat));
  dfloat *f_resQ = (dfloat*) calloc(mesh->K*p_Nfields*p_Np, sizeof(dfloat));
  dfloat *f_rhsQ = (dfloat*) calloc(mesh->K*p_Nfields*p_Np, sizeof(dfloat));
  c_Q    = device.malloc(sizeof(dfloat)*mesh->K*p_Np*p_Nfields, f_Q);
  c_resQ = device.malloc(sizeof(dfloat)*mesh->K*p_Np*p_Nfields, f_resQ);
  c_rhsQ = device.malloc(sizeof(dfloat)*mesh->K*p_Np*p_Nfields, f_rhsQ);

  c_vgeo = device.malloc(mesh->K*nvgeo*sizeof(dfloat), vgeo);
  c_fgeo = device.malloc(mesh->K*nfgeo*p_Nfaces*sizeof(dfloat), fgeo);
  c_vmapP  = device.malloc(p_Nfp*p_Nfaces*mesh->K*sizeof(int),h_vmapP);

  // build kernels
  if (sizeof(dfloat)==8){
    dgInfo.addDefine("USE_DOUBLE", 1);
  }else{
    dgInfo.addDefine("USE_DOUBLE", 0);
  }

  dgInfo.addDefine("p_EEL_size",mesh->EEL_val_vec.rows());
  dgInfo.addDefine("p_EEL_nnz",mesh->EEL_nnz);
  dgInfo.addDefine("p_L0_nnz",min(p_Nfp,7)); // max 7 nnz with L0 matrix

  printf("p_Nfields = %d\n",p_Nfields);
  dgInfo.addDefine("p_Nfields",      p_Nfields); // wave equation

  dgInfo.addDefine("p_N",      p_N);
  dgInfo.addDefine("p_KblkV",  KblkV);
  dgInfo.addDefine("p_KblkS",  KblkS);
  dgInfo.addDefine("p_KblkU",  KblkU);
  dgInfo.addDefine("p_KblkQ",  KblkQ);
  dgInfo.addDefine("p_KblkQf",  KblkQf);

  dgInfo.addDefine("p_Np",      p_Np);
  dgInfo.addDefine("p_Nfp",     p_Nfp);
  dgInfo.addDefine("p_Nfaces",  p_Nfaces);
  dgInfo.addDefine("p_NfpNfaces",     p_Nfp*p_Nfaces);

  // [JC] max threads
  dgInfo.addDefine("p_ceilNq",min(512,mesh->Nq));
  int T = max(p_Np,p_Nfp*p_Nfaces);
  dgInfo.addDefine("p_T",T);
  int Tq = max(p_Np,mesh->Nfq*p_Nfaces);
  dgInfo.addDefine("p_Tq",Tq);

  dgInfo.addDefine("p_Nvgeo",nvgeo);
  dgInfo.addDefine("p_Nfgeo",nfgeo);


  std::string src = "okl/WaveKernels.okl";
  std::cout << "using src = " << src.c_str() << std::endl;

#if USE_BERN
  printf("Building Bernstein kernels from %s\n",src.c_str());
  // bernstein kernels
  rk_volume_bern  = device.buildKernelFromSource(src.c_str(), "rk_volume_bern", dgInfo);

  printf("building rk_surface_bern from %s\n",src.c_str());

#if USE_SLICE_LIFT
  rk_surface_bern = device.buildKernelFromSource(src.c_str(), "rk_surface_bern_slice", dgInfo);
  //rk_surface_bern = device.buildKernelFromSource(src.c_str(), "rk_surface_bern_slice_loads", dgInfo);
  //rk_surface_bern = device.buildKernelFromSource(src.c_str(), "rk_surface_bern_slice_square", dgInfo);
  printf("using slice-by-slice bern surface kernel; more efficient for N > 6\n");
#else
  printf("using non-optimal bern surface kernel; more efficient for N < 6\n");
  rk_surface_bern = device.buildKernelFromSource(src.c_str(), "rk_surface_bern", dgInfo);
#endif
#endif
  // nodal kernels
  rk_volume  = device.buildKernelFromSource(src.c_str(), "rk_volume", dgInfo);
  rk_surface = device.buildKernelFromSource(src.c_str(), "rk_surface", dgInfo);
  rk_update  = device.buildKernelFromSource(src.c_str(), "rk_update", dgInfo);

  // estimate dt. may wish to replace with trace inequality constant
  dfloat CN = (p_N+1)*(p_N+3)/3.0;
  dfloat dt = .25/(CN*FscaleMax);

  return (dfloat) dt;
}

// compute quadrature nodes + error
void compute_error(Mesh *mesh, double time, dfloat *Q,
		   double(*uexptr)(double,double,double,double),
		   double &L2err, double &relL2err){

  // compute error
  L2err = 0.0;
  double L2norm = 0.0;
  int kChunk = max(mesh->K/10,1);
  for(int k=0;k<mesh->K;++k){

    //double J = mesh->J(0,k); // assuming J = constant (planar tet)
    for(int i = 0; i < mesh->Nq; ++i){

      double J;
      if (mesh->Jq.rows()>0){ // if wadg quadrature initialized
	J= mesh->Jq(i,k);
      }else{
	J = mesh->J(0,k);
      }

      // interp to cubature nodes
      double x = 0.0; double y = 0.0; double z = 0.0; double uq = 0.0;
      for(int j=0;j<p_Np;++j){

	double Vq = mesh->Vq(i,j);
        x += Vq*mesh->x(j,k);
        y += Vq*mesh->y(j,k);
        z += Vq*mesh->z(j,k);
        uq += Vq*Q[j+p_Np*p_Nfields*k]; // get field value for pressure
      }

      double uex = (*uexptr)(x,y,z,(double)time);
      double err = uq-uex;

      L2err += err*err*mesh->wq(i)*J;
      L2norm += uex*uex*mesh->wq(i)*J;
    }
  }
  L2err = sqrt(L2err);
  relL2err = sqrt(L2err)/sqrt(L2norm);

  return;

}


// times planar kernels
void time_kernels(Mesh *mesh){

  double gflops = 0.0;
  double bw = 0.0;

  FILE *timingFile = fopen ("blockTimings.txt","a");
  occa::initTimer(device);

  // nodal kernels
  for (int step = 0; step < 10; ++step){
    dfloat fdt = 1.f, rka = 1.f, rkb = 1.f;

    occa::tic("volume (nodal)");
    rk_volume(mesh->K, c_vgeo, c_Dr, c_Ds, c_Dt, c_Q, c_rhsQ);
    device.finish();
    dfloat elapsedV = occa::toc("volume (nodal)",rk_volume, gflops, bw * sizeof(dfloat));

    occa::tic("surface (nodal))");
    rk_surface(mesh->K, c_fgeo, c_Fmask, c_vmapP, c_LIFT, c_Q, c_rhsQ);
    device.finish();
    dfloat elapsedS = occa::toc("surface (nodal)",rk_surface, gflops, bw * sizeof(dfloat));

    occa::tic("update (nodal)");
    rk_update(mesh->K, rka, rkb, fdt, c_rhsQ, c_resQ, c_Q);
    device.finish();
    dfloat elapsedU = occa::toc("update (nodal)",rk_update, gflops, bw * sizeof(dfloat));

    timeV+=elapsedV;
    timeS+=elapsedS;
    timeU+=elapsedU;
  }
  occa::printTimer();
  printf("Nodal kernels: elapsed time per timestep: V = %g, S = %g, U = %g\n", timeV/10,timeS/10,timeU/10);
  double denom = (double) (10 * mesh->K * p_Np * p_Nfields); // 10 steps
  timeV /= denom;
  timeS /= denom;
  timeU /= denom;
  printf("Nodal kernels: elapsed time per dof per timestep: V = %g, S = %g, U = %g, Total = %g\n", timeV,timeS,timeU,timeV+timeS+timeU);
  fprintf(timingFile,"%%Nodal kernels for N = %d\n KblkV(%d,%d) = %4.4g; KblkS(%d,%d) = %4.4g; KblkU(%d,%d) = %4.4g;\n",
          p_N,p_N,KblkV,timeV,p_N,KblkS,timeS,p_N,KblkU,timeU);

  // now do Bernstein kernels
  timeV = 0.0;
  timeS = 0.0;
  for (int step = 0; step < 10; ++step){
    occa::tic("volume (bern)");
    rk_volume_bern(mesh->K, c_vgeo,
                   c_D_ids1, c_D_ids2, c_D_ids3, c_D_ids4, c_Dvals4,
                   c_Q, c_rhsQ);

    device.finish();
    dfloat elapsedV = occa::toc("volume (bern)",rk_volume_bern, gflops, bw * sizeof(dfloat));

    occa::tic("surface (bern)");
    rk_surface_bern(mesh->K, c_fgeo, c_Fmask, c_vmapP,
                    c_slice_ids,c_EEL_ids, c_EEL_vals, c_L0_ids, c_L0_vals, c_cEL,
                    c_Q, c_rhsQ);
    device.finish();
    dfloat elapsedS = occa::toc("surface (bern)",rk_surface_bern, gflops, bw * sizeof(dfloat));

    timeV+=elapsedV;
    timeS+=elapsedS;

  }
  occa::printTimer();

  printf("Bern kernels: elapsed time per timestep: V = %g, S = %g\n", timeV/10,timeS/10);
  timeV /= denom;
  timeS /= denom;

  printf("Bern kernels: elapsed time per dof per timestep: V = %g, S = %g, U = %g, Total = %g\n", timeV,timeS,timeU,timeV+timeS+timeU);
#if USE_SLICE_LIFT
  fprintf(timingFile,"%%Bern kernels for N = %d\n KblkVB(%d,%d) = %4.4g; KblkSB_slice(%d,%d) = %4.4g;\n", p_N,p_N,KblkV,timeV,p_N,KblkS,timeS);
#else
  fprintf(timingFile,"%%Bern kernels for N = %d\n KblkVB(%d,%d) = %4.4g; KblkSB(%d,%d) = %4.4g;\n", p_N,p_N,KblkV,timeV,p_N,KblkS,timeS);
#endif
  fclose(timingFile);

}

// time WADG kernels
void time_curved_kernels(Mesh *mesh,int nsteps){
  int K = mesh->K;
  double gflops = 0.0;
  double bw = 0.0;

  dfloat elapsedQ = 0.f;
  dfloat elapsedQf = 0.f;

  occa::initTimer(device);
  for (int step = 0; step < nsteps; ++step){
    dfloat fdt = 1.f;
    dfloat rka = 1.f;
    dfloat rkb = 1.f;

#if USE_SKEW

    occa::tic("volume (nodal WADG skew combined)");
    rk_volume_WADG_skew_combine(mesh->K, c_KlistAll,
                                c_vgeoq, c_Jq,
                                c_Vq, c_Vrq, c_Vsq, c_Vtq,
                                c_Pq, c_Prq, c_Psq, c_Ptq,
                                c_Vskew,c_Pskew,
                                c_Q, c_rhsQ);

    device.finish();
    dfloat elapsedV = occa::toc("volume (nodal WADG skew)",rk_volume_WADG_skew, gflops, bw * sizeof(dfloat));

    // dfloat elapsedV = occa::toc("volume (nodal WADG skew combined)",rk_volume_WADG_skew_combine, gflops, bw * sizeof(dfloat));


    occa::tic("write_surf_quad (nodal WADG skew)");
    kernel_write_surf_quad(mesh->K, c_KlistAll,
                           c_Fmask, c_VfqFace, c_Vfq,
                           c_Q, c_Qf);
    device.finish();
    elapsedQf = occa::toc("write_surf_quad (nodal WADG skew)", kernel_write_surf_quad, gflops, bw * sizeof(dfloat));


    occa::tic("surface (nodal WADG skew)");
    rk_surface_WADG_skew(mesh->K, c_KlistAll,
                         c_fgeoq, c_Pfq, c_mapPq,
                         c_Qf, c_rhsQ);
    dfloat elapsedS = occa::toc("surface (nodal WADG skew)",rk_surface_WADG_skew, gflops, bw * sizeof(dfloat));


#else // if use strong form

    occa::tic("volume (nodal WADG)");
    rk_volume_WADG(mesh->K, c_KlistAll,
                   c_vgeoq, c_Jq, c_Vrq, c_Vsq, c_Vtq, c_Pq,
                   c_Q, c_rhsQ);
    device.finish();
    dfloat elapsedV = occa::toc("volume (nodal)",rk_volume_WADG, gflops, bw * sizeof(dfloat));

    dfloat elapsedS = 0.f;
    occa::tic("surface (nodal WADG)");
    rk_surface_WADG(mesh->K, c_KlistAll,
                    c_fgeo,c_fgeoq,c_VfqFace,c_Pfq,c_Fmask,c_vmapP,
                    c_Q, c_rhsQ);
    device.finish();
    elapsedS = occa::toc("surface (nodal)",rk_surface_WADG, gflops, bw * sizeof(dfloat));

#endif

    occa::tic("update (nodal WADG)");
    rk_update_WADG(mesh->K,
                   c_Vq_reduced, c_Pq_reduced, c_Jq_reduced,
                   rka, rkb, fdt,
                   c_rhsQ, c_resQ, c_Q);
    device.finish();
    dfloat elapsedU = occa::toc("update",rk_update_WADG, gflops,bw*sizeof(dfloat));

    timeV+=elapsedV;
    timeS+=elapsedS;
    timeU+=elapsedU;
    timeQ+=elapsedQ;
    timeQf+=elapsedQf;
  }

  occa::printTimer();
  double denom = (double) (nsteps * mesh->K * p_Np * p_Nfields); // 10 steps
  timeV /= denom;
  timeS /= denom;
  timeU /= denom;
  timeQ /= denom;
  timeQf /= denom;

#if USE_SKEW
  FILE *timingFile = fopen ("blockTimingsSkew.txt","a");
  printf("Strong-weak form: elapsed time per dof per timestep: V = %g, S = %g, U = %g, Q = %g, Qf = %g, Total = %g\n",
         timeV,timeS,timeU,timeQ,timeQf,timeV+timeS+timeU+timeQ+timeQf);
  //  fprintf(timingFile,"Strong-weak form: elapsed time with N = %d, Kblk = %2d: V = %4.4g, S = %4.4g, U = %4.4g, Q = %4.4g, Qf = %4.4g\n",
  //          p_N,KblkV,timeV,timeS,timeU,timeQ, timeQf);
  fprintf(timingFile,"%%Strong-weak kernels for N = %d\n",p_N);
  fprintf(timingFile,"KblkV(%d,%d) = %4.4g; KblkS(%d,%d) = %4.4g; ",p_N,KblkV,timeV,p_N,KblkS,timeS);
  fprintf(timingFile,"KblkU(%d,%d) = %4.4g; KblkQ(%d,%d) = %4.4g; KblkQf(%d,%d) = %4.4g;\n",p_N,KblkV,timeU,p_N,KblkV,timeQ,p_N,KblkV,timeQf);
#else
  FILE *timingFile = fopen ("blockTimingsStrong.txt","a");
  printf("Strong form: elapsed time per dof per timestep: V = %g, S = %g, U = %g, Q = %g, Qf = %g, Total = %g\n",
         timeV,timeS,timeU,timeQ,timeQf,timeV+timeS+timeU+timeQ+timeQf);
  //  fprintf(timingFile,"Strong form: elapsed time with N = %d, Kblk = %2d: V = %4.4g, S = %4.4g, U = %4.4g, Q = %4.4g, Qf = %4.4g\n",
  //          p_N,KblkV,timeV,timeS,timeU,timeQ,timeQf);
  fprintf(timingFile,"%%Strong kernels for N = %d\n",p_N);
  fprintf(timingFile,"KblkV(%d,%d) = %4.4g; KblkS(%d,%d) = %4.4g; ",p_N,KblkV,timeV,p_N,KblkS,timeS);
  fprintf(timingFile,"KblkU(%d,%d) = %4.4g; KblkQ(%d,%d) = %4.4g; KblkQf(%d,%d) = %4.4g;\n",p_N,KblkV,timeU,p_N,KblkV,timeQ,p_N,KblkV,timeQf);
#endif

  fclose(timingFile);

}

// run RK
void Wave_RK_sample_error(Mesh *mesh, dfloat FinalTime, dfloat dt,
                          double(*uexptr)(double,double,double,double)){

  double time = 0;
  int    INTRK, tstep=0;

  int totalSteps = (int)floor(FinalTime/dt);
  int tchunk = max(totalSteps/10,1);

  int tsample = 2*(p_N+1)*(p_N+1); // sample at every (*) timesteps
  int num_samples = totalSteps/tsample + 1;

  double *L2err = (double*) calloc(num_samples,sizeof(double));
  double *tvec = (double*) calloc(num_samples,sizeof(double));
  int tstep_sample = 0;

  // for sampling L2 error in time
  FILE *L2errFile = fopen ("longTimeL2err.txt","w");
  dfloat *Q = (dfloat*) calloc(p_Nfields*mesh->K*p_Np, sizeof(dfloat));   // 4 fields

  /* outer time step loop  */
  while (time<FinalTime){
#if 1
    if (tstep%tsample==0){
      ++tstep_sample;
      double L2err, relL2err;
      WaveGetData3d(mesh, Q);
      compute_error(mesh, time, Q, uexptr, L2err, relL2err);
      fprintf(L2errFile,"t(%d) = %g; L2e(%d) = %g;\n",tstep_sample,time,tstep_sample,L2err);
    }
#endif

    if (tstep%tchunk==0){
      printf("on timestep %d/%d\n",tstep, totalSteps);
    }

    /* adjust final step to end exactly at FinalTime */
    if (time+dt > FinalTime) { dt = FinalTime-time; }

    for (INTRK=1; INTRK<=5; ++INTRK) {

      // compute DG rhs
      const dfloat fdt = dt;
      const dfloat fa = (float)mesh->rk4a[INTRK-1];
      const dfloat fb = (float)mesh->rk4b[INTRK-1];

      if (tstep==0 && INTRK==1){
        printf("running regular kernel\n\n");
      }
      RK_step(mesh, fa, fb, fdt);

    }

    time += dt;     /* increment current time */
    tstep++;        /* increment timestep */

  }

  fclose(L2errFile);
}

// run RK
void Wave_RK(Mesh *mesh, dfloat FinalTime, dfloat dt, int useWADG){

  double time = 0;
  int    INTRK, tstep=0;

  int totalSteps = (int)floor(FinalTime/dt);
  int tchunk = max(totalSteps/10,1);

  int tsample = 2*(p_N+1)*(p_N+1); // sample at every (*) timesteps
  int num_samples = totalSteps/tsample + 1;

  double *L2err = (double*) calloc(num_samples,sizeof(double));
  double *tvec = (double*) calloc(num_samples,sizeof(double));
  int tstep_sample = 0;

  /* outer time step loop  */
  while (time<FinalTime){

    if (tstep%tchunk==0){
      printf("on timestep %d/%d\n",tstep, totalSteps);

#if 0
      dfloat *Q = (dfloat*) calloc(p_Nfields*mesh->K*p_Np, sizeof(dfloat));
      WaveGetData3d(mesh, Q);
      ostringstream visname;
      visname << "p" << tstep_sample << ".msh";
      ++tstep_sample;
      printf("writing gmsh file to %s\n",visname.str().c_str());
      writeVisToGMSH(visname.str(),mesh,Q,0,p_Nfields);
#endif

    }

    /* adjust final step to end exactly at FinalTime */
    if (time+dt > FinalTime) { dt = FinalTime-time; }

    for (INTRK=1; INTRK<=5; ++INTRK) {

      // compute DG rhs
      const dfloat fdt = dt;
      const dfloat fa = (float)mesh->rk4a[INTRK-1];
      const dfloat fb = (float)mesh->rk4b[INTRK-1];

      if (useWADG==0){
	if (tstep==0 && INTRK==1){
	  printf("running regular kernel\n\n");
	}
	RK_step(mesh, fa, fb, fdt);
      }else if (useWADG==1){
	if (tstep==0 && INTRK==1){
	  printf("running planar + subelem WADG kernel\n");
	}
	RK_step_WADG_subelem(mesh, fa, fb, fdt, time);
      }else if (useWADG==2){
	if (tstep==0 && INTRK==1){
	  printf("running curvilinear WADG kernels\n");
	}
	RK_step_WADG(mesh, fa, fb, fdt);
      }

    }

    time += dt;     /* increment current time */
    tstep++;        /* increment timestep */

  }

}

// defaults to nodal!!
void RK_step_WADG_subelem(Mesh *mesh, dfloat rka, dfloat rkb, dfloat fdt, dfloat time){

  dfloat f0 = 10.0;
  dfloat tR = 1.0 / f0;
  dfloat at = M_PI*f0*(time-tR);
  dfloat ftime = 1e4*(1.0 - 2.0*at*at)*exp(-at*at); // ricker pulse
  //  if (time > tR){
  //    printf("ftime = %f\n",ftime);
  //  }
  rk_volume_elas(mesh->K, c_vgeo, c_Dr, c_Ds, c_Dt, c_Q, c_rhsQ);
  rk_surface_elas(mesh->K, c_fgeo, c_Fmask, c_vmapP, c_LIFT, c_Q, c_rhsQ);
  rk_update_elas(mesh->K, c_Vq_reduced, c_Pq_reduced,
  		 c_rhoq, c_lambdaq, c_muq, c_c11, c_c12,
		 ftime, c_fsrc,
		 rka, rkb, fdt,
  		 c_rhsQ, c_resQ, c_Q);
  device.finish();

}

void RK_step_WADG(Mesh *mesh, dfloat rka, dfloat rkb, dfloat fdt){

  int K = mesh->K;

#if USE_SKEW

  rk_volume_WADG_skew_combine(mesh->KCurved, c_KlistCurved,
                              c_vgeoq, c_Jq,
                              c_Vq, c_Vrq, c_Vsq, c_Vtq,
                              c_Pq, c_Prq, c_Psq, c_Ptq,
                              c_Vskew,c_Pskew,
                              c_Q, c_rhsQ);

  rk_surface_WADG_skew_combine(mesh->KCurved, c_KlistCurved,
			       c_fgeo,c_fgeoq,c_VfqFace,c_Pfq,c_Fmask,c_vmapP,
			       c_Q, c_rhsQ);


#else // if use strong form with increased quadrature strength

  rk_volume_WADG(mesh->KCurved, c_KlistCurved,
      		 c_vgeoq, c_Jq, c_Vrq, c_Vsq, c_Vtq, c_Pq,
      		 c_Q, c_rhsQ);

  rk_surface_WADG(mesh->KCurved, c_KlistCurved,
		  c_fgeo,c_fgeoq,c_VfqFace,c_Pfq,c_Fmask,c_vmapP,
		  c_Q, c_rhsQ);


#endif

  // == run planar elements separately
  rk_volume_planar(mesh->KPlanar, c_KlistPlanar,
      		   c_vgeo, c_Jplanar, c_Dr, c_Ds, c_Dt, c_Q, c_rhsQ);
  rk_surface_planar(mesh->KPlanar, c_KlistPlanar,
  		    c_fgeo, c_Jplanar, c_Fmask, c_vmapP, c_LIFT, c_Q, c_rhsQ);


  rk_update_WADG(mesh->K,
		 c_Vq_reduced, c_Pq_reduced,
                 c_Jq_reduced,
                 rka, rkb, fdt,
		 c_rhsQ, c_resQ, c_Q);

}


void RK_step(Mesh *mesh, dfloat rka, dfloat rkb, dfloat fdt){

  double gflops = 0.0;
  double bw = 0.0;
  int K = mesh->K;

#if USE_BERN
  //printf("using bernstein kernels\n");
  rk_volume_bern(mesh->K, c_vgeo,
		 c_D_ids1, c_D_ids2, c_D_ids3, c_D_ids4, c_Dvals4,
		 c_Q, c_rhsQ);

  rk_surface_bern(mesh->K, c_fgeo, c_Fmask, c_vmapP,
		  c_slice_ids,
		  c_EEL_ids, c_EEL_vals,
   		  c_L0_ids, c_L0_vals,
		  c_cEL,
      		  c_Q, c_rhsQ);
#else
  //printf("Using nodal kernels\n");
  rk_volume(mesh->K, c_vgeo, c_Dr, c_Ds, c_Dt, c_Q, c_rhsQ);

  rk_surface(mesh->K, c_fgeo, c_Fmask, c_vmapP, c_LIFT, c_Q, c_rhsQ);

#endif

  rk_update(mesh->K, rka, rkb, fdt, c_rhsQ, c_resQ, c_Q);
  //int Ntotal = p_Nfields*p_Np*mesh->K;
  //rk_update(Ntotal, rka, rkb, fdt, c_rhsQ, c_resQ, c_Q);

#if 1
  dfloat *f_Q = (dfloat*) calloc(p_Nfields*mesh->K*p_Np, sizeof(dfloat));
  c_Q.copyTo(f_Q);
  for(int fld = 0; fld < p_Nfields; ++fld){
    printf("Field solution %d: \n",fld);
    for(int i = 0; i < p_Np; ++i){
      //      for(int e = 0; e < mesh->K; ++e){
      for(int e = 0; e < 1; ++e){
	printf("%f ",f_Q[i + fld*p_Np + e*p_Nfields*p_Np]);
      }
      printf("\n");
    }
    printf("\n\n");
  }
#endif

}




// set initial condition
void WaveSetU0(Mesh *mesh, dfloat *Q, dfloat time, int field,
	       double(*uexptr)(double,double,double,double)){

  // write out field = fields 2-4 = 0 (velocity)
  for(int k = 0; k < mesh->K; ++k){

    // store local interpolant
    dfloat *Qloc = BuildVector(p_Np);
    for(int i = 0; i < p_Np; ++i){
      double x = mesh->x(i,k);
      double y = mesh->y(i,k);
      double z = mesh->z(i,k);
      Qloc[i] = (*uexptr)(x,y,z,0.0);
    }

#if USE_BERN // convert nodal values to bernstein coefficients

    dfloat *Qtmp = BuildVector(p_Np);
    for (int i = 0; i < p_Np; ++i){
      Qtmp[i] = 0.f;
      for (int j = 0; j < p_Np; ++j){
	//Qtmp[i] += mesh->invVB[i][j]*Qloc[j];
	Qtmp[i] += mesh->invVB(i,j)*Qloc[j];
      }
    }
    for (int i = 0; i < p_Np; ++i){
      Qloc[i] = Qtmp[i];
    }

#endif

    for (int i = 0; i < p_Np; ++i){
      int id = k*p_Np*p_Nfields + i + p_Np*field;
      Q[id] = Qloc[i];
    }

  }

}


// set initial condition
void WaveProjectU0(Mesh *mesh, dfloat *Q, dfloat time,int field,
		   double(*uexptr)(double,double,double,double)){

  int Nq = mesh->Nq;
  VectorXd wq = mesh->wq;
  MatrixXd Vq = mesh->Vq;
  printf("Nq = %d, size of xq = %d\n",Nq,mesh->xq.rows());
  printf("size of Vq = %d, %d\n",Vq.rows(),Vq.cols());

  MatrixXd Qloc(p_Np,1);

  // write out field = fields 2-4 = 0 (velocity)
  for(int k = 0; k < mesh->K; ++k){

    // compute mass matrix explicitly w/quadrature
    VectorXd Jq = mesh->Jq.col(k);
    MatrixXd Mloc = Vq.transpose() * wq.asDiagonal() * Jq.asDiagonal() * Vq;

    // compute fxn at quad nodes
    MatrixXd uq(Nq,1);
    for (int i = 0; i < Nq; ++i){
      double xq = mesh->xq(i,k);
      double yq = mesh->yq(i,k);
      double zq = mesh->zq(i,k);
      double uqi = (*uexptr)(xq,yq,zq,0.0);
      //printf("uq[%d] = %g\n",i,uqi);
      uq(i,0) = uqi*wq(i)*Jq(i);
    }
    MatrixXd b = mesh->Vq.transpose() * uq;
    Qloc = mldivide(Mloc,b);

#if USE_BERN // convert nodal values to bernstein coefficients

    Qloc =  mesh->invVB*Qloc;

#endif

    // set pressure
    for (int i = 0; i < p_Np; ++i){
      int id = k*p_Np*p_Nfields + i + field*p_Np;
      Q[id] = Qloc(i,0);
    }

    //cout << "Q on elem " << k << " initialized to " << endl << Qloc << endl;

  }


  //cout << "Vq = " << endl << mesh->Vq << endl;
  //cout << "Jq = " << endl << mesh->Jq << endl;
  //cout << "wq = " << endl << mesh->wq << endl;
  //cout << "xq = " << endl << mesh->xq << endl;
  //cout << "yq = " << endl << mesh->yq << endl;
  //cout << "zq = " << endl << mesh->zq << endl;

}

void WaveSetData3d(dfloat *Q){
  c_Q.copyFrom(Q);
}


void WaveGetData3d(Mesh *mesh, dfloat *Q){//dfloat *d_p, dfloat *d_fp, dfloat *d_fdpdn){
  c_Q.copyTo(Q);
#if USE_BERN // convert back to nodal representation for L2 error, etc
  printf("Converting back to nodal from Bern representation\n");
  dfloat *Qtmp = BuildVector(p_Np);
  for (int fld = 0; fld < p_Nfields; fld++){
    for (int e = 0; e < mesh->K; ++e){
      for (int i = 0; i < p_Np; ++i){
	Qtmp[i] = 0.f;
	for (int j = 0; j < p_Np; ++j){
	  int id = j + fld*p_Np + e*p_Nfields*p_Np;
	  Qtmp[i] += mesh->VB(i,j)*Q[id];
	}
      }
      for (int i = 0; i < p_Np; ++i){
	int id = i + fld*p_Np + e*p_Nfields*p_Np;
	Q[id] = Qtmp[i];
      }
    }
  }
#endif

}

void writeVisToGMSH(string fileName, Mesh *mesh,dfloat *Q, int iField, int Nfields){

  int timeStep = 0;
  double time = 0.0;
  int K = mesh->K;
  int Dim = 3;
  int N = p_N;


  MatrixXi monom(p_Np, Dim);
  MatrixXd vdm(p_Np, p_Np);
  for(int i=0, n=0; i<=N; i++){
    for(int j=0; j<=N; j++){
      for(int k=0; k<=N; k++){
	if(i+j+k <= N){
	  monom(n,0) = i;
	  monom(n,1) = j;
	  monom(n,2) = k;
	  n++;
	}
      }
    }
  }
  for(int m=0; m<p_Np; m++){
    for(int n=0; n<p_Np; n++){
      double r = mesh->r(n);
      double s = mesh->s(n);
      double t = mesh->t(n);
      vdm(m,n) = pow((r+1)/2.,monom(m,0)) *
	pow((s+1)/2.,monom(m,1)) * pow((t+1)/2.,monom(m,2));
    }
  }
  MatrixXd coeff = vdm.inverse();

  /// write the gmsh file
  ofstream *posFile;
  posFile = new ofstream(fileName.c_str());
  *posFile << "$MeshFormat" << endl;
  *posFile << "2.2 0 8" << endl;
  *posFile << "$EndMeshFormat" << endl;

  /// write the interpolation scheme
  *posFile << "$InterpolationScheme" << endl;
  *posFile << "\"MyInterpScheme\"" << endl;
  *posFile << "1" << endl;
  *posFile << "5 2" << endl;  // 5 2 = tets
  *posFile << p_Np << " " << p_Np << endl;  // size of matrix 'coeff'
  for(int m=0; m<p_Np; m++){
    for(int n=0; n<p_Np; n++)
      *posFile << coeff(m,n) << " ";
    *posFile << endl;
  }
  *posFile << p_Np << " " << Dim << endl;  // size of matrix 'monom'
  for(int n=0; n<p_Np; n++){
    for(int d=0; d<Dim; d++)
      *posFile << monom(n,d) << " ";
    *posFile << endl;
  }
  *posFile << "$EndInterpolationScheme" << endl;

  /// write element node data
  *posFile << "$ElementNodeData" << endl;
  *posFile << "2" << endl;
  *posFile << "\"" << "Field " << iField << "\"" << endl;  /// name of the view
  *posFile << "\"MyInterpScheme\"" << endl;
  *posFile << "1" << endl;
  *posFile << time << endl;
  *posFile << "3" << endl;
  *posFile << timeStep << endl;
  *posFile << "1" << endl;  /// ("numComp")
  *posFile << K << endl;  /// total number of elementNodeData in this file
  for(int k=0; k<K; k++){
    *posFile << mesh->EToGmshE(k) << " " << p_Np;
    for(int i=0; i<p_Np; i++)
      *posFile << " " << Q[i + iField*p_Np + k*p_Np*Nfields];
    *posFile << endl;
  }
  *posFile << "$EndElementNodeData" << endl;

  posFile->close();
  delete posFile;

}
