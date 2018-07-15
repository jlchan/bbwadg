#include "dfloat.h"
#include "fem.h"
//#include "types.h"

// may want to generalize later to multiple devices...

// sets up mesh-based occa parameters
void setupOccaMesh2d(Mesh *mesh, App *app){

  app->props = occa::getKernelProperties(); //props;

  app->props["defines/p_Np"] = mesh->Np; // number of dims
  app->props["defines/p_Np1"] = mesh->N + 1;  // (N+1)
  app->props["defines/p_Np2"] = (mesh->N + 1)*(mesh->N + 1);  // (N+1)  

  app->props["defines/p_Nq1"] = mesh->N + 1;  // (N+1)
  app->props["defines/p_Nq2"] = (mesh->N + 1)*(mesh->N + 1);  
  //  app->props["defines/p_Nq3"] = (mesh->N + 1)*(mesh->N + 1)*(mesh->N + 1);  
  
  app->props["defines/p_Nfaces"] = mesh->Nfaces;
  app->props["defines/p_Nfp"] = mesh->Nfp;
  app->props["defines/p_NfpNfaces"] = mesh->Nfp * mesh->Nfaces;
  app->props["defines/p_T"] = max(mesh->Nfp * mesh->Nfaces,mesh->Np);

  int Nvgeo = 5; 
  int Nfgeo = 3; 
  app->props["defines/p_Nvgeo"] = Nvgeo;
  app->props["defines/p_Nfgeo"] = Nfgeo;

  // switch dfloat type (double/float) in types.h
  if (sizeof(dfloat)==4){
    app->props["defines/USE_DOUBLE"] = 0;
  }else{
    app->props["defines/USE_DOUBLE"] = 1;
  }

  // allocate space for geofaces
  int K = mesh->K;
  int Np = mesh->Np;
  int NfpNfaces = (mesh->Nfp)*(mesh->Nfaces);  

  // interpolate to quad pts and store
  MatrixXd Vq   = mesh->Vq;
  MatrixXd rxJq = Vq*(mesh->rxJ);
  MatrixXd sxJq = Vq*(mesh->sxJ);
  MatrixXd ryJq = Vq*(mesh->ryJ);
  MatrixXd syJq = Vq*(mesh->syJ);
  MatrixXd Jq   = Vq*(mesh->J); // interp: may lose some accuracy

  MatrixXd vgeo(Nvgeo*(mesh->Vq.rows()),K);  
  vgeo << rxJq,
    ryJq,
    sxJq,
    syJq,
    Jq;	  

  // interp to face values 
  MatrixXd Vf   = mesh->Vf;
  MatrixXd rxJf = Vf*(mesh->rxJ);
  MatrixXd sxJf = Vf*(mesh->sxJ);
  MatrixXd ryJf = Vf*(mesh->ryJ);
  MatrixXd syJf = Vf*(mesh->syJ);
  MatrixXd Jf   = Vf*(mesh->J); // this may not be accurate. I don't think I use though.
  MatrixXd vfgeo(Nvgeo*(Vf.rows()),K);
  vfgeo << rxJf,
    ryJf,
    sxJf,
    syJf,
    Jf;	  
  
  // normals computed at face quad points already
  MatrixXd fgeo(Nfgeo*(mesh->Vf.rows()),K);
  fgeo << mesh->nxJ,
    mesh->nyJ,
    mesh->sJ;

  // init rhs, res, Qf to zero
  int Nfields = mesh->Nfields;
  MatrixXd rhs(Nfields*Np,K), res(Nfields*Np,K);
  rhs.fill(0.0);
  res.fill(0.0);
  MatrixXd Qf(Nfields*NfpNfaces,K),rhsf(Nfields*NfpNfaces,K);
  Qf.fill(0.0);
  rhsf.fill(0.0);  

  setOccaArray(app, rhs, app->o_rhs);
  setOccaArray(app, res, app->o_res);  

  setOccaArray(app, Qf,  app->o_Qf);
  setOccaArray(app, rhsf,  app->o_rhsf);
  
  setOccaArray(app, vgeo, app->o_vgeo);
  setOccaArray(app, vfgeo, app->o_vfgeo);  
  setOccaArray(app, fgeo, app->o_fgeo);

  // correct for Nfields > 1
  MatrixXi mapPqNfields(NfpNfaces,K);
  for (int i = 0; i < NfpNfaces*K; ++i){
    int idP = mesh->mapPq(i);
    int e = idP / NfpNfaces;
    int fidP = idP % NfpNfaces;
    mapPqNfields(i) = fidP + e*NfpNfaces*Nfields;
  }
  //cout << "mapPq = " << mesh->mapPq << endl;
  //cout << "mapPqNfields = " << mapPqNfields << endl;
  setOccaIntArray(app, mapPqNfields, app->o_mapPq);

  setOccaArray(app, mesh->D1D, app->o_D1D); 
  setOccaArray(app, mesh->Vf1D.row(0), app->o_Vf1D); // Vf1D rows 1,2 = mirror images
  setOccaArray(app, mesh->Lf1D.col(0), app->o_Lf1D); // Lf1D cols 1,2 = mirror images  
  //cout << "Vf1D = " << mesh->Vf1D.row(0) << endl;
  //cout << "Lf1D = " << mesh->Lf1D.col(0) << endl;
}


// sets up mesh-based occa parameters
void setupOccaMesh3d(Mesh *mesh, App *app){

  app->props = occa::getKernelProperties();

  app->props["defines/p_Np"] = mesh->Np; // number of dims

  app->props["defines/p_Nq1"] = mesh->N + 1;  // (N+1)
  app->props["defines/p_Nq2"] = (mesh->N + 1)*(mesh->N + 1);  
  app->props["defines/p_Nq3"] = (mesh->N + 1)*(mesh->N + 1)*(mesh->N + 1);  
    
  app->props["defines/p_Nfaces"] = mesh->Nfaces;
  app->props["defines/p_Nfp"] = mesh->Nfp;
  app->props["defines/p_NfpNfaces"] = mesh->Nfp * mesh->Nfaces;
  app->props["defines/p_T"] = max(mesh->Nfp * mesh->Nfaces,mesh->Np);

  int Nvgeo = 10;  // 3x3 geofacs + J
  int Nfgeo = 4; 
  app->props["defines/p_Nvgeo"] = Nvgeo;
  app->props["defines/p_Nfgeo"] = Nfgeo;

  // switch dfloat type (double/float) in types.h
  if (sizeof(dfloat)==4){
    app->props["defines/USE_DOUBLE"] = 0;
  }else{
    app->props["defines/USE_DOUBLE"] = 1;
  }

  // allocate space for geofaces
  int K = mesh->K;
  int Np = mesh->Np;
  int NfpNfaces = (mesh->Nfp)*(mesh->Nfaces);  

  // interpolate to quad pts and store
  MatrixXd Vq   = mesh->Vq;
  MatrixXd rxJq = Vq*(mesh->rxJ);
  MatrixXd sxJq = Vq*(mesh->sxJ);
  MatrixXd txJq = Vq*(mesh->txJ);  
  MatrixXd ryJq = Vq*(mesh->ryJ);
  MatrixXd syJq = Vq*(mesh->syJ);
  MatrixXd tyJq = Vq*(mesh->tyJ);  
  MatrixXd rzJq = Vq*(mesh->rzJ);
  MatrixXd szJq = Vq*(mesh->szJ);
  MatrixXd tzJq = Vq*(mesh->tzJ);  
  MatrixXd Jq   = Vq*(mesh->J); // interp: may lose some accuracy
  //cout << "Jq = " << Jq << endl;

  MatrixXd vgeo(Nvgeo*Vq.rows(),K);  
  vgeo << rxJq,
    ryJq,
    rzJq,
    sxJq,
    syJq,
    szJq,
    txJq,
    tyJq,
    tzJq,
    Jq;

  //printf("Nvgeo = %d, rows of Vq = %d\n",Nvgeo,Vq.rows());
  //  cout << "Vgeo = " << vgeo << endl;

  // interp to face values 
  MatrixXd Vf   = mesh->Vf;
  MatrixXd rxJf = Vf*(mesh->rxJ);
  MatrixXd sxJf = Vf*(mesh->sxJ);
  MatrixXd txJf = Vf*(mesh->txJ);  
  MatrixXd ryJf = Vf*(mesh->ryJ);
  MatrixXd syJf = Vf*(mesh->syJ);
  MatrixXd tyJf = Vf*(mesh->tyJ);  
  MatrixXd rzJf = Vf*(mesh->rzJ);
  MatrixXd szJf = Vf*(mesh->szJ);
  MatrixXd tzJf = Vf*(mesh->tzJ);  
  MatrixXd Jf   = Vf*(mesh->J); // this may not be accurate. I don't think I use though.
  MatrixXd vfgeo(Nvgeo*(Vf.rows()),K);
  vfgeo << rxJf,
    ryJf,
    rzJf,
    sxJf,
    syJf,
    szJf,
    txJf,
    tyJf,
    tzJf,
    Jf;	  

  // normals computed at face quad points already
  MatrixXd fgeo(Nfgeo*(Vf.rows()),K);
  fgeo << mesh->nxJ,
    mesh->nyJ,
    mesh->nzJ,    
    mesh->sJ;

  //printf("fgeo size = %d\n",fgeo.rows());
  //cout << "fgeo = " << fgeo << endl;

  // init rhs, res, Qf to zero
  int Nfields = mesh->Nfields;
  MatrixXd rhs(Nfields*Np,K), res(Nfields*Np,K);
  rhs.fill(0.0);
  res.fill(0.0);
  MatrixXd Qf(Nfields*NfpNfaces,K),rhsf(Nfields*NfpNfaces,K);
  Qf.fill(0.0);
  rhsf.fill(0.0);

#if 0
  printf("moving data to GPU: rhs arry is size %d, %d\n",rhs.rows(),rhs.cols());
  printf("moving data to GPU: rhsf arry is size %d, %d\n",rhsf.rows(),rhsf.cols());
  printf("moving data to GPU: vgeo arry is size %d, %d\n",vgeo.rows(),vgeo.cols());
  printf("moving data to GPU: vfgeo arry is size %d, %d\n",vfgeo.rows(),vfgeo.cols());
  printf("moving data to GPU: fgeo arry is size %d, %d\n",fgeo.rows(),fgeo.cols());
#endif  
  setOccaArray(app, rhs, app->o_rhs);
  setOccaArray(app, res, app->o_res);  
  setOccaArray(app, Qf,  app->o_Qf);
  setOccaArray(app, rhsf,  app->o_rhsf);
  setOccaArray(app, vgeo, app->o_vgeo);
  setOccaArray(app, vfgeo, app->o_vfgeo);
  setOccaArray(app, fgeo, app->o_fgeo);

  // correct for Nfields > 1
  MatrixXi mapPqNfields(NfpNfaces,K);
  for (long long int i = 0; i < NfpNfaces*K; ++i){
    int idP = mesh->mapPq(i);
    int e = idP / NfpNfaces;
    int fidP = idP % NfpNfaces;
    mapPqNfields(i) = fidP + e*NfpNfaces*Nfields;
  }
  setOccaIntArray(app, mapPqNfields, app->o_mapPq);

  setOccaArray(app, mesh->D1D, app->o_D1D); 
  setOccaArray(app, mesh->Vf1D.row(0), app->o_Vf1D); // Vf1D rows 1,2 = mirror images
  setOccaArray(app, mesh->Lf1D.col(0), app->o_Lf1D); // Lf1D cols 1,2 = mirror images    

}


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
		 

