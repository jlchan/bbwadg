#include "fem.h"
#include "Basis.h"

// make quads 
void InitRefData2d(Mesh *mesh, int N){

  int Np1 = (N+1);

  // nodal points (GCL and all)
  VectorXd r1D,w1D;
  JacobiGL(N, 0, 0, r1D,w1D);

  MatrixXd rmat,smat;
  meshgrid(r1D,r1D, smat,rmat);
  VectorXd r = Map<VectorXd>(rmat.data(), rmat.cols()*rmat.rows());
  VectorXd s = Map<VectorXd>(smat.data(), smat.cols()*smat.rows());
  //  cout << "r = " << r << endl;
  //  cout << "s = " << s << endl;  

  // quad points 
  VectorXd rq1D, wq1D;
  JacobiGL(N, 0, 0, rq1D, wq1D);
  //JacobiGQ(N, 0, 0, rq1D, wq1D);  
  
  //cout << "rq,wq for GQ = " << rq1D << ", " << wq1D << endl;
  
  meshgrid(rq1D,rq1D, rmat,smat);
  VectorXd rq = Map<VectorXd>(rmat.data(), rmat.cols()*rmat.rows());
  VectorXd sq = Map<VectorXd>(smat.data(), smat.cols()*smat.rows()); 
  meshgrid(wq1D,wq1D,rmat,smat);
  VectorXd wr = Map<VectorXd>(rmat.data(), rmat.cols()*rmat.rows());
  VectorXd ws = Map<VectorXd>(smat.data(), smat.cols()*smat.rows());
  VectorXd wq = wr.array()*ws.array();

  // face quad points
  VectorXd e(rq1D.size()); e.fill(1.0);
  int Nfaces = 4;
  int NfqNfaces = Nfaces * rq1D.size();  
  VectorXd rf(NfqNfaces); rf << rq1D, e, rq1D, -e;
  VectorXd sf(NfqNfaces); sf << -e, rq1D, e, rq1D;
  VectorXd wf(NfqNfaces); wf << wq1D, wq1D, wq1D, wq1D;
  
  // nodal operators
  MatrixXd V1D = Vandermonde1D(N,r1D);
  MatrixXd Vr1D = GradVandermonde1D(N,r1D);
  MatrixXd D1D = mrdivide(Vr1D,V1D);
  MatrixXd I1D = MatrixXd::Identity(Np1,Np1);

  MatrixXd V = Vandermonde2DQuad(N,r,s);
  MatrixXd Dr = kron(I1D,D1D);
  MatrixXd Ds = kron(D1D,I1D);
  //cout << "Dr = " << endl << Dr << endl;
  //cout << "Ds = " << endl << Ds << endl;  

  // interp to quadrature
  MatrixXd Vqtmp = Vandermonde2DQuad(N,rq,sq);
  MatrixXd Vftmp = Vandermonde2DQuad(N,rf,sf);
  MatrixXd Vq = mrdivide(Vqtmp,V);
  MatrixXd Vf = mrdivide(Vftmp,V);  
  
  // 1D quadrature operators
  MatrixXd Vq1D = Vandermonde1D(N,rq1D);
  MatrixXd Vrq1D = GradVandermonde1D(N,rq1D);
  MatrixXd Dq1D = mrdivide(Vrq1D,Vq1D);
  VectorXd rf1D(2); rf1D << -1.0,1.0;
  MatrixXd Vf1Dtmp = Vandermonde1D(N,rf1D);
  MatrixXd Vf1D = mrdivide(Vf1Dtmp,Vq1D);

  //  cout << "Dq1D = " << endl << Dq1D << endl;
  //  cout << "Vf1D = " << endl << Vf1D << endl;  

  mesh->N = N;
  mesh->Np = Np1*Np1;
  mesh->Nfp = Np1;
  mesh->Nfaces = Nfaces;
  mesh->Nverts = 4;

  // GLL nodes (for geofacs)
  mesh->r = r;
  mesh->s = s;  
  mesh->Dr = Dr;
  mesh->Ds = Ds;
  mesh->V = V;

  // interp to vol/face nodes
  mesh->rq = rq;
  mesh->sq = sq;
  mesh->wq = wq;  
  mesh->Vq = Vq;

  mesh->rf = rf;
  mesh->sf = sf;
  mesh->wf = wf;    
  mesh->Vf = Vf;

  //  cout << "rq = " << rq << endl;
  //  cout << "sq = " << sq << endl;  
  //cout << "rf = " << rf << endl;
  //cout << "sf = " << sf << endl;  
 
  // GQ nodes
  MatrixXd Q1D = wq1D.asDiagonal()*Dq1D;
  VectorXd invwq1D = wq1D.array().inverse();
  MatrixXd D1D_skew = invwq1D.asDiagonal()*(Q1D - Q1D.transpose());
  mesh->D1D = D1D_skew; // (Q - Q^T)
  //mesh->D1D = Dq1D;
  mesh->Vf1D = Vf1D;
  VectorXd wqInv = wq1D.array().inverse();
  mesh->Lf1D = wqInv.asDiagonal()*Vf1D.transpose(); // 1D lift - may not need

  //  cout << "D1D = " << mesh->D1D << endl;
  //  cout << "Vf1D = " << mesh->Vf1D << endl;
  //cout << "Lf1D = " << mesh->Lf1D << endl;   
  mesh->wq1D = wq1D;

  // LSRK-45 coefficients
  mesh->rk4a.resize(5);
  mesh->rk4a(0) =              0.0;
  mesh->rk4a(1) =  -567301805773.0 / 1357537059087.0;
  mesh->rk4a(2) = -2404267990393.0 / 2016746695238.0;
  mesh->rk4a(3) = -3550918686646.0 / 2091501179385.0;
  mesh->rk4a(4) = -1275806237668.0 /  842570457699.0;

  mesh->rk4b.resize(5);
  mesh->rk4b(0) =  1432997174477.0 /  9575080441755.0;
  mesh->rk4b(1) =  5161836677717.0 / 13612068292357.0;
  mesh->rk4b(2) =  1720146321549.0 /  2090206949498.0;
  mesh->rk4b(3) =  3134564353537.0 /  4481467310338.0;
  mesh->rk4b(4) =  2277821191437.0 / 14882151754819.0;

  mesh->rk4c.resize(6);
  mesh->rk4c(0) =              0.0;
  mesh->rk4c(1) =  1432997174477.0 / 9575080441755.0;
  mesh->rk4c(2) =  2526269341429.0 / 6820363962896.0;
  mesh->rk4c(3) =  2006345519317.0 / 3224310063776.0;
  mesh->rk4c(4) =  2802321613138.0 / 2924317926251.0;
  mesh->rk4c(5) =              1.0;
  
}

void QuadMesh2d(Mesh *mesh, int Nx, int Ny){
  // just make a uniform quadrilateral mesh on [-1,1]
  int Nxp = Nx+1;
  int Nyp = Ny+1;
  int K = Nx*Ny;
  int Nv = Nxp*Nyp;

  VectorXd x1D = VectorXd::LinSpaced(Nxp,-1,1);
  VectorXd y1D = VectorXd::LinSpaced(Nyp,-1,1);

  MatrixXd X,Y;
  meshgrid(x1D,y1D,X,Y);
  VectorXd x = Map<VectorXd>(X.data(), X.cols()*X.rows());
  VectorXd y = Map<VectorXd>(Y.data(), Y.cols()*Y.rows());

  // low:step:hi
  int low, hi, step, size;
  low = 0;  hi = Nx; step = 1; size = Nxp;
  VectorXi ix1D = VectorXi::LinSpaced(((hi-low)/step)+1,low,low+step*(size-1));

  low = 0;  hi = Ny; step = 1; size = Nyp;
  VectorXi iy1D = VectorXi::LinSpaced(((hi-low)/step)+1,low,low+step*(size-1));

  MatrixXi I,J;
  meshgrid(ix1D,iy1D,I,J);
  MatrixXi inds = I*Ny + (J+I);
  //  cout << "inds = " << endl << inds << endl;

  MatrixXi EToV(K,4);
  int k = 0;
  for (int i = 0; i < Ny; ++i){
    for (int j = 0; j < Nx; ++j){
      EToV(k,0) = inds(i,j);
      EToV(k,1) = inds(i,j+1);
      EToV(k,2) = inds(i+1,j+1);
      EToV(k,3) = inds(i+1,j);
      ++k;
    }
  }
  //  cout << EToV << endl;
  mesh->EToV = EToV;
  mesh->K = K;  
  mesh->Nv = Nv;
  mesh->VX = x;
  mesh->VY = y;      
  
}


// computes physical nodal position based on affine 
void MapNodes2d(Mesh *mesh){
  int Nfaces = mesh->Nfaces;
  int K = mesh->K;
  int Nv = mesh->Nv;
  MatrixXi EToV = mesh->EToV;

  VectorXd VX = mesh->VX;
  VectorXd VY = mesh->VY;  
 
  MatrixXd x1(EToV.cols(),K);
  MatrixXd y1(EToV.cols(),K);  
  for (int e = 0; e < K; ++e){
    for (int v = 0; v < EToV.cols(); ++v){
      int vid = EToV(e,v);
      x1(v,e) = VX(vid);
      y1(v,e) = VY(vid);      
    }
  }
  // interp to high order GLL nodes (matlab ordering)
  VectorXd r1(4),s1(4);
  r1 << -1.0, 1.0, 1.0, -1.0;
  s1 << -1.0, -1.0, 1.0, 1.0;

  VectorXd r = mesh->r;
  VectorXd s = mesh->s;
  
  MatrixXd Vdm1 = Vandermonde2DQuad(1,r1,s1);
  MatrixXd VdmN = Vandermonde2DQuad(1,r,s);
  //cout << "Vdm1 = " << Vdm1 << endl;
  MatrixXd V1 = mrdivide(VdmN,Vdm1);

  MatrixXd x = V1*x1;
  MatrixXd y = V1*y1;
  mesh->x = x;
  mesh->y = y;

}


void GeometricFactors2d(Mesh *mesh){

  int Nfaces = mesh->Nfaces;
  int K = mesh->K;
  int Nv = mesh->Nv;

  MatrixXd x = mesh->x;
  MatrixXd y = mesh->y;

  // compute geofacs
  MatrixXd Dr = mesh->Dr;
  MatrixXd Ds = mesh->Ds;
  MatrixXd xr = Dr*x;
  MatrixXd xs = Ds*x;
  MatrixXd yr = Dr*y;
  MatrixXd ys = Ds*y;

  MatrixXd rxJ = ys;
  MatrixXd sxJ = -yr;
  MatrixXd ryJ = -xs;
  MatrixXd syJ = xr;

  // for higher accuracy, interp xs, yr, sr xs to GQ points first    
  xr = mesh->Vq*xr;
  xs = mesh->Vq*xs;
  yr = mesh->Vq*yr;
  ys = mesh->Vq*ys;  
  /* */
  MatrixXd J = -xs.array()*yr.array() + xr.array()*ys.array();

  // interpolate J back to store at GLL points (we undo this later)
  J = mldivide(mesh->Vq,J);
 
  // all quantities stored at GLL points (must interp before using)
  mesh->rxJ = rxJ;
  mesh->sxJ = sxJ;
  mesh->ryJ = ryJ;
  mesh->syJ = syJ;
  mesh->J = J;
  
}

void Normals2d(Mesh *mesh){

  MatrixXd Vf = mesh->Vf;
  MatrixXd rxJf = Vf*mesh->rxJ;
  MatrixXd sxJf = Vf*mesh->sxJ;
  MatrixXd ryJf = Vf*mesh->ryJ;
  MatrixXd syJf = Vf*mesh->syJ;

  int Nfaces = mesh->Nfaces;
  int Nfp = mesh->Nfp;
  VectorXd e(Nfp); e.fill(1.0);
  VectorXd z(Nfp); z.fill(0.0);  
  VectorXd nrJ(Vf.rows()); nrJ << z, e, z, -e;
  VectorXd nsJ(Vf.rows()); nsJ << -e, z, e, z;

  MatrixXd nxJ = nrJ.asDiagonal()*rxJf + nsJ.asDiagonal()*sxJf;
  MatrixXd nyJ = nrJ.asDiagonal()*ryJf + nsJ.asDiagonal()*syJf;
  MatrixXd sJ = (nxJ.array().square() + nyJ.array().square()).sqrt();

  mesh->nxJ = nxJ;
  mesh->nyJ = nyJ;
  mesh->sJ = sJ;  
}

void MakeNodeMapsPeriodic2d(Mesh *mesh, MatrixXd xf, MatrixXd yf, double DX, double DY, MatrixXi &mapP){

  printf("Modifying node maps for periodicity\n");
  
  int Nfp = mesh->Nfp;
  int Nfaces = mesh->Nfaces;
  int K = mesh->K;

  VectorXi isBoundaryFace(Nfaces*K);
  isBoundaryFace.fill(0);
  for (int e = 0; e < K; ++e){
    for (int f = 0; f < Nfaces; ++f){
      if (e==mesh->EToE(e,f)){
	isBoundaryFace(f + e*Nfaces) = 1;
      }
    }
  }
 
  xf.resize(Nfp,Nfaces*K);
  yf.resize(Nfp,Nfaces*K);  
  vector<pair<int,int> > xfaces,yfaces;  
  for (int f1 = 0; f1 < Nfaces*K; ++f1){
    for (int f2 = 0; f2 < Nfaces*K; ++f2){
      if (isBoundaryFace(f1) && isBoundaryFace(f2)){
	// distance b/w faces = diff b/w min/max coeffs 
	double dx = .5*(fabs(xf.col(f1).maxCoeff()-xf.col(f2).maxCoeff())
			+ fabs(xf.col(f1).minCoeff()-xf.col(f2).minCoeff()));
	double dy = .5*(fabs(yf.col(f1).maxCoeff()-yf.col(f2).maxCoeff())
			+ fabs(yf.col(f1).minCoeff()-yf.col(f2).minCoeff()));
	
	if ((dx < NODETOL) & (fabs(dy-DY)<NODETOL)){
	  xfaces.push_back(make_pair(f1,f2)); 
	}else if ((dy < NODETOL) & (fabs(dx-DX)<NODETOL)){
	  yfaces.push_back(make_pair(f1,f2)); 
	}
      }
    }
  }
  printf("num xfaces = %d, num yfaces = %d\n",xfaces.size(),yfaces.size());
  //cout << "xf = " << xf << endl;
  //cout << "yf = " << yf << endl;  

  // find node maps in x
  for (int f = 0; f < xfaces.size(); ++f){
    int f1 = xfaces[f].first;
    int f2 = xfaces[f].second;
    for (int i = 0; i < Nfp; ++i){
      int idM = i + f1*Nfp;
      double xM = xf(idM);

      for (int j = 0; j < Nfp; ++j){
	int idP = j + f2*Nfp;
	double xP = xf(idP);
	
	if (fabs(xM-xP) < NODETOL){
	  mapP(idM) = idP;
	  break;
	}
      }      
    }
  }

  // find node maps in y
  for (int f = 0; f < yfaces.size(); ++f){
    int f1 = yfaces[f].first;
    int f2 = yfaces[f].second;
    for (int i = 0; i < Nfp; ++i){
      int idM = i + f1*Nfp;
      double yM = yf(idM);

      for (int j = 0; j < Nfp; ++j){
	int idP = j + f2*Nfp;
	double yP = yf(idP);
	
	if (fabs(yM-yP) < NODETOL){
	  mapP(idM) = idP;
	  break;
	}
      }      
    }
  }

  printf("Done making node maps periodic\n");
  
}

// ========= 3d routines ===============

void HexMesh3d(Mesh *mesh, int Nx, int Ny, int Nz){
  // just make a uniform quadrilateral mesh on [-1,1]
  int Nxp = Nx+1;
  int Nyp = Ny+1;
  int Nzp = Nz+1;  
  int K = Nx*Ny*Nz;
  int Nv = Nxp*Nyp*Nzp;

  VectorXd x1D = VectorXd::LinSpaced(Nxp,-1,1);
  VectorXd y1D = VectorXd::LinSpaced(Nyp,-1,1);
  VectorXd z1D = VectorXd::LinSpaced(Nzp,-1,1);  

  VectorXd x(Nv),y(Nv),z(Nv);
  int sk = 0;
  for (int k = 0; k < Nzp; ++k){
    for (int j = 0; j < Nyp; ++j){
      for (int i = 0; i < Nxp; ++i){	
	x(sk) = x1D(i);
	y(sk) = y1D(j);
	z(sk) = z1D(k);
	++sk;
      }
    }
  }

  MatrixXi EToV(K,8);
  for (int e = 0; e < K; ++e){
    int k = e/(Nx*Ny);
    int j = (e - k*Nx*Ny)/Nx;
    int i = e % Nx;

    //printf("e = %d, ijk = %d, %d, %d\n",e,i,j,k);
    
    EToV(e,0) = i     + Nxp*j     + Nxp*Nyp*k;
    EToV(e,1) = i     + Nxp*(j+1) + Nxp*Nyp*k;
    EToV(e,2) = (i+1) + Nxp*j     + Nxp*Nyp*k;
    EToV(e,3) = (i+1) + Nxp*(j+1) + Nxp*Nyp*k;
    EToV(e,4) = i     + Nxp*j     + Nxp*Nyp*(k+1);
    EToV(e,5) = i     + Nxp*(j+1) + Nxp*Nyp*(k+1);
    EToV(e,6) = (i+1) + Nxp*j     + Nxp*Nyp*(k+1);
    EToV(e,7) = (i+1) + Nxp*(j+1) + Nxp*Nyp*(k+1);
  }
  //cout << "VX = [" << x << "];" << endl;
  //cout << "VY = [" << y << "];" << endl;
  //cout << "VZ = [" << z << "];" << endl;  
  //cout << "EToV = [" << EToV << "];" << endl;
  
  mesh->EToV = EToV;
  mesh->K = K;  
  mesh->Nv = Nv;
  mesh->VX = x;
  mesh->VY = y;
  mesh->VZ = z;  
}

// make quads 
void InitRefData3d(Mesh *mesh, int N){

  int Np1 = (N+1);

  // nodal points (GCL and all)
  VectorXd r1D, w1D;
  JacobiGL(N, 0, 0, r1D, w1D);

  // quad points 
  VectorXd rq1D, wq1D;
  JacobiGQ(N, 0, 0, rq1D, wq1D); printf("using GQ nodes\n");
  //JacobiGL(N, 0, 0, rq1D, wq1D); printf("using GLL nodes\n");
  
  int Np2 = Np1*Np1;
  int Np3 = Np1*Np2;
  
  VectorXd r(Np3),s(Np3),t(Np3);
  VectorXd rq(Np3),sq(Np3),tq(Np3),wq(Np3);  
  int sk = 0;
  for (int k = 0; k < Np1; ++k){
    for (int i = 0; i < Np1; ++i){
      for (int j = 0; j < Np1; ++j){
	r(sk) = r1D(i);
	s(sk) = r1D(j);
	t(sk) = r1D(k);	

	rq(sk) = rq1D(i);
	sq(sk) = rq1D(j);
	tq(sk) = rq1D(k);

	wq(sk) = wq1D(i)*wq1D(j)*wq1D(k);	
	++sk;
      }
    }
  }
    
  // nodal operators
  MatrixXd V1D = Vandermonde1D(N,r1D);
  MatrixXd Vr1D = GradVandermonde1D(N,r1D);
  MatrixXd D1D = mrdivide(Vr1D,V1D);

  MatrixXd I1D = MatrixXd::Identity(Np1,Np1);
  MatrixXd I2D = MatrixXd::Identity(Np2,Np2);  

  MatrixXd V = Vandermonde3DHex(N,r,s,t);
  MatrixXd Dr = kron(I1D,kron(D1D,I1D));
  MatrixXd Ds = kron(I2D,D1D);
  MatrixXd Dt = kron(D1D,I2D);

  //cout << "I2D = " << I2D << endl;
  //  cout << "D1D = " << D1D << endl;
  //  cout << "Dr = " << endl << Dr << endl;
  //  cout << "Ds = " << endl << Ds << endl;
  //  cout << "Dt = " << endl << Dt << endl;    

  VectorXd rq2(Np2),sq2(Np2),wq2(Np2);
  sk = 0;
  for (int i = 0; i < Np1; ++i){
    for (int j = 0; j < Np1; ++j){
      rq2(sk) = rq1D(i);
      sq2(sk) = rq1D(j);      
      wq2(sk) = wq1D(i)*wq1D(j);
      ++sk;
    }
  }
  VectorXd e(Np2); e.fill(1.0);
 
  // face quad points: r -/+, s-/+, t-/+
  int Nfaces = 6;
  int NfqNfaces = Nfaces * Np2;
  VectorXd rf(NfqNfaces); rf << -e, e, rq2, rq2, rq2, rq2;
  VectorXd sf(NfqNfaces); sf << rq2, rq2, -e, e, sq2, sq2;
  VectorXd tf(NfqNfaces); tf << sq2, sq2, sq2, sq2, -e, e;
  VectorXd wf(NfqNfaces); wf << wq2, wq2, wq2, wq2, wq2, wq2;

#if 0
  cout << "rq2 = [" << rq2 << "];" << endl;
  cout << "sq2 = [" << sq2 << "];" << endl;  
  
  cout << "rc = [" << r << "];" << endl;
  cout << "sc = [" << s << "];" << endl;
  cout << "tc = [" << t << "];" << endl;
  cout << setprecision(15) << "wq = [" << wq << "];" << endl;

  cout << "rfc = [" << rf << "];" << endl;
  cout << "sfc = [" << sf << "];" << endl;
  cout << "tfc = [" << tf << "];" << endl;
  cout << setprecision(15) << "wf = [" << wf << "];" << endl;
#endif
  
  // interp to quadrature
  MatrixXd Vqtmp = Vandermonde3DHex(N,rq,sq,tq);
  MatrixXd Vftmp = Vandermonde3DHex(N,rf,sf,tf);
  MatrixXd Vq = mrdivide(Vqtmp,V);
  MatrixXd Vf = mrdivide(Vftmp,V);

#if 0
  cout << "Vf = " << Vf << endl;
#endif
  
  // 1D quadrature operators
  MatrixXd Vq1D = Vandermonde1D(N,rq1D);
  MatrixXd Vrq1D = GradVandermonde1D(N,rq1D);
  MatrixXd Dq1D = mrdivide(Vrq1D,Vq1D);
  VectorXd rf1D(2); rf1D << -1.0,1.0;
  MatrixXd Vf1Dtmp = Vandermonde1D(N,rf1D);
  MatrixXd Vf1D = mrdivide(Vf1Dtmp,Vq1D);

  mesh->N = N;
  mesh->Np = Np3;
  mesh->Nfp = Np2;
  mesh->Nfaces = Nfaces;
  mesh->Nverts = 8;

  // GLL nodes (for geofacs)
  mesh->r = r;
  mesh->s = s;
  mesh->t = t;    
  mesh->Dr = Dr;
  mesh->Ds = Ds;
  mesh->Dt = Dt;  
  mesh->V = V;

  // interp to vol/face nodes
  mesh->rq = rq;
  mesh->sq = sq;
  mesh->tq = tq;  
  mesh->wq = wq;  
  mesh->Vq = Vq;

  mesh->rf = rf;
  mesh->sf = sf;
  mesh->tf = tf;  
  mesh->wf = wf;    
  mesh->Vf = Vf;

  //  cout << "rq = " << rq << endl;
  //  cout << "sq = " << sq << endl;  
  //cout << "rf = " << rf << endl;
  //cout << "sf = " << sf << endl;  
 
  // GQ nodes
  MatrixXd Q1D = wq1D.asDiagonal()*Dq1D;
  VectorXd invwq1D = wq1D.array().inverse();
  MatrixXd D1D_skew = invwq1D.asDiagonal()*(Q1D - Q1D.transpose());
  mesh->D1D = D1D_skew; // (Q - Q^T)
  //mesh->D1D = Dq1D;
  mesh->Vf1D = Vf1D;
  VectorXd wqInv = wq1D.array().inverse();
  mesh->Lf1D = wqInv.asDiagonal()*Vf1D.transpose(); // 1D lift - may not need

#if 0
  cout << "D1D = " << mesh->D1D << endl;
  cout << "Vf1D = " << mesh->Vf1D << endl;
  cout << "Lf1D = " << mesh->Lf1D << endl;   
#endif

  // LSRK-45 coefficients
  mesh->rk4a.resize(5);
  mesh->rk4a(0) =              0.0;
  mesh->rk4a(1) =  -567301805773.0 / 1357537059087.0;
  mesh->rk4a(2) = -2404267990393.0 / 2016746695238.0;
  mesh->rk4a(3) = -3550918686646.0 / 2091501179385.0;
  mesh->rk4a(4) = -1275806237668.0 /  842570457699.0;

  mesh->rk4b.resize(5);
  mesh->rk4b(0) =  1432997174477.0 /  9575080441755.0;
  mesh->rk4b(1) =  5161836677717.0 / 13612068292357.0;
  mesh->rk4b(2) =  1720146321549.0 /  2090206949498.0;
  mesh->rk4b(3) =  3134564353537.0 /  4481467310338.0;
  mesh->rk4b(4) =  2277821191437.0 / 14882151754819.0;

  mesh->rk4c.resize(6);
  mesh->rk4c(0) =              0.0;
  mesh->rk4c(1) =  1432997174477.0 / 9575080441755.0;
  mesh->rk4c(2) =  2526269341429.0 / 6820363962896.0;
  mesh->rk4c(3) =  2006345519317.0 / 3224310063776.0;
  mesh->rk4c(4) =  2802321613138.0 / 2924317926251.0;
  mesh->rk4c(5) =              1.0;
  
}


// computes physical nodal position based on affine 
void MapNodes3d(Mesh *mesh){
  int Nfaces = mesh->Nfaces;
  int K = mesh->K;
  int Nv = mesh->Nv;
  MatrixXi EToV = mesh->EToV;

  VectorXd VX = mesh->VX;
  VectorXd VY = mesh->VY;
  VectorXd VZ = mesh->VZ;    
 
  MatrixXd x1(EToV.cols(),K);
  MatrixXd y1(EToV.cols(),K);
  MatrixXd z1(EToV.cols(),K);    
  for (int e = 0; e < K; ++e){
    for (int v = 0; v < EToV.cols(); ++v){
      int vid = EToV(e,v);
      x1(v,e) = VX(vid);
      y1(v,e) = VY(vid);
      z1(v,e) = VZ(vid);            
    }
  }
  // interp to high order GLL nodes (matlab ordering)
  VectorXd r1(8),s1(8),t1(8);
  r1 << -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0;
  s1 << -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0;  
  t1 << -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0;  

  VectorXd r = mesh->r;
  VectorXd s = mesh->s;
  VectorXd t = mesh->t;  
  
  MatrixXd Vdm1 = Vandermonde3DHex(1,r1,s1,t1);
  MatrixXd VdmN = Vandermonde3DHex(1,r,s,t);
  //cout << "Vdm1 = " << Vdm1 << endl;
  MatrixXd V1 = mrdivide(VdmN,Vdm1);

  MatrixXd x = V1*x1;
  MatrixXd y = V1*y1;
  MatrixXd z = V1*z1;  
  mesh->x = x;
  mesh->y = y;
  mesh->z = z;

  /*cout << "VX = [" << VX << "];" << endl;
  cout << "VY = [" << VY << "];" << endl;
  cout << "VZ = [" << VZ << "];" << endl;
  cout << "x1 = [" << x1 << "];" << endl;
  cout << "y1 = [" << y1 << "];" << endl;
  cout << "z1 = [" << z1 << "];" << endl;  
  cout << "x = [" << x << "];" << endl;
  cout << "y = [" << y << "];" << endl;
  cout << "z = [" << z << "];" << endl;  */

}

// use Kopriva interpolation trick at GLL points
// - guarantees GCL on uniform p meshes

void GeometricFactors3d(Mesh *mesh){

  int Nfaces = mesh->Nfaces;
  int K = mesh->K;

  MatrixXd x = mesh->x;
  MatrixXd y = mesh->y;
  MatrixXd z = mesh->z;  

  // compute geofacs
  MatrixXd Dr = mesh->Dr;
  MatrixXd Ds = mesh->Ds;
  MatrixXd Dt = mesh->Dt;

  MatrixXd xr = Dr*x;  MatrixXd xs = Ds*x;  MatrixXd xt = Dt*x;
  MatrixXd yr = Dr*y;  MatrixXd ys = Ds*y;  MatrixXd yt = Dt*y;  
  MatrixXd zr = Dr*z;  MatrixXd zs = Ds*z;  MatrixXd zt = Dt*z;  

  // use curl-conservative form to compute GCL-satisfying geofacs
  MatrixXd rxJ = Dt * (ys.array()*z.array()).matrix() - Ds * (yt.array()*z.array()).matrix();
  MatrixXd sxJ = Dr * (yt.array()*z.array()).matrix() - Dt * (yr.array()*z.array()).matrix();
  MatrixXd txJ = Ds * (yr.array()*z.array()).matrix() - Dr * (ys.array()*z.array()).matrix();  

  MatrixXd ryJ = -(Dt * (xs.array()*z.array()).matrix() - Ds * (xt.array()*z.array()).matrix());
  MatrixXd syJ = -(Dr * (xt.array()*z.array()).matrix() - Dt * (xr.array()*z.array()).matrix());
  MatrixXd tyJ = -(Ds * (xr.array()*z.array()).matrix() - Dr * (xs.array()*z.array()).matrix());  

  MatrixXd rzJ = -(Dt * (ys.array()*x.array()).matrix() - Ds * (yt.array()*x.array()).matrix());
  MatrixXd szJ = -(Dr * (yt.array()*x.array()).matrix() - Dt * (yr.array()*x.array()).matrix());
  MatrixXd tzJ = -(Ds * (yr.array()*x.array()).matrix() - Dr * (ys.array()*x.array()).matrix());

  MatrixXd J =
      xr.array()*(ys.array()*zt.array() - zs.array()*yt.array())
    - yr.array()*(xs.array()*zt.array() - zs.array()*xt.array())
    + zr.array()*(xs.array()*yt.array() - ys.array()*xt.array());
  
  // all quantities stored at GLL points (must interp before using)
  mesh->rxJ = rxJ;
  mesh->sxJ = sxJ;
  mesh->txJ = txJ;  
  mesh->ryJ = ryJ;
  mesh->syJ = syJ;
  mesh->tyJ = tyJ;  
  mesh->rzJ = rzJ;
  mesh->szJ = szJ;
  mesh->tzJ = tzJ;  
  mesh->J = J;

#if 0
  double err = (Dr*rxJ + Ds*sxJ + Dt*txJ).cwiseAbs().sum() +
    (Dr*ryJ + Ds*syJ + Dt*tyJ).cwiseAbs().sum() +
    (Dr*rzJ + Ds*szJ + Dt*tzJ).cwiseAbs().sum();
  printf("GCL err = %g\n",err);
#endif
  
#if 0
  printf("vgeo = [\n");
  for (int e = 0; e < K; ++e){
    for (int i = 0; i < rxJ.rows(); ++i){
      printf(" %f, %f, %f, %f, %f, %f, %f, %f, %f\n",
	     rxJ(i,e),sxJ(i,e),txJ(i,e),
	     ryJ(i,e),syJ(i,e),tyJ(i,e),
	     rzJ(i,e),szJ(i,e),tzJ(i,e));	     

    }
  }
  printf("];\n");  
#endif
  
}

// compute geofacs with degree Ngeo geometry mapping
void GeometricFactors3d_Ngeo(Mesh *mesh, int Ngeo){

  int Nfaces = mesh->Nfaces;
  int K = mesh->K;

  // degree Ngeo nodal points (GCL and all)
  int Np1 = (Ngeo+1);  
  VectorXd r1D, w1D;
  JacobiGL(Ngeo, 0, 0, r1D, w1D);  
  int Np2 = Np1*Np1;
  int Np3 = Np1*Np2;  
  VectorXd rNg(Np3),sNg(Np3),tNg(Np3);
  int sk = 0;
  for (int k = 0; k < Np1; ++k){
    for (int i = 0; i < Np1; ++i){
      for (int j = 0; j < Np1; ++j){
	rNg(sk) = r1D(i);
	sNg(sk) = r1D(j);
	tNg(sk) = r1D(k);	
	++sk;
      }
    }
  }
  // degree Ngeo operators
  MatrixXd V1D = Vandermonde1D(Ngeo,r1D);
  MatrixXd Vr1D = GradVandermonde1D(Ngeo,r1D);
  MatrixXd D1D = mrdivide(Vr1D,V1D);
  MatrixXd I1D = MatrixXd::Identity(Np1,Np1);
  MatrixXd I2D = MatrixXd::Identity(Np2,Np2);  
  MatrixXd Dr = kron(I1D,kron(D1D,I1D));
  MatrixXd Ds = kron(I2D,D1D);
  MatrixXd Dt = kron(D1D,I2D);

  // interp degree N nodes to Ngeo
  MatrixXd VdmNgeo = Vandermonde3DHex(mesh->N,rNg,sNg,tNg);
  MatrixXd VNNgeo = mrdivide(VdmNgeo,mesh->V); 
  
  // construct degree Ngeo nodes
  MatrixXd x = VNNgeo * mesh->x;
  MatrixXd y = VNNgeo * mesh->y;
  MatrixXd z = VNNgeo * mesh->z;
  
  // compute geofacs
  MatrixXd xr = Dr*x;  MatrixXd xs = Ds*x;  MatrixXd xt = Dt*x;
  MatrixXd yr = Dr*y;  MatrixXd ys = Ds*y;  MatrixXd yt = Dt*y;  
  MatrixXd zr = Dr*z;  MatrixXd zs = Ds*z;  MatrixXd zt = Dt*z;  

  // normal curl method: GCL satisfied if 2*Ngeo - 2 < N, or Ngeo < floor(N/2) + 1
  MatrixXd rxJ = (ys.array()*zt.array() - zs.array()*yt.array());
  MatrixXd ryJ = -(xs.array()*zt.array() - zs.array()*xt.array());
  MatrixXd rzJ = (xs.array()*yt.array() - ys.array()*xt.array());

  MatrixXd sxJ = -(yr.array()*zt.array() - zr.array()*yt.array());
  MatrixXd syJ =  (xr.array()*zt.array() - zr.array()*xt.array());
  MatrixXd szJ = -(xr.array()*yt.array() - yr.array()*xt.array());

  MatrixXd txJ =  (yr.array()*zs.array() - zr.array()*ys.array());
  MatrixXd tyJ = -(xr.array()*zs.array() - zr.array()*xs.array());
  MatrixXd tzJ = (xr.array()*ys.array() - yr.array()*xs.array());

  MatrixXd J =
      xr.array()*(ys.array()*zt.array() - zs.array()*yt.array())
    - yr.array()*(xs.array()*zt.array() - zs.array()*xt.array())
    + zr.array()*(xs.array()*yt.array() - ys.array()*xt.array());
  
  // all quantities stored at GLL points (must interp before using)
  MatrixXd VdmN = Vandermonde3DHex(Ngeo,mesh->r,mesh->s,mesh->t);
  MatrixXd VNgeo = Vandermonde3DHex(Ngeo,rNg,sNg,tNg);
  MatrixXd VNgeoN = mrdivide(VdmN,VNgeo); // interp degree N nodes to Ngeo  
  mesh->rxJ = VNgeoN*rxJ;
  mesh->sxJ = VNgeoN*sxJ;
  mesh->txJ = VNgeoN*txJ;  
  mesh->ryJ = VNgeoN*ryJ;
  mesh->syJ = VNgeoN*syJ;
  mesh->tyJ = VNgeoN*tyJ;  
  mesh->rzJ = VNgeoN*rzJ;
  mesh->szJ = VNgeoN*szJ;
  mesh->tzJ = VNgeoN*tzJ;  
  mesh->J = VNgeoN*J;
 
#if 0
  printf("vgeo = [\n");
  for (int e = 0; e < K; ++e){
    for (int i = 0; i < rxJ.rows(); ++i){
      printf(" %f, %f, %f, %f, %f, %f, %f, %f, %f\n",
	     rxJ(i,e),sxJ(i,e),txJ(i,e),
	     ryJ(i,e),syJ(i,e),tyJ(i,e),
	     rzJ(i,e),szJ(i,e),tzJ(i,e));	     

    }
  }
  printf("];\n");  
#endif
  
}


void Normals3d(Mesh *mesh){

  MatrixXd Vf = mesh->Vf;
  MatrixXd rxJf = Vf * mesh->rxJ;
  MatrixXd sxJf = Vf * mesh->sxJ;
  MatrixXd txJf = Vf * mesh->txJ;  

  MatrixXd ryJf = Vf * mesh->ryJ;
  MatrixXd syJf = Vf * mesh->syJ;
  MatrixXd tyJf = Vf * mesh->tyJ;
  
  MatrixXd rzJf = Vf * mesh->rzJ;
  MatrixXd szJf = Vf * mesh->szJ;
  MatrixXd tzJf = Vf * mesh->tzJ;

#if 0
  printf("vfgeo = [\n");
  for (int e = 0; e < mesh->K; ++e){
    for (int i = 0; i < rxJf.rows(); ++i){
      printf(" %f, %f, %f, %f, %f, %f, %f, %f, %f\n",
	     rxJf(i,e),sxJf(i,e),txJf(i,e),
	     ryJf(i,e),syJf(i,e),tyJf(i,e),
	     rzJf(i,e),szJf(i,e),tzJf(i,e));	     

    }
  }
  printf("];\n");  
#endif
  

  int Nfaces = mesh->Nfaces;
  int Nfp = mesh->Nfp;
  VectorXd e(Nfp); e.fill(1.0);
  VectorXd zz(Nfp); zz.fill(0.0);  
  VectorXd nrJ(Vf.rows()); nrJ << -e, e, zz, zz, zz, zz;
  VectorXd nsJ(Vf.rows()); nsJ << zz, zz, -e,e , zz, zz;
  VectorXd ntJ(Vf.rows()); ntJ << zz, zz, zz, zz, -e, e;

  //  for (int i = 0; i < Nfp*Nfaces; ++i){
  //    printf("nrstJ = %f, %f, %f\n",nrJ(i),nsJ(i),ntJ(i));
  //  }

  MatrixXd nxJ = nrJ.asDiagonal()*rxJf + nsJ.asDiagonal()*sxJf + ntJ.asDiagonal()*txJf;
  MatrixXd nyJ = nrJ.asDiagonal()*ryJf + nsJ.asDiagonal()*syJf + ntJ.asDiagonal()*tyJf;
  MatrixXd nzJ = nrJ.asDiagonal()*rzJf + nsJ.asDiagonal()*szJf + ntJ.asDiagonal()*tzJf;
  MatrixXd sJ = (nxJ.array().square() + nyJ.array().square() + nzJ.array().square()).sqrt();

  mesh->nxJ = nxJ;
  mesh->nyJ = nyJ;
  mesh->nzJ = nzJ;  
  mesh->sJ = sJ;

#if 0
  for (int e = 0; e < mesh->K; ++e){
    for (int i = 0; i < nxJ.rows(); ++i){
      printf("nxyzJ = %f, %f, %f, %f\n",nxJ(i,e),nyJ(i,e),nzJ(i,e),sJ(i,e));
    }
  }
#endif
  
}


void MakeNodeMapsPeriodic3d(Mesh *mesh, MatrixXd xf, MatrixXd yf, MatrixXd zf,
			    double DX, double DY, double DZ, MatrixXi &mapP){
  
  int Nfp = mesh->Nfp;
  int Nfaces = mesh->Nfaces;
  int K = mesh->K;

  VectorXi isBoundaryFace(Nfaces*K);
  isBoundaryFace.fill(0);
  for (int e = 0; e < K; ++e){
    for (int f = 0; f < Nfaces; ++f){
      if (e==mesh->EToE(e,f)){
	isBoundaryFace(f + e*Nfaces) = 1;
      }
    }
  }
  
  xf.resize(Nfp,Nfaces*K);
  yf.resize(Nfp,Nfaces*K);
  zf.resize(Nfp,Nfaces*K);    
  vector<pair<int,int> > xfaces,yfaces,zfaces;  
  for (int f1 = 0; f1 < Nfaces*K; ++f1){
    for (int f2 = 0; f2 < Nfaces*K; ++f2){
      if (isBoundaryFace(f1) && isBoundaryFace(f2)){
	// distance b/w faces = diff b/w min/max coeffs 
	double dx = .5*(fabs(xf.col(f1).maxCoeff()-xf.col(f2).maxCoeff())
			+ fabs(xf.col(f1).minCoeff()-xf.col(f2).minCoeff()));
	double dy = .5*(fabs(yf.col(f1).maxCoeff()-yf.col(f2).maxCoeff())
			+ fabs(yf.col(f1).minCoeff()-yf.col(f2).minCoeff()));
	double dz = .5*(fabs(zf.col(f1).maxCoeff()-zf.col(f2).maxCoeff())
			+ fabs(zf.col(f1).minCoeff()-zf.col(f2).minCoeff()));
	
	if ((fabs(dx-DX)<NODETOL) & (dy < NODETOL) & (dz < NODETOL)){
	  xfaces.push_back(make_pair(f1,f2)); 
	}else if ((fabs(dy-DY)<NODETOL) & (dx < NODETOL) & (dz < NODETOL)){
	  yfaces.push_back(make_pair(f1,f2)); 
	}else if ((fabs(dz-DZ)<NODETOL) & (dx < NODETOL) & (dy < NODETOL)){
	  zfaces.push_back(make_pair(f1,f2));
	}
      }
    }
  }
  printf("num xfaces = %d, num yfaces = %d, num zfaces = %d\n",xfaces.size(),yfaces.size(),zfaces.size());
  //  cout << "xf = " << xf << endl;
  //  cout << "yf = " << yf << endl;
  //  return;

  // find node maps in x  
  for (int f = 0; f < xfaces.size(); ++f){
    int f1 = xfaces[f].first;
    int f2 = xfaces[f].second;
    for (int i = 0; i < Nfp; ++i){
      int idM = i + f1*Nfp;
      double yM = yf(idM);
      double zM = zf(idM);

      bool match_found = false;	      

      for (int j = 0; j < Nfp; ++j){
	int idP = j + f2*Nfp;
	double yP = yf(idP);
	double zP = zf(idP);	

	double dist = fabs(yM-yP) + fabs(zM-zP);

	if (dist < NODETOL){
	  mapP(idM) = idP;
	  match_found = true;
	  break;	  
	}

      }
      if (match_found==false){
	printf("Match not found for xface\n");
      }
    }
  }

  // find node maps in y
  for (int f = 0; f < yfaces.size(); ++f){
    int f1 = yfaces[f].first;
    int f2 = yfaces[f].second;
    for (int i = 0; i < Nfp; ++i){
      int idM = i + f1*Nfp;
      double xM = xf(idM);
      double zM = zf(idM);      

      bool match_found = false;
      
      for (int j = 0; j < Nfp; ++j){
	int idP = j + f2*Nfp;
	double xP = xf(idP);
	double zP = zf(idP);      

	double dist = fabs(xM-xP) + fabs(zM-zP);	
	if (dist < NODETOL){
	  mapP(idM) = idP;
	  match_found=true;
	  break;
	}
      }
      if (match_found==false){
	printf("Match not found for yface\n");
      }

    }
  }

  // find node maps in z
  for (int f = 0; f < zfaces.size(); ++f){
    int f1 = zfaces[f].first;
    int f2 = zfaces[f].second;
    for (int i = 0; i < Nfp; ++i){
      int idM = i + f1*Nfp;
      double xM = xf(idM);
      double yM = yf(idM);      

      bool match_found = false;
      for (int j = 0; j < Nfp; ++j){
	int idP = j + f2*Nfp;
	double xP = xf(idP);
	double yP = yf(idP);      

	double dist = fabs(xM-xP) + fabs(yM-yP);	
	if (dist < NODETOL){
	  mapP(idM) = idP;
	  match_found = true;
	  break;
	}
      }
      if (match_found==false){
	printf("Match not found for zface\n");
      }

    }
  }
  
  printf("Done making node maps periodic\n");
  
}




// ======== dimension independent routines ============

// general: input nodes (ex: quadrature nodes), get map back
void BuildFaceNodeMaps(Mesh *mesh, MatrixXd xf, MatrixXd yf, MatrixXd zf, MatrixXi &mapP){

  int K       = mesh->K;
  int Nfaces  = mesh->Nfaces;

  mapP.resize(xf.rows(),K);
  mapP.fill(-1);

  int Nfpts = xf.rows() / Nfaces; // assume same # qpts per face (i.e. not hybrid meshes)

  double x1, y1, z1, x2, y2, z2, d12;

  //  cout << "EToF = " << endl << mesh->EToF << endl;
  //  cout << "EToE = " << endl << mesh->EToE << endl;  
  
  printf("Building face node maps\n");

  for(int k1=0;k1<K;++k1){

    for(int f1=0;f1<Nfaces;++f1){

      // find neighbor                                           
      int k2 = mesh->EToE(k1,f1);
      int f2 = mesh->EToF(k1,f1);

      //printf("k1 = %d, k2 = %d\n",k1,k2);

      if(k1==k2){
        for(int i = 0; i < Nfpts; ++i){
          mapP(i + f1*Nfpts,k1) = i + f1*Nfpts + xf.rows()*k1;
        }
      }else{

	//printf("Looking for matches on elem %d, face %d to elem %d, face %d\n",k1,f1,k2,f2);

        MatrixXd xM(Nfpts,3);
        MatrixXd xP(Nfpts,3);
        for(int i = 0; i < Nfpts; ++i){
          int id1 = i + Nfpts * f1;
          x1 = xf(id1,k1); y1 = yf(id1,k1); z1 = zf(id1,k1);
          int id2 = i + Nfpts * f2;
          x2 = xf(id2,k2); y2 = yf(id2,k2); z2 = zf(id2,k2);

          xM(i,0) = x1; xM(i,1) = y1; xM(i,2) = z1;
          xP(i,0) = x2; xP(i,1) = y2; xP(i,2) = z2;
        }

        for(int i = 0; i < Nfpts; ++i){

          double minDist = 1000.0;

          int id1 = i + Nfpts * f1;
          x1 = xf(id1,k1); y1 = yf(id1,k1); z1 = zf(id1,k1);
          bool node_found = false;
          for(int j = 0; j < Nfpts; ++j){

            int id2 = j + Nfpts * f2;
            x2 = xf(id2,k2); y2 = yf(id2,k2); z2 = zf(id2,k2);
            // find distance between these nodes                 
            d12 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
            minDist = min(minDist,d12);
            if (d12<NODETOL){
              mapP(id1,k1) = id2 + xf.rows()*k2;
              node_found = true;
              break;
            }
          }
          if(!node_found){
            printf("BuildFaceNodeMaps: lost node %d on elems %d, %d. min dist = %g\n",i,k1,k2,minDist);
            cout << "xM = " << endl << xM << endl;
            cout << "xP = " << endl << xP << endl;
          }
        }
      }

    }// faces
  }// k                                                                                                        

  return;
}


// connects quads (or hexes)
typedef unsigned long long ulint;
// check if "a < b" for faceInfo (i.e. all entries of a.first less than that of b.first)
static bool faceCompare(const std::pair<vector<int>,ulint> &a,const std::pair<vector<int>,ulint> &b){
  for (int i = 0; i < a.first.size(); ++i){
    if (a.first[i] != b.first[i]){
      return a.first[i] < b.first[i];
    }
  }
  return a.first.back() < b.first.back();
}

void ConnectElems(Mesh *mesh, int dim){

  int Nfaces = mesh->Nfaces;
  int K = mesh->K;
  int Nv = mesh->Nv;
  MatrixXi EToV = mesh->EToV;
  //cout << "EToV = " << endl << EToV << endl;
  
  // make face lists and fv scaling
  MatrixXi fvlist;
  if (dim==2){
    fvlist.resize(Nfaces,2); 
    fvlist << 0,1,
      1,2,
      2,3,
      3,0;
  }else if (dim==3){
    fvlist.resize(Nfaces,4);
    // todo: add 3D face list
    // faces ordered: r -/+, s -/+, t -/+
    // assume nodes = matlab meshgrid ordering
    /*
      .    5-----7
      .   /|    /|
      .  4-|---6 |
      .  | |   | |
      .  | 1---|-3
      .  |/    |/
      .  0-----2
     */
    // face ordering: right hand rule for outward normal (start shouldn't matter?)
    /*
    fvlist << 1,5,6,2,
      4,8,7,4,
      3,7,5,1,
      2,6,8,4,
      1,2,4,3,
      7,8,6,5;
    */
    fvlist << 0,4,5,1,
      3,7,6,2,
      2,6,4,0,
      1,5,7,3,
      0,1,3,2,
      6,7,5,4;
  }
  //cout << "fvlist = " << fvlist << endl;
  
  int Nfaceverts = fvlist.cols();

  // if K too large this can fail.
  vector<pair<vector<int>, int> > faceInfo;

  //cout << "K = " << K << endl;
  
  ulint i = 0;
  for (int f = 0; f < Nfaces; ++f){
    for (int e = 0; e < K; ++e){
      vector<int> fv; 
      for (int v = 0; v < Nfaceverts; ++v){
	//printf("e = %d, f = %d, v = %d: fv = %d\n",e,f,v,fvlist(f,v));
	fv.push_back(EToV(e,fvlist(f,v)));
      }
      std::sort(fv.begin(),fv.end()); // sort fv (nodes on face)
      //cout << fv.size() << endl; 
      //cout << fv[0] << ", " << fv[1] <<", " << fv[2] <<", " << fv[3] << endl;
      
      faceInfo.push_back(make_pair(fv,i));
      ++i;
    }
  }

  // sort by first elem (second goes along for ride)  
  /*
  for (int i = 0; i < faceInfo.size(); ++i){
    printf("faceInfo: %d %d, %d\n",faceInfo[i].first[0],faceInfo[i].first[1],faceInfo[i].second);
  }
  */
  std::sort(faceInfo.begin(),faceInfo.end(),faceCompare);
  
  // initialize EToE, EToF
  int low,hi;
  low = 0;  hi = K-1; 
  VectorXi eVec = VectorXi::LinSpaced((hi-low)+1,low,low+hi);
  low = 0;  hi = Nfaces-1;
  VectorXi fVec = VectorXi::LinSpaced((hi-low)+1,low,low+hi);  
  MatrixXi EToE = eVec.replicate(1,Nfaces);
  MatrixXi EToF = fVec.transpose().replicate(K,1);
  // cout << EToE << endl;
  // cout << EToF << endl;

  // check for consecutive matches
  for (int i = 0; i < Nfaces*K - 1; ++i){

    // check if all vertices match
    bool match = true;
    for (int v = 0; v < Nfaceverts; ++v){
      if (faceInfo[i].first[v]!=faceInfo[i+1].first[v]){
	match = false;
      }
    }

    if (match){  // if face[i] matches face[i+1]
      ulint id1 = faceInfo[i].second;
      ulint id2 = faceInfo[i+1].second;
      int e1 = EToE(id1);
      int f1 = EToF(id1);      

      int e2 = EToE(id2);
      int f2 = EToF(id2);            

      //printf("matched on i = %d: e = %d,%d, f = %d,%d\n",i,e1,e2,f1,f2);
      
      EToE(e1,f1) = e2;
      EToF(e1,f1) = f2;      
      EToE(e2,f2) = e1;
      EToF(e2,f2) = f1;      
    }
  }

  //  cout << "EToE = " << endl << EToE << endl;
  //  cout << "EToF = " << endl << EToF << endl;  
  
  mesh->EToE = EToE;
  mesh->EToF = EToF;  
}
