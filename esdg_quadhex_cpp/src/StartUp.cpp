#include "fem.h"
#include "Basis.h"

// make quads 
void InitRefData2d(Mesh *mesh, int N){

  int Np1 = (N+1);

  // nodal points (GCL and all)
  VectorXd r1D;
  JacobiGL(N, 0, 0, r1D);

  MatrixXd rmat,smat;
  meshgrid(r1D,r1D, smat,rmat);
  VectorXd r = Map<VectorXd>(rmat.data(), rmat.cols()*rmat.rows());
  VectorXd s = Map<VectorXd>(smat.data(), smat.cols()*smat.rows());
  //  cout << "r = " << r << endl;
  //  cout << "s = " << s << endl;  

  // quad points 
  VectorXd rq1D, wq1D;
  JacobiGQ(N, 0, 0, rq1D, wq1D);
  
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
  // cout << "Vf1D = " << endl << Vf1D << endl;  

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

  //cout << "D1D = " << mesh->D1D << endl;
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

  mesh->rk4c.resize(5);
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

// connects quads
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
  }
  int Nfaceverts = fvlist.cols();

  // if K too large this can fail.
  vector<pair<vector<int>, int> > faceInfo;

  ulint i = 0;
  for (int f = 0; f < Nfaces; ++f){
    for (int e = 0; e < K; ++e){
      vector<int> fv; 
      for (int v = 0; v < Nfaceverts; ++v){
	fv.push_back(EToV(e,fvlist(f,v)));
      }
      std::sort(fv.begin(),fv.end()); // sort fv (nodes on face)
      //cout << fv.size() << endl; 
      //cout << fv[0] << ", " << fv[1] << endl;
      
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

// computes physical nodal position
/*
void ComputePhysNodes(Mesh *mesh){
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
*/

void GeometricFactors2d(Mesh *mesh){

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
  //  cout << "x1 = " << endl << x1 << endl;
  //  cout << "y1 = " << endl << y1 << endl;  

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
  MatrixXd J = -xs.array()*yr.array() + xr.array()*ys.array();
  
  // cout << "VX = " << endl << x1 << endl;
  // cout << "VY = " << endl << y1 << endl;
  // cout << "VX = " << endl << VX << endl;
  // cout << "VY = " << endl << VY << endl;

  //  cout << "EToV = " << EToV << endl;
  //  cout << "x = " << endl << x << endl;
  //  cout << "y = " << endl << y << endl;
  
  /*
  cout << "rxJ = " << endl << rxJ << endl;
  cout << "ryJ = " << endl << ryJ << endl;
  cout << "sxJ = " << endl << sxJ << endl;
  cout << "syJ = " << endl << syJ << endl;    
  */

  // all quantities stored at GLL points (must interp before using)
  mesh->x = x;
  mesh->y = y;
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


void MakeNodeMapsPeriodic2D(Mesh *mesh, MatrixXd xf, MatrixXd yf, double DX, double DY, MatrixXi &mapP){
  
  int Nfp = mesh->Nfp;
  int Nfaces = mesh->Nfaces;
  int K = mesh->K;

  xf.resize(Nfp,Nfaces*K);
  yf.resize(Nfp,Nfaces*K);  
  vector<pair<int,int> > xfaces,yfaces;  
  for (int f1 = 0; f1 < Nfaces*K; ++f1){
    for (int f2 = 0; f2 < Nfaces*K; ++f2){

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
  //  printf("num xfaces = %d, num yfaces = %d\n",xfaces.size(),yfaces.size());
  //  cout << "xf = " << xf << endl;
  //  cout << "yf = " << yf << endl;  

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
