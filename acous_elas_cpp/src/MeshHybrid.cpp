#include "fem.h"

#include <iostream>
#include <fstream>


// reads meshes for other types of elements
int KTet = 0, KHex = 0, KPri = 0, KPyr = 0;
const int HEX=1,PRI=2,PYR=3,TET=4;
int firstElement (const int t){
 switch(t){
  case HEX : return 0;
  case PRI : return KHex;
  case PYR : return KHex + KPri;
  case TET : return KHex + KPri + KPyr;
  }
}

int getNumVerts (int e, MatrixXi EToV) {
  if (EToV(e,4) == -1) return 4;
  if (EToV(e,5) == -1) return 5;
  if (EToV(e,6) == -1) return 6;
  return 8;
}
int getNumFaces (int e, MatrixXi EToV) {
  if (EToV(e,4) == -1) return 4;
  if (EToV(e,5) == -1) return 5;
  if (EToV(e,6) == -1) return 5;
  return 6;
}
int getNumVertexOnFace (int e, int f, MatrixXi EToV) {
  if (EToV(e,4) == -1) return 3; // tet
  if (EToV(e,5) == -1) return f > 3 ? 4 : 3; // pyr
  if (EToV(e,6) == -1) return f > 1 ? 4 : 3; // prism
  return 4; // hex
}

// defines vertex orderings
int getVertexOnFace (int e, int f, int v, MatrixXi EToV) {
  int vid;

  if (EToV(e,4) == -1){   // tet
    static int fv[4][3] = {{1,3,2}, {1,2,4}, {3,1,4}, {2,3,4}};
    vid = fv[f][v]-1;
  }else if (EToV(e,5) == -1){  // pyr
    static const int fv[5][4] = {
      {1, 2, 5, 0},{4, 1, 5, 0},{2, 3, 5, 0},{3, 4, 5, 0},
      {1, 4, 3, 2}
    };
    vid = fv[f][v]-1;
  }else if (EToV(e,6) == -1){  // prism    
    static const int fv[5][4] = {
      {1, 3, 2, 0}, {4, 5, 6, 0},
      {1, 2, 5, 4}, {3, 1, 4, 6}, {2, 3, 6, 5}
    };
    vid = fv[f][v]-1;
  }else{  // hex
    static int fv[6][4] = {{1,4,3,2},{1,2,6,5},
			   {2,3,7,6},{3,4,8,7},
			   {4,1,5,8},{5,6,7,8}};
    vid = fv[f][v]-1;    
  }
  return EToV(e,vid);
    
}



Mesh *ReadGmshHybrid(char *filename){

  Mesh *mesh = (Mesh*) calloc(1, sizeof(Mesh));

  /// open the mesh file
  std::ifstream meshfile(filename);
  if(!meshfile.is_open()){
    std::cout << "Mesh reading ERROR: Mesh '" << filename << "' not opened" << std::endl;
    exit(1);
  }
  
  int dummy, allVertices, allElements;

  while(meshfile.good()){

    string line;
    getline(meshfile, line);

    /// scan vertices in gmsh file
    if(line.find("$Nodes") != string::npos){

      /// read number of vertices
      meshfile >> allVertices;
      mesh->Nv = allVertices;

      mesh->VX.resize(allVertices);
      mesh->VY.resize(allVertices);
      mesh->VZ.resize(allVertices);      

      /// resize vectors for the coordinates of vertices
      //     VX = (double*) calloc(mesh->Nv, sizeof(double));
      //     VY = (double*) calloc(mesh->Nv, sizeof(double));
      //     VZ = (double*) calloc(mesh->Nv, sizeof(double));     

      /// scan line of each vertex, ignoring gmsh index of vertex
      for (int v = 0; v < allVertices; ++v){
	int vx,vy,vz;
	//meshfile >> dummy >> VX >> VY >> VZ;
	//printf("VXYZ = [%.4f, %.4f, %.4f]\n",VX,VY,VZ);
	meshfile >> dummy >> vx >> vy >> vz;
	mesh->VX(v) = vx;
	mesh->VY(v) = vy;
	mesh->VZ(v) = vz;	
      }

    }

    /// scan elements in gmsh file
    if(line.find("$Elements") != string::npos){

      /// read number of elements
      meshfile >> allElements;
      getline(meshfile, line);

      /// record location in file
      int marker = meshfile.tellg();
      int eGmshNum, eGeoType, Ntags;
      VectorXi EToGmshE(allElements);

      /// count number of elements
      int KTri = 0, KQua = 0;
      //int KTet = 0, KHex = 0, KPri = 0, KPyr = 0;
      for (int k = 1; k <= allElements; ++k){
	meshfile >> eGmshNum >> eGeoType;
	switch(eGeoType){
	case 2: KTri++; break;
	case 3: KQua++; break;
	case 4: KTet++; break;
	case 5: KHex++; break;
	case 6: KPri++; break;
	case 7: KPyr++; break;
	case 1: break;  /// 2-node line
	case 15: break; /// 1-node point
	default:
	  cout << "HybriDG WARNING: unknow element type read in gmsh file ("
	       << eGeoType << ")" << endl;
	}
	getline(meshfile, line);
      }
      /// rewind to beginning of element list
      meshfile.seekg(marker);

      // number of element
      int KFace = KTri + KQua;
      mesh->K = KTet + KHex + KPri + KPyr;

      /// shape EToV to hold Elements triples
      int Nverts = 8; // max Nverts
      int NvertsPerFace = 4; // max verts per face

      //mesh->EToV = BuildIntMatrix(mesh->K,Nverts);
      mesh->EToV.resize(mesh->K,Nverts);
      mesh->EToV.fill(-1);
      VectorXi ETag(mesh->K); ETag.fill(0);
      MatrixXi physFToV(NvertsPerFace,KFace);
      VectorXi physFType(KFace,1); physFType.fill(0);
      VectorXi ENumGmsh(mesh->K,1); ENumGmsh.fill(0);
      VectorXi ENumGmsh_orig(mesh->K,1); ENumGmsh_orig.fill(0);

      /// save info on triangular/tetrahedral elements
      KFace = 1;
      int iHex = 0;
      int iPri = 0;
      int iTet = 0;
      int iPyr = 0;
      for (int k = 0; k < mesh->K; ++k) {
	meshfile >> eGmshNum >> eGeoType >> Ntags;
	/// Triangle
	if(eGeoType == 2){
	  meshfile >> physFType(KFace);           /// save the face group (physical gmsh label)
	  for(int tag=1; tag<Ntags; tag++)        /// dummy tags
	    meshfile >> dummy;
	  for(int v=0; v<3; ++v)
	    meshfile >> physFToV(v, KFace);       /// save the associated vertices
	  KFace++;
	}
	/// Quad
	else if(eGeoType == 3){
	  meshfile >> physFType(KFace);           /// save the face group (physical gmsh label)
	  for(int tag=1; tag<Ntags; tag++)        /// dummy tags
	    meshfile >> dummy;
	  for(int v=0; v<4; ++v)
	    meshfile >> physFToV(v, KFace);       /// save the associated vertices
	  KFace++;
	}
	
	/// Tetra
	else if(eGeoType == 4){
	  int posK = firstElement(TET) + iTet;
	  //printf("on elem %d, gmsh num = %d\n",posK,eGmshNum);
	  ENumGmsh(posK) = eGmshNum;              /// save the gmsh number
	  ENumGmsh_orig(posK) = eGmshNum;
	  meshfile >> ETag(posK);                /// save the element group (physical gmsh label)
	  for(int tag=1; tag<Ntags; tag++)        /// dummy tags
	    meshfile >> dummy;
	  for(int v=0; v<4; ++v)
	    meshfile >> mesh->EToV(posK,v);            /// save the associated vertices
	  iTet++;
	}
	/// Hex
	else if(eGeoType == 5){
	  int posK = firstElement(HEX) + iHex;
	  ENumGmsh(posK) = eGmshNum;              /// save the gmsh number
	  ENumGmsh_orig(posK,1) = eGmshNum;
	  meshfile >> ETag(posK);                /// save the element group (physical gmsh label)
	  for(int tag=1; tag<Ntags; tag++)        /// dummy tags
	    meshfile >> dummy;
	  for(int v=0; v<8; ++v)
	    meshfile >> mesh->EToV(posK,v);            /// save the associated vertices
	  iHex++;
	}
	/// Prism
	else if(eGeoType == 6){
	  int posK = firstElement(PRI) + iPri;
	  //printf("on elem %d, gmsh num = %d\n",posK,eGmshNum);
	  ENumGmsh(posK,1) = eGmshNum;              /// save the gmsh number
	  ENumGmsh_orig(posK,1) = eGmshNum;
	  meshfile >> ETag(posK);                /// save the element group (physical gmsh label)
	  for(int tag=1; tag<Ntags; tag++)        /// dummy tags
	    meshfile >> dummy;
	  for(int v=0; v<6; ++v)
	    meshfile >> mesh->EToV(posK,v);            /// save the associated vertices
	  iPri++;
	}
	/// Pyramid
	else if(eGeoType == 7){
	  int posK = firstElement(PYR) + iPyr;
	  ENumGmsh(posK) = eGmshNum;              /// save the gmsh number
	  ENumGmsh_orig(posK,1) = eGmshNum;
	  meshfile >> ETag(posK);                /// save the element group (physical gmsh label)
	  for(int tag=1; tag<Ntags; tag++)        /// dummy tags
	    meshfile >> dummy;
	  for(int v=0; v<5; ++v)
	    meshfile >> mesh->EToV(posK,v);            /// save the associated vertices
	  iPyr++;
	}
	getline(meshfile, line);
      }

      cout << "Hybrid mesh read in: " << allVertices << " vertices and " << mesh->K << " elements in GMSH file " << endl;
      cout << "      (" << KTri << " triangles, " << KQua << " quadrangles, " << KTet << " tetrahedron, " << KHex << " hexahedron, " << KPri << " prisms, " << KPyr << " pyramids)" <<endl;

      mesh->GX.resize(mesh->K,Nverts);
      mesh->GY.resize(mesh->K,Nverts);
      mesh->GZ.resize(mesh->K,Nverts);      

      int sk = 0, v;
      for (int k = 0; k < mesh->K; ++k){
	
        /* correct to 0-index */
        for(v=0; v<Nverts;++v){
	  if (mesh->EToV(k,v)>0){ // if not initialized to -1 already
	    --(mesh->EToV(k,v));
	  }
	}	
	
        for(v=0;v<Nverts;++v){
	  //          mesh->GX(k,v) = VX[mesh->EToV(k,v)];
	  //          mesh->GY(k,v) = VY[mesh->EToV(k,v)];
	  //	  mesh->GZ(k,v) = VZ[mesh->EToV(k,v)];
	  mesh->GX(k,v) = mesh->VX(mesh->EToV(k,v));
	  mesh->GY(k,v) = mesh->VY(mesh->EToV(k,v));
	  mesh->GZ(k,v) = mesh->VZ(mesh->EToV(k,v));	  
        }
	
      }
      
    }
  }

#if 1
  printf("VXYZ = \n");
  for (int i = 0; i < mesh->Nv; ++i){
    printf("%f %f %f\n",mesh->VX(i),mesh->VY(i),mesh->VZ(i));
  }
 
  printf("EToV = \n");
  for (int e = 0; e < mesh->K; ++e){
    for (int v = 0; v < 8; ++v){
      printf("%d ", mesh->EToV(e,v));
    }
    printf("\n");
  }
#endif

  return mesh;
}


typedef struct {
  int e, f;
  int v[4]; // max # faces = 4
}face4;

static int compare_faces_hybrid(const void *obj1, const void *obj2){
  face4 *e1 = (face4*) obj1;
  face4 *e2 = (face4*) obj2;

  int a1 = e1->v[0], a2 = e1->v[1], a3 = e1->v[2], a4 = e1->v[3];
  int b1 = e2->v[0], b2 = e2->v[1], b3 = e2->v[2], b4 = e2->v[3];

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

void FacePairHybrid(Mesh *mesh){
  
  /// build a face-to-vertex connectivity array "FToV"
  int NvertsPerFace = 4; // max # verts per face 
  int Nfaces = 6; // max number of faces
 
  MatrixXi FToV(NvertsPerFace+2, Nfaces*mesh->K);
  FToV.fill(-1);

  int counter = 0;
  for(int k=0; k < mesh->K; ++k){

    int Nf = getNumFaces(k,mesh->EToV);
    
    for(int f = 0; f < Nf; ++f){
      
      int Nv = getNumVertexOnFace(k, f, mesh->EToV);  /// number of vertices on the face

      /// raw sort of vertex numbers on this face
      VectorXi ns(NvertsPerFace); ns.fill(-1);      
      for(int i = 0; i < Nv; ++i){
	ns(i) = getVertexOnFace(k,f,i,mesh->EToV);
      }
      std::sort(ns.data(), ns.data()+ns.size()); // sort verts before using
      cout << "face " << counter << ", verts = " << endl << ns.transpose() << endl;
      for(int i=0; i< Nv; ++i){
	FToV(i, counter) = ns(i);  // vertex numbers
      }
      FToV(NvertsPerFace, counter) = k;  // element number
      FToV(NvertsPerFace+1, counter) = f;  /// face number

      //printf("counter = %d\n",counter);
      ++counter;
    }
  }
 
  /// sort by 1~4 row (forgot column major convention)
  int Ntotalfaces = counter;
  //  printf("Ntotalfaces = %d\n",Ntotalfaces);
  //  cout << "FToV = " << endl << FToV << endl;
  
  // sort to find matches
  face4 *myfaces = (face4*) calloc(Ntotalfaces, sizeof(face4));  
  for(int i = 0; i < Ntotalfaces; ++i){
    for (int v = 0; v < NvertsPerFace; ++v){
      myfaces[i].v[v] = FToV(v,i);
    }
    myfaces[i].e = FToV(NvertsPerFace,i);
    myfaces[i].f = FToV(NvertsPerFace+1,i);    
  }

  qsort(myfaces, Ntotalfaces, sizeof(face4), compare_faces_hybrid);

  /// build 'EToE' and 'EToF' connectivity arrays 
  mesh->EToE.resize(mesh->K,Nfaces); mesh->EToE.fill(-1);
  mesh->EToF.resize(mesh->K,Nfaces); mesh->EToF.fill(-1);

  for(int i = 0; i < Ntotalfaces; ++i){
    
    /// find neighbors
    if(myfaces[i].v[0]==myfaces[i+1].v[0] &&
       myfaces[i].v[1]==myfaces[i+1].v[1] &&
       myfaces[i].v[2]==myfaces[i+1].v[2] &&
       myfaces[i].v[3]==myfaces[i+1].v[3]){

      int v1 = myfaces[i].v[0];
      int v2 = myfaces[i].v[1];
      int v3 = myfaces[i].v[2];
      int v4 = myfaces[i].v[3];      

      int u1 = myfaces[i+1].v[0];
      int u2 = myfaces[i+1].v[1];
      int u3 = myfaces[i+1].v[2];
      int u4 = myfaces[i+1].v[3];      

      int e1 = myfaces[i].e;
      int f1 = myfaces[i].f;
      int e2 = myfaces[i+1].e;
      int f2 = myfaces[i+1].f;

#if 0
      printf("face %d: v1 = %d,%d,%d,%d, v2 = %d,%d,%d,%d,  e1 = %d, e2 = %d, f1 = %d, f2 = %d\n",
	     i,v1,v2,v3,v4,
	     u1,u2,u3,u4,
	     e1,e2,f1,f2);
#endif
      
      mesh->EToE(e1,f1) = e2;
      mesh->EToE(e2,f2) = e1;
      mesh->EToF(e1,f1) = f2;
      mesh->EToF(e2,f2) = f1;
    }
  }

  for(int e = 0; e < mesh->K; ++e){

    int Nf = getNumFaces(e,mesh->EToV);
    
    for(int f = 0; f < Nf; ++f){
      if (mesh->EToE(e,f)==-1){
	mesh->EToE(e,f) = e;
	mesh->EToF(e,f) = f;
      }
    }
  }
#if 1
  cout << "EToE = " << endl << mesh->EToE << endl;
  cout << "EToF = " << endl << mesh->EToF << endl;  
#endif
}


// start up for hybrid meshes
void StartUpWedge(Mesh *mesh){

  int N = p_N;
  int NpTri = (N+1)*(N+2)/2;
  int NpQuad = (N+1)*(N+1);
  int Np = (N+1)*NpTri;

  // get tri nodes from tet faces
  VectorXd r2,s2,t2;
  Nodes3D(N,r2,s2,t2);
  VectorXd rtri = r2.head(NpTri);
  VectorXd stri = s2.head(NpTri);

  // TP vol nodes
  VectorXd r1D, w1D;
  JacobiGL(N, 0, 0, r1D);
  //JacobiGQ(N, 0, 0, r1D,w1D);

  // for quad faces
  VectorXd r1Dq,w1Dq;
  //JacobiGL(N, 0, 0, r1Dq,w1Dq);
  JacobiGQ(N, 0, 0, r1Dq, w1Dq);

  // make wedge volume nodes
  mesh->r.resize(Np);
  mesh->s.resize(Np);
  mesh->t.resize(Np);  
  for (int i = 0; i < N+1; ++i){ 
    for (int j = 0; j < NpTri; ++j){
      mesh->r(j + i*NpTri) = rtri(j);
      mesh->s(j + i*NpTri) = stri(j);
      mesh->t(j + i*NpTri) = r1D(i);	
    }
  }

  // volume VDM matrix
  MatrixXd V,Vr,Vs,Vt;
  WedgeBasis(N,mesh->r,mesh->s,mesh->t,
	     V,Vr,Vs,Vt);
  
  int Nfaces = 5;
#if 0
  // get reference surface nodes
  mesh->Fmask.resize(NpQuad,Nfaces);
  mesh->Fmask.fill(-1); 
  int skf[Nfaces];
  for(int f=0; f<Nfaces; ++f){
    skf[f] = 0;
  }
  for (int i = 0; i < Np; ++i){
    double r = mesh->r(i);
    double s = mesh->s(i);
    double t = mesh->t(i);
    int f = -1;
    if (fabs(t+1.0)<NODETOL){
      mesh->Fmask(skf[0],0) = i;
      ++skf[0];
    }
    if (fabs(t-1.0)<NODETOL){
      mesh->Fmask(skf[1],1) = i;
      ++skf[1];
    }
    if (fabs(s+1.0)<NODETOL){
      mesh->Fmask(skf[2],2) = i;
      ++skf[2];
    }
    if (fabs(r+1.0)<NODETOL){      
      mesh->Fmask(skf[3],3) = i;
      ++skf[3];
    }
    if (fabs(r+s)<NODETOL){
      mesh->Fmask(skf[4],4) = i;
      ++skf[4];
    }
  }
#endif

  // construct face ndoes
  int NfpNfaces = 2*NpTri+3*NpQuad;
  VectorXd rf(NfpNfaces);
  VectorXd sf(NfpNfaces);
  VectorXd tf(NfpNfaces);
  
  // t = -1, use triangular nodes
  // (assume they're the first NpTri nodes of wedge)
  int sk = 0;
  for (int i = 0; i < NpTri; ++i){
    rf(sk) = mesh->r(i);
    sf(sk) = mesh->s(i);
    tf(sk) = -1.0;
    ++sk;
  }
  // t = 1
  for (int i = 0; i < NpTri; ++i){
    rf(sk) = mesh->r(i);
    sf(sk) = mesh->s(i);
    tf(sk) = 1.0;
    ++sk;    
  }
  // s = -1
  for (int i = 0; i < N+1; ++i){
    for (int j = 0; j < N+1; ++j){
      rf(sk) = r1Dq(i); 
      sf(sk) = -1.0;    
      tf(sk) = r1Dq(j); 
      ++sk;
    }
  }
  // r = -1
  for (int i = 0; i < N+1; ++i){
    for (int j = 0; j < N+1; ++j){
      rf(sk) = -1.0;
      sf(sk) = r1Dq(i);
      tf(sk) = r1Dq(j);
      ++sk;
    }
  }
  // r + s = 0
  for (int i = 0; i < N+1; ++i){
    for (int j = 0; j < N+1; ++j){
      rf(sk) = -r1Dq(j);
      sf(sk) = r1Dq(j);
      tf(sk) = r1Dq(i);
      ++sk;
    }
  }
    
  // map nodes to physical coords  
  VectorXd r1(6),s1(6),t1(6);
  r1 << -1., 1., -1., -1., 1., -1.;
  s1 << -1., -1., 1., -1., -1., 1.;
  t1 << -1., -1., -1., 1., 1., 1.;
  MatrixXd V1,V1r,V1s,V1t;
  WedgeBasis(1,r1,s1,t1,
	     V1,V1r,V1s,V1t);

  MatrixXd V1Ntmp,V1Nr,V1Ns,V1Nt;
  WedgeBasis(1,mesh->r,mesh->s,mesh->t,
	     V1Ntmp,V1Nr,V1Ns,V1Nt);
  MatrixXd V1N = mrdivide(V1Ntmp,V1);

  // make volume nodes
  mesh->x.resize(Np,mesh->K);
  mesh->y.resize(Np,mesh->K);
  mesh->z.resize(Np,mesh->K);  
  for(int e = 0; e < mesh->K; ++e){
    MatrixXd vxyz(6,3);
    for(int v = 0; v < 6; ++v){
      int vid = mesh->EToV(e,v);
      vxyz(v,0) = mesh->VX(vid);
      vxyz(v,1) = mesh->VY(vid);
      vxyz(v,2) = mesh->VZ(vid);      
    }
    MatrixXd xyz = V1N*vxyz;
    mesh->x.col(e) = xyz.col(0);
    mesh->y.col(e) = xyz.col(1);
    mesh->z.col(e) = xyz.col(2);    
  }

  // make surface nodes
  MatrixXd Vf,Vrf,Vsf,Vtf;
  WedgeBasis(N,rf,sf,tf,
	     Vf,Vrf,Vsf,Vtf);  
  mesh->Vfq = mrdivide(Vf,V);
  MatrixXd xf = mesh->Vfq * mesh->x;
  MatrixXd yf = mesh->Vfq * mesh->y;
  MatrixXd zf = mesh->Vfq * mesh->z;

#if 0
  cout << "x = [" << endl << mesh->x << "];"<<endl;
  cout << "y = [" << endl << mesh->y << "];"<<endl;
  cout << "z = [" << endl << mesh->z << "];"<<endl;  
  cout << "xf = [" << endl << xf << "];"<<endl;
  cout << "yf = [" << endl << yf << "];"<<endl;
  cout << "zf = [" << endl << zf << "];"<<endl;  
#endif
  
  // pair surface nodes
  MatrixXi mapM(NfpNfaces,mesh->K);
  MatrixXi mapP(NfpNfaces,mesh->K);

  for (int e = 0; e < mesh->K; ++e){

    // face node offsets
    VectorXi NfpSk(5);
    NfpSk << 0, NpTri, 2*NpTri, 2*NpTri+NpQuad,2*NpTri+2*NpQuad;
    
    int skf = 0;
    
    for (int f = 0; f < Nfaces; ++f){

      int Nfpts = (f < 2) ? NpTri : NpQuad;
      
      // find nbr
      int enbr = mesh->EToE(e,f);
      int fnbr = mesh->EToF(e,f);

      double x1, x2, y1,y2,z1,z2;

      for (int i = 0; i < Nfpts; ++i){
	mapM(i+NfpSk(f),e) = i + NfpSk(f) + NfpNfaces*e;
      }
      
      if (e==enbr){
	for(int i = 0; i < Nfpts; ++i){
	  mapP(i + NfpSk(f),e) = i + NfpSk(f) + NfpNfaces*e;
	  ++skf;
	}
      }else{

	MatrixXd xM(Nfpts,3);  MatrixXd xP(Nfpts,3);
        for(int i = 0; i < Nfpts; ++i){
	  int id1 = i + NfpSk(f);
	  x1 = xf(id1,e); y1 = yf(id1,e); z1 = zf(id1,e);
	  xM(i,0) = x1; xM(i,1) = y1; xM(i,2) = z1;
          int id2 = i + NfpSk(fnbr);
          x2 = xf(id2,enbr); y2 = yf(id2,enbr); z2 = zf(id2,enbr);  
          xP(i,0) = x2; xP(i,1) = y2; xP(i,2) = z2;
        }


	for(int i = 0; i < Nfpts; ++i){
	  
	  double minDist = 1000.0;
	  
	  int id1 = i + NfpSk(f);
	  x1 = xf(id1,e); y1 = yf(id1,e); z1 = zf(id1,e);
	  
          bool node_found = false;
	  for(int j = 0; j < Nfpts; ++j){
	    
	    int id2 = j + NfpSk(fnbr);
	    x2 = xf(id2,enbr); y2 = yf(id2,enbr); z2 = zf(id2,enbr);
	    
	    // find distance between these nodes
	    int d12 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
	    minDist = min(minDist,d12);
	    if (d12<NODETOL){
	      mapP(id1,e) = id2 + NfpNfaces*enbr;
	      node_found = true;
	      break;
	    }
	  }
          if(!node_found){
            printf("Lost node %d on elem %d, face %d. min dist = %g\n",i,e,f,minDist);
            cout << "xM = " << endl << xM << endl;
            cout << "xP = " << endl << xP << endl;
          }
	}
      }
    }
  }
  //  cout << "mapP = [" << endl << mapP << "];" << endl;


  // initialize tensor product wedge operators
  MatrixXd VTri = Vandermonde2D(N,rtri,stri);
  MatrixXd VrTri,VsTri;
  GradVandermonde2D(N,rtri,stri,VrTri,VsTri);
  mesh->Dr = mrdivide(VrTri,VTri);
  mesh->Ds = mrdivide(VsTri,VTri);

  MatrixXd V1D = Vandermonde1D(N,r1D);
  MatrixXd Vr1D;
  GradVandermonde1D(N,r1D,Vr1D);  
  mesh->Dt = mrdivide(Vr1D,V1D); // 1D op

  // quadrature operators
  VectorXd rqtri,sqtri,wqtri;
  tri_cubature(N, rqtri, sqtri, wqtri);  
  MatrixXd VqTri = Vandermonde2D(N,rqtri,sqtri);
  VqTri = mrdivide(VqTri,VTri); 
  MatrixXd VrqTri,VsqTri;
  GradVandermonde2D(N,rqtri,sqtri,VrqTri,VsqTri);
  
  MatrixXd V1Dq = Vandermonde1D(N,r1Dq);
  V1Dq = mrdivide(V1Dq,V1D);
  MatrixXd Vr1Dq;
  GradVandermonde1D(N,r1Dq,Vr1Dq);
  Vr1Dq = mrdivide(Vr1Dq,V1D);

  // cubature projection matrices
  MatrixXd invMtri = VTri*VTri.transpose();
  MatrixXd invM1D = V1D*V1D.transpose();
  MatrixXd PqTri = invMtri*VqTri.transpose()*wqtri.asDiagonal();
  MatrixXd Pq1D = invM1D*V1Dq.transpose()*w1Dq.asDiagonal();

  // for skew-sym wedges
  MatrixXd PrqTri = invMtri*VrqTri.transpose()*wqtri.asDiagonal();
  MatrixXd PsqTri = invMtri*VsqTri.transpose()*wqtri.asDiagonal();  
  MatrixXd Ptq1D = invM1D*Vr1Dq.transpose()*w1Dq.asDiagonal();    

  
  
}//end startup of wedge

#if 0
MatrixXd GeometricFactorsHybrid(Mesh *mesh, int e,
				MatrixXd Dr, MatrixXd Ds, MatrixXd Dt){
  
  double x1 = mesh->GX(e,0), y1 =  mesh->GY(e,0), z1 =  mesh->GZ(e,0);
  double x2 = mesh->GX(e,1), y2 =  mesh->GY(e,1), z2 =  mesh->GZ(e,1);
  double x3 = mesh->GX(e,2), y3 =  mesh->GY(e,2), z3 =  mesh->GZ(e,2);
  double x4 = mesh->GX(e,3), y4 =  mesh->GY(e,3), z4 =  mesh->GZ(e,3);

  int Nv = getNumVerts(e,mesh->EToV);
  MatrixXd xyz(Nv,3);
  for (int v = 0; v < Nv; ++v){
    xyz(v,0) = mesh->GX(e,v);
    xyz(v,2) = mesh->GY(e,v);
    xyz(v,3) = mesh->GZ(e,v);
  }

  MatrixXd Dxyzr = Dr*xyz;
  MatrixXd Dxyzs = Ds*xyz;
  MatrixXd Dxyzt = Dt*xyz;

  VectorXd dxdr=Dxyzr.col(0);  VectorXd dydr=Dxyzr.col(1);  VectorXd dzdr=Dxyzr.col(2);
  VectorXd dxds=Dxyzs.col(0);  VectorXd dyds=Dxyzs.col(1);  VectorXd dzds=Dxyzs.col(2);  
  VectorXd dxdt=Dxyzt.col(0);  VectorXd dydt=Dxyzt.col(1);  VectorXd dzdt=Dxyzt.col(2);  
  
  VectorXd J =
    dxdr.array()*(dyds.array()*dzdt.array()-dzds.array()*dydt.array())
    -dydr.array()*(dxds.array()*dzdt.array()-dzds.array()*dxdt.array())
    +dzdr.array()*(dxds.array()*dydt.array()-dyds.array()*dxdt.array());

  VectorXd drdx =  (dyds.array()*dzdt.array() - dzds.array()*dydt.array())/J.array();
  VectorXd drdy = -(dxds.array()*dzdt.array() - dzds.array()*dxdt.array())/J.array();
  VectorXd drdz =  (dxds.array()*dydt.array() - dyds.array()*dxdt.array())/J.array();

  VectorXd dsdx = -(dydr.array()*dzdt.array() - dzdr.array()*dydt.array())/J.array();
  VectorXd dsdy =  (dxdr.array()*dzdt.array() - dzdr.array()*dxdt.array())/J.array();
  VectorXd dsdz = -(dxdr.array()*dydt.array() - dydr.array()*dxdt.array())/J.array();

  VectorXd dtdx =  (dydr.array()*dzds.array() - dzdr.array()*dyds.array())/J.array();
  VectorXd dtdy = -(dxdr.array()*dzds.array() - dzdr.array()*dxds.array())/J.array();
  VectorXd dtdz =  (dxdr.array()*dyds.array() - dydr.array()*dxds.array())/J.array();

  MatrixXd vgeo(Dr.rows(),10);
  vgeo.col(0) = drdx;  vgeo.col(1) = drdy;  vgeo.col(2) = drdz;  
  vgeo.col(3) = dsdx;  vgeo.col(4) = dsdy;  vgeo.col(5) = dsdz;  
  vgeo.col(6) = dtdx;  vgeo.col(7) = dtdy;  vgeo.col(8) = dtdz;
  vgeo.col(9) = J;
  
}

MatrixXd sgeo NormalsHybrid(Mesh *mesh, int e,
			    MatrixXd Dr, MatrixXd Ds, MatrixXd Dt){

  MatrixXd vgeo;
  GeometricFactorsHybrid(mesh,e,Dr,Ds,Dt,vgeo);

  int Nv = getNumVerts(e,mesh->EToV);
  int Nf = getNumFaces(e,mesh->EToV);

  VectorXd drdx = vgeo.col(0); VectorXd drdy = vgeo.col(1); VectorXd drdz = vgeo.col(2);
  VectorXd dsdx = vgeo.col(3); VectorXd dsdy = vgeo.col(4); VectorXd dsdz = vgeo.col(5);
  VectorXd dtdx = vgeo.col(6); VectorXd dtdy = vgeo.col(7); VectorXd dtdz = vgeo.col(8);
  VectorXd J = vgeo.col(9);
 
  if (Nf==4){ // tet
    printf("Nf = %d\n",Nf);
  }else if (Nf==5 && Nv==5){ // pyramid
    printf("Nf = %d\n",Nf);
  }else if (Nf==5 && Nv==6){ // prism
    
    int fid = 0;
    

    
  }else{ // hex
    printf("Nf = %d\n",Nf);
  }
  
  MatrixXd sgeo(Dr.rows(),4);
  //  sgeo.col(0) = nx;
}

// general: input nodes, get map back
void BuildFaceMapsWedge(Mesh *mesh, MatrixXd xf, MatrixXd yf, MatrixXd zf,
			MatrixXi &mapP){

  int K       = mesh->K;
  int Nfaces  = mesh->Nfaces;
  //  printf("In buildfacenodemaps Nfaces = %d\n",Nfaces);

  mapP.resize(xf.rows(),K);
  mapP.fill(-1);

  int Nfpts = xf.rows() / Nfaces; // assume same # qpts per face

  int m;
  int k1,p1,n1,f1,k2,p2,n2,f2;

  double x1, y1, z1, x2, y2, z2, d12;

  printf("Hello %d\n", 1337);

  // first build local
  for(k1=0;k1<K;++k1){

    for(f1=0;f1<Nfaces;++f1){

      // find neighbor
      k2 = mesh->EToE(k1,f1);
      f2 = mesh->EToF(k1,f1);

      if(k1==k2){
	for(int i = 0; i < Nfpts; ++i){
	  mapP(i + f1*Nfpts,k1) = i + f1*Nfpts + xf.rows()*k1;
	}
      }else{

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

#endif

