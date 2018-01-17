%%
Globals3D
N = 5;

filename = 'Grid/cube1.msh';
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGmsh3D(filename);
VX = VX*2; VY = VY*2; VZ = VZ*2; % map -1,1
StartUp3D;
% 
% % make curvilinear
% a = .125;
% x = x + a*cos(.5*pi*x).*cos(.5*pi*y).*cos(.5*pi*z);
% y = y + a*cos(.5*pi*x).*cos(.5*pi*y).*cos(.5*pi*z);
% z = z + a*cos(.5*pi*x).*cos(.5*pi*y).*cos(.5*pi*z);

a = .05;
x = r - a*(r+s).^N;
y = s + a/2*(s+t).^N;
z = t + a/3*(r+t).^N;

[rq sq tq wq] = tet_cubature(2*N);
Vq = Vandermonde3D(N,rq,sq,tq)/V;
M =(Vq'*diag(wq)*Vq);
Pq = M\(Vq'*diag(wq));

Drq = Vq*Dr;
Dsq = Vq*Ds;
Dtq = Vq*Dt;

Nfp = (N+1)*(N+2)/2;
[rqtri sqtri wqtri] = Cubature2D(2*N);
Vfqf = Vandermonde2D(N,rqtri,sqtri)/Vandermonde2D(N,r(Fmask(:,1)),s(Fmask(:,1)));
Vfqf = kron(eye(4),Vfqf);
rfq = Vfqf*r(Fmask(:));
sfq = Vfqf*s(Fmask(:));
tfq = Vfqf*t(Fmask(:));
wfq = repmat(wqtri,4,1);
Vfq = Vandermonde3D(N,rfq,sfq,tfq)/V;
Lq = M\(Vfq'*diag(wfq));

Drfq = Vfq*Dr;
Dsfq = Vfq*Ds;
Dtfq = Vfq*Dt;

%%

% original cross product form
xrq = Drq*x; xsq = Dsq*x; xtq = Dtq*x;
yrq = Drq*y; ysq = Dsq*y; ytq = Dtq*y;
zrq = Drq*z; zsq = Dsq*z; ztq = Dtq*z;
J = xrq.*(ysq.*ztq-zsq.*ytq) - yrq.*(xsq.*ztq-zsq.*xtq) + zrq.*(xsq.*ytq-ysq.*xtq);

rxJ =  (ysq.*ztq - zsq.*ytq); 
sxJ = -(yrq.*ztq - zrq.*ytq); 
txJ =  (yrq.*zsq - zrq.*ysq); 

ryJ = -(xsq.*ztq - zsq.*xtq); 
syJ =  (xrq.*ztq - zrq.*xtq); 
tyJ = -(xrq.*zsq - zrq.*xsq); 

rzJ = (xsq.*ytq - ysq.*xtq);
szJ = -(xrq.*ytq - yrq.*xtq);
tzJ = (xrq.*ysq - yrq.*xsq);

%% make normals


xrq = Drfq*x; xsq = Dsfq*x; xtq = Dtfq*x;
yrq = Drfq*y; ysq = Dsfq*y; ytq = Dtfq*y;
zrq = Drfq*z; zsq = Dsfq*z; ztq = Dtfq*z;

Jf = xrq.*(ysq.*ztq-zsq.*ytq) - yrq.*(xsq.*ztq-zsq.*xtq) + zrq.*(xsq.*ytq-ysq.*xtq);

rxJf =  (ysq.*ztq - zsq.*ytq); 
sxJf = -(yrq.*ztq - zrq.*ytq); 
txJf =  (yrq.*zsq - zrq.*ysq); 

ryJf = -(xsq.*ztq - zsq.*xtq); 
syJf =  (xrq.*ztq - zrq.*xtq); 
tyJf = -(xrq.*zsq - zrq.*xsq); 

rzJf = (xsq.*ytq - ysq.*xtq);
szJf = -(xrq.*ytq - yrq.*xtq);
tzJf = (xrq.*ysq - yrq.*xsq);

Nfqf = length(rqtri);
e = ones(Nfqf,1); z = zeros(Nfqf,1);
nrJ = [z;z;e;-e];
nsJ = [z;-e;e;z];
ntJ = [-e;z;e;z];

nxJ = rxJf.*nrJ + sxJf.*nsJ + txJf.*ntJ;
nyJ = ryJf.*nrJ + syJf.*nsJ + tyJf.*ntJ;
nzJ = rzJf.*nrJ + szJf.*nsJ + tzJf.*ntJ;

nx = nxJ./Jf;
ny = nyJ./Jf;
nz = nzJ./Jf;

sJ = sqrt(nx.^2 + ny.^2 + nz.^2);
nx = nx./sJ;
ny = ny./sJ;
nz = nz./sJ;

sJ = sJ.*Jf;    

%% div-free projection
% min (1/2)*||G-Gh||^2_L2 = .5*(G,G) - (G,Gh) + .5*(Gh,Gh)
% st        (Gh,grad(v)) = L*(nx)
ids = 1:Np;

Div = [Dr'*M Ds'*M Dt'*M]; % weak divergence: linearly dependent cols
MM = blkdiag(M,M,M);

UDF = null(Div); % div-free V
[Uu, ~, ~] = svd(MM*UDF);
UD = Uu(:,size(UDF,2)+1:end); % L2 orthogonal component to UD
[Utest,~,~] = svd(M*ones(Np,1)); % remove constant mode from nodal basis
Utest = Utest(:,2:end);

fx = [rxJ;sxJ;txJ];
gx = nxJ;

UDFPq = UDF*((UDF'*MM*UDF)\(UDF'*kron(eye(3),Vq'*diag(wq))));
UDPq = UD*((Utest'*Div*UD)\(Utest'*Vfq'*diag(wfq)));
u_divfree = UDFPq*fx; % L2 projection
u_div = UDPq*gx; % weak lift with null space removed
ux = u_div + u_divfree;

norm(kron(eye(3),Vq)*ux - [rxJ;sxJ;txJ],'fro') % constrained projection
norm(kron(eye(3),Vq*Pq)*[rxJ;sxJ;txJ] - [rxJ;sxJ;txJ],'fro') % L2 projection

rxJ = ux(ids,:);
sxJ = ux(ids+Np,:);
txJ = ux(ids+2*Np,:);

Vq*[rxJ sxJ txJ]

%% apply derivatives of constant u
    
% norm(Dr'*M*rxJ + Ds'*M*sxJ + Dt'*M*txJ - M*Lq*nxJ,'fro')


