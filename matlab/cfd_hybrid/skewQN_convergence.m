%% 1D

clear
Globals1D

N = 2;
Kvec = [4 8 16 32 64 128 256];

a = 3;
u = @(x) exp(sin(a*x));
du = @(x) exp(sin(a*x)).*a.*cos(a*x);
v = @(x) x.^N + 0*cos(x);
iex = integral(@(x) du(x).*v(x),-1,1,'AbsTol',1e-10);
int_ex = integral(@(x) u(x),-1,1,'AbsTol',1e-14);

[rq wq] = JacobiGL(0,0,N); 
% [rq wq] = clenshaw_curtis(2*N+1,3); 

sk = 1;
for kk = 1:length(Kvec)
    K1D = Kvec(kk);
    
    [Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);
    StartUp1D;
        
    Vq = Vandermonde1D(N,rq)/V;
    M = (Vq'*diag(wq)*Vq);
    Pq = M\(Vq'*diag(wq));    
    Vf = Vandermonde1D(N,[-1;1])/V;
    xq = Vq*x;
    xf = Vf*x;            
    xN = [xq;xf];
    VN = [Vq;Vf];
    rxN = VN*rx;   
    wJq = diag(wq)*(Vq*J);
    
    Qr = diag(wq)*Vq*Dr*Pq ;
    B = diag([-1;1]);
    E = Vf*Pq;
    QN = [Qr - .5*E'*B*E .5*E'*B
        -.5*B*E .5*B];
    BN = QN+QN';
    
    PN = M\([Vq' Vf']);
    ustrong = PN*(rxN.*(QN*u(xN)));
    uweak = PN*(rxN.*((BN-QN')*u(xN)));
%     uweak = (M\(Vf'*diag([-1;1])*u(xf) - J.*rx.*(Dr'*M*Pq*u(xq))))./J;
    
    [rq2 wq2] = JacobiGR(0,0,N+4);    
    Vq2 = Vandermonde1D(N,rq2)/V;
    Pq2 = (Vq2'*diag(wq2)*Vq2)\(Vq2'*diag(wq2));
    xq2 = Vq2*x;
    wJq2 = diag(wq2)*(Vq2*J);
    errstrong(kk) = sqrt(sum(sum(wJq2.*(Vq2*ustrong - Vq2*Pq2*du(xq2)).^2)));
    errweak(kk) = sqrt(sum(sum(wJq2.*(Vq2*uweak - Vq2*Pq2*du(xq2)).^2)));
    
    qerr(kk) = abs(sum(sum(wJq.*(u(xq)))) - int_ex);
    
%     v = xN.^N;
%     wJq = diag(wq)*(Vq*J);
%     i1 = 0;
%     for e = 1:K
%         i1 = i1 + rx(1,e)*v(:,e)'*QN*u(xN(:,e));
%     end
    
%     errweakD(kk) = sqrt(sum(sum(wJq.*(
end

h = 2./Kvec;
h = h/h(1);
loglog(h,errstrong,'o--','linewidth',2,'markersize',16)
hold on
loglog(h,errweak,'x--','linewidth',2,'markersize',16)
loglog(h,qerr,'d--','linewidth',2,'markersize',16)

loglog(h,errstrong(1)*h.^N,'--','linewidth',2,'markersize',16)
loglog(h,errstrong(1)*h.^(N+1),'--','linewidth',2,'markersize',16)
loglog(h,errstrong(1)*h.^(N+2),'--','linewidth',2,'markersize',16)
text(h(end),.75*errstrong(1)*h(end)^N,'N','fontsize',16)
text(h(end),.75*errstrong(1)*h(end)^(N+1),'N+1','fontsize',16)
text(h(end),.75*errstrong(1)*h(end)^(N+2),'N+2','fontsize',16)

loglog(h,.1*qerr(1)*h.^(2*N),'--','linewidth',2,'markersize',16)
loglog(h,.1*qerr(1)*h.^(2*N+1),'--','linewidth',2,'markersize',16)
loglog(h,.01*qerr(1)*h.^(2*N+2),'--','linewidth',2,'markersize',16)
text(h(end),.1*qerr(1)*h(end)^(2*N),'2N','fontsize',16)
text(h(end),.1*qerr(1)*h(end)^(2*N+1),'2N+1','fontsize',16)
text(h(end),.01*qerr(1)*h(end)^(2*N+2),'2N+2','fontsize',16)
% return


%% tris

clear

useQuads = 0; mypath;

Globals2D

N = 2;
Kvec = [2 4 8 16 32 64];

sk = 1;
for K1D = Kvec
    [Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
    
%     a = .125;
%     VX = VX + a*sin(pi*VX).*sin(pi*VY);
%     VY = VY + a*sin(pi*VX).*sin(pi*VY);
    
    StartUp2D;
        
    [rq sq wq] = Cubature2D(2*N+2); % integrate u*v*c
    [rq1D wq1D] = JacobiGQ(0,0,N); % face quad 
%     [rq1D wq1D] = clenshaw_curtis(2*N-3,3); % face quad
    
    Vq = Vandermonde2D(N,rq,sq)/V;
    M = Vq'*diag(wq)*Vq;
    Pq = M\(Vq'*diag(wq)); % J's cancel out
    
    VqPq = Vq*Pq;
    
    xq = Vq*x; yq = Vq*y;
    Jq = Vq*J;
    
    rxJ = rx.*J; sxJ = sx.*J;
    ryJ = ry.*J; syJ = sy.*J;    
    
    rfq = [rq1D; -rq1D; -ones(size(rq1D))];
    sfq = [-ones(size(rq1D)); rq1D; -rq1D];
    wfq = [wq1D; wq1D; wq1D];
    Vq1D = Vandermonde1D(N,rq1D)/Vandermonde1D(N,JacobiGL(0,0,N));    
    
    Vfq = Vandermonde2D(N,rfq,sfq)/V;
    VfPq = Vfq*Pq;
    Vfqf = kron(eye(3),Vq1D);
    Mf = Vfq'*diag(wfq)*Vfq;
    Lq = M\(Vfq'*diag(wfq));
    
    Pq1D = (Vq1D'*diag(wq1D)*Vq1D) \ (Vq1D'*diag(wq1D));
    Pfqf = kron(eye(3),Pq1D);
        
%     nx = Vfqf*nx;
%     ny = Vfqf*ny;
%     sJ = Vfqf*sJ;
%     nxJ = (nx.*sJ);
%     nyJ = (ny.*sJ);
%     Fscale = Vfqf*Fscale;
    
    nrJ = [-zeros(size(rq1D)); ones(size(rq1D)); -ones(size(rq1D)); ];
    nsJ = [-ones(size(rq1D)); ones(size(rq1D)); -zeros(size(rq1D)); ];
    Nq = length(rq);
    nrJq = repmat(nrJ',Nq,1);
    nsJq = repmat(nsJ',Nq,1);
    
    WN = diag([wq;wfq]);
    Qr = diag(wq)*Vq*Dr*Pq;
    Qs = diag(wq)*Vq*Ds*Pq;
    
    Br = diag(wfq.*nrJ);
    Bs = diag(wfq.*nsJ);
    E = VfPq;
    QNr = [Qr-.5*E'*Br*E .5*E'*Br;
        -.5*Br*E .5*Br];
    QNs = [Qs-.5*E'*Bs*E .5*E'*Bs;
        -.5*Bs*E .5*Bs];
    
    QNrskew = .5*[Qr-Qr' VfPq'*diag(wfq.*nrJ);
        -diag(wfq.*nrJ)*VfPq diag(0*wfq)];
    QNsskew = .5*[Qs-Qs' VfPq'*diag(wfq.*nsJ);
        -diag(wfq.*nsJ)*VfPq diag(0*wfq)];
    QNrskew(abs(QNrskew)<1e-8) = 0; QNrskew = sparse(QNrskew);
    QNsskew(abs(QNsskew)<1e-8) = 0; QNsskew = sparse(QNsskew);
    
    % make skew symmetric diff matrices
    DNrskew = diag(1./[wq;wfq])*QNrskew;
    DNsskew = diag(1./[wq;wfq])*QNsskew;
    
    BNr = diag([zeros(size(Qr,1),1); nrJ]);
    BNs = diag([zeros(size(Qr,1),1); nsJ]);
    
    VqPN = Vq*(M\[Vq' Vfq']);
    
    f = @(x,y) exp(sin(x+y));
    dfex = @(x,y) exp(sin(x + y)).*cos(x + y);
    
    VN = [Vq;Vfq];
    rxJN = VN*rxJ;
    sxJN = VN*sxJ;
    %     fN = VN*Pq*f(xq,yq);
    fN = f(VN*x,VN*y);
    
    % skew deriv
    dfJ = VqPN*(rxJN.*((QNrskew + .5*WN*BNr)*fN) + sxJN.*((QNsskew + .5*WN*BNs)*fN));
    df = dfJ./(Vq*J);
    
    % strong/weak deriv
    dfJ = VqPN*(rxJN.*((QNr)*fN) + sxJN.*((QNs)*fN));
    dfstrong = dfJ./(Vq*J);
    dfJ = VqPN*(rxJN.*((WN*BNr-QNr')*fN) + sxJN.*((WN*BNs-QNs')*fN));
    dfweak = dfJ./(Vq*J);
    
    [rq2 sq2 wq2] = Cubature2D(3*N+4); % integrate u*v*c
    Vq2 = Vandermonde2D(N,rq2,sq2)/V;
    Pq2 = (Vq2'*diag(wq2)*Vq2)\(Vq2'*diag(wq2));
    xq2 = Vq2*x;
    yq2 = Vq2*y;    
    wJq = diag(wq)*(Vq*J);
    wJq2 = diag(wq2)*(Vq2*J);
    err(sk) = sqrt(sum(sum(wJq2.*(Vq2*Pq2*dfex(xq2,yq2)-Vq2*Pq*df).^2)));
    errstrong(sk) = sqrt(sum(sum(wJq2.*(Vq2*Pq2*dfex(xq2,yq2)-Vq2*Pq*dfstrong).^2)));
    errweak(sk) = sqrt(sum(sum(wJq2.*(Vq2*Pq2*dfex(xq2,yq2)-Vq2*Pq*dfweak).^2)));
    
%     qerr(sk) = abs(sum(sum(f(xq,yq).*wJq)) - sum(sum(f(xq2,yq2).*wJq2)));
    
    % test just surf quad err
    xf = Vfq*x;  yf = Vfq*y;
    
    [rq1D wq1D] = JacobiGQ(0,0,4*N+4); 
    rfq2 = [rq1D; -rq1D; -ones(size(rq1D))];
    sfq2 = [-ones(size(rq1D)); rq1D; -rq1D];
    wfq2 = [wq1D; wq1D; wq1D];    
    Vfqf2 = kron(eye(3),Vandermonde1D(N,rq1D)/Vandermonde1D(N,JacobiGL(0,0,N)));
    Vfq2 = Vandermonde2D(N,rfq2,sfq2)/V;
    xf2 = Vfq2*x;  yf2 = Vfq2*y;
    wsJ = diag(wfq)*(Vfqf*sJ);
    wsJ2 = diag(wfq2)*(Vfqf2*sJ);
    qerr(sk) = abs(sum(sum(f(xf,yf).*wsJ)) - sum(sum(f(xf2,yf2).*wsJ2)));
    
    sk = sk + 1;
end

h = 2./Kvec;
h = h/h(1);

% loglog(h,qerr,'o--','linewidth',2,'markersize',16)
% hold on
% loglog(h,qerr(1)*h.^(2*N-1),'--','linewidth',2)
% loglog(h,qerr(1)*h.^(2*N),'--','linewidth',2)
% loglog(h,qerr(1)*h.^(2*N+2),'--','linewidth',2)
% loglog(h,qerr(1)*h.^(2*N+4),'--','linewidth',2)
% text(h(end),.75*qerr(1)*h(end).^(2*N-1),'2N-1','fontsize',16)
% text(h(end),.75*qerr(1)*h(end).^(2*N),'2N','fontsize',16)
% text(h(end),.75*qerr(1)*h(end).^(2*N+2),'2N+2','fontsize',16)
% text(h(end),.75*qerr(1)*h(end).^(2*N+4),'2N+4','fontsize',16)
% return

loglog(h,errstrong,'o--','linewidth',2,'markersize',16)
hold on
loglog(h,errweak,'^--','linewidth',2,'markersize',16)
loglog(h,err,'x--','linewidth',2,'markersize',16)

legend('Strong','Weak','Skew','Location','Best')
% loglog(h,err(1)*h.^(N+3),'--','linewidth',2)
loglog(h,err(1)*h.^(N+2),'--','linewidth',2)
loglog(h,err(1)*h.^(N+1),'--','linewidth',2)
loglog(h,err(1)*h.^(N),':','linewidth',2)


% loglog(h,err(1)*h.^(N-1),'-.','linewidth',2)
% text(h(end),.75*err(1)*h(end).^(N+3),'N+3','fontsize',16)
text(h(end),.75*err(1)*h(end).^(N+2),'N+2','fontsize',16)
text(h(end),.75*err(1)*h(end).^(N+1),'N+1','fontsize',16)
text(h(end),.75*err(1)*h(end).^(N),'N','fontsize',16)
% text(h(end),.75*err(1)*h(end).^(N-1),'N-1','fontsize',16)

%% tets

clear

useQuads = 0; mypath;

Globals3D

N = 2;

sk = 1;
Kvec = 2:5;
for kk = Kvec
    [Nv, VX, VY, VZ, K, EToV] = MeshReaderGmsh3D(sprintf('cube%d.msh',kk));
    
%     a = .125;
%     VX = VX + a*sin(pi*VX).*sin(pi*VY).*sin(pi*VZ);
%     VY = VY + a*sin(pi*VX).*sin(pi*VY).*sin(pi*VZ);
%     VZ = VZ + a*sin(pi*VX).*sin(pi*VY).*sin(pi*VZ);
    
    StartUp3D;
        
    [rq sq tq wq] = tet_cubature(2*N+2); % integrate u*v*c
    [rqtri sqtri wqtri] = Cubature2D(2*N+3);
    
    e = ones(size(rqtri));
    rfq = [rqtri;  -e; rqtri; -(1+rqtri+sqtri)];
    sfq = [sqtri; rqtri;  -e; rqtri];
    tfq = [-e;  sqtri; sqtri; sqtri];
    wfq = [wqtri;wqtri;wqtri;wqtri];
    
    Vq = Vandermonde3D(N,rq,sq,tq)/V;
    M = Vq'*diag(wq)*Vq;
    Pq = M\(Vq'*diag(wq)); % J's cancel out
    
    VqPq = Vq*Pq;
    
    xq = Vq*x; yq = Vq*y; zq = Vq*z;
    Jq = Vq*J;
    
    rxJ = rx.*J; sxJ = sx.*J; txJ = tx.*J;
            
    
    Vfq = Vandermonde3D(N,rfq,sfq,tfq)/V;
    VfPq = Vfq*Pq;
    
    rtri = r(1:(N+1)*(N+2)/2);
    stri = s(1:(N+1)*(N+2)/2);
       
    zz = zeros(size(rqtri));
    nrJ = [zz; -e; zz;  e];
    nsJ = [zz;  zz; -e; e];
    ntJ = [-e; zz; zz ; e];    
   
    
    WN = diag([wq;wfq]);
    Qr = diag(wq)*Vq*Dr*Pq;
    Qs = diag(wq)*Vq*Ds*Pq;
    Qt = diag(wq)*Vq*Dt*Pq;
    
    Br = diag(wfq.*nrJ);
    Bs = diag(wfq.*nsJ);
    Bt = diag(wfq.*ntJ);
    E = VfPq;
    QNr = [Qr-.5*E'*Br*E .5*E'*Br;
        -.5*Br*E .5*Br];
    QNs = [Qs-.5*E'*Bs*E .5*E'*Bs;
        -.5*Bs*E .5*Bs];
    QNt = [Qt-.5*E'*Bt*E .5*E'*Bt;
        -.5*Bt*E .5*Bt];
    
    QNrskew = .5*[Qr-Qr' VfPq'*diag(wfq.*nrJ);
        -diag(wfq.*nrJ)*VfPq diag(0*wfq)];
    QNsskew = .5*[Qs-Qs' VfPq'*diag(wfq.*nsJ);
        -diag(wfq.*nsJ)*VfPq diag(0*wfq)];
    QNtskew = .5*[Qt-Qt' VfPq'*diag(wfq.*ntJ);
        -diag(wfq.*ntJ)*VfPq diag(0*wfq)];
    QNrskew(abs(QNrskew)<1e-8) = 0; QNrskew = sparse(QNrskew);
    QNsskew(abs(QNsskew)<1e-8) = 0; QNsskew = sparse(QNsskew);        
    
    BNr = diag([zeros(size(Qr,1),1); nrJ]);
    BNs = diag([zeros(size(Qr,1),1); nsJ]);
    BNt = diag([zeros(size(Qr,1),1); ntJ]);
    
    VqPN = Vq*(M\[Vq' Vfq']);
    
    f = @(x,y,z) exp(sin(x+y+z));
    dfex = @(x,y,z) exp(sin(x + y + z)).*cos(x + y + z);
    
    VN = [Vq;Vfq];
    rxJN = VN*rxJ;
    sxJN = VN*sxJ;
    txJN = VN*txJ;    
    fN = f(VN*x,VN*y,VN*z);
        
    % skew deriv
    dfJ = VqPN*(rxJN.*((QNrskew + .5*WN*BNr)*fN) + sxJN.*((QNsskew + .5*WN*BNs)*fN) + txJN.*((QNtskew + .5*WN*BNt)*fN));
    df = dfJ./(Vq*J);
    
    % strong/weak deriv
    dfJ = VqPN*(rxJN.*((QNr)*fN) + sxJN.*((QNs)*fN) + txJN.*((QNt)*fN));
    dfstrong = dfJ./(Vq*J);
    dfJ = VqPN*(rxJN.*((WN*BNr-QNr')*fN) + sxJN.*((WN*BNs-QNs')*fN) + txJN.*((WN*BNt-QNt')*fN));
    dfweak = dfJ./(Vq*J);
    
    [rq2 sq2 tq2 wq2] = tet_cubature(3*N+2); % integrate u*v*c
    Vq2 = Vandermonde3D(N,rq2,sq2,tq2)/V;
    Pq2 = (Vq2'*diag(wq2)*Vq2)\(Vq2'*diag(wq2));
    xq2 = Vq2*x;
    yq2 = Vq2*y;    
    zq2 = Vq2*z;    
    wJq2 = diag(wq2)*(Vq2*J);
    err(sk) = sqrt(sum(sum(wJq2.*(Vq2*Pq2*dfex(xq2,yq2,zq2)-Vq2*Pq*df).^2)));
    errstrong(sk) = sqrt(sum(sum(wJq2.*(Vq2*Pq2*dfex(xq2,yq2,zq2)-Vq2*Pq*dfstrong).^2)));
    errweak(sk) = sqrt(sum(sum(wJq2.*(Vq2*Pq2*dfex(xq2,yq2,zq2)-Vq2*Pq*dfweak).^2)));
    
    sk = sk + 1;
end

h = 2.^(-Kvec);
h = h/h(1);
loglog(h,errstrong,'o--','linewidth',2,'markersize',16)
hold on
loglog(h,errweak,'^--','linewidth',2,'markersize',16)
loglog(h,err,'x--','linewidth',2,'markersize',16)
legend('Strong','Weak','Skew','Location','Best')

a = 10;
loglog(h,a*err(1)*h.^(N+3),'--','linewidth',2)
loglog(h,a*err(1)*h.^(N+2),'--','linewidth',2)
loglog(h,a*err(1)*h.^(N+1),'--','linewidth',2)
loglog(h,a*err(1)*h.^(N),':','linewidth',2)
loglog(h,a*err(1)*h.^(N-1),'-.','linewidth',2)
text(h(end),a*err(1)*h(end).^(N+3),'N+3','fontsize',16)
text(h(end),a*err(1)*h(end).^(N+2),'N+2','fontsize',16)
text(h(end),a*err(1)*h(end).^(N+1),'N+1','fontsize',16)
text(h(end),a*err(1)*h(end).^(N),'N','fontsize',16)
text(h(end),a*err(1)*h(end).^(N-1),'N-1','fontsize',16)


%% quads

clear
useQuads = 1; mypath;
Globals2D

N = 3;
Kvec = [4 8 16 32 64 128];

[rq1D_vol wq1D_vol] = JacobiGL(0,0,N);
[rq1D_face wq1D_face] = JacobiGQ(0,0,N); 
% [rq1D_vol wq1D_vol] = clenshaw_curtis(2*N+2,2);
% [rq1D_face wq1D_face] = clenshaw_curtis(2*N+1,3);

% volume nodes
rq1D = rq1D_vol; wq1D = wq1D_vol;
[rq sq] = meshgrid(rq1D); 
rq = rq(:); sq = sq(:);
[wrq wsq] = meshgrid(wq1D);
wq = wrq(:).*wsq(:);

% face nodes
rq1D = rq1D_face;
wq1D = wq1D_face;
e = ones(size(rq1D));
rfq = [rq1D; e; rq1D; -e];
sfq = [-e; rq1D; e; rq1D];
wfq = [wq1D; wq1D; wq1D; wq1D];

rfq = [rq1D; e; rq1D; -e];
sfq = [-e; rq1D; e; rq1D];
wfq = [wq1D; wq1D; wq1D; wq1D];

sk = 1;
for K1D = Kvec
    
    [Nv, VX, VY, K, EToV] = QuadMesh2D(K1D);
    a = 0*.25;
    VX = VX + a*sin(pi*VX).*sin(pi*VY);
    VY = VY + a*sin(pi*VX).*sin(pi*VY);
    
%     if K1D==8
%         clf
%         plot(VX,VY,'o')
%         return
%     end
    
    % interp nodes
    [r1D w1D] = JacobiGL(0,0,N);
    [r s] = meshgrid(r1D); r = r(:); s = s(:);
    
    % map nodes
    r1 = [-1 1 1 -1]';
    s1 = [-1 -1 1 1]';
    V1 = Vandermonde2DQuad(1,r,s)/Vandermonde2DQuad(1,r1,s1);
    x = V1*VX(EToV)';
    y = V1*VY(EToV)';
    
    V1D = Vandermonde1D(N,r1D);
    D1D = GradVandermonde1D(N,r1D)/V1D;
    V = Vandermonde2DQuad(N,r,s);
    Dr = kron(D1D,eye(N+1));
    Ds = kron(eye(N+1),D1D);        
    
    Vq = Vandermonde2DQuad(N,rq,sq)/V;
    xq = Vq*x;
    yq = Vq*y;
    
    [rx,sx,ry,sy,J] = GeometricFactors2D(x,y,Dr,Ds);
    rxJ = rx.*J;
    sxJ = sx.*J;
    ryJ = ry.*J;
    syJ = sy.*J;
        
    V1D = Vandermonde1D(N,r1D);
    Vf = Vandermonde2DQuad(N,rfq,sfq)/V;
    
    xf = Vf*x;
    yf = Vf*y;
    
    nrJ = [0*e; e; 0*e; -e]; % sJ = 2 for all faces,
    nsJ = [-e; 0*e; e; 0*e];
    
    nxJ = diag(nrJ)*(Vf*rxJ) + diag(nsJ)*(Vf*sxJ);
    nyJ = diag(nrJ)*(Vf*ryJ) + diag(nsJ)*(Vf*syJ);
    sJ = sqrt(nxJ.^2 + nyJ.^2);
    
    % quadrature operators
    M = Vq'*diag(wq)*Vq;
    Pq = M\(Vq'*diag(wq));
    Lq = M\(Vf'*diag(wfq));
    VfPq = (Vf*Pq);
    VqPq = Vq*Pq;
    VqLq = Vq*Lq;
    
    WN = diag([wq;wfq]);
    
    Qr = diag(wq)*Vq*Dr*Pq;
    Qs = diag(wq)*Vq*Ds*Pq;
    QNrskew = .5*[Qr-Qr' VfPq'*diag(wfq.*nrJ);
        -diag(wfq.*nrJ)*VfPq diag(0*wfq)];
    QNsskew = .5*[Qs-Qs' VfPq'*diag(wfq.*nsJ);
        -diag(wfq.*nsJ)*VfPq diag(0*wfq)];
    
    E = VfPq;
    Br = diag(wfq.*nrJ);
    Bs = diag(wfq.*nsJ);
    QNr = [Qr-.5*E'*Br*E .5*E'*Br;
        -.5*Br*E .5*Br];
    QNs = [Qs-.5*E'*Bs*E .5*E'*Bs;
        -.5*Bs*E .5*Bs];
    
    QNrskew(abs(QNrskew)<1e-8) = 0; QNrskew = sparse(QNrskew);
    QNsskew(abs(QNsskew)<1e-8) = 0; QNsskew = sparse(QNsskew);
    
    % make skew symmetric diff matrices
    DNr = diag(1./[wq;wfq])*QNrskew;
    DNs = diag(1./[wq;wfq])*QNsskew;
    
    BNr = diag([zeros(size(Qr,1),1); nrJ]);
    BNs = diag([zeros(size(Qr,1),1); nsJ]);
        
    VqPN = Vq*(M\[Vq' Vf']);
    PN = (M\[Vq' Vf']);
    
    f = @(x,y) exp(sin(2*x+y)).*sin(1.5*pi*x).*sin(2*y);
    dfex = @(x,y) (3.*pi.*sin(2.*y).*exp(sin(2.*x + y)).*cos((3.*pi.*x)/2))/2 + 2.*cos(2.*x + y).*sin(2.*y).*exp(sin(2.*x + y)).*sin((3.*pi.*x)/2);
    
    VN = [Vq;Vf];
    rxJN = VN*rxJ;
    sxJN = VN*sxJ;
    %     fN = VN*Pq*f(xq,yq);
    fN = f(VN*x,VN*y);
    
    % for later testing
    rN = VN*r; sN = VN*s;
    WBNr = diag([zeros(size(Qr,1),1); nrJ.*wfq]);
    
    
    % skew deriv
    dfJ = VqPN*(rxJN.*((QNrskew + .5*WN*BNr)*fN) + sxJN.*((QNsskew + .5*WN*BNs)*fN));
    df = dfJ./(Vq*J);
    
    % strong/weak deriv
    dfJ = VqPN*(rxJN.*((QNr)*fN) + sxJN.*((QNs)*fN));
    dfstrong = dfJ./(Vq*J);
    dfJ = VqPN*(rxJN.*((WN*BNr-QNr')*fN) + sxJN.*((WN*BNs-QNs')*fN));
    dfweak = dfJ./(Vq*J);
    
    [rq1D wq1D] = JacobiGQ(0,0,3*N+3);
%     [rq1D wq1D] = JacobiGL(0,0,N);
    [rq2 sq2] = meshgrid(rq1D); rq2 = rq2(:); sq2 = sq2(:);
    [wrq wsq] = meshgrid(wq1D);
    wq2 = wrq(:).*wsq(:);
    Vq2 = Vandermonde2DQuad(N,rq2,sq2)/V;
    xq2 = Vq2*x;
    yq2 = Vq2*y;
    wJq2 = diag(wq2)*(Vq2*J);
    err(sk) = sqrt(sum(sum(wJq2.*(dfex(xq2,yq2)-Vq2*Pq*df).^2)));
    %errstrong(sk) = sqrt(sum(sum(wJq2.*(dfex(xq2,yq2)-Vq2*Pq*dfstrong).^2)));
    %errweak(sk) = sqrt(sum(sum(wJq2.*(dfex(xq2,yq2)-Vq2*Pq*dfweak).^2)));
    errstrong(sk) = sqrt(sum(sum(wJq2.*(Vq2*Pq*dfex(xq,yq)-Vq2*Pq*dfstrong).^2)));
    errweak(sk) = sqrt(sum(sum(wJq2.*(Vq2*Pq*dfex(xq,yq)-Vq2*Pq*dfweak).^2)));
        
    % test quad err
    int_ex = integral(@(x) 0*f(x,-1),-1,1,'AbsTol',1e-14) ...
        + integral(@(x) 0*f(x,1),-1,1,'AbsTol',1e-14) ...
        + integral(@(y) -1*f(-1,y),-1,1,'AbsTol',1e-14) ...
        + integral(@(y) 1*f(1,y),-1,1,'AbsTol',1e-14);
    wsJ = diag(wfq)*sJ;
    qerr(sk) = abs(sum(sum(diag(wfq)*(nxJ.*(f(xf,yf))))) - int_ex);
    
    sk = sk + 1;
end

h = 2./Kvec;
h = h/h(1);

% clf
% loglog(h,qerr,'*--','linewidth',2,'markersize',16); 
% hold on; 
% loglog(h,qerr(1)*h.^(2*N-1),'--')
% text(h(end),.75*qerr(1)*h(end).^(2*N-1),'2N-1','fontsize',16)
% loglog(h,qerr(1)*h.^(2*N+1),'--')
% text(h(end),.75*qerr(1)*h(end).^(2*N+1),'2N+1','fontsize',16)
% loglog(h,qerr(1)*h.^(2*N+2),'--')
% text(h(end),.75*qerr(1)*h(end).^(2*N+2),'2N+2','fontsize',16)
% return

loglog(h,errstrong,'o--','linewidth',2,'markersize',16)
hold on
loglog(h,errweak,'^--','linewidth',2,'markersize',16)
% loglog(h,err,'x--','linewidth',2,'markersize',16)
legend('Strong','Weak','Skew','Location','Best')
% loglog(h,errweak(1)*h.^(N+2),'-.','linewidth',2)
loglog(h,errweak(1)*h.^(N+1),'--','linewidth',2)
loglog(h,errweak(1)*h.^(N),':','linewidth',2)
loglog(h,errweak(1)*h.^(N-1),'-.','linewidth',2)

text(h(end),.75*errweak(1)*h(end).^(N-1),'N-1','fontsize',16)
text(h(end),.75*errweak(1)*h(end).^(N),'N','fontsize',16)
text(h(end),.75*errweak(1)*h(end).^(N+1),'N+1','fontsize',16)
% text(h(end),.75*errweak(1)*h(end).^(N+2),'N+2','fontsize',16)