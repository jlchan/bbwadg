function Wave3D_IGA_DG

Globals3D

NB = 3;
Ksub = 8;
smoothKnots = 25;

computeEigs = 0;
Testing = 0;

N = NB+Ksub-1;

Np = (N+1)^3;
Nfp = (N+1)^2;

FinalTime = 2;

%% reference elem stuff

global rxJ sxJ txJ ryJ syJ tyJ rzJ szJ tzJ
global nx ny nz Jq sJ
global DVq1D Vq1D Pq1D Prq1D D1D D1Dt LIFT Fmask
global mapM mapP mapB Nfaces

Nfaces = 6;
if Ksub==1
    r1D = JacobiGL(0,0,N);
else
    r1D = JacobiGL(0,0,N);
    [~, ~, ~, ~, ~, ~, ~, ~, VX] = bsplineVDM(NB,Ksub,r1D,smoothKnots); % VDM for interp, mass, M\S
    t = [VX(1)*ones(1,NB) VX VX(end)*ones(1,NB)]; % open knot vec
    for i = 1:N+1
        r1D(i,1) = mean(t((i+1):(i+NB))); % greville
    end
    %     plot(r1D,r1D*0,'o')
end

[r s t] = meshgrid(r1D); 
r = r(:); s = s(:); t = t(:); 

Fmask1 = find(abs(r+1)<1e-8); 
Fmask2 = find(abs(r-1)<1e-8); 
Fmask3 = find(abs(s+1)<1e-8); 
Fmask4 = find(abs(s-1)<1e-8); 
Fmask5 = find(abs(t+1)<1e-8); 
Fmask6 = find(abs(t-1)<1e-8); 
Fmask = [Fmask1(:) Fmask2(:) Fmask3(:) Fmask4(:) Fmask5(:) Fmask6(:)];

%rp1D = linspace(-1,1,50);
rp1D = linspace(-1,1,50);
[rp sp tp] = meshgrid(rp1D);
rp = rp(:); sp = sp(:); tp = tp(:);

% plotting fmask
Fmask1 = find(abs(rp+1)<1e-8); 
Fmask2 = find(abs(rp-1)<1e-8); 
Fmask3 = find(abs(sp+1)<1e-8); 
Fmask4 = find(abs(sp-1)<1e-8); 
Fmask5 = find(abs(tp+1)<1e-8); 
Fmask6 = find(abs(tp-1)<1e-8); 
Fmaskp = [Fmask1(:) Fmask2(:) Fmask3(:) Fmask4(:) Fmask5(:) Fmask6(:)];

% initialize TP operators
[BVDM M1D D1D R rq1D wq1D Vq1D DVq1D] = bsplineVDM(NB,Ksub,r1D,smoothKnots); % VDM for interp, mass, M\S
Pq1D = M1D\(Vq1D'*diag(wq1D));
Prq1D = M1D\(DVq1D'*diag(wq1D));

if Ksub >= NB
    Vq1D = sparse(Vq1D);
    Pq1D = sparse(Pq1D);
    Prq1D = sparse(Prq1D);
end

(N+1)
Nq = size(Vq1D,1)

D1Dt = M1D\(D1D'*M1D); % for skew-symmetric formulation

[rq sq tq] = meshgrid(rq1D); rq = rq(:); sq = sq(:); tq = tq(:);
[wrq wsq wtq] = meshgrid(wq1D); wq = wrq(:).*wsq(:).*wtq(:);

% face
[rq2 sq2] = meshgrid(rq1D); rq2 = rq2(:); sq2 = sq2(:); 
[wr2 ws2] = meshgrid(wq1D); wfq = wr2(:).*ws2(:); wfq = repmat(wfq,Nfaces,1);

% make face quadrature nodes
rfq = zeros(length(rq1D)^2,Nfaces);
sfq = zeros(length(rq1D)^2,Nfaces);
tfq = zeros(length(rq1D)^2,Nfaces);

% r = -/+ 1q
rfq(:,1) = -1; sfq(:,1) = rq2; tfq(:,1) = sq2; 
rfq(:,2) = 1;  sfq(:,2) = rq2; tfq(:,2) = sq2; 
% s = -/+ 1
rfq(:,3) = rq2; sfq(:,3) = -1; tfq(:,3) = sq2; 
rfq(:,4) = rq2; sfq(:,4) = 1;  tfq(:,4) = sq2; 
% t = -/+ 1
rfq(:,5) = rq2; sfq(:,5) = sq2; tfq(:,5) = -1; 
rfq(:,6) = rq2; sfq(:,6) = sq2; tfq(:,6) =  1; 

rfq = rfq(:); sfq = sfq(:); tfq = tfq(:);

% plot3(rfq,sfq,tfq,'o'); return 

% invM1D = inv(M1D);
Mf1D = zeros(N+1,2);
Mf1D(1,1) = 1;
Mf1D(N+1,2) = 1;
LIFT1D = M1D\Mf1D;

% make lift operator
ids = 1:(N+1)^3;
ids = reshape(ids,N+1,N+1,N+1);

LIFT = sparse((N+1)^3,6*(N+1)^2); % initialize as sparse matrix

% first/second face (r = +/- 1)
sk = 0;
for i = 1:(N+1)
    iid = squeeze(ids(i,:,:));
    for j = 1:(N+1)
        LIFT(iid(:,j),sk + (j-1)*(N+1) + i) = LIFT1D(:,1); % left-to-right
        LIFT(iid(:,j),sk + (N+1)^2 + (j-1)*(N+1) + i) = LIFT1D(:,2); % right-to-left
    end
end

% third/fourth face (s = +/- 1)
sk = sk + 2*(N+1)^2;
for i = 1:(N+1)
    iid = squeeze(ids(:,:,i));
    for j = 1:(N+1)
        LIFT(iid(:,j),sk + (i-1)*(N+1) + j) = LIFT1D(:,1); % left-to-right
        LIFT(iid(:,j),sk + (N+1)^2 + (i-1)*(N+1) + j) = LIFT1D(:,2); % right-to-left
    end
end

% fifth/sixth face (t = +/- 1)
sk = sk + 2*(N+1)^2;
for i = 1:(N+1)
    iid = squeeze(ids(:,i,:));
    for j = 1:(N+1)
        LIFT(iid(j,:),sk + (i-1)*(N+1) + j) = LIFT1D(:,1); % left-to-right
        LIFT(iid(j,:),sk + (N+1)^2 + (i-1)*(N+1) + j) = LIFT1D(:,2); % right-to-left
    end
end

% M = kron(M1D,kron(M1D,M1D));
% imagesc(M*LIFT)
% keyboard

%% make mesh and connectivity

addpath IGA_DG
addpath IGA_DG/Elbow3D

%% recover vertices and connectivity

[re se te] = meshgrid([-1;1]); 
re = re(:); se = se(:); te = te(:);

%     [vx vy vz] = Elbow3Dgeo(re,se,te);
%     [x y z geo] = Elbow3Dgeo(r,s,t);
%     [xq yq zq geo] = Elbow3Dgeo(rq,sq,tq);
%     [xfq yfq zfq fgeo] = Elbow3Dgeo(rfq,sfq,tfq);

[vx vy vz] = WarpedCube3Dgeo(re,se,te);
[x y z geo] = WarpedCube3Dgeo(r,s,t);
[xq yq zq geo] = WarpedCube3Dgeo(rq,sq,tq);
[xfq yfq zfq fgeo] = WarpedCube3Dgeo(rfq,sfq,tfq);

% plot3(x,y,z,'o'); axis equal; return

% make unique vertex list
[VXYZ idr idc] = uniquetol([vx(:) vy(:) vz(:)],1e-7,'ByRows',true);
VX = VXYZ(:,1); VY = VXYZ(:,2); VZ = VXYZ(:,3);
K = size(vx,2);
EToV = reshape(idc,8,K)';

% f1 = [1 5 7 3];
% f2 = [2 6 8 4];
% 
% f3 = [1 2 4 3];
% f4 = [5 6 8 7];
% 
% f5 = [1 2 6 5];
% f6 = [3 4 8 7];

% e = 1; f = f1;
% f = 1:8;
% plot3(VX(EToV(e,f)),VY(EToV(e,f)),VZ(EToV(e,f)),'o')
% text(VX(EToV(e,f))+.1,VY(EToV(e,f)),VZ(EToV(e,f)),num2str((1:4)'))

[EToE,EToF] = tiConnect3DHex(EToV);
% keyboard


%% get greville nodes for connectivity

K = size(x,2);

xf = reshape(x(Fmask(:),:),Nfp*Nfaces,K);
yf = reshape(y(Fmask(:),:),Nfp*Nfaces,K);
zf = reshape(z(Fmask(:),:),Nfp*Nfaces,K);

mapM = reshape(1:Nfp*Nfaces*K,Nfp*Nfaces,K);
mapP = reshape(1:Nfp*Nfaces*K,Nfp*Nfaces,K);

NODETOL = 1e-6;
for e = 1:K
    for f = 1:Nfaces
        enbr = EToE(e,f);
        fnbr = EToF(e,f);
        idM = (1:Nfp) + (f-1)*Nfp;
        idP = (1:Nfp) + (fnbr-1)*Nfp;
        
        if e==enbr
            mapP(idM,e) = mapM(idM,e);
        else
            xM = xf(idM,e); yM = yf(idM,e); zM = zf(idM,e);
            xP = xf(idP,enbr); yP = yf(idP,enbr); zP = zf(idP,enbr);
            
%             keyboard
            
            % find find volume node numbers of left and right nodes
            x1 = xf(idM,e); y1 = yf(idM,e); z1 = zf(idM,e);
            x2 = xf(idP,enbr); y2 = yf(idP,enbr); z2 = zf(idP,enbr);
            
            one = ones(1, Nfp);
            x1 = x1*one;  y1 = y1*one;  z1 = z1*one;
            x2 = x2*one;  y2 = y2*one;  z2 = z2*one;
            
            % Compute distance matrix
            D = (x1 -x2').^2 + (y1-y2').^2 + (z1-z2').^2;
            [idM, idP] = find(sqrt(abs(D)) < NODETOL);
            
            mapP(idM + (f-1)*Nfp,e) = idP + (fnbr-1)*Nfp + (enbr-1)*Nfaces*Nfp;
        end
    end
end

mapB = find(mapM==mapP);

vmapM = reshape(1:Np*K,Np,K);
vmapM = vmapM(Fmask(:),:);
vmapP = vmapM(mapP);
%     keyboard

err1 = norm(xf(mapM) - xf(mapP),'fro')+norm(yf(mapM) - yf(mapP),'fro')+norm(zf(mapM) - zf(mapP),'fro');
err2 = norm(x(vmapM) - x(vmapP),'fro')+norm(y(vmapM) - y(vmapP),'fro')+norm(z(vmapM) - z(vmapP),'fro');
if err1 + err2 > 1e-10
    keyboard
end
% keyboard


%% get quadrature nodes

rxJ = geo.rxJ; sxJ = geo.sxJ; txJ = geo.txJ;
ryJ = geo.ryJ; syJ = geo.syJ; tyJ = geo.tyJ;
rzJ = geo.rzJ; szJ = geo.szJ; tzJ = geo.tzJ;
Jq  = geo.J;

%% surface quadrature nodes and normals

Nfq = length(rfq)/Nfaces;

%     plot3(xfq,yfq,zfq,'o');  axis equal;   return
rxJf = fgeo.rxJ; sxJf = fgeo.sxJ; txJf = fgeo.txJ;
ryJf = fgeo.ryJ; syJf = fgeo.syJ; tyJf = fgeo.tyJ;
rzJf = fgeo.rzJ; szJf = fgeo.szJ; tzJf = fgeo.tzJ;
Jf   = fgeo.J;

fids = @(f) (1:Nfq) + (f-1)*Nfq;
nr = zeros(Nfq*Nfaces,1); 
ns = zeros(Nfq*Nfaces,1); 
nt = zeros(Nfq*Nfaces,1);
f = 1; nr(fids(f)) = -1;
f = 2; nr(fids(f)) = 1;
f = 3; ns(fids(f)) = -1;
f = 4; ns(fids(f)) = 1;
f = 5; nt(fids(f)) = -1;
f = 6; nt(fids(f)) = 1;

nr = repmat(nr,1,K);  
ns = repmat(ns,1,K);  
nt = repmat(nt,1,K);

nxJ = rxJf.*nr + sxJf.*ns + txJf.*nt;
nyJ = ryJf.*nr + syJf.*ns + tyJf.*nt;
nzJ = rzJf.*nr + szJf.*ns + tzJf.*nt;

nx = nxJ./Jf; ny = nyJ./Jf; nz = nzJ./Jf;
sJ = sqrt(nx.^2 + ny.^2 + nz.^2);
nx = nx./sJ; ny = ny./sJ; nz = nz./sJ;
sJ = sJ.*Jf;



% keyboard

if 0
    fids = 1:Nfq*Nfaces;
    for e = 1:K;
        clf
        plot3(xfq(fids,e),yfq(fids,e),zfq(fids,e),'o');
        hold on
        quiver3(xfq(fids,e),yfq(fids,e),zfq(fids,e),nx(fids,e),ny(fids,e),nz(fids,e))
        axis equal
        pause
    end
    keyboard
end
if min(Jq(:)) < 1e-8
    keyboard
end

wJq = spdiag(wq)*Jq;


%% for testing - assemble 3D global operators

if Testing
    global Vq Vrq Vsq Vtq Prq Psq Ptq Pq Vfq Pfq
    
%     V1D = Vandermonde1D(N,r1D);    
%     V = kron(kron(V1D,V1D),V1D);                
    
    Vq = kron(kron(Vq1D,Vq1D),Vq1D);
    Dr = kron(kron(eye(N+1),D1D),eye(N+1));
    Ds = kron(kron(eye(N+1),eye(N+1)),D1D);
    Dt = kron(kron(D1D,eye(N+1)),eye(N+1));
    
    Vrq = Vq*Dr; Vsq = Vq*Ds; Vtq = Vq*Dt;
    invM1D = inv(Vq1D'*diag(wq1D)*Vq1D);
    invM = kron(kron(invM1D,invM1D),invM1D);
    
    Pq = invM*(Vq'*diag(wq));
    Prq = invM*(Vrq'*diag(wq));
    Psq = invM*(Vsq'*diag(wq));
    Ptq = invM*(Vtq'*diag(wq));
    
    Vfq = kron(eye(6),kron(Vq1D,Vq1D));
    Pfq = kron(eye(6),kron(Pq1D,Pq1D));
    
    % check metric identities 
    norm(Vrq*(Pq*rxJ) + Vsq*(Pq*sxJ) + Vtq*(Pq*txJ),'fro')
    norm(Vrq*(Pq*ryJ) + Vsq*(Pq*syJ) + Vtq*(Pq*tyJ),'fro')
    norm(Vrq*(Pq*rzJ) + Vsq*(Pq*szJ) + Vtq*(Pq*tzJ),'fro')
            
    cc = Pq*xq; 
    norm((rxJ.*(Vrq*cc) + sxJ.*(Vsq*cc) + txJ.*(Vtq*cc))./Jq-1,'fro')
    norm((ryJ.*(Vrq*cc) + syJ.*(Vsq*cc) + tyJ.*(Vtq*cc))./Jq,'fro')
    norm((rzJ.*(Vrq*cc) + szJ.*(Vsq*cc) + tzJ.*(Vtq*cc))./Jq,'fro')
    
    cc = Pq*yq; 
    norm((rxJ.*(Vrq*cc) + sxJ.*(Vsq*cc) + txJ.*(Vtq*cc))./Jq,'fro')
    norm((ryJ.*(Vrq*cc) + syJ.*(Vsq*cc) + tyJ.*(Vtq*cc))./Jq-1,'fro')
    norm((rzJ.*(Vrq*cc) + szJ.*(Vsq*cc) + tzJ.*(Vtq*cc))./Jq,'fro')
    
    cc = Pq*zq; 
    norm((rxJ.*(Vrq*cc) + sxJ.*(Vsq*cc) + txJ.*(Vtq*cc))./Jq,'fro')
    norm((ryJ.*(Vrq*cc) + syJ.*(Vsq*cc) + tyJ.*(Vtq*cc))./Jq,'fro')
    norm((rzJ.*(Vrq*cc) + szJ.*(Vsq*cc) + tzJ.*(Vtq*cc))./Jq-1,'fro')
    
%     cc = Pq*(xq.^N+yq+zq);
%     norm(rxJ.*(Vrq*cc) + sxJ.*(Vsq*cc) + txJ.*(Vtq*cc) - N*xq.^(N-1).*Jq,'fro')
%     keyboard
    
    % check IBP    
    if 0
        e = 1;
        Vqf = zeros(6*Nq^2,Np);
        uu = zeros(Np,1);
        for i = 1:Np
            uu(i) = 1;
            ufq = zeros(Nq^2,6);
            for f = 1:6
                %tmp = Vq1D*reshape(uu(Fmask(:,f),:),N+1,N+1)*Vq1D';
                tmp = kron(Vq1D,Vq1D)*uu(Fmask(:,f),:);
                ufq(:,f) = tmp(:);
            end
            Vqf(:,i) = ufq(:);
            uu(i) = 0;
        end
        Vqf(abs(Vqf)<1e-8) = 0;
        
        Vq1Dnodal = Vandermonde1D(N,rq1D)/Vandermonde1D(N,r1D);
        fids = @(f) (1:Nfp) + (f-1)*Nfp;
        xfq = zeros(Nq^2,6);
        yfq = zeros(Nq^2,6);
        zfq = zeros(Nq^2,6);
        for f = 1:6
            xfq(:,f) = kron(Vq1Dnodal,Vq1Dnodal)*xf(fids(f),e);
            yfq(:,f) = kron(Vq1Dnodal,Vq1Dnodal)*yf(fids(f),e);
            zfq(:,f) = kron(Vq1Dnodal,Vq1Dnodal)*zf(fids(f),e);
        end
        xfq = xfq(:);
        yfq = yfq(:);
        zfq = zfq(:);
        norm(xfq(:)-Vqf*(Pq*xq(:,e)),'fro')
        norm(yfq(:)-Vqf*(Pq*yq(:,e)),'fro')
        norm(zfq(:)-Vqf*(Pq*zq(:,e)),'fro')
        %     keyboard
        
        u = Pq*xq(:,e);
        uwq = wq.*(Vq*u);
        uwfq = wfq.*(Vqf*u);
        dudxJ = rxJ.*(Vrq*u) + sxJ.*(Vsq*u) + txJ.*(Vtq*u);
        r1 = Vq'*(wq.*dudxJ);
        r2 = Vqf'*(uwfq.*nxJ) - (Vrq'*(rxJ.*uwq) + Vsq'*(sxJ.*uwq) + Vtq'*(txJ.*uwq));
        
        fids = @(f) (1:Nfq) + (f-1)*Nfq;
        e = 1; f = 2; plot3(xfq(fids(f),e),yfq(fids(f),e),zfq(fids(f),e),'o');
        hold on;
        quiver3(xfq(fids(f),e),yfq(fids(f),e),zfq(fids(f),e),nx(fids(f),e),ny(fids(f),e),nz(fids(f),e))
        [xfq(fids(f)),yfq(fids(f)),zfq(fids(f)) nx(fids(f)),ny(fids(f)),nz(fids(f))]
    end

        keyboard

end
% % check quadrature
% dx2q = rxJ.*(Vrq*(x.^2/2)) + sxJ.*(Vsq*(x.^2/2)) + txJ.*(Vtq*(x.^2/2));
% sum(sum(spdiag(wq)*dx2q.*yq.*zq))

% for e = 1:K
%     for f=1:Nfaces
%         if e==EToF(e,f) % if boundary
%             
%         end
%     end
% end
% keyboard

%% plotting

Vp1D_interp = Vandermonde1D(N,rp1D)/Vandermonde1D(N,r1D);
for e = 1:K
    xp(:,e) = matvec(Vp1D_interp,x(:,e),'all');
    yp(:,e) = matvec(Vp1D_interp,y(:,e),'all');
    zp(:,e) = matvec(Vp1D_interp,z(:,e),'all');        
end

% plot3(xp,yp,zp,'.'); return

Vp1D = bsplineVDM(NB,Ksub,rp1D,smoothKnots);


%% initial cond

k = 1;
pex = @(x,y,z,t) cos(k*pi*x/2).*cos(k*pi*y/2).*cos(k*pi*z/2).*cos(sqrt(3)*.5*k*pi*t);

% x0 = -1.4; y0 = -.5; z0 = 0;
% pex = @(x,y,z,t) exp(-5*((x-x0).^2 + (y-y0).^2 + (z-z0).^2));

global mapI
Nfp = size(xf,1)/Nfaces;
fids = @(f) (1:Nfp)+(f-1)*Nfp;
% vids = reshape(1:(N+1)^3*K,(N+1)^3,K);
mapI = [];
x0b = 1;
for e = 1:K
    for f = 1:6
        if e==EToE(e,f)
            xb = xf(fids(f),e);
            if norm(xb-x0b,'fro')<1e-8
                bids = fids(f) + (e-1)*Nfp*Nfaces;
                mapI = [mapI; bids(:)];                
            end            
        end
    end
end
% keyboard
% mapI = mapB(ids);

% plot3(xb,yb,zb,'o')
% hold on
% plot3(xb(mapI),yb(mapI),zb(mapI),'x')
% return

% p = matvec(inv(BVDM),pex(x,y,z,0),'all');
p = zeros(Np,K);
for e = 1:K
    p(:,e) = matvec(Pq1D,pex(xq(:,e),yq(:,e),zq(:,e),0),'all');
end
u = zeros(Np, K);
v = zeros(Np, K);
w = zeros(Np, K);

for e = 1:K
    pq(:,e) = matvec(Vq1D,p(:,e),'all');
end
err = wJq.*(pq - pex(xq,yq,zq,0)).^2;
init_cond_L2err = sqrt(sum(err(:)))

%% estimate max timestep

if 0
    U = randn((N+1)^3*K,4);
    for i = 1:10
        Uprev = U;
        [rhsp, rhsu, rhsv, rhsw] = acousticsRHS3Dq(reshape(U(:,1),Np,K),reshape(U(:,2),Np,K),reshape(U(:,3),Np,K),reshape(U(:,4),Np,K));
        U(:,1) = rhsp(:);
        U(:,2) = rhsu(:);
        U(:,3) = rhsv(:);
        U(:,4) = rhsw(:);
        
        lam = Uprev(:)'*U(:) / norm(Uprev(:));
        U = U/norm(U(:));
    end
    dt = .5/abs(lam)
else
% keyboard
    h = min(min(Jq)./max(sJ));
    if Ksub > 1
        dt = 3*h/((NB+1)*Ksub);
    else
        dt = .33*h/((NB+1)^2);
    end
end
% return


if computeEigs
    A = sparse((N+1)^3*4*K,(N+1)^3*4*K);
    U = zeros((N+1)^3*K,4);
    for i = 1:(N+1)^3*4*K
        U(i) = 1;
        if Testing
            [rhsp, rhsu, rhsv, rhsw] = acousticsRHS3Dtest(reshape(U(:,1),(N+1)^3,K),reshape(U(:,2),(N+1)^3,K),reshape(U(:,3),(N+1)^3,K),reshape(U(:,4),(N+1)^3,K));
        else
            [rhsp, rhsu, rhsv, rhsw] = acousticsRHS3Dq(reshape(U(:,1),(N+1)^3,K),reshape(U(:,2),(N+1)^3,K),reshape(U(:,3),(N+1)^3,K),reshape(U(:,4),(N+1)^3,K),0);
        end
        A(:,i) = sparse([rhsp(:);rhsu(:);rhsv(:);rhsw(:)]);
        U(i) = 0;
        if (mod(i,100)==0)
            disp(sprintf('on i = %d out of %d\n',i,(N+1)^3*4*K))
        end
    end
    lam = eig(full(A));
    plot(lam,'o')
    title(sprintf('max real lam = %g\n',max(real(lam))))
    keyboard
end



%%

% Runge-Kutta residual storage
resp = zeros(Np,K);
resu = zeros(Np,K);
resv = zeros(Np,K);
resw = zeros(Np,K);

Nsteps = ceil(FinalTime/dt)
dt = FinalTime/Nsteps;

% outer time step loop
vids = reshape(1:length(rp)*K,length(rp),K);
ids = vids(Fmaskp(:),:);
for tstep = 1:Nsteps
    for INTRK = 1:5
        
        timelocal = tstep*dt + rk4c(INTRK)*dt;
        
        if Testing            
            [rhsp,rhsu,rhsv,rhsw] = acousticsRHS3Dtest(p,u,v,w);
        else
            [rhsp,rhsu,rhsv,rhsw] = acousticsRHS3Dq(p,u,v,w,timelocal);
        end
        
        % initiate and increment Runge-Kutta residuals
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        resw = rk4a(INTRK)*resw + dt*rhsw;
        
        % update fields
        u = u+rk4b(INTRK)*resu;
        v = v+rk4b(INTRK)*resv;
        w = w+rk4b(INTRK)*resw;
        p = p+rk4b(INTRK)*resp;
        
    end;
    
    if nargin==0 && (mod(tstep,10)==0 || tstep==Nsteps)
        for e = 1:K
            vv(:,e) = matvec(Vp1D,p(:,e),'all');
        end
        clf
        color_line3(xp(ids),yp(ids),zp(ids),vv(ids),'.');
        axis equal; axis tight
        view(3)
        colorbar
%         caxis([-1.1 1.2])
        title(sprintf('time = %f',tstep*dt))
        drawnow

    end
    
    if mod(tstep,25)==0
        disp(sprintf('on tstep %d out of %d\n',tstep,Nsteps))
    end
end

% keyboard
% return

% J = matvec(Vq1D,J,'all');
for e = 1:K
    pq(:,e) = matvec(Vq1D,p(:,e),'all');
end
wJq = spdiag(wq)*Jq;
pterr = pq - pex(xq,yq,zq,FinalTime);
err = wJq.*(pterr).^2;
L2err = sqrt(sum(err(:)));
L2err

return
% keyboard

ids = find(abs(yq) < .5);
figure
color_line3(xq(ids),yq(ids),zq(ids),pterr(ids),'.');
axis equal


% Nsteps
% L2err * Nsteps % error * steps (assume same overall work)

% work = ceil(FinalTime/dt)*(N+1)^4; % total work
% L2err * work


function [rhsp rhsu rhsv rhsw] = acousticsRHS3Dq(p,u,v,w,time)

global rxJ sxJ txJ ryJ syJ tyJ rzJ szJ tzJ
global nx ny nz Jq sJ
global DVq1D Vq1D Pq1D Prq1D D1D D1Dt LIFT Fmask
global mapM mapP mapB Nfaces

N = size(Vq1D,2)-1;
Nq = size(Vq1D,1);
Nfp = (N+1)^2; 
K = size(p,2);

pf = p(Fmask(:),:);
uf = u(Fmask(:),:);
vf = v(Fmask(:),:);
wf = w(Fmask(:),:);

dp   = pf(mapP)-pf(mapM);
uavg = uf(mapP)+uf(mapM);
vavg = vf(mapP)+vf(mapM);
wavg = wf(mapP)+wf(mapM);
du   = uf(mapP)-uf(mapM);
dv   = vf(mapP)-vf(mapM);
dw   = wf(mapP)-wf(mapM);

opt=1;
if opt==1
    % neumann in pipe
    dp(mapB) = 0;
    uavg(mapB) = 0;
    vavg(mapB) = 0;
    wavg(mapB) = 0;
    du(mapB) = -2*uf(mapB);
    dv(mapB) = -2*vf(mapB);
    dw(mapB) = -2*wf(mapB);

    % velocity pulse @ inflow
    global mapI
    t0 = .5;
    if time < t0
        % uP-uM = 2*(uB-uM),  uP+uM = 2*uB
        ut = -(1-cos(2*pi*time/(t0))).*(time < t0);
        uavg(mapI) = ut;
        du(mapI) = 2*(ut-uf(mapI));        
    end
elseif opt==2  % dirichlet
    dp(mapB) = 2*(-pf(mapB));
else % neumann
    dp(mapB) = 0;
    uavg(mapB) = 0;
    vavg(mapB) = 0;
    wavg(mapB) = 0;
    du(mapB) = -2*uf(mapB);
    dv(mapB) = -2*vf(mapB);
    dw(mapB) = -2*wf(mapB);
end

for e = 1:K       
    dpK = reshape(dp(:,e),Nfp,Nfaces);
    uavgK = reshape(uavg(:,e),Nfp,Nfaces);
    vavgK = reshape(vavg(:,e),Nfp,Nfaces);
    wavgK = reshape(wavg(:,e),Nfp,Nfaces);
    duK = reshape(du(:,e),Nfp,Nfaces);
    dvK = reshape(dv(:,e),Nfp,Nfaces);
    dwK = reshape(dw(:,e),Nfp,Nfaces);
    
    rxJK = rxJ(:,e); ryJK = ryJ(:,e); rzJK = rzJ(:,e);
    sxJK = sxJ(:,e); syJK = syJ(:,e); szJK = szJ(:,e);
    txJK = txJ(:,e); tyJK = tyJ(:,e); tzJK = tzJ(:,e);        
    
    Nfq = Nq^2;
    nxK = reshape(nx(:,e),Nfq,Nfaces);
    nyK = reshape(ny(:,e),Nfq,Nfaces);
    nzK = reshape(nz(:,e),Nfq,Nfaces);
    sJK = reshape(sJ(:,e),Nfq,Nfaces);    
    
    pe = p(:,e);
    pr = matvec(Vq1D,matvec(D1D,pe,1),'all');
    ps = matvec(Vq1D,matvec(D1D,pe,2),'all');
    pt = matvec(Vq1D,matvec(D1D,pe,3),'all');
    
    uq = matvec(Vq1D,u(:,e),'all');
    vq = matvec(Vq1D,v(:,e),'all');
    wq = matvec(Vq1D,w(:,e),'all');
    
    dpdx = matvec(Pq1D,rxJK.*pr + sxJK.*ps + txJK.*pt,'all');
    dpdy = matvec(Pq1D,ryJK.*pr + syJK.*ps + tyJK.*pt,'all');
    dpdz = matvec(Pq1D,rzJK.*pr + szJK.*ps + tzJK.*pt,'all');
    
    ur = matvec(Pq1D,rxJK.*uq + ryJK.*vq + rzJK.*wq,'all');
    ur = matvec(D1Dt,ur,1);
    us = matvec(Pq1D,sxJK.*uq + syJK.*vq + szJK.*wq,'all');
    us = matvec(D1Dt,us,2);
    ut = matvec(Pq1D,txJK.*uq + tyJK.*vq + tzJK.*wq,'all');
    ut = matvec(D1Dt,ut,3);   

    divU = -(ur + us + ut);
    
    % Define field differences at faces
    tau = 1;
    fluxp = zeros((N+1)^2,Nfaces);
    fluxu = zeros((N+1)^2,Nfaces);
    fluxv = zeros((N+1)^2,Nfaces);
    fluxw = zeros((N+1)^2,Nfaces);
    for f = 1:Nfaces
        dpfq = (Vq1D*reshape(dpK(:,f),N+1,N+1))*Vq1D';
        ufq  = (Vq1D*reshape(uavgK(:,f),N+1,N+1))*Vq1D';
        vfq  = (Vq1D*reshape(vavgK(:,f),N+1,N+1))*Vq1D';
        wfq  = (Vq1D*reshape(wavgK(:,f),N+1,N+1))*Vq1D';        
        ndotUavg = nxK(:,f).*ufq(:) + nyK(:,f).*vfq(:) + nzK(:,f).*wfq(:);
        
        ufq = (Vq1D*reshape(duK(:,f),N+1,N+1))*Vq1D';
        vfq = (Vq1D*reshape(dvK(:,f),N+1,N+1))*Vq1D';
        wfq = (Vq1D*reshape(dwK(:,f),N+1,N+1))*Vq1D';        
        ndotdU = nxK(:,f).*ufq(:) + nyK(:,f).*vfq(:) + nzK(:,f).*wfq(:);
        
        fp = (tau*dpfq(:) - ndotUavg).*sJK(:,f);
        fU = (tau*ndotdU(:) - dpfq(:)).*sJK(:,f);
        
        fluxp(:,f) = reshape((Pq1D*reshape(fp,Nq,Nq))*Pq1D',(N+1)^2,1);
        fluxu(:,f) = reshape((Pq1D*reshape(fU.*nxK(:,f),Nq,Nq))*Pq1D',(N+1)^2,1);
        fluxv(:,f) = reshape((Pq1D*reshape(fU.*nyK(:,f),Nq,Nq))*Pq1D',(N+1)^2,1);
        fluxw(:,f) = reshape((Pq1D*reshape(fU.*nzK(:,f),Nq,Nq))*Pq1D',(N+1)^2,1);               
    end
    
    rp =  -divU + .5*LIFT*(fluxp(:));
    ru =  -dpdx + .5*LIFT*(fluxu(:));
    rv =  -dpdy + .5*LIFT*(fluxv(:));
    rw =  -dpdz + .5*LIFT*(fluxw(:));

%     rp =  -divU;    ru =  -dpdx;
%     rv =  -dpdy;    rw =  -dpdz;
%     
%     rp = .5*LIFT*(fluxp(:));  ru = .5*LIFT*(fluxu(:));
%     rv = .5*LIFT*(fluxv(:));  rw = .5*LIFT*(fluxw(:));
    
    rhsp(:,e) = matvec(Pq1D,matvec(Vq1D,rp,'all')./Jq(:,e),'all');
    rhsu(:,e) = matvec(Pq1D,matvec(Vq1D,ru,'all')./Jq(:,e),'all');
    rhsv(:,e) = matvec(Pq1D,matvec(Vq1D,rv,'all')./Jq(:,e),'all');
    rhsw(:,e) = matvec(Pq1D,matvec(Vq1D,rw,'all')./Jq(:,e),'all');

end

function [rhsp rhsu rhsv rhsw] = acousticsRHS3Dtest(p,u,v,w)

global rxJ sxJ txJ ryJ syJ tyJ rzJ szJ tzJ
global nx ny nz Jq sJ
global DVq1D Vq1D Pq1D Prq1D D1D D1Dt LIFT Fmask
global mapM mapP mapB Nfaces

N = size(Vq1D,2)-1;
Nq = size(Vq1D,1);
Nfp = (N+1)^2; 
K = size(p,2);

pf = p(Fmask(:),:);
uf = u(Fmask(:),:);
vf = v(Fmask(:),:);
wf = w(Fmask(:),:);

dp   = pf(mapP)-pf(mapM);
uavg = uf(mapP)+uf(mapM);
vavg = vf(mapP)+vf(mapM);
wavg = wf(mapP)+wf(mapM);
du = uf(mapP)-uf(mapM);
dv = vf(mapP)-vf(mapM);
dw = wf(mapP)-wf(mapM);

if 1  % dirichlet
    dp(mapB) = -2*pf(mapB);
else % neumann
    dp(mapB) = 0;
    uavg(mapB) = 0;
    vavg(mapB) = 0;
    wavg(mapB) = 0;    
end

global Vq Vrq Vsq Vtq Prq Psq Ptq Pq Vfq Pfq 

pr = Vrq*p; ps = Vsq*p; pt = Vtq*p;
dpdx = Pq*(rxJ.*pr + sxJ.*ps + txJ.*pt);
dpdy = Pq*(ryJ.*pr + syJ.*ps + tyJ.*pt);
dpdz = Pq*(rzJ.*pr + szJ.*ps + tzJ.*pt);

uq = Vq*u; vq = Vq*v; wq = Vq*w;
ur = Prq*(rxJ.*uq + ryJ.*vq + rzJ.*wq);
us = Psq*(sxJ.*uq + syJ.*vq + szJ.*wq);
ut = Ptq*(txJ.*uq + tyJ.*vq + tzJ.*wq);
divU = -(ur + us + ut);

dp = Vfq*dp;

uavg = Vfq*uavg; vavg = Vfq*vavg; wavg = Vfq*wavg;
ndotUavg = nx.*uavg + ny.*vavg + nz.*wavg;

rp =  -divU - .5*LIFT*(Pfq*(ndotUavg.*sJ));
ru =  -dpdx - .5*LIFT*(Pfq*(dp.*nx.*sJ));
rv =  -dpdy - .5*LIFT*(Pfq*(dp.*ny.*sJ));
rw =  -dpdz - .5*LIFT*(Pfq*(dp.*nz.*sJ));

rhsp = Pq*((Vq*rp)./Jq);
rhsu = Pq*((Vq*ru)./Jq);
rhsv = Pq*((Vq*rv)./Jq);
rhsw = Pq*((Vq*rw)./Jq);

% apply 1D tensor product matvec
function b = matvec(A,x,dim)

N = round((length(x(:)))^(1/3));
M = size(A,1);
x = reshape(x,N,N,N);

if strcmp(dim,'all')
    % apply in all dimensions    
    b1 = zeros(M,M,N); % apply xy
    for i = 1:N
        xi = reshape(x(:,:,i),N,N);
        b1i = A*xi*A';
        b1(:,:,i) = b1i;
    end
    b = reshape(b1,M^2,N)*A';
elseif dim==1  % apply only to one dimension
    b = zeros(M,N,N); % apply z to xy result
    for i = 1:N
        b(:,:,i) = reshape(x(:,:,i),N,N)*A';
    end
elseif dim==2
    b = zeros(N,M,N);
    for i = 1:N
        b(:,:,i) = A*reshape(x(:,:,i),N,N);
    end
elseif dim==3
    b = zeros(N,N,M);
    for i = 1:N
        b(i,:,:) = reshape(x(i,:,:),N,N)*A';
    end
end

b = b(:);

return

