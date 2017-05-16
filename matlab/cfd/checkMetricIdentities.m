clear
Globals2D

N = 4;
K1D = 8;

[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);

% Initialize solver and construct grid and metric
StartUp2D; 

vmapM = reshape(vmapM,Nfp*Nfaces,K);
vmapP = reshape(vmapP,Nfp*Nfaces,K);
mapM = reshape(mapM,Nfp*Nfaces,K);
mapP = reshape(mapP,Nfp*Nfaces,K);

a = 1/8;
x = x + a*cos(.5*3*pi*y).*cos(pi/2*x);
y = y + a*sin(.5*3*pi*x).*cos(pi/2*y);

[rp sp] = EquiNodes2D(50); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

[rx,sx,ry,sy,J] = GeometricFactors2D(x,y,Dr,Ds);
[nx, ny, sJ] = Normals2D();

Nq = 3*N; % in 3D, this is 4N-2. 
[rq sq wq] = Cubature2D(Nq);
Vq = Vandermonde2D(N,rq,sq)/V;
xq = Vq*x; yq = Vq*y;
Pq = (V*V')*Vq'*diag(wq);
[Vrq Vsq] = GradVandermonde2D(N,rq,sq);
Vrq = Vrq/V; Vsq = Vsq/V;

[rq1D wq1D] = JacobiGQ(0,0,Nq);
wfq = [wq1D;wq1D;wq1D];
rfq = [rq1D; -rq1D; -ones(size(rq1D))];
sfq = [-ones(size(rq1D));rq1D; -rq1D];
% plot(rfq,sfq,'o');return

Vfq = Vandermonde2D(N,rfq,sfq)/V;

V1D = Vandermonde1D(N,JacobiGL(0,0,N));
Vfqf = Vandermonde1D(N,rq1D)/V1D;
Pfqf = V1D*V1D'*Vfqf'*diag(wq1D);
Vfqf = blkdiag(Vfqf,Vfqf,Vfqf);
Pfqf = blkdiag(Pfqf,Pfqf,Pfqf);

% DivHat( J*[ rx sx; ry sy]) = DivHat(J*G^T) = 0
norm(Dr*(rx.*J) + Ds*(sx.*J),'fro')
norm(Dr*(ry.*J) + Ds*(sy.*J),'fro')

rxJ = rx.*J; sxJ = sx.*J;
ryJ = ry.*J; syJ = sy.*J;

rxJ = rxJ(Fmask(:),:); sxJ = sxJ(Fmask(:),:);
ryJ = ryJ(Fmask(:),:); syJ = syJ(Fmask(:),:);

nxJ = nx.*sJ;
nyJ = ny.*sJ;

fids = 1:N+1;    

% reference normals
f = 1; ffids = fids + (f-1)*Nfp;
nr(ffids,1) = 0; ns(ffids,1) = -1;
f = 2; ffids = fids + (f-1)*Nfp;
nr(ffids) = 1; ns(ffids) = 1;
f = 3; ffids = fids + (f-1)*Nfp;
nr(ffids) = -1; ns(ffids) = 0;
nr = repmat(nr,1,K);
ns = repmat(ns,1,K);

norm(nxJ - (rxJ.*nr + sxJ.*ns),'fro')
norm(nyJ - (ryJ.*nr + syJ.*ns),'fro')
% return

% quadrature arrays
rxJ = Vq*(rx.*J); sxJ = Vq*(sx.*J);
ryJ = Vq*(ry.*J); syJ = Vq*(sy.*J);
Jq = Vq*J; wJq = spdiag(wq)*Jq;
nxJ = Vfqf*(nx.*sJ); nyJ = Vfqf*(ny.*sJ);

if 0
    % test IBP - shouldn't be satisfied unless Nq = 3*N
    e = 1;
    u = x(:,e).^N + y(:,e).^N;
    dudxJ = rxJ(:,e).*(Vrq*u) + sxJ(:,e).*(Vsq*u);
    
    r1 = Vq'*diag(wq)*dudxJ;
    r2 = Vrq'*diag(rxJ(:,e).*wq)*Vq*u + Vsq'*diag(sxJ(:,e).*wq)*Vq*u - Vfq'*diag(wfq.*nxJ(:,e))*Vfq*u;
    norm(r1 + r2)
end

% compute physical derivative and check accuracy
f = @(x,y) sin(pi*x).*sin(pi*y);
fx = @(x,y) pi*cos(pi*x).*sin(pi*y);
fy = @(x,y) pi*sin(pi*x).*cos(pi*y);

% f = @(x,y) exp(x).*exp(y);
% fx = @(x,y) f(x,y);
% fy = @(x,y) f(x,y);

u = Pq*f(xq,yq);
sum(sum(wJq.*(f(xq,yq)-Vq*u).^2))

if 0   % discrete gradient = J*G*Dhat(u)
    du = u(vmapP)-u(vmapM);
    ur = Dr*u - .5*LIFT*(nr.*du);
    us = Ds*u - .5*LIFT*(ns.*du);
    ux = (Vq*(Pq*(rxJ.*(Vq*ur) + sxJ.*(Vq*us))))./Jq;
    uy = (Vq*(Pq*(ryJ.*(Vq*ur) + syJ.*(Vq*us))))./Jq;
    sum(sum(wJq.*((ux-fx(xq,yq)).^2 + (uy-fy(xq,yq)).^2)))
end

% discrete divergence = Dhat(J*G^T*u)
fx = @(x,y) sin(pi*x);
fy = @(x,y) sin(pi*y);
divf = @(x,y) pi*(cos(pi*x)+cos(pi*y));

% project (J*G^T*u)
ux = Pq*fx(xq,yq); uy = Pq*fy(xq,yq);
uxq = Vq*ux; uyq = Vq*uy;
ur = Pq*(rxJ.*uxq + ryJ.*uyq);
us = Pq*(sxJ.*uxq + syJ.*uyq);

% compute jump differently
du1 = ur(vmapP).*nr(mapP) + ur(vmapM).*nr(mapM);
du2 = us(vmapP).*ns(mapP) + us(vmapM).*ns(mapM);
dUnhat = du1 + du2;
divu = Dr*ur + Ds*us - .5*LIFT*(dUnhat);

sum(sum(wJq.*(divf(xq,yq) - (Vq*divu)./Jq).^2))


