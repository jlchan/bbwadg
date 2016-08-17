clear

% Driver script for solving the 3D IPDG acoustic wave equation
Globals3D;

% N = inN;
% Order of polymomials used for approximation
N = 4;

% % % single element
[VX VY VZ] = Nodes3D(1);
[VX VY VZ] = xyztorst(VX,VY,VZ);
K = 1; EToV = 1:length(VX);

% Initialize solver and construct grid and metric
StartUp3D;

[rq sq tq wq] = tet_cubature(2*N+1);
Vq = Vandermonde3D(N,rq,sq,tq)/V;
Pq = Vq'*diag(wq);

[rqt sqt wqtri] = Cubature2D(2*N);
VfqFace = Vandermonde2D(N,rqt,sqt)/Vandermonde2D(N,r(Fmask(:,1)),s(Fmask(:,1)));

f = 1;
% hold on;
% plot3(r(Fmask(:,f)),s(Fmask(:,f)),t(Fmask(:,f)),'o');
% text(r(Fmask(:,f))+.1,s(Fmask(:,f)),t(Fmask(:,f)),num2str((1:size(Fmask,1))'))

[rfq sfq tfq wfq] = tet_surface_cubature(2*N);
plot3(rfq,sfq,tfq,'o')
% text(rfq+.1,sfq,tfq,num2str((1:length(rfq))'))
wfq = [wqtri; wqtri; wqtri*sqrt(2); wqtri];
Vfq = Vandermonde3D(N,rfq,sfq,tfq)/V;

% rf = Vfq*r;
% sf = Vfq*s;
% tf = Vfq*t;
rf = VfqFace*r(Fmask);
sf = VfqFace*s(Fmask);
tf = VfqFace*t(Fmask);

plot3(rf,sf,tf,'o')
