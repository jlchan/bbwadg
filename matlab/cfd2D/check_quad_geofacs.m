clear -globals
Globals2D

N = 2;
K1D = 1;
FinalTime = 5;

wadgProjEntropyVars = abs(a)>1e-8;
CFL = .5;
global tau
tau = 1;

[Nv, VX, VY, K, EToV] = QuadMesh2D(K1D,K1D);

VX = VX + .25*randn(size(VX));
VY = VY + .25*randn(size(VX));
StartUp2D;
BuildPeriodicMaps2D(max(VX)-min(VX),max(VY)-min(VY));

V1D = Vandermonde1D(N,JacobiGL(0,0,N));

[rq1D wq1D] = JacobiGQ(0,0,N);
[rq sq] = meshgrid(rq1D);
[wr ws] = meshgrid(wq1D);
rq = rq(:);
sq = sq(:);
wq = wr(:).*ws(:);

Vq = Vandermonde2DQuad(N,rq,sq)/V;
[rx,sx,ry,sy,J] = GeometricFactors2D(Vq*x,Vq*y,Vq*Dr/Vq,Vq*Ds/Vq);
rxJ = reshape(rx.*J,N+1,N+1);
sxJ = reshape(sx.*J,N+1,N+1);
ryJ = reshape(ry.*J,N+1,N+1);
syJ = reshape(sy.*J,N+1,N+1);

vv = rxJ; color_line3(rq,sq,vv,vv,'.')

Qr = diag(wq)*(Vq*Dr/Vq);
Qs = diag(wq)*(Vq*Ds/Vq);
Qx = diag(rxJ(:))*Qr + diag(sxJ(:))*Qs;
Qy = diag(ryJ(:))*Qr + diag(syJ(:))*Qs;

% PlotMesh2D;axis on;return

