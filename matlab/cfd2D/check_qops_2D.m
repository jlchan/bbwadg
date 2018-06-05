Globals2D

N = 2;
K1D = 1;
FinalTime = 10;
CFL = .25;
global tau
tau = 1;
a = 0/16; % curv warping

Lx = 10; Ly = 5;
if 0   
    [Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D,round(K1D*Ly/Lx));
    VX = VX/max(abs(VX));  VY = VY/max(abs(VY));
    VX = VX*Lx + Lx; VY = VY*Ly;
else
    [Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
end

StartUp2D;
BuildPeriodicMaps2D(max(VX)-min(VX),max(VY)-min(VY));

% plotting nodes
[rp sp] = EquiNodes2D(15); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;
% plot(xp,yp,'o'); return

global Vq Pq Lq Lqf Vff Vf Pfqf VqPq
global rxJ sxJ ryJ syJ rxJf sxJf ryJf syJf
global nxJ nyJ nrJ nsJ nrJq nsJq
global mapPq
Nq = 2*N;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c

[rq1D wq1D] = JacobiGQ(0,0,ceil(Nq/2));
rfq = [rq1D; -rq1D; -ones(size(rq1D))];
sfq = [-ones(size(rq1D)); rq1D; -rq1D];
wfq = [wq1D; wq1D; wq1D];
Vq1D = Vandermonde1D(N,rq1D)/Vandermonde1D(N,JacobiGL(0,0,N));
% plot(rfq,sfq,'o')
Nfq = length(rq1D);

opt=2;
if opt==1
    % compute shu vol quadrature
    r0 = -1/3;
    s0 = -1/3;
    rq = [r0;rfq];
    sq = [s0;sfq];    
else
    % compute hicken quadrature
    r0 = -1/3;
    s0 = -1/3;
    r1D = JacobiGL(0,0,N);
    rf = [r1D; -r1D; -ones(size(r1D))];
    sf = [-ones(size(r1D)); r1D; -r1D];
    rq = [r0;rf];
    sq = [s0;sf];    
    [rs] = uniquetol([rq sq],'ByRows',true);
    rq = rs(:,1); sq = rs(:,2);
end

Vqtmp = Vandermonde2D(2*N-1,rq,sq);
b = zeros(size(Vqtmp,2),1); b(1) = sqrt(2);
wq = (Vqtmp')\b;
norm(Vqtmp'*wq-b)

[Vqmodal] = Vandermonde2D(N,rq,sq);
Vq = Vqmodal/V;
Vrq = GradVandermonde2D(N,rq,sq);
M = Vq'*diag(wq)*Vq;
Pq = M\(Vq'*diag(wq)); % J's cancel out
VqPq = Vq*Pq;

Vf = Vandermonde2D(N,rfq,sfq)/V;
Mf = Vf'*diag(wfq)*Vf;
Lq = M\(Vf'*diag(wfq));

Pq1D = (Vq1D'*diag(wq1D)*Vq1D) \ (Vq1D'*diag(wq1D));

nrJ = [-zeros(size(rq1D)); ones(size(rq1D)); -ones(size(rq1D)); ];
nsJ = [-ones(size(rq1D)); ones(size(rq1D)); -zeros(size(rq1D)); ];
Nq = length(rq);
nrJq = repmat(nrJ',Nq,1);
nsJq = repmat(nsJ',Nq,1);

% flux differencing operators
global Drq Dsq VfPq VqLq
Drq = (Vq*Dr*Pq - .5*Vq*Lq*diag(nrJ)*Vf*Pq);
Dsq = (Vq*Ds*Pq - .5*Vq*Lq*diag(nsJ)*Vf*Pq);
VfPq = (Vf*Pq);
VqLq = Vq*Lq;


DNr = [Drq .5*Vq*Lq*diag(nrJ);-.5*diag(nrJ)*Vf*Pq .5*diag(nrJ)];
DNs = [Dsq .5*Vq*Lq*diag(nsJ);-.5*diag(nsJ)*Vf*Pq .5*diag(nsJ)];
WN = diag([wq; wfq]);

QNr = WN*DNr;
% QNs = WN*DNs;
% norm(QNr+QNr' - diag([wq*0; nrJ.*wfq]),'fro')
% norm(QNs+QNs' - diag([wq*0; nsJ.*wfq]),'fro')

if opt==1
    Ef = [zeros(length(rfq),1) eye(length(rfq))]; % shu operators - diagonal bdry norms
else
    Ef = Vf*Pq; % hicken operators - dense bdry norms
end
Iq = eye(length(rq));

Qr = [Iq Ef']*QNr*[Iq;Ef];
Dr_sbp = diag(1./wq)*Qr;

if opt==1
    norm(Qr+Qr' - Ef'*diag(wfq.*nrJ)*Ef,'fro')
else
    Dr_sbp - Vq*Dr*Pq 

    norm(Qr+Qr' - (Vf*Pq)'*diag(wfq.*nrJ)*(Vf*Pq),'fro')
end
Er = (Vf*Pq)'*diag(wfq.*nrJ)*(Vf*Pq);

%% build hicken sbp ops

[U S V] = svd(Vqmodal);
W = U(:,size(Vqmodal,2)+1:size(Vqmodal,1));
Vtilde = [Vqmodal W];
m = size(Vqmodal,1)-size(Vqmodal,2);
H = diag(wq);

opt=2;
if opt==1
    Wx = zeros(size(Vqmodal,1),m);
    Gx = zeros(m,m);
    Fx = W'*H*Vrq;
elseif opt==2
    Fx = W'*H*Vrq*0;
    Wx = pinv(Vqmodal'*H)*(-Vrq'*H*W);
    Gx = W'*H*Wx + Wx'*H*W;
end

Etilde = Vqmodal'*H*Vrq + Vrq'*H*Vqmodal;
Eblk = [Etilde Fx';
    Fx Gx];
    
Ex = inv(Vtilde)'*Eblk*inv(Vtilde);

Vxtilde = [Vrq Wx];
Sx = H*Vxtilde*inv(Vtilde)-.5*Ex;
Qx = Sx + .5*Ex;

Dx = H\Qx;


