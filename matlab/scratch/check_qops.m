clear
Globals2D

N = 2;
Nq = 2*N;

[r s] = Nodes2D(N); [r s] = xytors(r,s);
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
% [r s w] = QNodes2D(N,Nq); [r s] = xytors(r,s); % integrate u*v*c
% r = rq; s = sq; wq = w;

[rq1D wq1D] = JacobiGQ(0,0,ceil(Nq/2));
rfq = [rq1D; -rq1D; -ones(size(rq1D))];
sfq = [-ones(size(rq1D)); rq1D; -rq1D];
wfq = [wq1D; wq1D; wq1D];


opt = 1;
if opt==1
    V = Vandermonde2D(N,r,s);
    [Vr Vs] = GradVandermonde2D(N,r,s); 
    Dr = Vr/V; Ds = Vs/V;
    Vq = Vandermonde2D(N,rq,sq)/V;
    Vfq = Vandermonde2D(N,rfq,sfq)/V;
elseif opt==2
    [Vq Vrq Vsq] = bern_basis_tri(N,rq,sq);
    Dr = Vq\Vrq; Ds = Vq\Vsq;
    Vfq = bern_basis_tri(N,rfq,sfq);
    
elseif opt==3
    V = Vandermonde2D(N,r,s);
    [Vr Vs] = GradVandermonde2D(N,r,s); 
    Dr = V\Vr; Ds = V\Vs;
    Vq = Vandermonde2D(N,rq,sq);
    Vfq = Vandermonde2D(N,rfq,sfq);
end

M = Vq'*diag(wq)*Vq;
Pq = M\(Vq'*diag(wq)); % J's cancel out
VqPq = Vq*Pq;


Vq1D = Vandermonde1D(N,rq1D)/Vandermonde1D(N,JacobiGL(0,0,N));
% plot(rfq,sfq,'o')
Nfq = length(rq1D);
Nq = length(rq);
Np = length(r);
Nfaces = 3;

Vfqf = kron(eye(Nfaces),Vq1D);
Mf = Vfq'*diag(wfq)*Vfq;
Lq = M\(Vfq'*diag(wfq));

Pq1D = (Vq1D'*diag(wq1D)*Vq1D) \ (Vq1D'*diag(wq1D));
Pfqf = kron(eye(Nfaces),Pq1D);

nrJ = [-zeros(size(rq1D)); ones(size(rq1D)); -ones(size(rq1D)); ];
nsJ = [-ones(size(rq1D)); ones(size(rq1D)); -zeros(size(rq1D)); ];
Nq = length(rq);
nrJq = repmat(nrJ',Nq,1);
nsJq = repmat(nsJ',Nq,1);

% flux differencing operators
global Drq Dsq VfPq VqLq
Drq = (Vq*Dr*Pq - .5*Vq*Lq*diag(nrJ)*Vfq*Pq);
Dsq = (Vq*Ds*Pq - .5*Vq*Lq*diag(nsJ)*Vfq*Pq);
VfPq = (Vfq*Pq);
VqLq = Vq*Lq;

DNr = [Drq .5*VqLq*diag(nrJ); -.5*diag(nrJ)*VfPq .5*diag(nrJ)*eye(Nfq*3)];

W = diag([wq;wfq]);
Q = W*DNr;
norm(Q + Q' - diag([0*wq;wfq.*nrJ]),'fro')
size(DNr)

% VV'*W*DNr*VV


%% build shu matrices

N = 2;
Nq = 2*N-1;

[r s] = Nodes2D(N); [r s] = xytors(r,s);
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
% [r s w] = QNodes2D(N,Nq); [r s] = xytors(r,s); % integrate u*v*c
% r = rq; s = sq; wq = w;
rq = -1/3; sq = -1/3; % for N = 2;

[rq1D wq1D] = JacobiGQ(0,0,ceil(Nq/2));
rfq = [rq1D; -rq1D; -ones(size(rq1D))];
sfq = [-ones(size(rq1D)); rq1D; -rq1D];
wfq = [wq1D; wq1D; wq1D];
nrJ = [-zeros(size(rq1D)); ones(size(rq1D)); -ones(size(rq1D)); ];
nsJ = [-ones(size(rq1D)); ones(size(rq1D)); -zeros(size(rq1D)); ];

rq = [rq;rfq];
sq = [sq;sfq];
plot(rq,sq,'o')

Vq = Vandermonde2D(Nq,rq,sq);
b = zeros(size(Vq,2),1); b(1) = sqrt(2);
wq = pinv(Vq')*b;

V = Vandermonde2D(N,rq,sq);
[Vr Vs] = GradVandermonde2D(N,rq,sq);
Dr = V\Vr;
Ds = V\Vs;
I = eye(Nfq*Nfaces) ;

Dr = .5*(diag(1./wq) + V*V')*diag([zeros(length(rq)-length(rfq),1);wfq.*nrJ])*(eye(length(wq)) - V*V'*diag(wq)) + V*Dr*V'*diag(wq);
Ds = .5*(diag(1./wq) + V*V')*diag([zeros(length(rq)-length(rfq),1);wfq.*nsJ])*(eye(length(wq)) - V*V'*diag(wq)) + V*Ds*V'*diag(wq);
% 
% Qr = diag(wq)*Dr;
return
% replicate using non-orthonormal mats
VV = Vandermonde2D(N,r,s);
[Vr Vs] = GradVandermonde2D(N,r,s);
Dr = Vr/VV; Ds = Vs/VV;

Vq = Vandermonde2D(N,rq,sq);
M = Vq'*diag(wq)*Vq;
Pq = M\(Vq'*diag(wq));

% .5*(diag(1./wq) + V*V')*diag([zeros(length(rq)-length(rfq),1);wfq.*nrJ])*(eye(length(wq)) - V*V'*diag(wq)) + V*Dr*V'*diag(wq)

% Dr = .5*diag(1./wq)*(diag(wq) + Pq'*Vq')*diag([zeros(length(rq)-length(rfq),1);wfq.*nrJ])*(eye(length(wq)) - Vq*Pq) + Vq*Dr*Pq



