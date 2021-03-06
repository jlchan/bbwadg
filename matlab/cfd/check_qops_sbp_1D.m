clear

[H,Q,rr]=CSBPp2(8);

N = 2;

rq = rr;
dx = rr(2)-rr(1);
wq = diag(H);

% N = 4;
% r = JacobiGL(0,0,N);
% r = JacobiGQ(0,0,N);
% 
% [rq wq] = JacobiGL(0,0,N); rq = .9999999999*rq; % N+1 -> deg 2(N+1)-3
% [rq wq] = JacobiGL(0,0,N+1); % 
% 
% % rq = [-1;rq;1];
% % wq = [0;wq;0];
% return

Nq = length(rq);

r = JacobiGL(0,0,N);
V = Vandermonde1D(N,r);
Dr = GradVandermonde1D(N,r)/V;
Vq = Vandermonde1D(N,rq)/V;
Vf = Vandermonde1D(N,[-1 1])/V;

% [Vq, ~, Dr] = bsplineVDM(1,14,rr,0,N); % VDM for interp, mass, M\S
% Vf = zeros(2,size(Vq,2));Vf(1) = 1; Vf(end) = 1;

M = Vq'*diag(wq)*Vq;
LIFT = M\Vf';
Lq = LIFT;

Pq = M\(Vq'*diag(wq));
Pq(abs(Pq)<1e-8) = 0;

VqPq = Vq*Pq;
tf = zeros(2,length(rq)); tf(1) = 1; tf(end) = 1;
Drq = Vq*Dr*Pq;
Dfq = .5*Vq*LIFT*diag([-1,1])*(tf - tf*VqPq);

VqLq = Vq*Lq;
VfPq = Vf*Pq;

nrJ = [-1;1];
DNr = [Drq-.5*Vq*Lq*diag(nrJ)*VfPq .5*VqLq*diag(nrJ); -.5*diag(nrJ)*VfPq .5*diag(nrJ)*eye(2)];
W = diag(wq);

% face extraction
E = zeros(2,Nq); E(1,1) = 1; E(end,end) = 1;
wfq = [1;1];
W = diag([wq;wfq]);
Q2 = [eye(Nq);E]'*diag([wq;wfq])*DNr*[eye(Nq);E]


return

%% build shu operators
V = Vandermonde1D(N,rq);
Vr = GradVandermonde1D(N,rq);
D1 = Vr*V'*diag(wq) + .5*diag(1./wq)* (eye(length(wq)) + diag(wq)*V*V')*diag([-1 zeros(1,length(rq)-2) 1])*(eye(length(wq))-V*V'*diag(wq));

% face extraction
E = zeros(2,Nq); E(1,1) = 1; E(end,end) = 1;

VQ = Vandermonde1D(2*N,rq);
b = zeros(size(VQ,2),1); b(1) = sqrt(2);
wq = pinv(VQ')*b;

% build shu operator but my way
V = Vandermonde1D(N,r);
Vq = Vandermonde1D(N,rq)/V;
Vf = Vandermonde1D(N,[-1;1])/V;
Dr = GradVandermonde1D(N,r)/V;
M = Vq'*diag(wq)*Vq;
Pq = M\(Vq'*diag(wq));
Lq = M\(Vf');
VfPq = Vf*Pq;
VqLq = Vq*Lq;

wqinv = 1./wq; wqinv(abs(wqinv) > 1e10) = 1; % pinv of W
D1 = Vq*Dr*Pq + .5*diag(wqinv)*(eye(length(wq)) + Vq*Pq)'*E'*diag([-1,1])*E*(eye(length(wq))-Vq*Pq);
W = diag(wq);

DNr = [Drq-.5*Vq*Lq*diag(nrJ)*VfPq .5*VqLq*diag(nrJ); -.5*diag(nrJ)*VfPq .5*diag(nrJ)*eye(2)];

D1 = Vq*Dr*Pq + .5*diag(wqinv)*(E' + Pq'*Vf')*diag([-1,1])*(E-Vf*Pq);

wfq = [1;1];
W = diag([wq;wfq]);
D2 = diag(wqinv)*[eye(Nq);E]'*diag([wq;wfq])*DNr*[eye(Nq);E]
Q2 = W*D2




