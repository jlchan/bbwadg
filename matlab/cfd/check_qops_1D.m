clear

N = 4;


r = JacobiGL(0,0,N);
% r = JacobiGQ(0,0,N);

[rq wq] = JacobiGL(0,0,N); rq = .9999999999*rq; % N+1 -> deg 2(N+1)-3
[rq wq] = JacobiGL(0,0,N+1); % 

% rq = [-1;rq;1];
% wq = [0;wq;0];

Nq = length(rq);

V = Vandermonde1D(N,r);
Dr = GradVandermonde1D(N,r)/V;
Vq = Vandermonde1D(N,rq)/V;
Vf = Vandermonde1D(N,[-1 1])/V;

M = Vq'*diag(wq)*Vq;
LIFT = M\Vf';
Lq = LIFT;

Pq = M\(Vq'*diag(wq));
Pq(abs(Pq)<1e-8) = 0;

rp = linspace(-1,1,50); F = eye(N+1);
Vp = Vandermonde1D(N,rp)/V;

VqPq = Vq*Pq;
Drq = Vq*Dr*Pq;

VqLq = Vq*Lq;
VfPq = Vf*Pq;

nrJ = [-1;1];
DNr = [Drq-.5*Vq*Lq*diag(nrJ)*VfPq .5*VqLq*diag(nrJ); -.5*diag(nrJ)*VfPq .5*diag(nrJ)*eye(2)];

W = diag([wq;1;1]);

Ef = (Vf*Pq)'*diag([1;1]);
% Pq = eye(length(rq));

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




