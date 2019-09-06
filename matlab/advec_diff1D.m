clear
% close all;

N = 3; % polynomial degree
K = 8;
CFL = .25;
FinalTime = 2;

global Fmask mapP tau rxJ nxJ J M Dr L 
tau = 0; % penalty parameter - try varying between [0,1]

%% set up reference element

% nodes + quadrature
r = JacobiGL(0,0,N); % interpolation points
[rq wq] = JacobiGQ(0,0,N+1); % compute quadrature points

% VDM matrices
V = Vandermonde1D(N,r); % modal VDM
Vrq = GradVandermonde1D(N,rq)/V;
Vq = Vandermonde1D(N,rq)/V;

% mass and convection matrices
M = Vq'*diag(wq)*Vq;
Q = Vq'*diag(wq)*Vrq;

B = zeros(N+1,2);
B(1,1) = 1;
B(end,end) = 1;
Fmask = [1;N+1]; % face node indices

Dr = M\Q;
L = M\B;


%% new decoupled operators
global QNr VN

Pq = M\(Vq'*diag(wq));
Vf = Vandermonde1D(N,[-1,1])/V;

Q = Pq'*Q*Pq;

E = Vf*Pq;
B = diag([-1 1]);
QNr = [Q - .5*E'*B*E .5*E'*B;
    -.5*B*E .5*B];

VN = [Vq;Vf];

%% physical elements

VX = linspace(-1,1,K+1);

x = zeros(N+1,K);
for e = 1:K
    h = VX(e+1)-VX(e);
    x(:,e) = h*(r+1)/2 + VX(e);
end

rp = linspace(-1,1,25)';
Vp = Vandermonde1D(N,rp)/V;
xp = Vp*x;

Nfp = 1;
Nfaces = 2;
mapP = reshape(1:Nfp*Nfaces*K,Nfp*Nfaces,K);
for e = 1:K
    mapP(:,e) = [(Nfp*Nfaces)*e-2, Nfp*Nfaces*(e+1)-1];
end
mapP(1) = 1; mapP(end) = Nfp*Nfaces*K;
mapP(end) = 1; mapP(1) = Nfp*Nfaces*K;
mapP = reshape(mapP,Nfp*Nfaces,K);

% geometric terms - these may be pedantic, but will set the stage for 2D/3D
nxJ = repmat([-1;1],1,K); % "scaled" normals
J = h/2*ones(N+1,K); % determinant of Jacobian of mapping
rx = 2/h*ones(N+1,K); % geometric mapping term for derivative
rxJ = rx.*J; % scaled geometric mapping

%% compute matrix

u = zeros(N+1,K);
A = zeros((N+1)*K);
for i = 1:(N+1)*K
   u(i) = 1;
   rhs = diffRHS(u);
   A(:,i) = rhs(:);
   u(i) = 0;
end

keyboard
%%

function rhs = diffRHS(u)

global Fmask mapP tau rxJ nxJ J M Dr L

uf = u(Fmask(:),:);
uavg = .5*(uf(mapP)+uf);

% computes viscous flux
uflux = uavg - uf;
q = (rxJ.*(Dr*u) + L*(uflux.*nxJ))./J;

% computes "2nd derivative"
qf = q(Fmask(:),:);
qavg = .5*(qf(mapP)+qf);

qflux = (qavg-qf).*nxJ; % BR1
tau = 1; % must scale as O(1/h)
qflux = (qavg-qf).*nxJ + tau*(uf(mapP)-uf); % LDG with C11 = tau, C12 = 0

% Du = (rxJ.*(Dr*u))./J;
% Duf = Du(Fmask(:),:);
% Duavg = .5*(Duf(mapP)+Duf);
% tau = 25;
% qflux = (Duavg - qf).*nxJ + tau*(uf(mapP)-uf); % IPDG

rhs = (rxJ.*(Dr*q) + L*(qflux))./J; % dq/dx DG style

end