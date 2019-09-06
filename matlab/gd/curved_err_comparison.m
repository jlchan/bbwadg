clear

N = 5;
K = 8;

VX = (-K/2:K/2);

% piecewise quadrature
rq1D = []; wq = [];
h = 1; % for constructing GD basis
for e = 1:K
    [rqe wqe] = JacobiGQ(0,0,N);
    D1D = GradVandermonde1D(N,rqe)/Vandermonde1D(N,rqe);
    rq1D = [rq1D; h*(rqe+1)/2 + VX(e)];
    wq = [wq; h/2*wqe];
end

% map to [-1,1]
rq = rq1D*2/K;
wq = wq*2/K;

%% form VDMs

if K > ceil(N/2)
    % GD (uses [-K/2,K/2])
    VqGD = GDVDM(N,K,rq*K/2);
else
    disp('Increase K or decrease N to use GD')
end

% B-spline 
smoothKnots = 1;
VqBS = bsplineVDM(N,K,rq,smoothKnots);

% polynomial 
r = JacobiGL(0,0,N); 
V = Vandermonde1D(N,r);
VqP = Vandermonde1D(N,rq)/V;
D1D = GradVandermonde1D(N,r)/V;

%% make 2D ops

[x y] = meshgrid(r);
x = x(:); y = y(:);

[wr ws] = meshgrid(wq);
wq = wr(:).*ws(:);

Dr = kron(D1D,speye(size(D1D,1)));
Ds = kron(speye(size(D1D,1)),D1D);
Vq = kron(VqP,VqP);

%% curve domain

xq = Vq*x;
yq = Vq*y;

% non-affine mappings
a = 1/8;
x = x + a*cos(.5*3*pi*y).*cos(pi/2*x);
y = y + a*sin(.5*3*pi*x).*cos(pi/2*y);

xqJ = Vq*x;
yqJ = Vq*y;
[~,~,~,~,Jq] = GeometricFactors2D(x,y,Vq*Dr,Vq*Ds);

%% check error

% Vq1D = VqGD;
Vq1D = VqBS;
% Vq1D = VqP;
Vq = kron(Vq1D,Vq1D);
PqJ = (Vq'*diag(wq.*Jq)*Vq)\(Vq'*diag(wq.*Jq));
Pq = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq));

% affine mesh
u = cos(pi/2*xq).*cos(pi/2*yq);
L2errAffine = sqrt(sum(sum(wq.*(u-Vq*Pq*u).^2)));

% curved mesh
uJ = cos(pi/2*xqJ).*cos(pi/2*yqJ);
L2errCurved = sqrt(sum(sum(wq.*Jq.*(uJ-Vq*PqJ*uJ).^2)));

% err = h^(-(N+1)), h = 
ndofs = size(Vq,2);
fprintf('Affine err = %g, curved err = %g, increase = %g\n',L2errAffine,L2errCurved,L2errCurved/L2errAffine)

