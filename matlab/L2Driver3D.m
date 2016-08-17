% Driver script for solving the 3D vacuum Maxwell's equations

Globals3D;

% Polynomial order of approximation 
N = 2;

% Read in Mesh
% [Nv, VX, VY, VZ, K, EToV] = MeshReaderGambit3D('cubeK6.neu');

% Initialize solver and construct grid and metric
StartUp3D;

% Build cubature
[za,wa] = JacobiGQ(0, 0, N);
[zb,wb] = JacobiGQ(1, 0, N);
[zc,wc] = JacobiGQ(2, 0, N);

qa = zeros(N+1, N+1, N+1);
qb = zeros(N+1, N+1, N+1);
qc = zeros(N+1, N+1, N+1);
qw = zeros(N+1,N+1,N+1);

for i=1:N+1
  for j=1:N+1
    for k=1:N+1

      qa(i,j,k) = za(i);
      qb(i,j,k) = zb(j);
      qc(i,j,k) = zc(k);

      qw(i,j,k) = wa(i)*wb(j)*wc(k);
    end
  end
end

qr = 0.25*(1+qa).*(1-qb).*(1-qc)-1;
qs = 0.50*(1+qb).*(1-qc)-1;
qt = qc;

qr = qr(:);
qs = qs(:);
qt = qt(:);
qw = qw(:);

qI = Vandermonde3D(N, qr, qs, qt)/V;
qx = qI*x;
qy = qI*y;
qz = qI*z;

uex = @(xo,yo,zo) sin(pi*xo);


quex = feval(uex, qx, qy, qz);

quh = 0*quex;

qerr = 0;
for k=1:K
  quh(:,k) = qI*((qI'*diag(qw(:))*qI)\(qI'*diag(qw(:))*quex(:,k)));
  
  qerr = qerr + J(1,k)*qw'*(( quh(:,k) - quex(:,k) ).^2);
end

err = sqrt(sum(qerr))
  





