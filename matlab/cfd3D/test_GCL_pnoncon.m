clear
N = 2;
N2 = 3;
r1D1 = JacobiGL(0,0,N);
[rr s1 t1] = meshgrid(r1D1);
rr = rr(:); s1 = s1(:); t1 = t1(:);

[r1D2 w1D2] = JacobiGL(0,0,N2);
[rs s2 t2] = meshgrid(r1D2);
rs = rs(:); s2 = s2(:); t2 = t2(:);

% interface at x = 1
x1 = rr;
y1 = s1;
z1 = t1;
x2 = rs+2;
y2 = s2;
z2 = t2;

ids1 = find(abs(x1-1)<1e-8);
ids2 = find(abs(x2-1)<1e-8);

a = .25/N^2;
x1 = x1 + a*randn(size(rr));
y1 = y1 + a*randn(size(rr));
z1 = z1 + a*randn(size(rr));

x2 = x2 + 1*a*randn(size(rs));
y2 = y2 + 1*a*randn(size(rs));
z2 = z2 + 1*a*randn(size(rs));

V1D1 = Vandermonde1D(N,r1D1);
D1D1 = GradVandermonde1D(N,r1D1)/V1D1;
I1 = eye(length(r1D1));
Dr1 = kron(kron(I1,D1D1),I1);
Ds1 = kron(kron(I1,I1),D1D1);
Dt1 = kron(kron(D1D1,I1),I1);

V1D2 = Vandermonde1D(N2,r1D2);
D1D2 = GradVandermonde1D(N2,r1D2)/V1D2;
I2 = eye(length(r1D2));
Dr2 = kron(kron(I2,D1D2),I2);
Ds2 = kron(kron(I2,I2),D1D2);
Dt2 = kron(kron(D1D2,I2),I2);

F1D = eye(N+1); %(Vandermonde1D(N-1,r1D1)/Vandermonde1D(N-1,JacobiGL(0,0,N-1))) * (Vandermonde1D(N,JacobiGL(0,0,N-1))/V1D1);
Fr1 = kron(kron(I1,F1D),I1);
Fs1 = kron(kron(I1,I1),F1D);
Ft1 = kron(kron(F1D,I1),I1);

F1D = eye(N2+1); %(Vandermonde1D(N2-1,r1D2)/Vandermonde1D(N2-1,JacobiGL(0,0,N2-1))) * (Vandermonde1D(N2,JacobiGL(0,0,N2-1))/V1D2);
Fr2 = kron(kron(I2,F1D),I2);
Fs2 = kron(kron(I2,I2),F1D);
Ft2 = kron(kron(F1D,I2),I2);

V1D_12 = Vandermonde1D(N2,r1D1); % interp from N2 to N nodes
V12 = kron(V1D_12,V1D_12) / kron(V1D2,V1D2);

% enforce continuity: N1 = N nodes set by N2 = N-1 nodes.
x1(ids1) = V12*x2(ids2);
y1(ids1) = V12*y2(ids2);
z1(ids1) = V12*z2(ids2);

rr = Fr1*(Dr1*y1).*z1;
rs = Fs1*(Ds1*y1).*z1;
rt = Ft1*(Dt1*y1).*z1;
rx1 = Dt1*(rs) - Ds1*(rt); % this is the problematic one 
sx1 = Dr1*(rt) - Dt1*(rr);
tx1 = Ds1*(rr) - Dr1*(rs);

rr = Fr1*(Dr1*x1).*z1;
rs = Fs1*(Ds1*x1).*z1;
rt = Ft1*(Dt1*x1).*z1;
ry1 = -(Dt1*(rs) - Ds1*(rt));
sy1 = -(Dr1*(rt) - Dt1*(rr));
ty1 = -(Ds1*(rr) - Dr1*(rs));

rr = Fr1*(Dr1*y1).*x1;
rs = Fs1*(Ds1*y1).*x1;
rt = Ft1*(Dt1*y1).*x1;
rz1 = -(Dt1*(rs) - Ds1*(rt));
sz1 = -(Dr1*(rt) - Dt1*(rr));
tz1 = -(Ds1*(rr) - Dr1*(rs));

% ==========

rr = Fr2*(Dr2*y2).*z2;
rs = Fs2*(Ds2*y2).*z2;
rt = Ft2*(Dt2*y2).*z2;
rx2 = Dt2*(rs) - Ds2*(rt); 
sx2 = Dr2*(rt) - Dt2*(rr); 
tx2 = Ds2*(rr) - Dr2*(rs);

rr = Fr2*(Dr2*x2).*z2;
rs = Fs2*(Ds2*x2).*z2;
rt = Ft2*(Dt2*x2).*z2;
ry2 = -(Dt2*(rs) - Ds2*(rt));
sy2 = -(Dr2*(rt) - Dt2*(rr));
ty2 = -(Ds2*(rr) - Dr2*(rs));

rr = Fr2*(Dr2*y2).*x2;
rs = Fs2*(Ds2*y2).*x2;
rt = Ft2*(Dt2*y2).*x2;
rz2 = -(Dt2*(rs) - Ds2*(rt));
sz2 = -(Dr2*(rt) - Dt2*(rr));
tz2 = -(Ds2*(rr) - Dr2*(rs));

% [~,rids] = find(abs(Dr2(ids2,:))>1e-8);
% [~,sids] = find(abs(Ds2(ids2,:))>1e-8);
% [~,tids] = find(abs(Dt2(ids2,:))>1e-8);
% plot3(r2,s2,t2,'o')
% hold on;
% rstids = rids;
% plot3(r2(rstids),s2(rstids),t2(rstids),'x')
% return

% plot3(r2(sids),s2(sids),t2(sids),'o')
% plot3(r2(tids),s2(tids),t2(tids),'o')

fprintf('GCL for elem 1 = %f, GCL for elem 2 = %f\n',...
    norm(Dr1*rx1 + Ds1*sx1 + Dt1*tx1,'fro'),norm(Dr2*rx2 + Ds2*sx2 + Dt2*tx2,'fro'))

e = ones((N+1)^2,1); zz = 0*e;
nrJ1 = e; nsJ1 = zz; ntJ1 = zz;

e = ones((N2+1)^2,1); zz = 0*e;
nrJ2 = -e; nsJ2 = zz; ntJ2 = zz;

nxJ1 = rx1(ids1).*nrJ1 + sx1(ids1).*nsJ1 + tx1(ids1).*ntJ1;
nyJ1 = ry1(ids1).*nrJ1 + sy1(ids1).*nsJ1 + ty1(ids1).*ntJ1;
nzJ1 = rz1(ids1).*nrJ1 + sz1(ids1).*nsJ1 + tz1(ids1).*ntJ1;

nxJ2 = rx2(ids2).*nrJ2 + sx2(ids2).*nsJ2 + tx2(ids2).*ntJ2;
nyJ2 = ry2(ids2).*nrJ2 + sy2(ids2).*nsJ2 + ty2(ids2).*ntJ2;
nzJ2 = rz2(ids2).*nrJ2 + sz2(ids2).*nsJ2 + tz2(ids2).*ntJ2;


% e = 1;plot3(x{e},y{e},z{e},'o')
% hold on
% e = 2;plot3(x{e},y{e},z{e},'x')

plot3(x1(ids1),y1(ids1),z1(ids1),'o','markersize',12)
hold on
plot3(V12*x2(ids2),V12*y2(ids2),V12*z2(ids2),'x','markersize',12)

quiver3(x1(ids1),y1(ids1),z1(ids1),nxJ1,nyJ1,nzJ1)
hold on;
quiver3(V12*x2(ids2),V12*y2(ids2),V12*z2(ids2),V12*nxJ2,V12*nyJ2,V12*nzJ2)
% quiver3(x2(ids2),y2(ids2),z2(ids2),nxJ2,nyJ2,nzJ2)

fprintf('normal errors: nx = %f, ny = %f, nz = %f\n',...
    norm(nxJ1+V12*nxJ2,'fro'),norm(nyJ1+V12*nyJ2,'fro'),norm(nzJ1+V12*nzJ2,'fro'))

% ry = -(Dt*((Ds*x).*z) - Ds*((Dt*x).*z));
norm(ry1(ids1) - V12*ry2(ids2),'fro')

% check quadrature for nxJ
% [rq1D wq1D] = JacobiGQ(0,0,N2);
% [rq1D wq1D] = JacobiGQ(0,0,N2);
% [wr ws] = meshgrid(wq1D); wq = wr(:).*ws(:);
% Iq = Vandermonde1D(N2,JacobiGQ(0,0,N2))/Vandermonde1D(N2,JacobiGL(0,0,N2));
% [wr ws] = meshgrid(w1D2); wGLLq = wr(:).*ws(:);
% wGLLq'*(nxJ2).^2
% wq'*(kron(Iq,Iq)*nxJ2).^2

% Dt1(ids1,ids1)
% V12*Dt2(ids2,ids2)
