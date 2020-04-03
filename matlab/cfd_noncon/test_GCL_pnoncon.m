clear
N = 5;
N2 = 6;
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

V21 = kron(Vandermonde1D(N,r1D2),Vandermonde1D(N,r1D2)) / kron(V1D1,V1D1);

a = .1/N^2;
x1 = x1 + a*randn(size(rr));
y1 = y1 + a*randn(size(rr));
z1 = z1 + a*randn(size(rr));
x2 = x2 + a*randn(size(rs));
y2 = y2 + a*randn(size(rs));
z2 = z2 + a*randn(size(rs));

x2(ids2) = V21*x1(ids1);
y2(ids2) = V21*y1(ids1);
z2(ids2) = V21*z1(ids1);

% enforce continuity: N1 = N nodes set by N2 = N-1 nodes.
rr1 = (Dr1*y1).*z1;
rs1 = (Ds1*y1).*z1;
rt1 = (Dt1*y1).*z1;
rr2 = (Dr2*y2).*z2;
rs2 = (Ds2*y2).*z2;
rt2 = (Dt2*y2).*z2;
rr2(ids2) = V21*rr1(ids1);
rs2(ids2) = V21*rs1(ids1);
rt2(ids2) = V21*rt1(ids1);
rx1 = Dt1*(rs1) - Ds1*(rt1); 
sx1 = Dr1*(rt1) - Dt1*(rr1);
tx1 = Ds1*(rr1) - Dr1*(rs1);
rx2 = Dt2*(rs2) - Ds2*(rt2); 
sx2 = Dr2*(rt2) - Dt2*(rr2);
tx2 = Ds2*(rr2) - Dr2*(rs2);

rr1 = (Dr1*x1).*z1;
rs1 = (Ds1*x1).*z1;
rt1 = (Dt1*x1).*z1;
rr2 = (Dr2*x2).*z2;
rs2 = (Ds2*x2).*z2;
rt2 = (Dt2*x2).*z2;
rr2(ids2) = V21*rr1(ids1);
rs2(ids2) = V21*rs1(ids1);
rt2(ids2) = V21*rt1(ids1);
ry1 = Dt1*(rs1) - Ds1*(rt1); 
sy1 = Dr1*(rt1) - Dt1*(rr1);
ty1 = Ds1*(rr1) - Dr1*(rs1);
ry2 = Dt2*(rs2) - Ds2*(rt2); 
sy2 = Dr2*(rt2) - Dt2*(rr2);
ty2 = Ds2*(rr2) - Dr2*(rs2);

rr1 = (Dr1*y1).*x1;
rs1 = (Ds1*y1).*x1;
rt1 = (Dt1*y1).*x1;
rr2 = (Dr2*y2).*x2;
rs2 = (Ds2*y2).*x2;
rt2 = (Dt2*y2).*x2;
rr2(ids2) = V21*rr1(ids1);
rs2(ids2) = V21*rs1(ids1);
rt2(ids2) = V21*rt1(ids1);
rz1 = Dt1*(rs1) - Ds1*(rt1); 
sz1 = Dr1*(rt1) - Dt1*(rr1);
tz1 = Ds1*(rr1) - Dr1*(rs1);
rz2 = Dt2*(rs2) - Ds2*(rt2); 
sz2 = Dr2*(rt2) - Dt2*(rr2);
tz2 = Ds2*(rr2) - Dr2*(rs2);

fprintf('GCL for elem 1 = %g, GCL for elem 2 = %g\n',...
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

plot3(V21*x1(ids1),V21*y1(ids1),V21*z1(ids1),'o','markersize',12)
hold on
plot3(x2(ids2),y2(ids2),z2(ids2),'x','markersize',12)

quiver3(V21*x1(ids1),V21*y1(ids1),V21*z1(ids1),V21*nxJ1,V21*nyJ1,V21*nzJ1)
hold on;
quiver3(x2(ids2),y2(ids2),z2(ids2),nxJ2,nyJ2,nzJ2)
% quiver3(x2(ids2),y2(ids2),z2(ids2),nxJ2,nyJ2,nzJ2)

fprintf('normal errors: nx = %g, ny = %g, nz = %g\n',...
    norm(V21*nxJ1+nxJ2,'fro'),norm(V21*nyJ1+nyJ2,'fro'),norm(V21*nzJ1+nzJ2,'fro'))

