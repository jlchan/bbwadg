% two-element test of GCL for non-conforming interfaces

clear
N = 4;
r1D = JacobiGL(0,0,N);
D = GradVandermonde1D(N,r1D)/Vandermonde1D(N,r1D);
I = eye(N+1);
Dr = kron(I,kron(D,I));
Ds = kron(I,kron(I,D));
Dt = kron(D,kron(I,I));
[r,s,t] = meshgrid(r1D);
r = r(:); s = s(:); t = t(:);

e = ones((N+1)^2,1);
z = 0*e;
nrJ = [-e;e;z;z;z;z];
nsJ = [z;z;-e;e;z;z];
ntJ = [z;z;z;z;-e;e];

% fine nodes
rp1D = linspace(-1,1,25);
[rp,sp,tp] = meshgrid(rp1D);
rp = rp(:); sp = sp(:); tp = tp(:);
Vp1D = Vandermonde1D(N,rp1D)/Vandermonde1D(N,r1D);
Vp = kron(Vp1D,kron(Vp1D,Vp1D));

% find face nodes on left/right of each hex
fid1 = find(abs(1-r)<1e-10);
fid2 = find(abs(1+r)<1e-10);

rf = r(fid1); sf = s(fid1); tf = t(fid1);
Vf{1} = hex_basis(N,.5*(1+rf),.5*(1+sf),.5*(1+tf)-1)/hex_basis(N,r,s,t);
% Vf{1} = hex_basis(N,.5*(1+rf),.5*(1+sf),tf)/hex_basis(N,r,s,t);
rf = r(fid2); sf = s(fid2); tf = t(fid2);
Vf{2} = hex_basis(N,rf,sf,tf)/hex_basis(N,r,s,t);
    
% make curved mapping
a = .25;
x1 = [r 1+.5*(1+r)];
y1 = [s .5*(1+s)];
z1 = [t .5*(1+t)-1];
% z1 = [t t]; 
dx = sin(x1).*sin(y1).*sin(z1);
x = x1 + a*dx;
y = y1 + a*dx;
z = z1 + a*dx;

% compute geometric terms
Fr = (Dr*y).*z;
Fs = (Ds*y).*z;
Ft = (Dt*y).*z;
% .5 = scaling factor for each direction under isotropic refinement
Fr(fid2,2) = .5*Vf{1}*Fr(:,1);
Fs(fid2,2) = .5*Vf{1}*Fs(:,1);
Ft(fid2,2) = .5*Vf{1}*Ft(:,1);
rxJ = Dt*(Fs) - Ds*(Ft);
sxJ = Dr*(Ft) - Dt*(Fr);
txJ = Ds*(Fr) - Dr*(Fs);

Fr = (Dr*x).*z;
Fs = (Ds*x).*z;
Ft = (Dt*x).*z;
Fr(fid2,2) = .5*Vf{1}*Fr(:,1);
Fs(fid2,2) = .5*Vf{1}*Fs(:,1);
Ft(fid2,2) = .5*Vf{1}*Ft(:,1);
ryJ = -(Dt*(Fs) - Ds*(Ft));
syJ = -(Dr*(Ft) - Dt*(Fr));
tyJ = -(Ds*(Fr) - Dr*(Fs));

Fr = (Dr*y).*x;
Fs = (Ds*y).*x;
Ft = (Dt*y).*x;
Fr(fid2,2) = .5*Vf{1}*Fr(:,1);
Fs(fid2,2) = .5*Vf{1}*Fs(:,1);
Ft(fid2,2) = .5*Vf{1}*Ft(:,1);
rzJ = -(Dt*(Fs) - Ds*(Ft));
szJ = -(Dr*(Ft) - Dt*(Fr));
tzJ = -(Ds*(Fr) - Dr*(Fs));

plot3(x,y,z,'.')
hold on
for e = 1:2
    scl = (e==1) - (e==2); % scaling factor = -1 or 1
    nxJ(:,e) = scl*(Vf{e}*rxJ(:,e)); 
    nyJ(:,e) = scl*(Vf{e}*ryJ(:,e)); 
    nzJ(:,e) = scl*(Vf{e}*rzJ(:,e)); 
    plot3(Vf{e}*x(:,e),Vf{e}*y(:,e),Vf{e}*z(:,e),'x')
    quiver3(Vf{e}*x(:,e),Vf{e}*y(:,e),Vf{e}*z(:,e),nxJ(:,e),nyJ(:,e),nzJ(:,e))
end

% factor of 4 = isotropic refinement results in a mortar face that's 4x smaller 
nx_err = norm(sum(nxJ*diag([1 4]),2))
ny_err = norm(sum(nyJ*diag([1 4]),2))
nz_err = norm(sum(nzJ*diag([1 4]),2))
gcl_err = norm(Dr*rxJ+Ds*sxJ+Dt*txJ,'fro')+ norm(Dr*ryJ+Ds*syJ+Dt*tyJ,'fro')+ norm(Dr*rzJ+Ds*szJ+Dt*tzJ,'fro')
