function Wave3D_IGA

Globals3D

NB = 4;
Ksub = 16;
K1D = 1;
smoothKnots = 0;
useQuadrature = 0;

if Ksub==1
    CN = (NB+1)^2/2;
    dt = 2/CN; % size of domain
else
    CN = (NB+2)*Ksub;
    dt = 2/CN;
end

N = NB+Ksub-1;

Np = (N+1)^3;
dofs = Np*K1D^3;

FinalTime = .53;

% return

%% mesh stuff

global rxJ sxJ txJ ryJ syJ tyJ rzJ szJ tzJ
global nx ny nz Jq sJ
global DVq1D Vq1D Pq1D D1D D1Dt LIFT Fmask

if Ksub==1
    [r1D w1D] = JacobiGL(0,0,N);
else
    VX = linspace(-1,1,Ksub+1);
    t = [VX(1)*ones(1,NB) VX VX(end)*ones(1,NB)]; % open knot vec
    for i = 1:N+1
        r1D(i,1) = mean(t((i+1):(i+NB))); % greville
    end
%     plot(r1D,r1D*0,'o')
end

[r s t] = meshgrid(r1D);
r = r(:); s = s(:); t = t(:);

Nfaces = 6;
Fmask1 = find(abs(r+1)<1e-8);
Fmask2 = find(abs(r-1)<1e-8);
Fmask3 = find(abs(s+1)<1e-8);
Fmask4 = find(abs(s-1)<1e-8);
Fmask5 = find(abs(t+1)<1e-8);
Fmask6 = find(abs(t-1)<1e-8);
Fmask = [Fmask1(:) Fmask2(:) Fmask3(:) Fmask4(:) Fmask5(:) Fmask6(:)];

rp1D = linspace(-1,1,100);
[rp sp tp] = meshgrid(rp1D);
rp = rp(:); sp = sp(:); tp = tp(:);

% initialize TP operators
[BVDM M1D D1D R rq1D wq1D Vq1D DVq1D] = bsplineVDM(NB,Ksub,r1D,smoothKnots); % VDM for interp, mass, M\S
Pq1D = M1D\(Vq1D'*diag(wq1D));
if Ksub >= NB
    Vq1D = sparse(Vq1D);
    Pq1D = sparse(Pq1D);
end

(N+1)
Nq = size(Vq1D,1)

D1Dt = M1D\(D1D'*M1D); % for skew-symmetric formulation

[rq sq tq] = meshgrid(rq1D); rq = rq(:); sq = sq(:); tq = tq(:);
[wrq wsq wtq] = meshgrid(wq1D); wq = wrq(:).*wsq(:).*wtq(:);

% invM1D = inv(M1D);

Mf1D = zeros(N+1,2);
Mf1D(1,1) = 1;
Mf1D(N+1,2) = 1;
LIFT1D = M1D\Mf1D;

% make lift operator
ids = 1:(N+1)^3;
ids = reshape(ids,N+1,N+1,N+1);

LIFT = sparse((N+1)^3,6*(N+1)^2); % initialize as spares matrix

% first/second face (r = +/- 1)
sk = 0;
for i = 1:(N+1)
    iid = squeeze(ids(i,:,:));
    for j = 1:(N+1)
        LIFT(iid(:,j),sk + (j-1)*(N+1) + i) = LIFT1D(:,1); % left-to-right
        LIFT(iid(:,j),sk + (N+1)^2 + (j-1)*(N+1) + i) = LIFT1D(:,2); % right-to-left
    end
end

% third/fourth face (s = +/- 1)
sk = sk + 2*(N+1)^2;
for i = 1:(N+1)
    iid = squeeze(ids(:,:,i));
    for j = 1:(N+1)
        LIFT(iid(:,j),sk + (i-1)*(N+1) + j) = LIFT1D(:,1); % left-to-right
        LIFT(iid(:,j),sk + (N+1)^2 + (i-1)*(N+1) + j) = LIFT1D(:,2); % right-to-left
    end
end

% third/fourth face (s = +/- 1)
sk = sk + 2*(N+1)^2;
for i = 1:(N+1)
    iid = squeeze(ids(:,i,:));
    for j = 1:(N+1)
        LIFT(iid(j,:),sk + (i-1)*(N+1) + j) = LIFT1D(:,1); % left-to-right
        LIFT(iid(j,:),sk + (N+1)^2 + (i-1)*(N+1) + j) = LIFT1D(:,2); % right-to-left
    end
end

%% non-affine mappings

K = 1; % assume one element

% no mapping yet
x = r; y = s; z = t;
% plot3(x,y,z,'o'); return

% change of variables factors
rxJ = zeros((N+1)^3,1); sxJ = zeros((N+1)^3,1); txJ = zeros((N+1)^3,1);
ryJ = zeros((N+1)^3,1); syJ = zeros((N+1)^3,1); tyJ = zeros((N+1)^3,1);
rzJ = zeros((N+1)^3,1); szJ = zeros((N+1)^3,1); tzJ = zeros((N+1)^3,1);

rxJ(:) = 1; sxJ(:) = 0; txJ(:) = 0;
ryJ(:) = 0; syJ(:) = 1; tyJ(:) = 0;
rzJ(:) = 0; szJ(:) = 0; tzJ(:) = 1;

J = ones(size(rxJ));

% normals
nx = zeros(size(Fmask(:)));
ny = zeros(size(Fmask(:)));
nz = zeros(size(Fmask(:)));

fids = 1:(N+1)^2;
nx(fids) = -1; ny(fids) = 0; nz(fids) = 0;
fids = fids + (N+1)^2;
nx(fids) = 1; ny(fids) = 0; nz(fids) = 0;
fids = fids + (N+1)^2;
nx(fids) = 0; ny(fids) = -1; nz(fids) = 0;
fids = fids + (N+1)^2;
nx(fids) = 0; ny(fids) = 1; nz(fids) = 0;
fids = fids + (N+1)^2;
nx(fids) = 0; ny(fids) = 0; nz(fids) = -1;
fids = fids + (N+1)^2;
nx(fids) = 0; ny(fids) = 0; nz(fids) = 1;

sJ = nx.^2 + ny.^2 + nz.^2;
nx = nx./sJ; ny = ny./sJ; nz = nz./sJ;

%% 
% plotting
Vp1D = bsplineVDM(NB,Ksub,rp1D,smoothKnots)/bsplineVDM(NB,Ksub,r1D,smoothKnots);
xp = matvec(Vp1D,x,'all');
yp = matvec(Vp1D,y,'all');
zp = matvec(Vp1D,z,'all');
Vp1D = bsplineVDM(NB,Ksub,rp1D,smoothKnots);

% make quadrature-based geofacs
Vq_nodal_1D = bsplineVDM(NB,Ksub,rq1D,smoothKnots)/bsplineVDM(NB,Ksub,r1D,smoothKnots);
xq = matvec(Vq_nodal_1D,x,'all');
yq = matvec(Vq_nodal_1D,y,'all');
zq = matvec(Vq_nodal_1D,z,'all');

if useQuadrature
    rxJ = matvec(Vq1D,rxJ,'all'); sxJ = matvec(Vq1D,sxJ,'all'); txJ = matvec(Vq1D,txJ,'all');
    ryJ = matvec(Vq1D,ryJ,'all'); syJ = matvec(Vq1D,syJ,'all'); tyJ = matvec(Vq1D,tyJ,'all');
    rzJ = matvec(Vq1D,rzJ,'all'); szJ = matvec(Vq1D,szJ,'all'); tzJ = matvec(Vq1D,tzJ,'all');
    Jq = matvec(Vq1D,J,'all');

    nx = reshape(nx,(N+1)^2,6); ny = reshape(ny,(N+1)^2,6); nz = reshape(nz,(N+1)^2,6);
    sJ = reshape(sJ,(N+1)^2,6);
    Nfq = size(Vq1D,1)^2;
    nxq = zeros(Nfq,6); nyq = zeros(Nfq,6); nzq = zeros(Nfq,6);
    for f = 1:6
        tmp = (Vq1D * reshape(nx(:,f),N+1,N+1))*Vq1D';
        nxq(:,f) = tmp(:);
        tmp = (Vq1D * reshape(ny(:,f),N+1,N+1))*Vq1D';
        nyq(:,f) = tmp(:);
        tmp = (Vq1D * reshape(nz(:,f),N+1,N+1))*Vq1D';
        nzq(:,f) = tmp(:);
        tmp = (Vq1D * reshape(sJ(:,f),N+1,N+1))*Vq1D';
        sJq(:,f) = tmp(:);
    end
    nx = nxq; ny = nyq; nz = nzq;
    sJ = sJq;
end

%% initial cond

k = 1;
pex = @(x,y,z,t) cos(k*pi*x/2).*cos(k*pi*y/2).*cos(k*pi*z/2).*cos(sqrt(3)*.5*k*pi*t);

% x0 = 0; y0 = 0; z0 = 0;
% pex = @(x,y,z,t) exp(-25*((x-x0).^2 + (y-y0).^2 + (z-z0).^2));

% p = matvec(inv(BVDM),pex(x,y,z,0),'all');
p = matvec(Pq1D,pex(xq,yq,zq,0),'all');
u = zeros(Np, 1);
v = zeros(Np, 1);
w = zeros(Np, 1);

% compute init cond error
[rq2 wq2] = JacobiGQ(0,0,2*N); % "fine" quadrature for error
[wq2r wq2s wq2t] = meshgrid(wq2);
wq2 = wq2r(:).*wq2s(:).*wq2t(:);
Vq2 = bsplineVDM(NB,Ksub,rq2,smoothKnots);
J2 = matvec(Vq2,J,'all');

xq2 = matvec(Vandermonde1D(N,rq2)/Vandermonde1D(N,r1D),x,'all');
yq2 = matvec(Vandermonde1D(N,rq2)/Vandermonde1D(N,r1D),y,'all');
zq2 = matvec(Vandermonde1D(N,rq2)/Vandermonde1D(N,r1D),z,'all');

init_cond_err = sqrt(sum(sum(wq2.*J2.*(matvec(Vq2,p,'all') - pex(xq2,yq2,zq2,0)).^2)))

% keyboard
% return
% pp = matvec(Vp1D,p,'all'); ids = yp > 0; h = color_line3(xp(ids),yp(ids),zp(ids),pp(ids),'.'); set(h,'markersize',32); return

%% test operators

if 0
    r2 = VDM\r; s2 = VDM\s; t2 = VDM\t;
    Dr = kron(kron(eye(N+1),D1D),eye(N+1));
    norm(Dr*r2 - matvec(D1D,r2,1))
    Ds = kron(kron(eye(N+1),eye(N+1)),D1D);
    norm(Ds*s2 - matvec(D1D,s2,2))
    Dt = kron(kron(D1D,eye(N+1)),eye(N+1));
    norm(Dt*t2 - matvec(D1D,t2,3))
    
    g = VDM\(r + 2*s + 3*t + r.*s + s.*t + r.*s.*t);
    norm(VDM*(Dr*g) - (1 + s + s.*t))
    norm(VDM*(Ds*g) - (2 + r + t + r.*t))
    norm(VDM*(Dt*g) - (3 + s + r.*s))
end


%% estimate max timestep

U = randn((N+1)^3,4);
for i = 1:10
    Uprev = U;
    if useQuadrature
        [rhsp, rhsu, rhsv, rhsw] = acousticsRHS3Dq(U(:,1),U(:,2),U(:,3),U(:,4));
    else
        [rhsp, rhsu, rhsv, rhsw] = acousticsRHS3D(U(:,1),U(:,2),U(:,3),U(:,4));
    end
    U(:,1) = rhsp;
    U(:,2) = rhsu;
    U(:,3) = rhsv;
    U(:,4) = rhsw;    
    
    lam = Uprev(:)'*U(:) / norm(Uprev(:));
    U = U/norm(U(:));    
end
dt = .95/abs(lam)
% return

if 0
    A = sparse((N+1)^3*4,(N+1)^3*4);
    U = zeros((N+1)^3,4);
    for i = 1:(N+1)^3*4
        U(i) = 1;
        [rhsp, rhsu, rhsv, rhsw] = acousticsRHS3D(U(:,1),U(:,2),U(:,3),U(:,4));
        A(:,i) = sparse([rhsp(:);rhsu(:);rhsv(:);rhsw(:)]);
        U(i) = 0;
        if (mod(i,100)==0)
            disp(sprintf('on i = %d out of %d\n',i,(N+1)^3*4))
        end
    end
    lam = eig(full(A));
    plot(lam,'o')
    title(sprintf('max real lam = %g\n',max(real(lam))))
    keyboard
end

%%

time = 0;

% Runge-Kutta residual storage
resp = zeros(Np,K);
resu = zeros(Np,K); 
resv = zeros(Np,K); 
resw = zeros(Np,K);

Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

% outer time step loop
ids = find(abs(yp)<2.5e-2);
for tstep = 1:Nsteps
    for INTRK = 1:5
        
        timelocal = tstep*dt + rk4c(INTRK)*dt;
        
        if useQuadrature
            [rhsp,rhsu,rhsv,rhsw] = acousticsRHS3Dq(p,u,v,w);
        else
            [rhsp,rhsu,rhsv,rhsw] = acousticsRHS3D(p,u,v,w);
        end
        
        % initiate and increment Runge-Kutta residuals
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        resw = rk4a(INTRK)*resw + dt*rhsw;
        
        % update fields
        u = u+rk4b(INTRK)*resu;
        v = v+rk4b(INTRK)*resv;
        w = w+rk4b(INTRK)*resw;
        p = p+rk4b(INTRK)*resp;
        
    end;
    
    if 1 && nargin==0 && (mod(tstep,10)==0 || tstep==Nsteps)
        clf
        vv = matvec(Vp1D,p,'all');
        color_line3(xp(ids),yp(ids),zp(ids),vv(ids),'.');
        axis equal; axis tight
        view(3)
        colorbar
        title(sprintf('time = %f',tstep*dt))
        drawnow
    end
    
    if mod(tstep,25)==0
        disp(sprintf('on tstep %d out of %d\n',tstep,Nsteps))
    end
end

% J = matvec(Vq1D,J,'all');
err = wq2.*J2.*(matvec(Vq2,p,'all') - pex(xq2,yq2,zq2,FinalTime)).^2;
L2err = sqrt(sum(err(:)));

L2err
% Nsteps
% L2err * Nsteps % error * steps (assume same overall work)

% work = ceil(FinalTime/dt)*(N+1)^4; % total work
% L2err * work


function [rhsp, rhsu, rhsv, rhsw] = acousticsRHS3D(p,u,v,w)

global rxJ sxJ txJ ryJ syJ tyJ rzJ szJ tzJ
global nx ny nz J sJ
global DVq1D Vq1D Pq1D D1D LIFT Fmask

% Define field differences at faces
dp = -2*p(Fmask(:)); 
ndotdU = zeros(size(dp)); %u(Fmask(:)).*nx + v(Fmask(:)).*ny + w(Fmask(:)).*nz; 

tau = 1;
fluxp =  tau*dp - ndotdU;
fluxu =  (tau*ndotdU - dp);

pr = matvec(D1D,p,1); 
ps = matvec(D1D,p,2); 
pt = matvec(D1D,p,3); 
dpdx = rxJ.*pr + sxJ.*ps + txJ.*pt;
dpdy = ryJ.*pr + syJ.*ps + tyJ.*pt;
dpdz = rzJ.*pr + szJ.*ps + tzJ.*pt;

ur = matvec(D1D,rxJ.*u + ryJ.*v + rzJ.*w,1); 
us = matvec(D1D,sxJ.*u + syJ.*v + szJ.*w,2); 
ut = matvec(D1D,txJ.*u + tyJ.*v + tzJ.*w,3); 
divU = ur + us + ut;

rhsp =  -divU + .5*LIFT*(sJ.*fluxp);
rhsu =  -dpdx + .5*LIFT*(sJ.*fluxu.*nx);
rhsv =  -dpdy + .5*LIFT*(sJ.*fluxu.*ny);
rhsw =  -dpdz + .5*LIFT*(sJ.*fluxu.*nz);

rhsp = rhsp./J;
rhsu = rhsu./J;
rhsv = rhsv./J;
rhsw = rhsw./J;

function [rhsp rhsu rhsv rhsw] = acousticsRHS3Dq(p,u,v,w)

global rxJ sxJ txJ ryJ syJ tyJ rzJ szJ tzJ
global nx ny nz Jq sJ
global DVq1D Vq1D Pq1D PDq1D D1D D1Dt LIFT Fmask

N = round((length(p))^(1/3))-1;
Nq = size(Vq1D,1);

pr = matvec(Vq1D,matvec(D1D,p,1),'all');
ps = matvec(Vq1D,matvec(D1D,p,2),'all');
pt = matvec(Vq1D,matvec(D1D,p,3),'all');
dpdx = matvec(Pq1D,rxJ.*pr + sxJ.*ps + txJ.*pt,'all');
dpdy = matvec(Pq1D,ryJ.*pr + syJ.*ps + tyJ.*pt,'all');
dpdz = matvec(Pq1D,rzJ.*pr + szJ.*ps + tzJ.*pt,'all');

uq = matvec(Vq1D,u,'all'); 
vq = matvec(Vq1D,v,'all'); 
wq = matvec(Vq1D,w,'all'); 
ur = matvec(Pq1D,rxJ.*uq + ryJ.*vq + rzJ.*wq,'all');
us = matvec(Pq1D,sxJ.*uq + syJ.*vq + szJ.*wq,'all');
ut = matvec(Pq1D,txJ.*uq + tyJ.*vq + tzJ.*wq,'all');
ur = matvec(D1Dt,ur,1);
us = matvec(D1Dt,us,2);
ut = matvec(D1Dt,ut,3);
divU = -(ur + us + ut);

% Define field differences at faces
tau = 1;
dp = -2*p(Fmask); 
uf = reshape(u(Fmask(:)),(N+1)^2,6);
vf = reshape(v(Fmask(:)),(N+1)^2,6);
wf = reshape(w(Fmask(:)),(N+1)^2,6);
fluxp = zeros((N+1)^2,6);
fluxu = zeros((N+1)^2,6);
fluxv = zeros((N+1)^2,6);
fluxw = zeros((N+1)^2,6);
for f = 1:6
    dpfq = (Vq1D*reshape(dp(:,f),N+1,N+1))*Vq1D';
    ndotdU = zeros(size(dpfq)); 
    ufq = (Vq1D*reshape(uf(:,f),N+1,N+1))*Vq1D';
    vfq = (Vq1D*reshape(vf(:,f),N+1,N+1))*Vq1D';
    wfq = (Vq1D*reshape(wf(:,f),N+1,N+1))*Vq1D';      
    ndotUavg = nx(:,f).*ufq(:) + ny(:,f).*vfq(:) + nz(:,f).*wfq(:);    
%     ndotUavg = 0*ndotUavg; % fix when using skew
    fluxp(:,f) = reshape((Pq1D*reshape((tau*dpfq(:) - 2*ndotUavg).*sJ(:,f),Nq,Nq))*Pq1D',(N+1)^2,1);
    fU = (tau*ndotdU(:) - dpfq(:)).*sJ(:,f);
    fluxu(:,f) = reshape((Pq1D*reshape(fU.*nx(:,f),Nq,Nq))*Pq1D',(N+1)^2,1);
    fluxv(:,f) = reshape((Pq1D*reshape(fU.*ny(:,f),Nq,Nq))*Pq1D',(N+1)^2,1);
    fluxw(:,f) = reshape((Pq1D*reshape(fU.*nz(:,f),Nq,Nq))*Pq1D',(N+1)^2,1);    
end

rhsp =  -divU + .5*LIFT*(fluxp(:));
rhsu =  -dpdx + .5*LIFT*(fluxu(:));
rhsv =  -dpdy + .5*LIFT*(fluxv(:));
rhsw =  -dpdz + .5*LIFT*(fluxw(:));

rhsp = matvec(Pq1D,matvec(Vq1D,rhsp,'all')./Jq,'all');
rhsu = matvec(Pq1D,matvec(Vq1D,rhsu,'all')./Jq,'all');
rhsv = matvec(Pq1D,matvec(Vq1D,rhsv,'all')./Jq,'all');
rhsw = matvec(Pq1D,matvec(Vq1D,rhsw,'all')./Jq,'all');


% apply 1D tensor product matvec
function b = matvec(A,x,dim)

N = round((length(x(:)))^(1/3));
M = size(A,1);
x = reshape(x,N,N,N);

if strcmp(dim,'all')
    % apply in all dimensions
    
    b1 = zeros(M,M,N); % apply xy
    b = zeros(M,M,M); % apply z to xy result
    for i = 1:N
        b1(:,:,i) = (A*reshape(x(:,:,i),N,N))*A';
    end
    for i = 1:M
        b(i,:,:) = reshape(b1(i,:,:),M,N)*A';
    end
    
elseif dim==1  % apply only to one dimension
    b = zeros(M,N,N); % apply z to xy result
    for i = 1:N
        b(:,:,i) = reshape(x(:,:,i),N,N)*A';
    end
elseif dim==2
    b = zeros(N,M,N);
    for i = 1:N
        b(:,:,i) = A*reshape(x(:,:,i),N,N);
    end
elseif dim==3
    b = zeros(N,N,M);
    for i = 1:N
        b(i,:,:) = reshape(x(i,:,:),N,N)*A';
    end   
end

b = b(:);

return

