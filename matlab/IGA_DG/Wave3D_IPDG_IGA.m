function Wave3D_IGA

Globals3D

NB = 1;
Ksub = 2;
K1D = 1;
smoothKnots = 0;
useQuadrature = 1;

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
global DPq1D DVq1D Vq1D Pq1D D1D D1Dt LIFT Fmask M1D

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
DPq1D = M1D\(DVq1D'*diag(wq1D));
if Ksub >= NB
    Vq1D = sparse(Vq1D);
    Pq1D = sparse(Pq1D);
    DVq1D = sparse(DVq1D);
    DPq1D = sparse(DPq1D);
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

global Grr Grs Grt Gss Gst Gtt
Grr = (rxJ.^2 + ryJ.^2 + rzJ.^2)./(Jq.^2);
Grs = (rxJ.*sxJ + ryJ.*syJ + rzJ.*szJ)./(Jq.^2);
Grt = (rxJ.*txJ + ryJ.*tyJ + rzJ.*tzJ)./(Jq.^2);
Gss = (sxJ.^2 + syJ.^2 + szJ.^2)./(Jq.^2);
Gst = (sxJ.*txJ + syJ.*tyJ + szJ.*tzJ)./(Jq.^2);
Gtt = (txJ.^2 + tyJ.^2 + tzJ.^2)./(Jq.^2);

global S1D
S1D = M1D\(DVq1D'*diag(wq1D)*DVq1D);
% keyboard

%% initial cond

k = 1;
%pex = @(x,y,z,t) cos(k*pi*x/2).*cos(k*pi*y/2).*cos(k*pi*z/2).*cos(sqrt(3)*.5*k*pi*t);
pex = @(x,y,z,t) sin(k*pi*x/2).*sin(k*pi*y/2).*sin(k*pi*z/2).*cos(sqrt(3)*.5*k*pi*t);

% x0 = 0; y0 = 0; z0 = 0;
% pex = @(x,y,z,t) exp(-25*((x-x0).^2 + (y-y0).^2 + (z-z0).^2));

% pex = @(x,y,z,t) x + y + z;

% p = matvec(inv(BVDM),pex(x,y,z,0),'all');
p = matvec(Pq1D,pex(xq,yq,zq,0),'all');
u = zeros(Np, 1);

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
% pp = matvec(Vp1D,p,'all'); ids = yp > 0; color_line3(xp(ids),yp(ids),zp(ids),pp(ids),'.'); return


%% estimate max timestep

if 0
    U = randn((N+1)^3,2);
    for i = 1:20
        Uprev = U;
        if useQuadrature
            [rhsp, rhsu] = acousticsRHS3Dq(U(:,1),U(:,2));
        else
            [rhsp, rhsu] = acousticsRHS3D(U(:,1),U(:,2));
        end
        U(:,1) = rhsp;
        U(:,2) = rhsu;
        
        lam = Uprev(:)'*U(:) / norm(Uprev(:));
        U = U/norm(U(:));
%         abs(lam)
    end
    dt = .95/abs(lam)
    % return
end
dt = .25/(NB*max(NB,Ksub)*K1D);

if 0 && ((N+1)^3*2)<1000
    A = sparse((N+1)^3*2,(N+1)^3*2);
    U = zeros((N+1)^3,2);
    for i = 1:(N+1)^3*2
        U(i) = 1;
        [rhsp, rhsu] = acousticsRHS3Dq(U(:,1),U(:,2));
        A(:,i) = sparse([rhsp(:);rhsu(:)]);
        U(i) = 0;
        if (mod(i,100)==0)
            disp(sprintf('on i = %d out of %d\n',i,(N+1)^3*2))
        end
    end
    A = full(A((N+1)^3 + (1:(N+1)^3),1:(N+1)^3));
    lam = eig(A,kron(kron(M1D,M1D),M1D));
    plot(lam,'o')
    title(sprintf('max real lam = %g\n',max(real(lam))))
    keyboard
    
%     A = full(A((N+1)^3 + (1:(N+1)^3),1:(N+1)^3))
    
end

%%

time = 0;

% Runge-Kutta residual storage
resp = zeros(Np,K);
resu = zeros(Np,K); 

Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

% outer time step loop
ids = find(abs(yp) < 2.5e-2);

for tstep = 1:Nsteps
    for INTRK = 1:5
        
        timelocal = tstep*dt + rk4c(INTRK)*dt;
        
        [rhsp,rhsu] = acousticsRHS3Dq(p,u);
        
        % initiate and increment Runge-Kutta residuals
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        
        % update fields
        p = p+rk4b(INTRK)*resp;
        u = u+rk4b(INTRK)*resu;
        
    end;
    
    if 1 && nargin==0 && (mod(tstep,10)==0 || tstep==Nsteps)
        clf
        vv = matvec(Vp1D,p,'all');
%         vv = pex(xp,yp,zp,tstep*dt);
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



function [rhsp rhsu] = acousticsRHS3Dq(p,u)

global rxJ sxJ txJ ryJ syJ tyJ rzJ szJ tzJ
global nx ny nz Jq sJ
global PDq1D DVq1D Vq1D Pq1D D1D D1Dt LIFT Fmask M1D

N = round((length(p))^(1/3))-1;

if 0
    pr = matvec(Vq1D,matvec(D1D,p,1),'all');
    ps = matvec(Vq1D,matvec(D1D,p,2),'all');
    pt = matvec(Vq1D,matvec(D1D,p,3),'all');
    
    global Grr Grs Grt Gss Gst Gtt
    dp1 = matvec(Pq1D,(Grr.*pr + Grs.*ps + Grt.*pt).*Jq,'all');
    dp2 = matvec(Pq1D,(Grs.*pr + Gss.*ps + Gst.*pt).*Jq,'all');
    dp3 = matvec(Pq1D,(Grt.*pr + Gst.*ps + Gtt.*pt).*Jq,'all');
    
    rhsu = -(matvec(D1Dt,dp1,1) + matvec(D1Dt,dp2,2) + matvec(D1Dt,dp3,3));
    % rhsu = rhsu + ; % add symmetric boundary contribution 
    rhsu = matvec(Pq1D,matvec(Vq1D,rhsu,'all')./Jq,'all');

    rhsp = u;           
    rhsp(Fmask(:),:) = 0;    
    
else    
    
    global S1D
%     p(Fmask(:),:) = 0;
    rhsu = -(matvec(S1D,p,1) + matvec(S1D,p,2) + matvec(S1D,p,3));
%     rhsu(Fmask(:),:) = 0;
    rhsp = u;
end



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

