function Wave2D_spectra

Globals2D
global bx by Vq Pq Vfq Pfq Vrq Vsq Prq Psq


% Polynomial order used for approximation
N = 3;

% Read in Mesh
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Grid/Other/squarereg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Grid/Other/squareireg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Grid/Other/block2.neu');

K1D = 2;
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);

% Initialize solver and construct grid and metric
StartUp2D;

% % rebuild maps for periodic
% BuildPeriodicMaps2D(2,2);

[re se] = EquiNodes2D(N); [re se] = xytors(re,se);
Ve = Vandermonde2D(N,re,se)/V;
xe = Ve*x; ye = Ve*y;

[rp sp] = EquiNodes2D(50); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

[rq sq wq] = Cubature2D(3*N);
% [rq sq wq] = QNodes2D(N); [rq sq] = xytors(rq,sq);
Vq = Vandermonde2D(N,rq,sq)/V;
[Vrq Vsq] = GradVandermonde2D(N,rq,sq);
Vrq = Vrq/V; Vsq = Vsq/V;
xq = Vq*x; yq = Vq*y;
M = Vq'*diag(wq)*Vq;
Pq = M\(Vq'*diag(wq));
Prq = M\(Vrq'*diag(wq));
Psq = M\(Vsq'*diag(wq));


[r1D, w1D] = JacobiGQ(0,0,2*N);
% rfq = [r1D; r1D; -ones(size(r1D))];
% sfq = [-ones(size(r1D)); -r1D; r1D];
wfq = [w1D; w1D; w1D];

Vfq1 = Vandermonde1D(N,r1D)/Vandermonde1D(N,r(Fmask(:,1)));
Vfq2 = Vandermonde1D(N,r1D)/Vandermonde1D(N,r(Fmask(:,2)));
Vfq3 = Vandermonde1D(N,r1D)/Vandermonde1D(N,s(Fmask(:,3)));
Vfq = blkdiag(Vfq1,Vfq2,Vfq3);
Pfq = (Vfq'*diag(wfq)*Vfq)\(Vfq'*diag(wfq));

%% build projection to CG

R = getCGRestriction()';
P = R*((R'*kron(diag(J(1,:)),M)*R)\(R'*kron(diag(J(1,:)),M)));
P = blkdiag(P,eye(Np*K*2));

%% eigs

rLGL = JacobiGQ(0,0,N); rmin = abs(rLGL(1)-rLGL(2));
dtscale = dtscale2D; dt = min(dtscale)*rmin*2/3

if 1
    e = zeros(3*Np*K,1);
    A = zeros(3*Np*K);
    for tau = [10]
        for i = 1:3*Np*K
            e(i) = 1;
            ids = 1:Np*K;
            p = reshape(e(ids),Np,K);
            u = reshape(e(ids + Np*K),Np,K);
            v = reshape(e(ids + 2*Np*K),Np,K);
            
            % dudt + dudx= 0
            [rhsp rhsu rhsv] = WaveRHS(p,u,v,tau);
            A(:,i) = [rhsp(:); rhsu(:); rhsv(:)];
            
            e(i) = 0;
            if mod(i,100)==0
                disp(sprintf('on column %d out of %d\n',i,Np*K))
            end
        end
        [W D] = eig(A);
        d = diag(D);
        plot(real(d),imag(d),'o')
        axis equal
%         xlim([-100 1])
        %axis([-100 1 -50 50])
        hold on
        title(sprintf('tau = %d\n',tau))
        drawnow
    end
    keyboard
    ids = find(abs(real(d)) < 1);
    d = d(ids);
    W = W(:,ids);
        
    for i = 1:length(ids)
        u = reshape(W(1:Np*K,i),Np,K);        
        v = reshape(W((1:Np*K) + (Np*K),i),Np,K);        
        w = reshape(W((1:Np*K) + 2*(Np*K),i),Np,K);        
        figure(1)
        clf
        vv = Vp*abs(u);
        color_line3(xp,yp,vv,vv,'.')
        view(3)
        
        figure(2)
        clf
        subplot(1,2,1)
        vv = Vp*abs(v);
        color_line3(xp,yp,vv,vv,'.')
        view(3)
        subplot(1,2,2)
        vv = Vp*abs(w);
        color_line3(xp,yp,vv,vv,'.')
        view(3)
        pause
    end       
    
%     [W2 D2] = eig(P*A);
%     d2 = diag(D2);
%     plot(d2,'x')
%     return
end

%%
FinalTime = 5.0;

time = 0;

% p = exp(-10*(x.^2 + y.^2));
k = 1; % frequency of solution
W = (2*k-1)/2*pi;
p = cos(W*x).*cos(W*y);

u = zeros(Np, K); v = zeros(Np, K);


% Runge-Kutta residual storage
resu = zeros(Np,K); resv = zeros(Np,K); resp = zeros(Np,K);

% compute time step size
CN = (N+1)^2/2; % guessing...
%dt = 1/(CN*max(abs(sJ(:)))*max(abs(1./J(:))));
dt = .1/(CN*max(Fscale(:)));

% outer time step loop
tstep = 0;
% figure
ids = 1:Np*K;
while (time<FinalTime)
    if(time+dt>FinalTime), dt = FinalTime-time; end
    
    for INTRK = 1:5
                
        [rhsp, rhsu, rhsv] = WaveRHS(p,u,v,tau);
        rhs = P*[rhsp(:); rhsu(:); rhsv(:)];
        a = 4/50;
        rhsp = (a)*rhsp + (1-a)*reshape(rhs(ids),Np,K);
        rhsu = (a)*rhsu + (1-a)*reshape(rhs(ids + Np*K),Np,K);
        rhsv = (a)*rhsv + (1-a)*reshape(rhs(ids + 2*Np*K),Np,K);
        
        % initiate and increment Runge-Kutta residuals
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        
        % update fields
        u = u+rk4b(INTRK)*resu;
        v = v+rk4b(INTRK)*resv;
        p = p+rk4b(INTRK)*resp;
        
    end;
    
    if 1 && nargin==0 && mod(tstep,10)==0
        clf
        vv = Vp*p;
        color_line3(xp,yp,vv,vv,'.');        
        axis([-1 1 -1 1 -1.25 1.25]);
%         PlotField2D(N+1, x, y, p); view(2)
        
        drawnow
    end
    
    % Increment time
    time = time+dt; tstep = tstep+1;
end

function rhsu = AdvecRHS(u,tau)

Globals2D
global bx by Vq Pq Vfq Pfq Vrq Vsq Prq Psq

du = zeros(Nfp*Nfaces,K);
bn = bx(vmapM).*nx(:) + by(vmapM).*ny(:);
bnf = reshape((bn - tau*abs(bn(:))),Nfp*Nfaces,K);
ujump = reshape(u(vmapP)-u(vmapM),Nfp*Nfaces,K);

du(:) = .5*Pfq*((Vfq*ujump).*(Vfq*bnf)); % project ujump/bnf
uq = Vq*u;
fx = Pq*((Vq*bx).*uq);
fy = Pq*((Vq*by).*uq);

% fx = bx.*u; fy = by.*u;
% du(:) = 0.5*(ujump).*bnf;

ux = rx.*(Dr*fx) + sx.*(Ds*fx);
uy = ry.*(Dr*fy) + sy.*(Ds*fy);
rhsu = -((ux + uy) + LIFT*(Fscale.*(du)));

function [rhsp, rhsu, rhsv] = WaveRHS(p,u,v,tau)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;

% Define field differences at faces
dp = zeros(Nfp*Nfaces,K); dp(:) = p(vmapP)-p(vmapM);
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapP)-u(vmapM);
dv = zeros(Nfp*Nfaces,K); dv(:) = v(vmapP)-v(vmapM);

% Impose reflective boundary conditions (p+ = -p-)
du(mapB) = 0; dv(mapB) = 0;
dp(mapB) = -2*p(vmapB);

% evaluate upwind fluxes
ndotdU = nx.*du + ny.*dv;
fluxp =  tau*dp - ndotdU;

fluxu =  (tau*ndotdU - dp).*nx;
fluxv =  (tau*ndotdU - dp).*ny;

pr = Dr*p; ps = Ds*p;
dpdx = rx.*pr + sx.*ps;
dpdy = ry.*pr + sy.*ps;
divU = Dr*(u.*rx + v.*ry) + Ds*(u.*sx + v.*sy);

% compute right hand sides of the PDE's
rhsp =  -divU + LIFT*(Fscale.*fluxp)/2.0;
rhsu =  -dpdx + LIFT*(Fscale.*fluxu)/2.0;
rhsv =  -dpdy + LIFT*(Fscale.*fluxv)/2.0;
