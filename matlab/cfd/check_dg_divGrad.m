function Advec2D_cfl

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

% rebuild maps for periodic
BuildPeriodicMaps2D(2,2);

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

% % Set initial conditions
cx = .1; cy = .1;
D = (x-cx).^2 + (y-cy).^2;
u = exp(-D*5^2).*(1-x.^2).*(1-y.^2);
u = sin(pi*x).*sin(pi*y);

bx = ones(size(x));
by = .5*ones(size(x));
bx = -y;
by = x;


% plot3(x,y,-pi*cos(pi*x).*sin(pi*y),'o');return

% divB = rx.*(Dr*bx) + sx.*(Ds*bx) + ry.*(Dr*by) + sy.*(Ds*by);

FinalTime = 2;

%% eigs

rLGL = JacobiGQ(0,0,N); rmin = abs(rLGL(1)-rLGL(2));
dtscale = dtscale2D; dt = min(dtscale)*rmin*2/3

if 1
    e = zeros(Np*K,1);
    Grad = zeros(2*Np*K,Np*K);
    for i = 1:Np*K
        e(i) = 1;
        ids = 1:Np*K;
        u = reshape(e(ids),Np,K);
                
        [dudx dudy] = grad(u);
        
        Grad(:,i) = [dudx(:); dudy(:)];
        e(i) = 0;
        if mod(i,100)==0
            disp(sprintf('on column %d out of %d\n',i,Np*K))
        end
    end
    
    Div = zeros(Np*K,2*Np*K);
    e = zeros(Np*K,2);
    for i = 1:2*Np*K
        e(i) = 1;                                
        divU = div(reshape(e(:,1),Np,K),reshape(e(:,2),Np,K));
        Div(:,i) = divU(:);
        e(i) = 0;
        if mod(i,100)==0
            disp(sprintf('on column %d out of %d\n',i,Np*K))
        end
    end
        
    M = kron(diag(J(1,:)),inv(V*V'));
    MM = blkdiag(M,M);
    B = [eye(Np*K);eye(Np*K)];
    
    norm(M*Div + Grad'*blkdiag(M,M),'fro')
    
    keyboard
    
    return
    
end

%%
time = 0;

% Runge-Kutta residual storage
resu = zeros(Np,K);

% compute time step size
rLGL = JacobiGQ(0,0,N); rmin = abs(rLGL(1)-rLGL(2));
dtscale = dtscale2D; dt = min(dtscale)*rmin*2/3

% outer time step loop
while (time<FinalTime)
    
    if(time+dt>FinalTime), dt = FinalTime-time; end
    
    for INTRK = 1:5
        
        rhsu = AdvecRHS(u);
        
        % initiate and increment Runge-Kutta residuals
        resu = rk4a(INTRK)*resu + dt*rhsu;
        
        % update fields
        u = u+rk4b(INTRK)*resu;
    end
    
    clf;
    up = Vp*u;
    color_line3(xp,yp,up,up,'.');drawnow
    
    % Increment time
    time = time+dt;
end
return

function [dudx dudy] = grad(u)

Globals2D

du = zeros(Nfp*Nfaces,K); 
du(:) = u(vmapP)-u(vmapM);

pr = Dr*u; ps = Ds*u;
dpdx = rx.*pr + sx.*ps;
dpdy = ry.*pr + sy.*ps;

dudx =  dpdx + LIFT*(Fscale.*du.*nx)/2.0;
dudy =  dpdy + LIFT*(Fscale.*du.*ny)/2.0;

function divU = div(u,v)

Globals2D

du = reshape(u(vmapP)-u(vmapM),Nfp*Nfaces,K);
dv = reshape(v(vmapP)-v(vmapM),Nfp*Nfaces,K);

dUn = nx.*du + ny.*dv;

dudx = rx.*(Dr*u) + sx.*(Ds*u);
dudy = ry.*(Dr*v) + sy.*(Ds*v);
divU =  dudx + dudy + .5*LIFT*(Fscale.*dUn);

