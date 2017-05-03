function Advec2D_cfl

Globals2D
global bx by Vq Pq Vfq Pfq Vrq Vsq Prq Psq M


% Polynomial order used for approximation
N = 4;

% Read in Mesh
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Grid/Other/squarereg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Grid/Other/squareireg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Grid/Other/block2.neu');

K1D = 8;
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
% VX = (VX+1)/2;
% VY = (VY+1)/2;

FinalTime = 2*pi;

% Initialize solver and construct grid and metric
StartUp2D;

% rebuild maps for periodic
BuildPeriodicMaps2D(2,2);

[rp sp] = EquiNodes2D(25); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

Nq = 2*N+1;
[rq sq wq] = Cubature2D(Nq);
Vq = Vandermonde2D(N,rq,sq)/V;

[Vrq Vsq] = GradVandermonde2D(N,rq,sq);
Vrq = Vrq/V; Vsq = Vsq/V;
xq = Vq*x; yq = Vq*y;
M = Vq'*diag(wq)*Vq;
Pq = M\(Vq'*diag(wq));
Prq = M\(Vrq'*diag(wq));
Psq = M\(Vsq'*diag(wq));

[r1D, w1D] = JacobiGQ(0,0,ceil(Nq/2));
% rfq = [r1D; r1D; -ones(size(r1D))];
% sfq = [-ones(size(r1D)); -r1D; r1D];
wfq = [w1D; w1D; w1D];

Vfq1 = Vandermonde1D(N,r1D)/Vandermonde1D(N,r(Fmask(:,1)));
Vfq2 = Vandermonde1D(N,r1D)/Vandermonde1D(N,r(Fmask(:,2)));
Vfq3 = Vandermonde1D(N,r1D)/Vandermonde1D(N,s(Fmask(:,3)));
Vfq = blkdiag(Vfq1,Vfq2,Vfq3);
Pfq = (Vfq'*diag(wfq)*Vfq)\(Vfq'*diag(wfq));

% % Set initial conditions
u = sin(pi*x).*sin(pi*y);
a = 100;
x0 = 0; y0 = -.25;
u0 = @(x,y) exp(-a*((x-x0).^2+(y-y0).^2));

r2 = @(x,y) x.^2 + y.^2;
rt = @(x,y) 2.5*sqrt(r2(x,y-1/3));
u0 = @(x,y) ((1+ cos(pi*rt(x,y)))/2).^2 .* (rt(x,y)<=1);
u0 = @(x,y) (rt(x,y)<=.5);

global bxex byex xq yq
xq = Vq*x; yq = Vq*y;
bxex = @(x,y,t) ones(size(x)); byex = @(x,y,t) .5*ones(size(x));
bxex = @(x,y,t) -y; byex = @(x,y,t) x;
% bxex = @(x,y,t) -y*cos(t); byex = @(x,y,t) x*cos(t);

% bxex = @(x,y,t) 1 + .75*sin(4*pi*y); byex = @(x,y) 1 + .75*sin(4*pi*x);

% bx = rand(Np,K); by = rand(Np,K);
% bx = bxex(x,y); by = byex(x,y); % continuous

% divB = rx.*(Dr*bx) + sx.*(Ds*bx) + ry.*(Dr*by) + sy.*(Ds*by);

u = Pq*u0(xq,yq);

% color_line3(xp,yp,Vp*bx,Vp*bx,'.');return


%% eigs

if 0
    e = zeros(Np*K,1);
    A = zeros(Np*K);
    for i = 1:Np*K
        e(i) = 1;
        ids = 1:Np*K;
        u = reshape(e(ids),Np,K);
        
        % dudt + dudx= 0
        rhsu = AdvecRHS(u);
        A(:,i) = rhsu(:);
        
        e(i) = 0;
        if mod(i,100)==0
            disp(sprintf('on column %d out of %d\n',i,Np*K))
        end
    end
    [W D] = eig(A);
    lam = diag(D);
    clf
    plot(lam,'x')    
    %         axis([-100 1 -50 50])
    %         hold on
    M = kron(diag(J(1,:)),inv(V*V'));
    MA = M*A;
    norm(MA+MA','fro') % check skew
    max(real(lam))
    keyboard
    
end

%%
time = 0;

% Runge-Kutta residual storage
resu = zeros(Np,K);

% compute time step size
CN = (N+1)*(N+2)/2;
dt = 2./(CN*max(Fscale(:)));

Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

for tstep = 1:Nsteps
        
    for INTRK = 1:5        
        timelocal = tstep*dt + rk4c(INTRK)*dt;
        
        rhsu = AdvecRHS(u,timelocal);        
        resu = rk4a(INTRK)*resu + dt*rhsu;
        u = u+rk4b(INTRK)*resu;
    end
        
    if mod(tstep,10)==0
        clf;
        up = Vp*u;
        color_line3(xp,yp,up,up,'.');drawnow
    end
end

wJq = diag(wq)*(Vq*J);
err = wJq.*(Vq*u - u0(xq,yq)).^2;
L2err = sqrt(sum(err(:)))
return

function rhsu = AdvecRHS(u,t)

Globals2D
global bx by Vq Pq Vfq Pfq Vrq Vsq Prq Psq M

global bxex byex xq yq

bx = Pq*(bxex(xq,yq,t));
by = Pq*(byex(xq,yq,t));

opt = 1;
if opt==1
 
    % compute Div_h(P_N(beta*u))
    uq = Vq*u;
    
    fx = Pq*((Vq*bx).*uq);
    fy = Pq*((Vq*by).*uq);
    dfdx = rx.*(Dr*fx) + sx.*(Ds*fx);
    dfdy = ry.*(Dr*fy) + sy.*(Ds*fy);
    du = reshape(u(vmapP)-u(vmapM),Nfp*Nfaces,K);    
    dfx = reshape(fx(vmapP)-fx(vmapM),Nfp*Nfaces,K);
    dfy = reshape(fy(vmapP)-fy(vmapM),Nfp*Nfaces,K);
    dfn = dfx.*nx + dfy.*ny;

    % upwinding
    bn = reshape(bx(vmapM).*nx(:) + by(vmapM).*ny(:),Nfp*Nfaces,K);    
    dfn = dfn + reshape(Pfq*((Vfq*du).*abs(Vfq*bn)),Nfp*Nfaces,K);    
    rhs1 = dfdx + dfdy + .5*LIFT*(Fscale.*dfn); 
    
    % compute P_N(beta * Grad_h(u))
    ur = Dr*u; us = Ds*u;
    ux = rx.*ur + sx.*us;
    uy = ry.*ur + sy.*us;    
    uxq = ux + .5*LIFT*(Fscale.*du.*nx);
    uyq = uy + .5*LIFT*(Fscale.*du.*ny);
    
    % local 
    rhs2 = Pq*((Vq*bx).*(Vq*uxq)) + Pq*((Vq*by).*(Vq*uyq));
    
    rhsu = .5*(rhs1 + rhs2);
%     rhsu = rhs1;

elseif opt==2 % kopriva + gassner split form
    
%     Vq = eye(Np);
%     Pq = eye(Np);
%     Vfq = eye(Nfp*Nfaces);
%     Pfq = Vfq;

    Prq = V*V' * Dr'* M;
    Psq = V*V' * Ds'* M;
    
    % -.5*((u, beta*Grad(v)) + (u, Div(P_N(beta*v))))    
    % -.5*(P_N(beta*u), Grad(v) + .5*(P_N(beta*Grad(u)),v)))  - .5*<b_n*u,v>
    
    uq = Vq*u;
    bxu = Pq*((Vq*bx).*uq);
    byu = Pq*((Vq*by).*uq);
    
    vol1 = Prq*(rx.*bxu) + Psq*(sx.*bxu) + Prq*(ry.*byu) + Psq*(sy.*byu);    
    ux = rx.*(Dr*u) + sx.*(Ds*u);
    uy = ry.*(Dr*u) + sy.*(Ds*u);
    vol2 = Pq*((Vq*bx).*(Vq*ux) + (Vq*by).*(Vq*uy));
    vol = .5*(-vol1 + vol2);
    
    uavg = reshape(u(vmapP) + u(vmapM),Nfp*Nfaces,K);
    bxavg = .5*Vfq*reshape(bx(vmapP) + bx(vmapM),Nfp*Nfaces,K);
    byavg = .5*Vfq*reshape(by(vmapP) + by(vmapM),Nfp*Nfaces,K);
    afn = Pfq*((Vfq*nx).*(bxavg).*(Vfq*uavg) + (Vfq*ny).*(byavg).*(Vfq*uavg));
    afn = afn - Pfq*(((Vfq*nx).*(bxavg) + (Vfq*ny).*(byavg)).*(Vfq*u(Fmask(:),:)));
    rhsu = vol + .5*LIFT*(Fscale.*afn);
    
    % upwinding 
    bn = reshape(bx(vmapM).*nx(:) + by(vmapM).*ny(:),Nfp*Nfaces,K);
    ujump = reshape(u(vmapP)-u(vmapM),Nfp*Nfaces,K);
    du = -reshape(Pfq*((Vfq*ujump).*abs(Vfq*bn)),Nfp*Nfaces,K);
    rhsu = rhsu - .5*LIFT*(Fscale.*du);   

elseif opt==3 % split form but 
        
    

end 

