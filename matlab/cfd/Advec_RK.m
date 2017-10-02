% nodal wave RK driver

function Advec_RK

% Driver script for solving the 3D IPDGtion equations
Globals3D;

% Order of polymomials used for approximation
N = 9;
FinalTime = 1;

% % % % single element
% [VX VY VZ] = Nodes3D(1); [VX VY VZ] = xyztorst(VX,VY,VZ); K = 1; EToV = 1:length(VX);

% GMSH meshes
% cubeTetra1
filename = 'Grid/cube1.msh';
% filename = 'Grid/sphere112.msh';
%filename = 'Grid/sphere400.msh';
% filename = 'Grid/sphere1202.msh';
% filename = 'Grid/sphere2736.msh';
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGmsh3D(filename);

% VX = 2*VX; VY = 2*VY; VZ = 2*VZ; % biunit cube

% Initialize solver and construct grid and metric
StartUp3D;

%% cubature and plotting

JK = J(1,:);

[rq sq tq wq] = tet_cubature(2*N+1);
Vq = Vandermonde3D(N, rq, sq, tq)/V;
xq = Vq*x; yq = Vq*y; zq = Vq*z; % phys cubature nodes

% for plotting - build coordinates of all the nodes
[rp sp tp] = EquiNodes3D(25);
Vp = Vandermonde3D(N,rp,sp,tp)*invV;
xp = Vp*x; yp = Vp*y; zp = Vp*z;


%% 

rad = @(x,y,z) sqrt((x+.0).^2 + y.^2 + z.^2);
uex = @(x,y,z,t) exp(-25*rad(x,y,z));

% compute time step size
dt = .125/max(((N+1)^2)*Fscale(:));

% set initial conditions
u = uex(x,y,z,0);

% up = Vp*u; ids = yp > 0; color_line3(xp(ids),yp(ids),zp(ids),up(ids),'.');return
%% get matrix

u = zeros(Np,K);
A = zeros(Np*K);
for i = 1:Np*K
    u(i) = 1;
    rhsu = AdvecRHS(u);
    A(:,i) = rhsu(:);
    u(i) = 0;    
end

VB = bern_basis_tet(N,r,s,t);
T1 = kron(speye(K),inv(VB));
T2 = kron(speye(K),VB);
AB = T1*A*T2;
AB(abs(AB)<1e-8) = 0;

P = {};
for e = 1:K
    ids = (1:Np) + (e-1)*Np;
    P{e} = AB(ids,ids); % block Jacobi
end
P = blkdiag(P{:});
P(abs(P)<1e-6) = 0;

(Np*Np*K)/nnz(P)

keyboard
%%

% outer time step loop
time = 0;
tstep = 0;

totalSteps = floor(FinalTime/dt);

resu = zeros(Np,K); 

% while(tstep<1)
while (time<FinalTime)
    
    if(time+dt>FinalTime), dt = (FinalTime-time);  end;
        
    % low storage RK
    for INTRK = 1:5
        
        rhsu = AdvecRHS(u);
        
        resu = rk4a(INTRK)*resu + dt*rhsu;
        u = u + rk4b(INTRK)*resu;
        
    end;

    time = time+dt;
    tstep = tstep+1; % Increment time
         
    if 1 && (mod(tstep,10)==0)
        clf
        up = Vp*reshape(u,Np,K);
        
        ids = yp>0;
        color_line3(xp(ids),yp(ids),zp(ids),up(ids),'.');
        title(sprintf('time t = %f',time))
%         caxis([0,1])
        colorbar
        axis([min(xp(:)) max(xp(:)) min(yp(:)) max(yp(:)) min(zp(:)) max(zp(:))])
        view(3);
                
        drawnow
    end
    if mod(tstep,floor(totalSteps/10))==0
        disp(sprintf('on timestep %d / %d\n',tstep,totalSteps));
    end
    
        
end
% keyboard
uexq = uex(xq,yq,zq,time);
err = (Vq*p-uexq).^2;
for e = 1:K
    err(:,e) = err(:,e).*wq*JK(e);
end
err = sqrt(sum(err(:)));

% err


function [rhsu] = AdvecRHS(u)
% Purpose  : Evaluate RHS flux in 3D DG

Globals3D;

Dru = Dr*u; Dsu = Ds*u; Dtu = Dt*u;
dudx = rx.*Dru + sx.*Dsu + tx.*Dtu;

divF = dudx;
flux_u = zeros(Nfp*Nfaces,K);
uP = u(vmapP);
ids = abs(x(vmapP)+.5)<1e-8;
uP(ids) = 0;
flux_u(:) = .5*(uP-u(vmapM)).*(nx(:) - abs(nx(:)));

rhsu = -(divF + LIFT*(Fscale.*flux_u));  % simplified version


