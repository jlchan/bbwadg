% nodal wave RK driver

function Elasticity3D

% Driver script for solving the 3D IPDGtion equations
Globals3D;

% Order of polymomials used for approximation
N = 2;
FinalTime = 2;

% % % % single element
% [VX VY VZ] = Nodes3D(1); [VX VY VZ] = xyztorst(VX,VY,VZ); K = 1; EToV = 1:length(VX);

% GMSH meshes
% cubeTetra1
filename = 'Grid/cube1.msh';
% filename = 'Grid/sphere48.msh';
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGmsh3D(filename);
% VX = 2*VX; VY = 2*VY; VZ = 2*VZ; % biunit cube

% Initialize solver and construct grid and metric
StartUp3D;

%% cubature and plotting

global Nfld mu lambda Vq Pq tau useWADG

JK = J(1,:);

[rq sq tq wq] = tet_cubature(2*N+1);
Vq = Vandermonde3D(N, rq, sq, tq)/V;
Pq = (V*V')*Vq' * diag(wq);
xq = Vq*x; yq = Vq*y; zq = Vq*z; % phys cubature nodes

% for plotting - build coordinates of all the nodes
[rp sp tp] = EquiNodes3D(40);
Vp = Vandermonde3D(N,rp,sp,tp)/V;
xp = Vp*x; yp = Vp*y; zp = Vp*z;


%% problem params

Nfld = 9;  %(u1,u2,u3, 6 stresses)

mu = ones(size(x));
lambda = ones(size(x));

tau0 = 1;
for fld = 1:Nfld
    tau{fld} = tau0;
    if fld > 2
        %tau{fld} = tau0./max(2*mu(:)+lambda(:));
        tau{fld} = tau0;
    end
end

% compute time step size
CN = (N+1)*(N+3)/3;
dt = 2/(max(mu(:)+lambda(:))*CN*max(Fscale(:)));

% set initial conditions
for fld = 1:Nfld
    U{fld} = zeros(size(x));
end
r2 = @(x,y,z) x.^2 + y.^2 + z.^2;
U{1} = Pq*exp(-100*r2(xq,yq,zq));
% U{1} = exp(-16^2*r2(x,y,z));
% keyboard
% ids = yp>0; p = Vp*U{1}; plot3(VX,VY,VZ,'k.','markersize',32); hold on;color_line3(xp(ids),yp(ids),zp(ids),p(ids),'.'); view(3); colorbar; return

%% eigs

if 0
    tauvec = 0;
    for ii = 1:length(tauvec)
        for fld = 1:Nfld
            tau{fld} = tauvec(ii);
            if fld > 2
                tau{fld} = tauvec(ii)./max(mu(:)+lambda(:));
            end
        end
        u = zeros(Nfld*Np*K,1);
        rhs = zeros(Nfld*Np*K,1);
        A = zeros(Nfld*Np*K);
        ids = 1:Np*K;
        for i = 1:Nfld*Np*K
            u(i) = 1;
            for fld = 1:Nfld
                U{fld} = reshape(u(ids + (fld-1)*Np*K),Np,K);
            end
            rU = ElasRHS3D(U);
            u(i) = 0;
            for fld = 1:Nfld
                rhs(ids + (fld-1)*Np*K) = rU{fld};
            end
            A(:,i) = rhs(:);
            if (mod(i,100)==0)
                disp(sprintf('on col %d out of %d\n',i,Np*K*Nfld))
            end
        end
        lam = eig(A);        
        points{ii} = [real(lam) imag(lam)];
        plot(lam,'.','markersize',32)
        hold on
        title(sprintf('Largest real part = %g\n',max(real(lam))))
        axis equal
%         drawnow
%         max(abs(lam))
    end
    keyboard
end

%%  time step loop
time = 0;
tstep = 0;

totalSteps = floor(FinalTime/dt);

res = cell(Nfld,1);
for fld = 1:Nfld
    res{fld} = zeros(size(x));
end

ids = yp > 0;
ax = [min(xp(:)) max(xp(:)) min(yp(:)) max(yp(:)) min(zp(:)) max(zp(:))];

% while(tstep<1)
while (time<FinalTime)
    
    if(time+dt>FinalTime), dt = (FinalTime-time);  end;
        
    % low storage RK
    for INTRK = 1:5
        
        [rhs] = ElasRHS3D(U);
        
        for fld = 1:Nfld
            res{fld} = rk4a(INTRK)*res{fld} + dt*rhs{fld};
            U{fld} = U{fld} + rk4b(INTRK)*res{fld};
        end
        
    end;

    time = time+dt;
    tstep = tstep+1; % Increment time
         
    if 1 && (mod(tstep,5)==0)
        clf        

        p = U{1};
%         p = U{4}+U{5}+U{6};
        up = Vp*reshape(p,Np,K);
               
        color_line3(xp(ids),yp(ids),zp(ids),up(ids),'.');        
        title(sprintf('time t = %f',time))
%         caxis([-1,1]*.075)
        colorbar
        axis(ax)
        view(0,0)
                
        drawnow
    end
    if mod(tstep,floor(totalSteps/10))==0
        disp(sprintf('on timestep %d / %d\n',tstep,totalSteps));
    end
    
        
end


function [rhs] = ElasRHS3D(U)
% Purpose  : Evaluate RHS flux in 3D DG

Globals3D;

global Nfld mu lambda Vq Pq tau useWADG

% Define field differences at faces
for fld = 1:Nfld
    u = U{fld};
    
    % compute jumps
    dU{fld} = zeros(Nfp*Nfaces,K);
    dU{fld}(:) = u(vmapP)-u(vmapM);
    
    ur = Dr*u; 
    us = Ds*u; 
    ut = Dt*u; 
    Ux{fld} = rx.*ur + sx.*us + tx.*ut;
    Uy{fld} = ry.*ur + sy.*us + ty.*ut;
    Uz{fld} = rz.*ur + sz.*us + tz.*ut;
end

divSx = Ux{4} + Uy{9} + Uz{8}; 
divSy = Ux{9} + Uy{5} + Uz{7}; 
divSz = Ux{8} + Uy{7} + Uz{6}; 
du1dx = Ux{1}; 
du2dy = Uy{2}; 
du3dz = Uz{3}; 
du1 = Uy{3} + Uz{2};
du2 = Ux{3} + Uz{1};
du3 = Ux{2} + Uy{1};

% velocity fluxes
nSx = nx.*dU{4} + ny.*dU{9} + nz.*dU{8};
nSy = nx.*dU{9} + ny.*dU{5} + nz.*dU{7};
nSz = nx.*dU{8} + ny.*dU{7} + nz.*dU{6};

if 0 % traction BCs
    nSx(mapB) = -2*(nx(mapB).*U{4}(vmapB) + ny(mapB).*U{9}(vmapB) + nz(mapB).*U{8}(vmapB));
    nSy(mapB) = -2*(nx(mapB).*U{9}(vmapB) + ny(mapB).*U{5}(vmapB) + nz(mapB).*U{7}(vmapB));
    nSz(mapB) = -2*(nx(mapB).*U{8}(vmapB) + ny(mapB).*U{7}(vmapB) + nz(mapB).*U{6}(vmapB));
else % basic abcs
    nSx(mapB) = -(nx(mapB).*U{4}(vmapB) + ny(mapB).*U{9}(vmapB) + nz(mapB).*U{8}(vmapB));
    nSy(mapB) = -(nx(mapB).*U{9}(vmapB) + ny(mapB).*U{5}(vmapB) + nz(mapB).*U{7}(vmapB));
    nSz(mapB) = -(nx(mapB).*U{8}(vmapB) + ny(mapB).*U{7}(vmapB) + nz(mapB).*U{6}(vmapB));
    dU{1}(mapB) = -U{1}(vmapB);
    dU{2}(mapB) = -U{2}(vmapB);
    dU{3}(mapB) = -U{3}(vmapB);
end

% evaluate central fluxes
fc{1} = nSx;
fc{2} = nSy;
fc{3} = nSz;
fc{4} = dU{1}.*nx;
fc{5} = dU{2}.*ny;
fc{6} = dU{3}.*nz;
fc{7} = dU{3}.*ny + dU{2}.*nz;
fc{8} = dU{3}.*nx + dU{1}.*nz;
fc{9} = dU{2}.*nx + dU{1}.*ny;

% penalization terms - reapply An
fp{1} = nx.*fc{4} + ny.*fc{9} + nz.*fc{8};
fp{2} = nx.*fc{9} + ny.*fc{5} + nz.*fc{7};
fp{3} = nx.*fc{8} + ny.*fc{7} + nz.*fc{6};
fp{4} = fc{1}.*nx;
fp{5} = fc{2}.*ny;
fp{6} = fc{3}.*nz;
fp{7} = fc{3}.*ny + fc{2}.*nz;
fp{8} = fc{3}.*nx + fc{1}.*nz;
fp{9} = fc{2}.*nx + fc{1}.*ny;

flux = cell(Nfld,1);
for fld = 1:Nfld   
    flux{fld} = zeros(Nfp*Nfaces,K);    
    flux{fld}(:) = tau{fld}*fp{fld}(:) + fc{fld}(:);
end

% compute right hand sides of the PDE's
rr{1} =  divSx   +  LIFT*(Fscale.*flux{1})/2.0;
rr{2} =  divSy   +  LIFT*(Fscale.*flux{2})/2.0;
rr{3} =  divSz   +  LIFT*(Fscale.*flux{3})/2.0;
rr{4} =  du1dx   +  LIFT*(Fscale.*flux{4})/2.0;
rr{5} =  du2dy   +  LIFT*(Fscale.*flux{5})/2.0;
rr{6} =  du3dz   +  LIFT*(Fscale.*flux{6})/2.0;
rr{7} =  du1     +  LIFT*(Fscale.*flux{7})/2.0;
rr{8} =  du2     +  LIFT*(Fscale.*flux{8})/2.0;
rr{9} =  du3     +  LIFT*(Fscale.*flux{9})/2.0;

for fld = 1:Nfld
    Lf{fld} = LIFT*(Fscale.*flux{fld})/2.0;
end
[divSx(:,1) divSy(:,1) divSz(:,1) du1dx(:,1) du2dy(:,1) du3dz(:,1)]
[flux{1}(:,1) flux{2}(:,1) flux{3}(:,1) flux{4}(:,1) flux{5}(:,1) flux{6}(:,1)]
[Lf{1}(:,1) Lf{2}(:,1) Lf{3}(:,1) Lf{4}(:,1) Lf{5}(:,1) Lf{6}(:,1)]
[rr{1}(:,1) rr{2}(:,1) rr{3}(:,1) rr{4}(:,1) rr{5}(:,1) rr{6}(:,1)]
keyboard

% C = [2*mu+lambda       lambda       0
%      lambda       2*mu+lambda       0
%           0       0                mu/2]
rhs{1} = rr{1};
rhs{2} = rr{2};
rhs{3} = rr{3};
rhs{4} = (2*mu+lambda).*rr{4} + lambda.*(rr{5} + rr{6});
rhs{5} = (2*mu+lambda).*rr{5} + lambda.*(rr{4} + rr{6});
rhs{6} = (2*mu+lambda).*rr{6} + lambda.*(rr{4} + rr{5});
rhs{7} = (mu) .* rr{7};
rhs{8} = (mu) .* rr{8};
rhs{9} = (mu) .* rr{9};



