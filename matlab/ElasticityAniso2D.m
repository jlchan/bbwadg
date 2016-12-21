function ElasticityAniso2D

% clear all, clear
clear -global *

Globals2D

K1D = 16;
N = 5;
c_flag = 0;
FinalTime = .5; % milliseconds

% filename = 'Grid/Other/block2.neu';
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(filename);
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);

StartUp2D;

% BuildPeriodicMaps2D(1,1);

% BCType = zeros(size(EToV));
% for e = 1:K
%     BCType(e,EToE(e,:)==e) = 6;
% end
% 
% BuildBCMaps2D;
% nref = 1;
% for ref = 1:nref
%     Refine2D(ones(size(EToV)));
%     StartUp2D;
%     BuildBCMaps2D;
% end


[rp sp] = EquiNodes2D(25); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

Nq = 2*N+1;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
Pq = V*V'*Vq'*diag(wq); % J's cancel out
Mref = inv(V*V');
xq = Vq*x; yq = Vq*y;
Jq = Vq*J;

M = inv(V*V');
wqJ = diag(wq)*(Vq*J);


%%
global mapBL vmapBL t0

t0 = .0;
% t0 = .075;
% t0 = .05;

% find x = 0 faces 
fbids = reshape(find(vmapP==vmapM),Nfp,nnz(vmapP==vmapM)/Nfp);
xfb = x(vmapP(fbids));
bfaces = sum(abs(xfb),1)<1e-8;

fbids = fbids(:,bfaces); 
fbids = fbids(:);

mapBL = fbids;
vmapBL = vmapP(mapBL);

%%
global Nfld mu lambda rho Vq Pq tau useWADG
Nfld = 5; %(u1,u2,sxx,syy,sxy)

rho = ones(size(xq));

global C11 C12 C13 C22 C23 C33

useWADG = 0;

a = .5;
k = 2;
mu = 1+a*cos(k*pi*xq).*cos(k*pi*yq);
lambda = 1+a*cos(k*pi*xq).*cos(k*pi*yq);

% anisotropic
C11a = (2*mu+lambda)/4;
C12a = lambda/2;
C13a = 0;
C22a = (2*mu+lambda);
C23a = 0;
C33a = mu;

% isotropic media
C11b = 2*mu+lambda;
C12b = lambda;
C13b = 0;
C22b = 2*mu+lambda;
C23b = 0;
C33b = mu;

C11 = C11a.*(yq < 0) + C11b.*(yq > 0);
C12 = C12a.*(yq < 0) + C12b.*(yq > 0);
C13 = C13a.*(yq < 0) + C13b.*(yq > 0);
C22 = C22a.*(yq < 0) + C22b.*(yq > 0);
C23 = C23a.*(yq < 0) + C23b.*(yq > 0);
C33 = C33a.*(yq < 0) + C33b.*(yq > 0);

% vv = C11; color_line3(xq,yq,vv,vv,'.'); return

Cnorm = zeros(size(xq));
Cmax = 0;
for i = 1:length(xq(:))
    C = [C11(i) C12(i) C13(i);
        C12(i) C22(i) C23(i);
        C13(i) C23(i) C33(i)];
    Cnorm(i) = norm(C);
    Cmax = max(Cmax,norm(C));
end
if useWADG==0
    rho = Pq*rho;
    C11 = Pq*C11;
    C12 = Pq*C12;
    C13 = Pq*C13;
    C22 = Pq*C22;
    C23 = Pq*C23;
    C33 = Pq*C33;
end


tau0 = 1;
for fld = 1:5
    tau{fld} = tau0*ones(size(x));
    if fld > 2        
        tau{fld} = tau0*ones(size(x));        
%         tau{fld} = tau0./max(Cnorm(:))*ones(size(x));
    end
end


%% params setup

x0 = 0;
y0 = 0.025;
p = @(x,y) exp(-25^2*((x-x0).^2 + (y-y0).^2));
% p = exp(-50^2*((x-x0).^2));
u = zeros(Np, K);

U{1} = p(x,y);
U{2} = u;
U{3} = u;
U{4} = u;
U{5} = u;

% global rick ptsrc 
% f0 = 170;
% tR = 1/f0;
% rick = @(t) 1e6*(1 - 2*(pi*f0*(t-tR)).^2).*exp(-(pi*f0*(t-tR)).^2);
% % tt = 0:.0001:.1; plot(tt,rick(tt));return
% ptsrc = exp(-(250)^2*((xq-x0).^2 + (yq-y0).^2));
% % r2 = (yq - x0).^2 + (yq-y0).^2; delta = 1e-5; ptsrc = cos(pi*r2*.5./delta).*(abs(r2)<delta);
% ptsrc = Pq*ptsrc;
% ptsrc = ptsrc/max(abs(ptsrc(:)));

err = sqrt(sum(sum(wqJ.*(p(xq,yq) - Vq*U{1}).^2)))
vv = Vp*U{1};
% color_line3(xp,yp,vv,vv,'.');return

%%

time = 0;

% Runge-Kutta residual storage
for fld = 1:Nfld
    res{fld} = zeros(Np,K);
end

% compute time step size
CN = (N+1)^2/2; % guessing...
dt = 2/(Cmax*CN*max(Fscale(:)))
ceil(FinalTime/dt)

% outer time step loop
tstep = 0;

% figure
% colormap(gray)
% colormap(hot)
while (time<FinalTime)   
            
    if(time+dt>FinalTime), dt = FinalTime-time; end
    
    for INTRK = 1:5
        
        timeloc = time + rk4c(INTRK)*dt;
        rhs = ElasRHS2D(U,timeloc);
                        
        % initiate and increment Runge-Kutta residuals
        for fld = 1:Nfld
            res{fld} = rk4a(INTRK)*res{fld} + dt*rhs{fld};
            U{fld} = U{fld} + rk4b(INTRK)*res{fld};
        end
        
    end;
    
    if 1 && mod(tstep,10)==0
        clf
        
%         p = U{3} + U{4}; % trace(S)        
        p = U{1};

        vv = Vp*p;
        vv = abs(vv);
        color_line3(xp,yp,vv,vv,'.');
        axis tight        
        axis equal
        title(sprintf('time = %f',time));
        colorbar;
%         view(3); 
        drawnow
        
        
    end
    
    % Increment time
    time = time+dt; tstep = tstep+1;        
    
    if mod(tstep,100)==0
        disp(sprintf('On timestep %d out of %d\n',tstep,round(FinalTime/dt)))
    end
end

keyboard



function [rhs] = ElasRHS2D(U,time)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;

global Nfld mu lambda rho Vq Pq tau useWADG
global C11 C12 C13 C22 C23 C33

% Define field differences at faces
for fld = 1:Nfld
    u = U{fld};
    
    % compute jumps
    dU{fld} = zeros(Nfp*Nfaces,K);
    dU{fld}(:) = u(vmapP)-u(vmapM);
    
    ur = Dr*u;
    us = Ds*u;
    Ux{fld} = rx.*ur + sx.*us;
    Uy{fld} = ry.*ur + sy.*us;
end

divSx = Ux{3} + Uy{5}; % d(Sxx)dx + d(Sxy)dy
divSy = Ux{5} + Uy{4}; % d(Sxy)dx + d(Syy)dy
du1dx = Ux{1}; % du1dx
du2dy = Uy{2}; % du2dy
du12dxy = Ux{2} + Uy{1}; % du2dx + du1dy

% velocity fluxes
nSx = nx.*dU{3} + ny.*dU{5};
nSy = nx.*dU{5} + ny.*dU{4};
 
opt=2;
if opt==1 % traction BCs
    %     global mapBL vmapBL t0
    %     if time < t0
    %         dU{1}(mapBL) = 2*(sin(pi*time/t0)*(time<t0) - U{1}(vmapBL));
    %         dU{2}(mapBL) = -2*U{2}(vmapBL);
    %     else
    nSx(mapB) = -2*(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
    nSy(mapB) = -2*(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));
    %     end    
elseif opt==2 % basic ABCs
    nSx(mapB) = -(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
    nSy(mapB) = -(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));
    dU{1}(mapB) = -U{1}(vmapB);
    dU{2}(mapB) = -U{2}(vmapB);    
elseif opt==3 % zero velocity
    dU{1}(mapB) = -2*U{1}(vmapB);
    dU{2}(mapB) = -2*U{2}(vmapB);    
    
end

% stress fluxes
nUx = dU{1}.*nx;
nUy = dU{2}.*ny;
nUxy = dU{2}.*nx + dU{1}.*ny;

% evaluate central fluxes
fc{1} = nSx;
fc{2} = nSy;
fc{3} = nUx;
fc{4} = nUy;
fc{5} = nUxy;

% penalization terms - reapply An
fp{1} = nx.*fc{3} + ny.*fc{5};
fp{2} = nx.*fc{5} + ny.*fc{4};
fp{3} = fc{1}.*nx;
fp{4} = fc{2}.*ny;
fp{5} = fc{2}.*nx + fc{1}.*ny;


flux = cell(5,1);
for fld = 1:Nfld   
    flux{fld} = zeros(Nfp*Nfaces,K);
    flux{fld}(:) = fc{fld}(:) + tau{fld}(vmapM).*fp{fld}(:);
end

% compute right hand sides of the PDE's
global rick ptsrc 
rr{1} =  divSx   +  LIFT*(Fscale.*flux{1})/2.0;
rr{2} =  divSy   +  LIFT*(Fscale.*flux{2})/2.0; % + rick(time)*ptsrc;
rr{3} =  du1dx   +  LIFT*(Fscale.*flux{3})/2.0;
rr{4} =  du2dy   +  LIFT*(Fscale.*flux{4})/2.0;
rr{5} =  du12dxy +  LIFT*(Fscale.*flux{5})/2.0;

% C = [2*mu+lambda       lambda       0
%      lambda       2*mu+lambda       0
%           0       0                mu]
if useWADG
    for fld = 1:Nfld
        rr{fld} = Vq*rr{fld};
    end
    rhs{1} = Pq*(rr{1}./rho);
    rhs{2} = Pq*(rr{2}./rho);
    rhs{3} = Pq*(C11.*rr{3} + C12.*rr{4} + C13.*rr{5});
    rhs{4} = Pq*(C12.*rr{3} + C22.*rr{4} + C23.*rr{5});
    rhs{5} = Pq*(C13.*rr{3} + C23.*rr{4} + C33.*rr{5});
else
    rhs{1} = rr{1}./rho;
    rhs{2} = rr{2}./rho;
    rhs{3} = C11.*rr{3} + C12.*rr{4} + C13.*rr{5};
    rhs{4} = C12.*rr{3} + C22.*rr{4} + C23.*rr{5};
    rhs{5} = C13.*rr{3} + C23.*rr{4} + C33.*rr{5};
end

return;

