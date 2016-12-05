function Elasticity2D_disc

% clear all, clear
clear -global *

Globals2D

% K1D = 8;
N = 2;
c_flag = 0;
FinalTime = .25;

% filename = 'Grid/Other/block2.neu';
filename = 'Grid/Other/circA01.neu';
% filename = 'Grid/CNS2D/couette.neu';
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(filename);
[Nv, VX, VY, K, EToV] = MeshReaderGambitBC2D(filename);
% aa = 1;
% % [Nv, VX, VY, K, EToV] = unif_tri_mesh(aa*K1D,K1D);
% [Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D,aa*K1D+1);
% VY = aa*VY;

addpath('./distmesh');
r1 = 1.0;
r2 = .5;
h = 2*.0625;
fd=@(p) ddiff(dcircle(p,0,0,r1),dcircle(p,0,0,r2));
% fd = @(p) dcircle(p,0,0,1);
[p,t]=distmesh2d(fd,@huniform,h,1.1*[-1,-1;1,1],1.1*[-1,-1;-1,1;1,-1;1,1]);
close all
VX = p(:,1)'; VY = p(:,2)'; EToV = t; K = size(EToV,1);

StartUp2D;

BCType = zeros(size(EToV));
fids = [1 2;2 3; 3 1];
for e = 1:K
    for f = 1:3
        fid = fids(f,:);
        rr = sqrt(VX(EToV(e,fid)).^2 + VY(EToV(e,fid)).^2);
        if all(abs(rr-r1)<1e-8)
            BCType(e,f) = 1;
        elseif all(abs(rr-r2)<1e-8)
            BCType(e,f) = 2;
        end
    end
end
[k1,f1] = find(BCType==1);
[k2,f2] = find(BCType==2);

% Push all boundary faces to unit cylinder
Nq = 2*N+1;
% MakeCylinder2D([k1,f1], r1, 0, 0);
% MakeCylinder2D([k2,f2], r2, 0, 0);

% PlotMesh2D; axis on; return

% for i = 1:length(mapB)
%     clf
%     plot(x,y,'.')
%     hold on
%     id = vmapM(mapB(i));
%    plot(x(id), y(id),'o','markersize',32)
%    id = vmapP(mapB(i));
%    plot(x(id), y(id),'x','markersize',32)
%    pause
% end
% return

%%
[rp sp] = EquiNodes2D(50); [rp sp] = xytors(rp,sp);
rp = .99*rp; sp = .99*sp;
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;


global Vq wq Vrq Vsq
global Vfq wfq VfqFace
global rxq sxq ryq syq Jq Jfq nxq nyq sJq
global Pq Pfq


Nq = 2*N+1;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
Pq = V*V'*Vq'*diag(wq); % J's cancel out
Mref = inv(V*V');
xq = Vq*x; yq = Vq*y;
Jq = Vq*J;

[Vrq Vsq] = GradVandermonde2D(N,rq,sq);
Vrq = Vrq/V; Vsq = Vsq/V;

% interp nodal to quadrature
[rq1D, wq1D] = JacobiGQ(0, 0, Nq);
rfq = [rq1D, -rq1D, -ones(size(rq1D))];
sfq = [-ones(size(rq1D)), rq1D, -rq1D];
rfq = rfq(:); sfq = sfq(:); 
wfq = repmat(wq1D,Nfaces,1);
Vfq = Vandermonde2D(N,rfq,sfq)/V;
xfq = Vfq*x; yfq = Vfq*y;

[Vrfq Vsfq] = GradVandermonde2D(N,rfq,sfq);
Vrfq = Vrfq/V; Vsfq = Vsfq/V; 

V1D = Vandermonde1D(N,JacobiGL(0,0,N));
VfqFace = Vandermonde1D(N,rq1D)/V1D;
VfqFace = blkdiag(VfqFace,VfqFace,VfqFace); % repeat for 3 faces

% projection matrices
Pq = (V*V')*(Vq'*diag(wq));
Pfq = (V*V')*(Vfq'*diag(wfq));

[rxq,sxq,ryq,syq,Jq] = GeometricFactors2D(x,y,Vrq,Vsq);
[rxf,sxf,ryf,syf,Jfq] = GeometricFactors2D(x,y,Vrfq,Vsfq);

xrf = Vrfq*x; yrf = Vrfq*y; 
xsf = Vsfq*x; ysf = Vsfq*y; 
Jfq = xrf.*ysf-xsf.*yrf;

fid = (1:length(rq1D));
nxq(fid, :) =  yrf(fid, :); nyq(fid, :) = -xrf(fid, :);
fid = fid + (Nq+1);
nxq(fid, :) =  ysf(fid, :)-yrf(fid, :); nyq(fid, :) = -xsf(fid, :)+xrf(fid, :);
fid = fid + (Nq+1);
nxq(fid, :) = -ysf(fid, :); nyq(fid, :) =  xsf(fid, :);
sJq = sqrt(nxq.*nxq+nyq.*nyq); 
nxq = nxq./sJq; nyq = nyq./sJq;

% quiver(xfq(:),yfq(:),nxq(:),nyq(:))
% return


%% BC stuff

global mapBd vmapBd mapBt vmapBt

% find x = 0 faces
fbids = reshape(find(vmapP==vmapM),Nfp,nnz(vmapP==vmapM)/Nfp);
% xfb = x(vmapM(fbids));
% bfaces = sum(abs(abs(xfb)),1)<1e-8;
yfb = y(vmapM(fbids));
bfaces = sum(abs(abs(yfb)-.5),1)<1e-8;

% yfb = y(vmapM(fbids));
% bfaces = sum(abs(ptsfb),1)<1e-8;

fbids = fbids(:,bfaces);
fbids = fbids(:);

mapBt = fbids;
vmapBt = vmapP(fbids);

% remove Dirichlet BCs for traction
mapBd = setdiff(mapB,mapBt);
vmapBd = vmapP(mapBd);

if 0
    hold on
    plot(x(vmapB),y(vmapB),'bo')
%     plot(x(vmapBt),y(vmapBt),'bo')
    
%     plot(x(vmapBd),y(vmapBd),'rd')
    plot(x,y,'k.')
    return
end
%%

global Nfld mu lambda rho tau useWADG
Nfld = 5; %(u1,u2,sxx,syy,sxy)

rho = ones(size(x));
lambda = ones(size(x));
mu = ones(size(x));

tau0 = 1;
rhoF = max(rho(vmapM),rho(vmapP));
muF = max(mu(vmapM),mu(vmapP));
lambdaF = max(lambda(vmapM),lambda(vmapP));

rhoF = reshape(rhoF,Nfp*Nfaces,K);
muF = reshape(muF,Nfp*Nfaces,K);
lambdaF = reshape(lambdaF,Nfp*Nfaces,K);
tau{1} = tau0./rhoF;
tau{2} = tau0./rhoF;
tau{3} = tau0./(2*muF + lambdaF);
tau{4} = tau0./(2*muF + lambdaF);
tau{5} = tau0./(2*muF + lambdaF);

useWADG = 0;
if useWADG
    rho = Vq*rho;
    mu = Vq*mu;
    lambda = Vq*lambda;
end

%% params setup

% cp = sqrt((2*mu + lambda)/rho);
% cs = sqrt(mu/rho);

% assume constant
x0 = 0; y0 = .75;
p = exp(-16^2*((x-x0).^2 + (y-y0).^2));
u = zeros(Np,K);

rad = @(x,y) sqrt(x.^2 + y.^2);

alpha = 28.28892419267416;
B = 1.069885069032360;
w = alpha*sqrt((2*mu(1) + lambda(1))/rho(1));
pex = @(x,y,t) (besselj(1,alpha*rad(x,y)) + B*bessely(1,alpha*rad(x,y)))*cos(w*t);
dudr = @(x,y,t) cos(t*w)*(B*(alpha*bessely(0, alpha*rad(x,y)) - bessely(1, alpha*rad(x,y))./rad(x,y)) + alpha*besselj(0, alpha*rad(x,y)) - besselj(1, alpha*rad(x,y))./rad(x,y));
u1x = @(x,y,t) dudr(x,y,t).*x./rad(x,y); % fix - u1 = ur*cos(theta), u1x = du/dr * dr/dx * dcos(theta)/dtheta * dtheta/dx
u2y = @(x,y,t) dudr(x,y,t).*y./rad(x,y);
u12xy = @(x,y,t) dudr(x,y,t).*(x + y)./rad(x,y);

vr = @(x,y,t) -w*sin(t*w)*(B*bessely(1, alpha*rad(x,y)) + besselj(1, alpha*rad(x,y)));
v1 = @(x,y,t) vr(x,y,t).*cos(atan2(y,x));
v2 = @(x,y,t) vr(x,y,t).*sin(atan2(y,x));

sxx = @(x,y,t) ((2*mu+lambda) .* u1x(x,y,t)) + (lambda.*u2y(x,y,t));
syy = @(x,y,t) lambda.*u1x(x,y,t) + ((2*mu+lambda) .* u2y(x,y,t));
sxy = @(x,y,t) mu .* u12xy(x,y,t);

if 0
    for t = 0:.01:2
        clf
        vv = v1(xp,yp,t);
        color_line3(xp,yp,vv,vv,'.')
        %     view(3)
        axis([-1 1 -1 1 -10 10])
        drawnow
    end
    return
end

U{1} = v1(x,y,0);
U{2} = v2(x,y,0);
U{3} = sxx(x,y,0);
U{4} = syy(x,y,0);
U{5} = sxy(x,y,0);

% vv = Vp*U{3};color_line3(xp,yp,vv,vv,'.');return

%%

if 0
    tau0 = 1;
    
    rhoF = max(rho(vmapM),rho(vmapP));
    muF = max(mu(vmapM),mu(vmapP));
    lambdaF = max(lambda(vmapM),lambda(vmapP));
    
    rhoF = reshape(rhoF,Nfp*Nfaces,K);
    muF = reshape(muF,Nfp*Nfaces,K);
    lambdaF = reshape(lambdaF,Nfp*Nfaces,K);
    for fld = 1:5
        tau{fld} = tau0./rhoF;
        if fld > 2
            tau{fld} = tau0./(2*muF+lambdaF);
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
        rU = ElasRHS2D(U,0);
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
    plot(lam,'.','markersize',32)
    hold on
    title(sprintf('Largest real part = %g\n',max(real(lam))))
    axis equal
    %         drawnow
    max(abs(lam))

    keyboard
end
%%

time = 0;

% Runge-Kutta residual storage
for fld = 1:Nfld
    res{fld} = zeros(Np,K);
end

% compute time step size
CN = (N+1)^2/2; % guessing...
dt = 2/(sqrt(max((2*mu(:)+lambda(:))./rho(:)))*CN*max(Fscale(:)));

% outer time step loop
tstep = 0;


figure
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
        p = U{3}+U{4}; % trace(S)
%         p = U{5};
        vp1 = Vp*p;        
        color_line3(xp,yp,vp1,vp1,'.');
        axis equal
        axis tight
%         axis([-1 1 -1 1 -2 2])
        view(2)
        
        title(sprintf('time = %f',time))
        colorbar        
        drawnow
        
        
    end
    
    % Increment time
    time = time+dt; tstep = tstep+1;
    if mod(tstep,100)==0
        disp(sprintf('On timestep %d out of %d\n',tstep,round(FinalTime/dt)))
    end
end

keyboard

Nq = 3*N;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
xq = Vq*x; yq = Vq*y;
Jq = Vq*J;

wqJ = diag(wq)*(Vq*J);
err = (Vq*U{1} - v1(xq,yq,FinalTime)).^2 + (Vq*U{2} - v2(xq,yq,FinalTime)).^2;
uq = v1(xq,yq,FinalTime).^2 + v2(xq,yq,FinalTime).^2;
errU = sqrt(sum(wqJ(:).*err(:)))/sqrt(sum(wqJ(:).*uq(:)));

errU

return
serr = (Vq*U{3} - sxx(xq,yq,FinalTime)).^2 + (Vq*U{4} - syy(xq,yq,FinalTime)).^2 + (Vq*U{5} - sxy(xq,yq,FinalTime)).^2;
sq = sxx(xq,yq,FinalTime).^2 + syy(xq,yq,FinalTime).^2 + sxy(xq,yq,FinalTime).^2;
errS = sqrt(sum(wqJ(:).*serr(:)))/sqrt(sum(wqJ(:).*sq(:)));

figure
pterr = abs(Vp*U{1} - v1(xp,yp,FinalTime)) + abs(Vp*U{2} - v2(xp,yp,FinalTime));
color_line3(xp,yp,pterr,pterr,'.')
title(sprintf('velocity L2 err = %g\n',errU))

figure
pterr = abs(Vp*U{3} - sxx(xp,yp,FinalTime)) + abs(Vp*U{4} - syy(xp,yp,FinalTime)) + abs(Vp*U{5} - sxy(xp,yp,FinalTime));
color_line3(xp,yp,pterr,pterr,'.')
title(sprintf('sigma L2 err = %g\n',errS))

L2err = sqrt(sum(wqJ(:).*(err(:)+serr(:))))/sqrt(sum(wqJ(:).*(uq(:)+sq(:))))
errU
errS



function [rhs] = ElasRHS2D(U,time)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;

global Nfld mu lambda rho Vq Pq tau useWADG
global mapBd vmapBd mapBt vmapBt

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

global v1 v2
opt = 1;
if opt==1 % traction BCs
    nSx(mapB) = -2*(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
    nSy(mapB) = -2*(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));
    
elseif opt==2 % basic ABCs
    nSx(mapB) = -(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
    nSy(mapB) = -(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));
    dU{1}(mapB) = -U{1}(vmapB);
    dU{2}(mapB) = -U{2}(vmapB);
    
elseif opt==3
    
    dU{1}(mapB) = 2*(v1(x(vmapB),y(vmapB),time)-U{1}(vmapB));
    dU{2}(mapB) = 2*(v2(x(vmapB),y(vmapB),time)-U{2}(vmapB));
    
else
    
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
    flux{fld}(:) = fc{fld}(:) + tau{fld}(mapM).*fp{fld}(:);
end

% compute right hand sides of the PDE's
rr{1} =  divSx   +  LIFT*(Fscale.*flux{1})/2.0;
rr{2} =  divSy   +  LIFT*(Fscale.*flux{2})/2.0;
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
    rhs{3} = Pq*((2*mu+lambda).*rr{3} + lambda.*rr{4});
    rhs{4} = Pq*(lambda.*rr{3} + (2*mu+lambda).*rr{4});
    rhs{5} = Pq*((mu) .* rr{5});
else
    rhs{1} = rr{1}./rho;
    rhs{2} = rr{2}./rho;
    rhs{3} = (2*mu+lambda).*rr{3} + lambda.*rr{4};
    rhs{4} = lambda.*rr{3} + (2*mu+lambda).*rr{4};
    rhs{5} = (mu) .* rr{5};
end

return;



% copied in b/c matlab is funny
function MakeCylinder2D(faces, ra,xo,yo)

% Function: MakeCylinder2D(faces, ra, xo, yo)
% Purpose:  Use Gordon-Hall blending with an isoparametric map to modify a list
%           of faces so they conform to a cylinder of radius r centered at (xo,yo)
Globals2D;

NCurveFaces = size(faces,1);
vflag = zeros(size(VX));
for n=1:NCurveFaces
    
    % move vertices of faces to be curved onto circle
    k = faces(n,1); f = faces(n,2);
    v1 = EToV(k, f); v2 = EToV(k, mod(f,Nfaces)+1);
    
    % compute polar angles of start and end face vertices relative to circle center
    theta1 = atan2(VY(v1)-yo,VX(v1)-xo);
    theta2 = atan2(VY(v2)-yo,VX(v2)-xo);
    
    % move vertices onto circle
    newx1 = xo + ra*cos(theta1); newy1 = yo + ra*sin(theta1);
    newx2 = xo + ra*cos(theta2); newy2 = yo + ra*sin(theta2);
    
    % update mesh vertex locations
    VX(v1) = newx1; VX(v2) = newx2; VY(v1) = newy1; VY(v2) = newy2;
    
    % store modified vertex numbers
    vflag(v1) = 1;  vflag(v2) = 1;
end

% map modified vertex flag to each element
vflag = vflag(EToV);

% locate elements with at least one modified vertex
ks = find(sum(vflag,2)>0);

% build coordinates of all the corrected nodes
va = EToV(ks,1)'; vb = EToV(ks,2)'; vc = EToV(ks,3)';
x(:,ks) = 0.5*(-(r+s)*VX(va)+(1+r)*VX(vb)+(1+s)*VX(vc));
y(:,ks) = 0.5*(-(r+s)*VY(va)+(1+r)*VY(vb)+(1+s)*VY(vc));

for n=1:NCurveFaces  % deform specified faces
    k = faces(n,1); f = faces(n,2);
    
    % find vertex locations for this face and tangential coordinate
    if(f==1) v1 = EToV(k,1); v2 = EToV(k,2); vr = r; end
    if(f==2) v1 = EToV(k,2); v2 = EToV(k,3); vr = s; end
    if(f==3) v1 = EToV(k,1); v2 = EToV(k,3); vr = s; end
    fr = vr(Fmask(:,f));
    x1 = VX(v1); y1 = VY(v1); x2 = VX(v2); y2 = VY(v2);    
    
    % move vertices at end points of this face to the cylinder
    theta1 = atan2(y1-yo, x1-xo); theta2 = atan2(y2-yo, x2-xo);
    
    % check to make sure they are in the same quadrant
    if ((theta2 > 0) && (theta1 < 0)), theta1 = theta1 + 2*pi; end;
    if ((theta1 > 0) && (theta2 < 0)), theta2 = theta2 + 2*pi; end;
    
    % in case they cross the theta = 0 line
    if abs(theta2-theta1) > pi
        if theta2 > theta1
            theta2 = theta2 - 2*pi;
        else
            theta1 = theta1 - 2*pi;
        end
    end
    
    % distribute N+1 nodes by arc-length along edge
    theta = 0.5*theta1*(1-fr) + 0.5*theta2*(1+fr);

    % evaluate deformation of coordinates
    fdx = xo + ra*cos(theta)-x(Fmask(:,f),k);
    fdy = yo + ra*sin(theta)-y(Fmask(:,f),k);

    % build 1D Vandermonde matrix for face nodes and volume nodes
    Vface = Vandermonde1D(N, fr);  Vvol  = Vandermonde1D(N, vr);
    % compute unblended volume deformations
    vdx = Vvol*(Vface\fdx); vdy = Vvol*(Vface\fdy);

    % blend deformation and increment node coordinates
    ids = find(abs(1-vr)>1e-7); % warp and blend
    if(f==1) blend = -(r(ids)+s(ids))./(1-vr(ids)); end;
    if(f==2) blend =      +(r(ids)+1)./(1-vr(ids)); end;
    if(f==3) blend = -(r(ids)+s(ids))./(1-vr(ids)); end;
    
%     if (n==13)
%         keyboard
%     end
    
    x(ids,k) = x(ids,k)+blend.*vdx(ids);
    y(ids,k) = y(ids,k)+blend.*vdy(ids);

end

% repair other coordinate dependent information
Fx = x(Fmask(:), :); Fy = y(Fmask(:), :);
[rx,sx,ry,sy,J] = GeometricFactors2D(x, y,Dr,Ds);
[nx, ny, sJ] = Normals2D(); Fscale = sJ./(J(Fmask,:));
return

