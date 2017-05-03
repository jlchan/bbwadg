% nodal wave RK driver

function Wave_RK_quadrature

% Driver script for solving the 3D IPDGtion equations
Globals3D;

% Order of polymomials used for approximation
N = 2;
FinalTime = .01;

% % % % single element
% [VX VY VZ] = Nodes3D(1); [VX VY VZ] = xyztorst(VX,VY,VZ); K = 1; EToV = 1:length(VX);

% GMSH meshes
% cubeTetra1
filename = 'Grid/cubeTetra1.msh';
% filename = 'Grid/oneTet.msh';
% filename = 'Grid/sphere112.msh';
%filename = 'Grid/sphere400.msh';
% filename = 'Grid/sphere1202.msh';
% filename = 'Grid/sphere2736.msh';
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGmsh3D(filename);

% VX = 2*VX; VY = 2*VY; VZ = 2*VZ; % biunit cube

% Initialize solver and construct grid and metric
StartUp3D;

%% cubature and plotting

global Vq Vrq Vsq Vtq Pq Prq Psq Ptq 
global rxq sxq txq ryq syq tyq rzq szq tzq Jq 
global Vfq Vfqsurf Pfq nxq nyq nzq sJq 

[rq sq tq wq] = tet_cubature(2*N+1);
Vq = Vandermonde3D(N, rq, sq, tq)/V;
xq = Vq*x; yq = Vq*y; zq = Vq*z; % phys cubature nodes

rxq = Vq*rx; sxq = Vq*sx; txq = Vq*tx;
ryq = Vq*ry; syq = Vq*sy; tyq = Vq*ty;
rzq = Vq*rz; szq = Vq*sz; tzq = Vq*tz;
Jq  = Vq*J;

[Vrq Vsq Vtq] = GradVandermonde3D(N,rq,sq,tq);
Vrq = Vrq/V; Vsq = Vsq/V; Vtq = Vtq/V;

Pq  = V*V' * Vq'*diag(wq);
Prq = V*V' * Vrq'*diag(wq);
Psq = V*V' * Vsq'*diag(wq);
Ptq = V*V' * Vtq'*diag(wq);

[rfq sfq tfq wfq] = tet_surface_cubature(2*N+1);
wfq = repmat(wfq(1:length(wfq)/4),4,1); % remove scaling by face area 
Vfq = Vandermonde3D(N,rfq,sfq,tfq)/V;
Pfq = V*V' * Vfq'*diag(wfq);

rqtri = rfq(1:length(rfq)/4);
sqtri = sfq(1:length(rfq)/4);
Vtri = Vandermonde2D(N,r(Fmask(:,1)),s(Fmask(:,1)));
Vqtri = Vandermonde2D(N,rqtri,sqtri)/Vtri;
Vfqsurf = blkdiag(Vqtri,Vqtri,Vqtri,Vqtri);
nxq = Vfqsurf*nx;
nyq = Vfqsurf*ny;
nzq = Vfqsurf*nz;
sJq = Vfqsurf*sJ;

% rf = Vfqsurf*r(Fmask(:));
% sf = Vfqsurf*s(Fmask(:));
% tf = Vfqsurf*t(Fmask(:));
% plot3(rf,sf,tf,'o')
% return

% for plotting - build coordinates of all the nodes
[rp sp tp] = EquiNodes3D(25);
Vp = Vandermonde3D(N,rp,sp,tp)*invV;
xp = Vp*x; yp = Vp*y; zp = Vp*z;

%% 

uex = @(x,y,z,t) cos(pi*x).*cos(pi*y).*cos(pi*z)*cos(sqrt(3)*pi*t);
% uex = @(x,y,z,t) sqrt(pi/2)*besselj(.5,pi*sqrt(x.^2 + y.^2 + z.^2))./sqrt(pi*sqrt(x.^2 + y.^2 + z.^2)).*cos(pi*t);
% uex = @(x,y,z,t) sin(pi*sqrt(x.^2 + y.^2 + z.^2))./(pi*sqrt(x.^2 + y.^2 + z.^2)).*cos(pi*t);

% compute time step size
dt = .125/max( ((N+1)^2)*Fscale(:))

% set initial conditions
%p = V\uex(x,y,z,0);
p = uex(x,y,z,0);
u = zeros(Np,K); v = zeros(Np,K); w = zeros(Np,K); % zero velocities

Q = [
 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
 2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
 3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
 4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
 5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28
 6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29
 7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30
 8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0];

% p = Q(1:Np,:);


% p(:) = 0:Np*K-1; % testing

% outer time step loop
time = 0;
tstep = 0;

totalSteps = floor(FinalTime/dt);

resp = zeros(Np,K);
resu = zeros(Np,K); 
resv = zeros(Np,K);
resw = zeros(Np,K);

% while(tstep<1)
while (time<FinalTime)
    
    if(time+dt>FinalTime), dt = (FinalTime-time);  end;
        
    % low storage RK
    for INTRK = 1:5
        
%         [rhsp rhsu rhsv rhsw] = WaveRHS3D(p,u,v,w);
        [rhsp rhsu rhsv rhsw] = WaveRHS_quad(p,u,v,w);
        
        resp = rk4a(INTRK)*resp + dt*rhsp;        
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        resw = rk4a(INTRK)*resw + dt*rhsw;        
        
        p = p + rk4b(INTRK)*resp;
        u = u + rk4b(INTRK)*resu;
        v = v + rk4b(INTRK)*resv;
        w = w + rk4b(INTRK)*resw;
       
        keyboard
    end;

    time = time+dt;
    tstep = tstep+1; % Increment time
         
    if 1 && (mod(tstep,10)==0)
        clf
        up = Vp*reshape(p,Np,K);
        
        ids = yp>0;
        color_line3(xp(ids),yp(ids),zp(ids),up(ids),'.');
        title(sprintf('time t = %f',time))
        caxis([-1,1])
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
JK = J(1,:);

for e = 1:K
    err(:,e) = err(:,e).*wq*JK(e);
end
err = sqrt(sum(err(:)));

err


function [rhsp rhsu rhsv rhsw] = WaveRHS_quad(p,u,v,w)

global Vq Vrq Vsq Vtq Pq Prq Psq Ptq 
global rxq sxq txq ryq syq tyq rzq szq tzq Jq 
global Vfq Vfqsurf Pfq nxq nyq nzq sJq 

Globals3D;

pr = Vrq*p; ps = Vsq*p; pt = Vtq*p;
ur = Vrq*u; us = Vsq*u; ut = Vtq*u;
vr = Vrq*v; vs = Vsq*v; vt = Vtq*v;
wr = Vrq*w; ws = Vsq*w; wt = Vtq*w;

dpdxq = rxq.*pr + sxq.*ps + txq.*pt;
dpdyq = ryq.*pr + syq.*ps + tyq.*pt;
dpdzq = rzq.*pr + szq.*ps + tzq.*pt;

uq = Vq*u; vq = Vq*v; wq = Vq*w;
Ur = rxq.*uq + ryq.*vq + rzq.*wq;
Us = sxq.*uq + syq.*vq + szq.*wq;
Ut = txq.*uq + tyq.*vq + tzq.*wq;

dpdx = Pq*dpdxq;
dpdy = Pq*dpdyq;
dpdz = Pq*dpdzq;
divU = -(Prq*Ur + Psq*Us + Ptq*Ut);

p_jump = zeros(Nfp*Nfaces,K);
nu_jump = zeros(Nfp*Nfaces,K);
nu_avg = zeros(Nfp*Nfaces,K);

p_jump(:) = p(vmapP)-p(vmapM);
p_jump(mapB) = -2*p(vmapB); % boundary conditions

uM = u(vmapM); vM = v(vmapM); wM = w(vmapM);
uP = u(vmapP); uP(mapB) = uM(mapB);
vP = v(vmapP); vP(mapB) = vM(mapB);
wP = w(vmapP); wP(mapB) = wM(mapB);

u_jump  = uP - uM; v_jump  = vP - vM; w_jump  = wP - wM;
u_avg   = uP + uM; v_avg   = vP + vM; w_avg   = wP + wM;
nu_jump(:) = nx(:).*u_jump + ny(:).*v_jump + nz(:).*w_jump;
nu_avg(:)  = .5*(nx(:).*u_avg + ny(:).*v_avg + nz(:).*w_avg);

% flux_p = (.5*p_jump - nu_avg).*Fscale;
% flux_u = .5*(nu_jump - p_jump).*Fscale;
% rhsp = -divU + LIFT*(flux_p);
% rhsu = -dpdx + LIFT*(nx.*flux_u);  % simplified version
% rhsv = -dpdy + LIFT*(ny.*flux_u);  % simplified version
% rhsw = -dpdz + LIFT*(nz.*flux_u);  % simplified version
% return

p_jump = Vfqsurf*p_jump;
nu_jump = Vfqsurf*nu_jump;
nu_avg = Vfqsurf*nu_avg;
flux_p = (.5*p_jump - nu_avg).*sJq;
flux_u = .5*(nu_jump - p_jump).*sJq;

rhsp = -divU + (Pfq*(flux_p))./J;
rhsu = -dpdx + (Pfq*(nxq.*flux_u))./J;
rhsv = -dpdy + (Pfq*(nyq.*flux_u))./J;  
rhsw = -dpdz + (Pfq*(nzq.*flux_u))./J;  

% rhsp = -divU + Pfq*(Vfqsurf*flux_p);
% keyboard


function [rhsp rhsu rhsv rhsw] = WaveRHS3D(p,u,v,w)

Globals3D;

Drp = Dr*p; Dsp = Ds*p; Dtp = Dt*p;
Dru = Dr*u; Dsu = Ds*u; Dtu = Dt*u;
Drv = Dr*v; Dsv = Ds*v; Dtv = Dt*v;
Drw = Dr*w; Dsw = Ds*w; Dtw = Dt*w;

dpdx = rx.*Drp + sx.*Dsp + tx.*Dtp;
dpdy = ry.*Drp + sy.*Dsp + ty.*Dtp;
dpdz = rz.*Drp + sz.*Dsp + tz.*Dtp;

dudx = rx.*Dru + sx.*Dsu + tx.*Dtu;
dvdy = ry.*Drv + sy.*Dsv + ty.*Dtv;
dwdz = rz.*Drw + sz.*Dsw + tz.*Dtw;
divU = dudx + dvdy + dwdz;

p_jump = zeros(Nfp*Nfaces,K);
nu_jump = zeros(Nfp*Nfaces,K);

p_jump(:) = p(vmapP)-p(vmapM);
p_jump(mapB) = -2*p(vmapB); % boundary conditions

u_jump  = u(vmapP) - u(vmapM);
v_jump  = v(vmapP) - v(vmapM);
w_jump  = w(vmapP) - w(vmapM);
nu_jump(:) = nx(:).*u_jump + ny(:).*v_jump + nz(:).*w_jump;
nu_jump(mapB) = 0;

flux_p = .5*(p_jump - nu_jump).*Fscale;
flux_u = .5*(nu_jump - p_jump).*Fscale;

rhsp = -divU + LIFT*(flux_p);
rhsu = -dpdx + LIFT*(nx.*flux_u);  % simplified version
rhsv = -dpdy + LIFT*(ny.*flux_u);  % simplified version
rhsw = -dpdz + LIFT*(nz.*flux_u);  % simplified version

% keyboard

