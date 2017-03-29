function Wave_RK_BB

% Driver script for solving the 3D IPDGtion equations
Globals3D;

% Order of polymomials used for approximation
N = 2;
FinalTime = 0;

% % % % single element
% [VX VY VZ] = Nodes3D(1); [VX VY VZ] = xyztorst(VX,VY,VZ); K = 1; EToV = 1:length(VX);

% GMSH meshes
cubeTetra1

% VX = 2*VX; VY = 2*VY; VZ = 2*VZ; % biunit cube

% Initialize solver and construct grid and metric
StartUp3D;
K

keyboard

% LIFTnodal = LIFT;

%% cubature and plotting

JK = J(1,:);

[rq sq tq wq] = tet_cubature(2*N+1);
Vq = Vandermonde3D(N, rq, sq, tq)/V;
xq = Vq*x; yq = Vq*y; zq = Vq*z; % phys cubature nodes

% for plotting - build coordinates of all the nodes
[rp sp tp] = EquiNodes3D(25);
Vp = Vandermonde3D(N,rp,sp,tp)*invV;
xp = Vp*x; yp = Vp*y; zp = Vp*z;

%% bernstein conversion

Vp = bern_basis_tet(N,rp,sp,tp); 
Vq = bern_basis_tet(N,rq,sq,tq); 

[V Vr Vs Vt V1 V2 V3 V4] = bern_basis_tet(N,r,s,t);
Dr = V\Vr; Dr(abs(Dr)<1e-8) = 0;
Ds = V\Vs; Ds(abs(Ds)<1e-8) = 0;
Dt = V\Vt; Dt(abs(Dt)<1e-8) = 0;

global D1 D2 D3 D4
D1 = V\V1; D1(abs(D1)<1e-8) = 0;
D2 = V\V2; D2(abs(D2)<1e-8) = 0;
D3 = V\V3; D3(abs(D3)<1e-8) = 0;
D4 = V\V4; D4(abs(D4)<1e-8) = 0;

%% BB conversion

% M = inv(V*V');  
% invM = V*V';
M = Vq'*diag(wq)*Vq;
invM = inv(M);

%% lift matrix
Emat = zeros(Np, Nfaces*Nfp);
Mf = zeros(Np, Nfp);

for face=1:Nfaces   
    % WARNING: ASSUMES ORDERING OF FACE NODES/FMASK
    [rt st] = Nodes2D(N); [rt st] = xytors(rt,st); 
    [rfq sfq wf] = Cubature2D(2*N); 

    Vfq = bern_basis_tri(N,rfq,sfq);
    massFace = Vfq'*diag(wf)*Vfq;
    
    idr = Fmask(:,face);
    idc = (face-1)*Nfp+1:face*Nfp;
    
    Emat(idr, idc) = Emat(idr, idc) + massFace;
    if face==1
        Mf(idr,idc)=massFace;
    end
end

% inv(mass matrix)*\I_n (L_i,L_j)_{edge_n}
LIFT = invM*Emat;
LIFT(abs(LIFT)<1e-8) = 0;

%% EEL matrix

LIFTf = LIFT(:,1:Nfp);

% get permutations of rows for 2nd-4th cols of LIFT
p = zeros(Np,Nfaces-1);
for f = 1:Nfaces-1
    ids = (1:Nfp) + f*Nfp;
    for i = 1:Np % match rows of first cols with second       
        diff = kron(ones(Np,1),LIFTf(i,:)) - LIFT(:,ids);
        [~,j] = min(sum(abs(diff),2));
        p(j,f) =i;
    end
end

[r2D s2D] = Nodes2D(N); [r2D s2D] = xytors(r2D,s2D);

% 2d degree ops
ED = {}; 
for i = 0:N     
    EDi = bern_basis_tri(N,r2D,s2D)\bern_basis_tri(N-i,r2D,s2D);
    EDi(abs(EDi)<1e-8) = 0;
    ED{i+1} = EDi;    
end
LIFTf = LIFT(:,1:Nfp);
EE = blkdiag(ED{:}); EE = EE';

% get l_i scalings
R = EE*kron(ones(N+1,1),LIFTf(1:Nfp,1:Nfp))./LIFTf;  % ratio
block_starts = []; id = 1;
for i = 0:N
    block_starts = [block_starts id];
    Ni = N-i; Nfpi = (Ni+1)*(Ni+2)/2; % decreasing block sizes
    id = id + Nfpi;
end
R = R(:,1);
l = 1./R(block_starts); %R = R(~isnan(R)); [~, perm] = uniquetol(R); l = 1./R(sort(perm)); % scaling constants
EDL = cell(N+1,1);
for i = 0:N
    EDL{i+1} = ED{i+1}*l(i+1);
end

% built in scalings
EL = [];
for i = 0:N
    EL = [EL;EDL{i+1}'];
end

global L0 EEL
EEL=([EL EL(p(:,1),:) EL(p(:,2),:) EL(p(:,3),:)]);
L0 = LIFTf(1:Nfp,1:Nfp);

%% 

uex = @(x,y,z,t) cos(pi*x).*cos(pi*y).*cos(pi*z)*cos(sqrt(3)*pi*t);

% compute time step size
dt = .125/max( ((N+1)^2)*Fscale(:));

% set initial conditions
p = V\uex(x,y,z,0);
u = zeros(Np,K); v = zeros(Np,K); w = zeros(Np,K); % zero velocities

% p(:) = 0:Np*K-1; p = V\p; % TESTING
keyboard

% outer time step loop
time = 0;
tstep = 0;

totalSteps = floor(FinalTime/dt);

resp = zeros(Np,K);
resu = zeros(Np,K); 
resv = zeros(Np,K);
resw = zeros(Np,K);

while (time<FinalTime)
    
    if(time+dt>FinalTime), dt = (FinalTime-time);  end;
        
    % low storage RK
    for INTRK = 1:5        
        
        [rhsp rhsu rhsv rhsw] = WaveRHS3D(p,u,v,w);
        
        resp = rk4a(INTRK)*resp + dt*rhsp;        
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        resw = rk4a(INTRK)*resw + dt*rhsw;
        
        p = p + rk4b(INTRK)*resp;
        u = u + rk4b(INTRK)*resu;
        v = v + rk4b(INTRK)*resv;
        w = w + rk4b(INTRK)*resw;
       
    end;

    time = time+dt;
    tstep = tstep+1; % Increment time
         
    if (mod(tstep,10)==0)
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
err = (Vq*p-uex(xq,yq,zq,time)).^2;
for e = 1:K
    err(:,e) = err(:,e).*wq*JK(e);
end
err = sqrt(sum(err(:)));

err


function [rhsp rhsu rhsv rhsw] = WaveRHS3D(p,u,v,w)
% Purpose  : Evaluate RHS flux in 3D DG

Globals3D;

% BB specifics
global D1 D2 D3 D4
global L0 EEL

D1p = D1*p; D2p = D2*p; D3p = D3*p; D4p = D4*p;
Drp = .5*(D2p-D1p); 
Dsp = .5*(D3p-D1p); 
Dtp = .5*(D4p-D1p);

D1u = D1*u; D2u = D2*u; D3u = D3*u; D4u = D4*u;
Dru = .5*(D2u-D1u); 
Dsu = .5*(D3u-D1u); 
Dtu = .5*(D4u-D1u);

D1v = D1*v; D2v = D2*v; D3v = D3*v; D4v = D4*v;
Drv = .5*(D2v-D1v); 
Dsv = .5*(D3v-D1v); 
Dtv = .5*(D4v-D1v);

D1w = D1*w; D2w = D2*w; D3w = D3*w; D4w = D4*w;
Drw = .5*(D2w-D1w); 
Dsw = .5*(D3w-D1w); 
Dtw = .5*(D4w-D1w);

% Drp = Dr*p; Dsp = Ds*p; Dtp = Dt*p;
% Dru = Dr*u; Dsu = Ds*u; Dtu = Dt*u;
% Drv = Dr*v; Dsv = Ds*v; Dtv = Dt*v;
% Drw = Dr*w; Dsw = Ds*w; Dtw = Dt*w;

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

%rhsp = -divU + LIFT*(flux_p);
%rhsu = -dpdx + LIFT*(nx.*flux_u);  % simplified version
%rhsv = -dpdy + LIFT*(ny.*flux_u);  % simplified version
%rhsw = -dpdz + LIFT*(nz.*flux_u);  % simplified version

L0p = kron(eye(Nfaces),L0)*flux_p; 
Lu = kron(eye(Nfaces),L0)*flux_u; 

rhsp = -divU + EEL*L0p;
rhsu = -dpdx + EEL*(nx.*Lu);  % simplified version
rhsv = -dpdy + EEL*(ny.*Lu);  % simplified version
rhsw = -dpdz + EEL*(nz.*Lu);  % simplified version
% keyboard
