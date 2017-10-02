% nodal wave RK driver

function Wave_RK

% Driver script for solving the 3D IPDGtion equations
Globals3D;

% Order of polymomials used for approximation
N = 3;
FinalTime = .5;

% % % % single element
% [VX VY VZ] = Nodes3D(1); [VX VY VZ] = xyztorst(VX,VY,VZ); K = 1; EToV = 1:length(VX);

% GMSH meshes
% cubeTetra1
filename = 'Grid/cubeTetra1.msh';
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

% %% bernstein conversion
% 
% Vp = bern_basis_tet(N,rp,sp,tp); 
% Vq = bern_basis_tet(N,rq,sq,tq); 
% 
% [V Vr Vs Vt V1 V2 V3 V4] = bern_basis_tet(N,r,s,t);
% Dr = V\Vr; Dr(abs(Dr)<1e-8) = 0;
% Ds = V\Vs; Ds(abs(Ds)<1e-8) = 0;
% Dt = V\Vt; Dt(abs(Dt)<1e-8) = 0;
% 
% global D1 D2 D3 D4
% D1 = V\V1; D1(abs(D1)<1e-8) = 0;
% D2 = V\V2; D2(abs(D2)<1e-8) = 0;
% D3 = V\V3; D3(abs(D3)<1e-8) = 0;
% D4 = V\V4; D4(abs(D4)<1e-8) = 0;
% 
% %% BB conversion
% 
% % M = inv(V*V');  
% % invM = V*V';
% M = Vq'*diag(wq)*Vq;
% invM = inv(M);
% 
% %% lift matrix
% Emat = zeros(Np, Nfaces*Nfp);
% Mf = zeros(Np, Nfp);
% 
% for face=1:Nfaces   
%     % WARNING: ASSUMES ORDERING OF FACE NODES/FMASK
%     [rt st] = Nodes2D(N); [rt st] = xytors(rt,st); 
%     [rfq sfq wf] = Cubature2D(2*N); 
% 
%     Vfq = bern_basis_tri(N,rfq,sfq);
%     massFace = Vfq'*diag(wf)*Vfq;
%     
%     idr = Fmask(:,face);
%     idc = (face-1)*Nfp+1:face*Nfp;
%     
%     Emat(idr, idc) = Emat(idr, idc) + massFace;
%     if face==1
%         Mf(idr,idc)=massFace;
%     end
% end
% 
% % inv(mass matrix)*\I_n (L_i,L_j)_{edge_n}
% LIFT = invM*Emat;
% LIFT(abs(LIFT)<1e-8) = 0;

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
        
        [rhsp rhsu rhsv rhsw] = WaveRHS3D(p,u,v,w);
        
        resp = rk4a(INTRK)*resp + dt*rhsp;        
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        resw = rk4a(INTRK)*resw + dt*rhsw;
        
        p = p + rk4b(INTRK)*resp;
        u = u + rk4b(INTRK)*resu;
        v = v + rk4b(INTRK)*resv;
        w = w + rk4b(INTRK)*resw;
       
%         keyboard
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
for e = 1:K
    err(:,e) = err(:,e).*wq*JK(e);
end
err = sqrt(sum(err(:)));

err


function [rhsp rhsu rhsv rhsw] = WaveRHS3D(p,u,v,w)
% Purpose  : Evaluate RHS flux in 3D DG

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

