function hybrid_wedge_driver

clearvars -global

% Driver script for solving the 3D wave equations
hybridgGlobals3D
hybridgGlobalFlags

% useSEM = 1;
useSEM = 0;
useLSC = 0;
nodalLIFT = 1;

global DWr DWs DWt rxW ryW rzW sxW syW szW txW tyW tzW
global wsJW nxW nyW nzW % quadrature
global Dx Dy Dz Ltri
global TESTING

% wedge_mesh1
% VZ(1) = -2;
% VZ(1) = VZ(1) + .25;
% VZ(6) = VZ(6) + .25;
wedge_mesh2
VZ(18) = VZ(18) + .25; % make wedge_mesh non-affine - match wedge_mesh2_perturb.msh

p = [4 5 6 1 2 3]; EToV(:,1:6) = EToV(:,p); % permute for wedges and jacobian > 0

% tet_wedge % FIX TOP BOUNDARY
%tet_wedge_flat
% gmsh_simple
% VY = VY + .25*VX;
% VX = VX + .5*VY;
% VZ = .5*VZ;

% % % % single elem
% u = [-1 1 -1 -1 1 -1]; v = [-1 -1 1 -1 -1 1]; w = [-1 -1 -1 1 1 1];
% VX = v(:); VY = w(:); VZ = u(:); % flipping coordinates for Gmsh
% tmp = VZ; VZ = VY(end:-1:1); VY = tmp(end:-1:1); % make vertical wedge, flip to preserve ordering
% VX = VX + .5*VY;  VY = VY + .25*VX;  VZ = VZ*.5; % rotate affinely
% VZ(1) = VZ(1) + .5; % make non-affine
% VZ(6) = VZ(6) + .5;
% EToV = 1:6; EToE = ones(1,5); EToF = 1:5;

K = size(EToV,1);
N = 4;

% hybridgStartUp
WedgeTetStartUp; % replaces LSC wedge with nodal wedge
keyboard
% drawMesh

% NqTri = length(rqW)/(N+1);
% ids = 1:NqTri;
% plot3(rqW(ids),tqW(ids),sz(ids,1),'o');
% [rt st] = Nodes2D(1); [rt st] = xytors(rt,st);
% Vt = Vandermonde2D(1,rqW(ids),tqW(ids));
% J(ids,1) - Vt*(Vt\J(ids,1))
% h = color_line3(rqW(ids),sqW(ids),tqW(ids),tz(ids,1),'.'); set(h,'markersize',32)

% build Ltri matrix for testing


%% check eigs
if 1
    Np = NpW; % assume all-wedge mesh
    A = zeros(4*Np*K);
    
    U = zeros(4*Np,K);
    ids = 1:Np*K;
    
    for i = 1:4*Np*K
        
        U(i) = 1;
        p = reshape(U(ids),Np,K);
        u = reshape(U(ids+Np*K),Np,K);
        v = reshape(U(ids+2*Np*K),Np,K);
        w = reshape(U(ids+3*Np*K),Np,K);
        U(i) = 0;
        
        [ps us vs ws] = NodalSurfaceInterp(p,u,v,w);
        %         [rhsp rhsu rhsv rhsw] = WaveRHS3D(p,u,v,w);
        [rhsp rhsu rhsv rhsw] = RHSall(p,u,v,w,ps,us,vs,ws);
        
        rhsAll = [rhsp(:); rhsu(:); rhsv(:); rhsw(:)];
        A(:,i) = rhsAll(:);
        if (mod(i,Np)==0)
            fprintf('computing column %d out of %d\n',i,4*Np*K)
        end
    end
    A(abs(A)<1e-8) = 0;
    [W D] = eig(A); lam = diag(D);
    for e = 1:K
        MW = VW'*diag(wJ(1:NcW,e))*VW;
        Mblk{e} = MW;
    end
    Mblk = blkdiag(Mblk{:});
    AM = blkdiag(Mblk,Mblk,Mblk,Mblk)*A; % undo inverse mass matrix
    
    hold on;plot(lam,'o','linewidth',2)
    hold on;plot(1i*linspace(min(imag(lam)),max(imag(lam)),100),'k--','linewidth',2)
    
    keyboard
end


%% set initial conditions

% ids = y > 0;
ids = y > -inf;

% exact sol
pexf = @(x,y,z,time) cos(pi/2*x).*cos(pi/2.*y).*cos(pi/2.*z)*cos(sqrt(3)*pi/2*time);
uexf = @(x,y,z,time) sin(pi/2*x).*cos(pi/2.*y).*cos(pi/2.*z)*sin(sqrt(3)*pi/2*time)*1/sqrt(3);
vexf = @(x,y,z,time) cos(pi/2*x).*sin(pi/2.*y).*cos(pi/2.*z)*sin(sqrt(3)*pi/2*time)*1/sqrt(3);
wexf = @(x,y,z,time) cos(pi/2*x).*cos(pi/2.*y).*sin(pi/2.*z)*sin(sqrt(3)*pi/2*time)*1/sqrt(3);

f = @(x,y,z) pexf(x,y,z,0);

p = zeros(NpMax,K);
p(1:NpW,wedgK) = f(xW,yW,zW); % interp for nodal wedges
p(1:NpT,tetK ) = f(xT,yT,zT); % interp for nodal tets

u = zeros(NpMax,K); v = zeros(NpMax,K); w = zeros(NpMax,K); % velocity

TESTING = 0;
if TESTING
    % for testing
    for ee = 1:length(tetK)
        p(1:NpT,tetK(ee)) = 1:NpT;
        u(1:NpT,tetK(ee)) = 1:NpT;
    end
    for ee = 1:length(wedgK)
        p(1:NpW,wedgK(ee)) = (1:NpW)    + ee;
        u(1:NpW,wedgK(ee)) = (1:NpW) +1 + ee;
        v(1:NpW,wedgK(ee)) = (1:NpW) +2 + ee;
        w(1:NpW,wedgK(ee)) = (1:NpW) +3 + ee;
    end
end

% get info of u at all cubature points
%[ps us vs ws] = SurfaceInterp(p,u,v,w);
[ps us vs ws] = NodalSurfaceInterp(p,u,v,w);


%% error in init cond

pexc = zeros(NcMax,K);
time = 0;
pexc(1:NcW,wedgK) = pexf(xqW, yqW, zqW, time);
pexc(1:NcT,tetK)  = pexf(xqT, yqT, zqT, time);

pc = zeros(NcMax,K);
pc(1:NcW,wedgK) = VW*p(1:NpW,wedgK);
pc(1:NcT,tetK ) = VT*p(1:NpT,tetK );

err2 = 0;
for e = 1:K
    diff = pexc(:,e) - pc(:,e);
    err2 = err2 + wJ(:,e)'*(diff.^2);
end
init_err = sqrt(err2);
disp(sprintf('Init cond L2 err = %f\n',init_err))


%%
resp = zeros(NpMax,K); resu = zeros(NpMax,K); resv = zeros(NpMax,K); resw = zeros(NpMax,K);

% % compute time step size using trace inequality of pyramid
% [C CK] = estimateCFL(); % CK = local CFLs
% dt = 1/C

% estimate just using Fscale - not quite right but seems to work ok in practice
CW = [ 9.926135933275306  18.563576705381948  29.030325215439685  42.988345972840094 58.802145509223905  78.006158337859318  99.271490513773074];
CT = [ 12.222306675423649  20.461046774606864  29.175757913298817  41.653800010846517 54.447306071124011  71.104042033768295  88.317462093087173];
FscaleT = CT(N)*Fscale(1:NfcT,tetK);
if isempty(FscaleT)
    FscaleT = 0;
end
FscaleW = CW(N)*Fscale(1:NfcW,wedgK);
dt = min(max(1./FscaleT(:)),max(1./FscaleW(:)))
% dt = dt/10
% dt = 0.00147547; % copied from c++, N=3, pri3_perturb
dt = 0.0025 % wedge_mesh2_perturb, N=1

% keyboard

FinalTime = .01;
% FinalTime = 5*dt;

useRK = 0
if useRK
    
    % outer time step loop
    time = 0; tstep = 0;
    up = zeros(NcMax,K);
    while time < FinalTime
        
        if(time+dt>FinalTime),
            dt = FinalTime-time;
        end
        
        % low storage RK
        for INTRK = 1:5
            
            [rhsp rhsu rhsv rhsw] = RHSall(p,u,v,w,ps,us,vs,ws);
            
            resp = rk4a(INTRK)*resp + dt*rhsp;
            resu = rk4a(INTRK)*resu + dt*rhsu;
            resv = rk4a(INTRK)*resv + dt*rhsv;
            resw = rk4a(INTRK)*resw + dt*rhsw;
            
            p = p + rk4b(INTRK)*resp;
            u = u + rk4b(INTRK)*resu;
            v = v + rk4b(INTRK)*resv;
            w = w + rk4b(INTRK)*resw;
            
            %         [ps us vs ws] = SurfaceInterp(p,u,v,w);
            [ps us vs ws] = NodalSurfaceInterp(p,u,v,w);
        end
        
        % plot nodal
        if (0 && mod(tstep,2)==0)
            %         up(1:NcW,wedgK) = abs(up(1:NcW,wedgK) - uexf(x(1:NcW,wedgK),y(1:NcW,wedgK),z(1:NcW,wedgK),time+dt));
            %         up(1:NcT,tetK)  = abs(up(1:NcT,tetK)  - uexf(x(1:NcT,tetK),y(1:NcT,tetK),z(1:NcT,tetK),time+dt));
            upW = abs(p(1:NpW,wedgK) - uexf(xW,yW,zW,time+dt));
            upT = abs(p(1:NpT,tetK) - uexf(xT,yT,zT,time+dt));
            clf
            %         drawMesh
            color_line3(xW,yW,zW,upW,'.');
            color_line3(xT,yT,zT,upT,'.');
            %         caxis([-1 1]);
            view(30,30);  xlabel('x'); ylabel('y'); colorbar
            title(sprintf('t = %f, max val = %f',time,max(abs(up(:)))))
            axis equal
            drawnow
        end
        if mod(tstep,50)==0
            fprintf('Timestep %d/%d\n',tstep,round(FinalTime/dt))
        end
        
        % Increment time
        time = time+dt; tstep = tstep+1;
        
    end
    
else
    % for AB timestepping: 2 previous histories. reuse v for rhsu2
    for i = 1:2
        rhsp_hist{i} = zeros(NpMax,K);
        rhsu_hist{i} = zeros(NpMax,K);
        rhsv_hist{i} = zeros(NpMax,K);
        rhsw_hist{i} = zeros(NpMax,K);
    end
    ab = [23/12 -16/12 5/12];
    
    % outer time step loop
    time = 0; tstep = 0;
    up = zeros(NcMax,K);
    while time < FinalTime
        
        if(time+dt>FinalTime), dt = (FinalTime-time);  end;
        
        if (tstep < 4)
            p(1:NpW,wedgK) = pexf(xW,yW,zW,time); % interp for nodal wedges
            u(1:NpW,wedgK) = uexf(xW,yW,zW,time); % interp for nodal wedges
            v(1:NpW,wedgK) = vexf(xW,yW,zW,time); % interp for nodal wedges
            w(1:NpW,wedgK) = wexf(xW,yW,zW,time); % interp for nodal wedges
            [ps us vs ws] = NodalSurfaceInterp(p,u,v,w);
            fprintf('on step %d and time %f, setting interpolant\n',tstep,time);
            e = 1;
            U = [p(1:NpW,e) u(1:NpW,e) v(1:NpW,e) w(1:NpW,e)]
        end
        
        [rhsp rhsu rhsv rhsw] = RHSall(p,u,v,w,ps,us,vs,ws);
        
%         U = [p(1:NpW,:);u(1:NpW,:);v(1:NpW,:);w(1:NpW,:)];
        rhs = [rhsp(1:NpW,:);rhsu(1:NpW,:);rhsv(1:NpW,:);rhsw(1:NpW,:)];
%         out = [repmat((0:NpW-1)',4,1) U(:,1) rhs(:,1)]; out(abs(out)<1e-8) = 0;        
%         [rhsp_hist{1}(1:NpW,:);rhsu_hist{1}(1:NpW,:);rhsv_hist{1}(1:NpW,:);rhsw_hist{1}(1:NpW,:)]
%         [rhsp_hist{2}(1:NpW,:);rhsu_hist{2}(1:NpW,:);rhsv_hist{2}(1:NpW,:);rhsw_hist{2}(1:NpW,:)]

        if tstep==3
            keyboard
        end


        p = p + dt*(ab(3)*rhsp_hist{1} + ab(2)*rhsp_hist{2} + ab(1)*rhsp);
        u = u + dt*(ab(3)*rhsu_hist{1} + ab(2)*rhsu_hist{2} + ab(1)*rhsu);
        v = v + dt*(ab(3)*rhsv_hist{1} + ab(2)*rhsv_hist{2} + ab(1)*rhsv);
        w = w + dt*(ab(3)*rhsw_hist{1} + ab(2)*rhsw_hist{2} + ab(1)*rhsw);        
        
        
        
        [ps us vs ws] = NodalSurfaceInterp(p,u,v,w);
        
        rhsp_hist{1} = rhsp_hist{2}; rhsp_hist{2} = rhsp;
        rhsu_hist{1} = rhsu_hist{2}; rhsu_hist{2} = rhsu;
        rhsv_hist{1} = rhsv_hist{2}; rhsv_hist{2} = rhsv;
        rhsw_hist{1} = rhsw_hist{2}; rhsw_hist{2} = rhsw;
        
        time = time+dt;
        tstep = tstep+1; % Increment time                
                
        
        % plot nodal
        if (0 && mod(tstep,2)==0)
            %         up(1:NcW,wedgK) = abs(up(1:NcW,wedgK) - uexf(x(1:NcW,wedgK),y(1:NcW,wedgK),z(1:NcW,wedgK),time+dt));
            %         up(1:NcT,tetK)  = abs(up(1:NcT,tetK)  - uexf(x(1:NcT,tetK),y(1:NcT,tetK),z(1:NcT,tetK),time+dt));
            upW = abs(p(1:NpW,wedgK) - uexf(xW,yW,zW,time+dt));
            upT = abs(p(1:NpT,tetK) - uexf(xT,yT,zT,time+dt));
            clf
            %         drawMesh
            color_line3(xW,yW,zW,upW,'.');
            color_line3(xT,yT,zT,upT,'.');
            %         caxis([-1 1]);
            view(30,30);  xlabel('x'); ylabel('y'); colorbar
            title(sprintf('t = %f, max val = %f',time,max(abs(up(:)))))
            axis equal
            drawnow
        end
        if mod(tstep,50)==0
            fprintf('Timestep %d/%d\n',tstep,round(FinalTime/dt))
        end
        
    end
    
    
end

%% 

% [rqW tqW sqW wW] = wedge_cub(N+1);
% VW = wedge_basis(N,rqW,sqW,tqW)/VWnodal;
% xqW = VW*xW; 
% yqW = VW*yW;
% zqW = VW*zW;
% NcW = length(rqW);
% wJ(1:NcW,wedgK) = diag(wW)*(VW*JW);

% compute error
pexc = zeros(NcMax,K);
pexc(1:NcW,wedgK) = pexf(xqW, yqW, zqW, time);
pexc(1:NcT,tetK)  = pexf(xqT, yqT, zqT, time);

pc = zeros(NcMax,K);
pc(1:NcW,wedgK) = VW*p(1:NpW,wedgK);
pc(1:NcT,tetK ) = VT*p(1:NpT,tetK );

err2 = 0;
err2K = zeros(K,1);
norm2 = 0;
for e = 1:K
    diff = pexc(:,e) - pc(:,e);
    err2K(e) = wJ(:,e)'*(diff.^2);
    norm2K = wJ(:,e)'*(pc(:,e).^2);
    err2 = err2 + err2K(e);
    norm2 = norm2 + wJ(:,e)'*(pc(:,e).^2);
end
err = sqrt(err2);
p_norm = sqrt(norm2);
rel_err = err/p_norm;
fprintf('Abs L2 err = %f, rel L2 err = %f\n',err,rel_err)

%% velocity error

uexc = zeros(NcMax,K);
% uexf = @(x,y,z,time) cos(pi/2*x).*cos(pi/2.*y).*cos(pi/2.*z)*cos(sqrt(3)*pi/2*time);

uexc = zeros(NcMax,K);
vexc = zeros(NcMax,K);
wexc = zeros(NcMax,K);
uexc(1:NcW,wedgK) = uexf(xqW, yqW, zqW, time);
vexc(1:NcW,wedgK) = vexf(xqW, yqW, zqW, time);
wexc(1:NcW,wedgK) = wexf(xqW, yqW, zqW, time);

uc = zeros(NcMax,K);
vc = zeros(NcMax,K);
wc = zeros(NcMax,K);
uc(1:NcW,wedgK) = VW*u(1:NpW,wedgK);
vc(1:NcW,wedgK) = VW*v(1:NpW,wedgK);
wc(1:NcW,wedgK) = VW*w(1:NpW,wedgK);

err2u = 0;
err2v = 0;
err2w = 0;
for e = 1:K
    diffu = uexc(:,e) - uc(:,e);
    err2u = err2u + wJ(:,e)'*(diffu.^2);
    diffv = vexc(:,e) - vc(:,e);
    err2v = err2v + wJ(:,e)'*(diffv.^2);
    diffw = wexc(:,e) - wc(:,e);
    err2w = err2w + wJ(:,e)'*(diffw.^2);
end
erru = sqrt(err2u);
errv = sqrt(err2v);
errw = sqrt(err2w);
fprintf('L2 err = %f, %f, %f\n',erru,errv,errw)

%%

keyboard


% computes RHS for all element types
function [rhsp rhsu rhsv rhsw] = RHSall(p,u,v,w,ps,us,vs,ws)

hybridgGlobals3D
hybridgGlobalFlags

rhsp = zeros(size(p));
rhsu = zeros(size(p));
rhsv = zeros(size(p));
rhsw = zeros(size(p));

% wedges
VWfnodal = zeros(length(SmaskW),NpW);
[rhspK rhsuK rhsvK rhswK] = Nodal_wedge_RHS(p,u,v,w,ps,us,vs,ws,VW,VWr,VWs,VWt,VWfnodal,wedgK);
rhsp(1:NpW,wedgK) = rhspK; rhsu(1:NpW,wedgK) = rhsuK;
rhsv(1:NpW,wedgK) = rhsvK; rhsw(1:NpW,wedgK) = rhswK;

% tet - switch w/nodal DG
[rhspK rhsuK rhsvK rhswK] = Nodal_tet_RHS(p,u,v,w,ps,us,vs,ws,VTr,VTs,VTt,VTf);
rhsp(1:NpT,tetK ) = rhspK; rhsu(1:NpT,tetK ) = rhsuK;
rhsv(1:NpT,tetK ) = rhsvK; rhsw(1:NpT,tetK ) = rhswK;

% nodal RHS for tets: "collocation" DG
function [rhsp rhsu rhsv rhsw] = Nodal_tet_RHS(p,u,v,w, ...
    p_surface,u_surface,v_surface,w_surface,...
    Dr,Ds,Dt,Vf)

hybridgGlobals3D;

if nnz(tetK)==0
    rhsp = []; rhsu = []; rhsv = []; rhsw = [];
    return
end

% element-type specifics
Np = NpT; Nfc = size(Vf,1);

% face cubature
% wsJK = wsJ(1:Nfc,tetK);
nxK = nx(1:Nfc,tetK); nyK = ny(1:Nfc,tetK); nzK = nz(1:Nfc,tetK);

% nodal dofs
pK = p(1:Np,tetK); uK = u(1:Np,tetK); vK = v(1:Np,tetK); wK = w(1:Np,tetK);

[p_jump nu_jump] = SurfaceFlux(p_surface,u_surface,v_surface,w_surface,Vf,tetK);
flux_p = .5*(p_jump - nu_jump);
flux_u = .5*(nu_jump - p_jump);

% volume derivative operators
dudx = rxT.*(Dr*uK) + sxT.*(Ds*uK) + txT.*(Dt*uK);
dvdy = ryT.*(Dr*vK) + syT.*(Ds*vK) + tyT.*(Dt*vK);
dwdz = rzT.*(Dr*wK) + szT.*(Ds*wK) + tzT.*(Dt*wK);
divU = dudx + dvdy + dwdz;

pr = Dr*pK; ps = Ds*pK; pt = Dt*pK;
dpdx = rxT.*pr + sxT.*ps + txT.*pt;
dpdy = ryT.*pr + syT.*ps + tyT.*pt;
dpdz = rzT.*pr + szT.*ps + tzT.*pt;

FscaleK = Fscale(1:NfcT,tetK);
flux_uc = FscaleK.*flux_u;
rhsp = -divU + LIFTT*(FscaleK.*flux_p);
rhsu = -dpdx + LIFTT*(nxK.*flux_uc);
rhsv = -dpdy + LIFTT*(nyK.*flux_uc);
rhsw = -dpdz + LIFTT*(nzK.*flux_uc);

% keyboard
% [rhsp;rhsu;rhsv;rhsw]
% [LIFTT*(FscaleK.*flux_p); LIFTT*(nxK.*flux_uc);LIFTT*(nyK.*flux_uc);LIFTT*(nzK.*flux_uc)]

return


function [rhsp rhsu rhsv rhsw] = Nodal_wedge_RHS(p,u,v,w, ...
    p_surface,u_surface,v_surface,w_surface,...
    V,Vr,Vs,Vt,Vf,typeK)

hybridgGlobals3D;
global DWr DWs DWt rxW ryW rzW sxW syW szW txW tyW tzW
global Dx Dy Dz Ltri
global wsJW nxW nyW nzW % quadrature

if nnz(typeK)==0
    rhsp = []; rhsu = []; rhsv = []; rhsw = [];
    return
end

% element-type specifics
Nc = size(V,1); Np = size(V,2); Nfc = size(Vf,1);

nxK = nx(1:Nfc,typeK); nyK = ny(1:Nfc,typeK); nzK = nz(1:Nfc,typeK);
pK = p(1:Np,typeK); uK = u(1:Np,typeK); vK = v(1:Np,typeK); wK = w(1:Np,typeK);

rxK = rx(1:Nc,typeK); sxK = sx(1:Nc,typeK); txK = tx(1:Nc,typeK);
ryK = ry(1:Nc,typeK); syK = sy(1:Nc,typeK); tyK = ty(1:Nc,typeK);
rzK = rz(1:Nc,typeK); szK = sz(1:Nc,typeK); tzK = tz(1:Nc,typeK);
wJK = wJ(1:Nc,typeK);

% opt = 'skew';
% opt = 'quadrature';
opt = 'weak derivs';
% opt = 'nodal';
switch opt
    case 'weak derivs'
        
        for e = 1:length(typeK)
            dpdx(:,e) = Dx{e}*pK(:,e);
            dpdy(:,e) = Dy{e}*pK(:,e);
            dpdz(:,e) = Dz{e}*pK(:,e);
            dudx(:,e) = Dx{e}*uK(:,e);
            dvdy(:,e) = Dy{e}*vK(:,e);
            dwdz(:,e) = Dz{e}*wK(:,e);
        end
        divU = dudx + dvdy + dwdz;
%         keyboard
        
    case 'quadrature'
        
        ur = VWr*uK; us = VWs*uK; ut = VWt*uK;
        vr = VWr*vK; vs = VWs*vK; vt = VWt*vK;
        wr = VWr*wK; ws = VWs*wK; wt = VWt*wK;
        ux = VW'*(wJK.*(rxK.*ur + sxK.*us + txK.*ut));
        vy = VW'*(wJK.*(ryK.*vr + syK.*vs + tyK.*vt));
        wz = VW'*(wJK.*(rzK.*wr + szK.*ws + tzK.*wt));
        
        pr = VWr*pK; ps = VWs*pK; pt = VWt*pK;
        px = VW'*(wJK.*(rxK.*pr + sxK.*ps + txK.*pt));
        py = VW'*(wJK.*(ryK.*pr + syK.*ps + tyK.*pt));
        pz = VW'*(wJK.*(rzK.*pr + szK.*ps + tzK.*pt));
        
        uvw = ux + vy + wz;
        for e = 1:length(typeK)
            dU = MW{e}\[px(:,e) py(:,e) pz(:,e) uvw(:,e)];
            dpdx(:,e) = dU(:,1);
            dpdy(:,e) = dU(:,2);
            dpdz(:,e) = dU(:,3);
            divU(:,e) = dU(:,4);
        end
        
    case 'nodal'
        
        % volume derivative operators
        pr = DWr*pK; ps = DWs*pK; pt = DWt*pK;
        dpdx = rxW.*pr + sxW.*ps + txW.*pt;
        dpdy = ryW.*pr + syW.*ps + tyW.*pt;
        dpdz = rzW.*pr + szW.*ps + tzW.*pt;
        
        dudx = rxW.*(DWr*uK) + sxW.*(DWs*uK) + txW.*(DWt*uK);
        dvdy = ryW.*(DWr*vK) + syW.*(DWs*vK) + tyW.*(DWt*vK);
        dwdz = rzW.*(DWr*wK) + szW.*(DWs*wK) + tzW.*(DWt*wK);
        divU = dudx + dvdy + dwdz;
        
end

rhsp = -divU;
rhsu = -dpdx;
rhsv = -dpdy;
rhsw = -dpdz;

[p_jump nu_jump nu_avg] = SurfaceFlux(p_surface,u_surface,v_surface,w_surface,Vf,typeK);
tau = 0;
flux_p = .5*(tau*p_jump - nu_jump);
flux_u = .5*(tau*nu_jump - p_jump);

for e = 1:length(typeK)
    rhsp(:,e) = rhsp(:,e) + LIFTW{e}*flux_p(:,e);
    rhsu(:,e) = rhsu(:,e) + LIFTW{e}*(nxK(:,e).*flux_u(:,e));
    rhsv(:,e) = rhsv(:,e) + LIFTW{e}*(nyK(:,e).*flux_u(:,e));
    rhsw(:,e) = rhsw(:,e) + LIFTW{e}*(nzK(:,e).*flux_u(:,e));
end

%% checking application of ops
if 0
    e = 2;
    
    J = J(1:NcW,:);
    NqTri = size(J,1)/(N+1);
    Jtri = reshape(J(:,e),NqTri,N+1);
    Jtri = Jtri(:,e); % constant in extruded direction
    [rqtri sqtri wtri] = Cubature2D(2*N+1);
    [rtri ttri] = Nodes2D(N); [rtri ttri] = xytors(rtri,ttri); % warning - only valid for N = 1,2
    Vtri = Vandermonde2D(N,rtri,ttri);
    Vqtri = Vandermonde2D(N,rqtri,sqtri)/Vtri;
    Mtri = Vqtri'*diag(wtri.*Jtri)*Vqtri;
    Mhat_tri = Vqtri'*diag(wtri)*Vqtri;
    Ltri = Mtri\Mhat_tri;
    
    sxJ = sxW.*JW;
    sxJ = reshape(sxJ(:,e),(N+1)*(N+2)/2,N+1);
    sxJ = sxJ(1,:);
    
    syJ = syW.*JW;
    syJ = reshape(syJ(:,e),(N+1)*(N+2)/2,N+1);
    syJ = syJ(1,:);
    
    szJ = szW.*JW;
    szJ = szJ(1,e); % constant everywhere in elem
    
    [rtri ttri] = Nodes2D(N); [rtri ttri] = xytors(rtri,ttri); % warning - only valid for N = 1,2
    Vtri = Vandermonde2D(N,rtri,ttri);
    [Vrtri Vttri] = GradVandermonde2D(N,rtri,ttri);
    Drtri = Vrtri/Vtri;
    Dttri = Vttri/Vtri;
    D1D = GradVandermonde1D(N,JacobiGL(0,0,N))/Vandermonde1D(N,JacobiGL(0,0,N));
    
    norm(rxW(1,e)*kron(eye(N+1),Drtri) + txW(1,e)*kron(eye(N+1),Dttri) + kron(diag(sxJ)*D1D,Ltri) - Dx{e},'fro')
    norm(kron(eye(N+1),rxW(1,e)*Drtri + txW(1,e)*Dttri) + kron(diag(sxJ)*D1D,Ltri) - Dx{e},'fro')
    norm(ryW(1,e)*kron(eye(N+1),Drtri) + tyW(1,e)*kron(eye(N+1),Dttri) + kron(diag(syJ)*D1D,Ltri) - Dy{e},'fro')
    norm(                                                            kron(szJ*D1D,Ltri) - Dz{e},'fro')
end

%%
global TESTING
if TESTING
    e = 2;
    keyboard
%     [-divU;-dpdx;-dpdy;-dpdz]
    [-divU(:,e) -dpdx(:,e) -dpdy(:,e) -dpdz(:,e)]
    [LIFTW{e}*flux_p(:,e) LIFTW{e}*(nxK(:,e).*flux_u(:,e)) LIFTW{e}*(nyK(:,e).*flux_u(:,e)) LIFTW{e}*(nzK(:,e).*flux_u(:,e))]
    [rhsp; rhsu; rhsv; rhsw]
end

return




% independent of regular DG or LSC-DG
function [p_jump nu_jump nu_avg] = SurfaceFlux(p_surface,u_surface,v_surface,w_surface,...
    Vf,typeK)

hybridgGlobals3D;

if nnz(typeK)==0
    p_jump = [];    nu_jump = [];    nu_avg = [];
    return
end

Nfc = size(Vf,1); Np = size(Vf,2);
nxK = nx(1:Nfc,typeK); nyK = ny(1:Nfc,typeK); nzK = nz(1:Nfc,typeK);

mapMK = mapM(1:Nfc,typeK); mapPK = mapP(1:Nfc,typeK);

% flux gather
p_jump  = p_surface(mapPK) - p_surface(mapMK);
u_jump  = u_surface(mapPK) - u_surface(mapMK);
v_jump  = v_surface(mapPK) - v_surface(mapMK);
w_jump  = w_surface(mapPK) - w_surface(mapMK);
nu_jump = nxK.*u_jump + nyK.*v_jump + nzK.*w_jump;

u_avg  = .5*(u_surface(mapPK) + u_surface(mapMK));
v_avg  = .5*(v_surface(mapPK) + v_surface(mapMK));
w_avg  = .5*(w_surface(mapPK) + w_surface(mapMK));
nu_avg = nxK.*u_avg + nyK.*v_avg + nzK.*w_avg;

% apply BCs, todo: try p^+ = - p ^ -, which is energy stable
mapBK = mapMK == mapPK;
p_jump(mapBK) = -2*p_surface(mapMK(mapBK));

% nu_jump(mapBK) = 0; % extra neumann BCs??

% function [ps us vs ws] = SurfaceInterp(p,u,v,w)
%
% hybridgGlobals3D
% hybridgGlobalFlags
%
% ps = zeros(NfcMax,K);
% us = zeros(NfcMax,K);
% vs = zeros(NfcMax,K);
% ws = zeros(NfcMax,K);
%
% ps(1:NfcW,wedgK) = VWf*p(1:NpW,wedgK);
% us(1:NfcW,wedgK) = VWf*u(1:NpW,wedgK);
% vs(1:NfcW,wedgK) = VWf*v(1:NpW,wedgK);
% ws(1:NfcW,wedgK) = VWf*w(1:NpW,wedgK);
%
% ps(1:NfcT,tetK ) = VTf*p(1:NpT,tetK );
% us(1:NfcT,tetK ) = VTf*u(1:NpT,tetK );
% vs(1:NfcT,tetK ) = VTf*v(1:NpT,tetK );
% ws(1:NfcT,tetK ) = VTf*w(1:NpT,tetK );


function [ps us vs ws] = NodalSurfaceInterp(p,u,v,w)

hybridgGlobals3D
hybridgGlobalFlags

ps = zeros(NfcMax,K);
us = zeros(NfcMax,K);
vs = zeros(NfcMax,K);
ws = zeros(NfcMax,K);

ps(1:NfcW,wedgK) = p(SmaskW,wedgK);
us(1:NfcW,wedgK) = u(SmaskW,wedgK);
vs(1:NfcW,wedgK) = v(SmaskW,wedgK);
ws(1:NfcW,wedgK) = w(SmaskW,wedgK);

ps(1:NfcT,tetK ) = p(SmaskT,tetK );
us(1:NfcT,tetK ) = u(SmaskT,tetK );
vs(1:NfcT,tetK ) = v(SmaskT,tetK );
ws(1:NfcT,tetK ) = w(SmaskT,tetK );

