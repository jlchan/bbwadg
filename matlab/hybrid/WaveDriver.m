function WaveDriver

% Driver script for solving the 3D wave equations
hybridgGlobals3D
hybridgGlobalFlags
useSEM = 1; useLSC = 1; useNodalTets = 1; useSkewRHS = 0;

% Cube01; VX = 2*VX-1; VY = 2*VY-1; VZ = 2*VZ-1;
prism_mesh
% prism_Twist2
% hybrid_mesh
% tet_mesh
hybrid_mesh2; EToV = EToV'; EToE = EToE'; EToF = EToF'; 
% hex_mesh

% VX = [    -1     1    1     -1    -1     1    1     -1]'; 
% VY = [    -1    -1     1     1    -1    -1     1     1]';  
% VZ = [    -1    -1    -1    -1     1     1     1     1]';
% EToV = 1:8; EToE = 0*EToV; EToF = 0*EToV;

% pyr_mesh
% prism_mesh2
wedge_mesh
% EToV(EToV==0) = -1; EToV = EToV'; EToE = EToE'; EToF = EToF'; 

K = size(EToV,1);

N = 4;

% % max 8 vertices
% [VX VY VZ EToVP EToEP EToFP] = makeHexPyrMesh(2);  K = size(EToVP,1);
% EToV = -ones(K,8); EToV(:,1:5) = EToVP;
% EToE = zeros(K,6); EToE(:,1:5) = EToEP;
% EToF = zeros(K,6); EToF(:,1:5) = EToFP;
% for e = 1:K
%     for f = 1:5 % # faces on pyr
%         if EToE(e,f) == e
%             EToE(e,f) = 0;  EToF(e,f) = 0;
%         end
%     end
% end

% a = .05;
% ids = abs(abs(VX)-1)>1e-4;
% dV = (1:length(VX(ids)))';
% VX(ids) = VX(ids) + a*dV;
% ids = abs(abs(VY)-1)>1e-4;
% VY(ids) = VY(ids) + a*dV;
% ids = abs(abs(VZ)-1)>1e-4;
% VZ(ids) = VZ(ids) + a*dV;

hybridgStartUp

% drawMesh; 

% for e = 1:K
%     v = EToV(e,:); v = v(v > 0); NvK = nnz(v);
%     clf
%     switch NvK
%         case 4
%             NfcK = NfcT;
%         case 5 
%             NfcK = NfcP;
%         case 6 
%             NfcK = NfcW;
%         case 8            
%             NfcK = NfcH;
%     end
%     xfK = xf(1:NfcK,e); yfK = yf(1:NfcK,e); zfK = zf(1:NfcK,e);    
%     nxK = nx(1:NfcK,e); nyK = ny(1:NfcK,e); nzK = nz(1:NfcK,e);
%     quiver3(xfK,yfK,zfK,nxK,nyK,nzK);hold on
%     plot3(xfK,yfK,zfK,'.')
%     drawMesh(e)
%     axis off
%     pause
% end

%% set initial conditions

ids = y > 0;

% exact sol
uexf = @(x,y,z,time) cos(pi/2*x).*cos(pi/2.*y).*cos(pi/2.*z)*cos(sqrt(3)*pi/2*time);
% uexf = @(x,y,z,time) sin(pi*x).*sin(pi.*y).*sin(pi.*z)*cos(sqrt(3)*pi*time);
f = @(x,y,z) uexf(x,y,z,0);

b = zeros(NpMax,K);
b(1:NpH,hexK ) = VH'*(wJ(1:NcH,hexK ).*f(xqH,yqH,zqH));
if useLSC 
    b(1:NpW,wedgK) = VW'*(sqJW.*wJ(1:NcW,wedgK).*f(xqW,yqW,zqW));
else
    b(1:NpW,wedgK) = VW'*(wJ(1:NcW,wedgK).*f(xqW,yqW,zqW));
end
b(1:NpP,pyrK ) = VP'*(wJ(1:NcP,pyrK ).*f(xqP,yqP,zqP));
b(1:NpT,tetK ) = VT'*(wJ(1:NcT,tetK ).*f(xqT,yqT,zqT));

p = invM.*b; % nodal/LSC inverse = Identity built into invM - extra storage but OK
if useNodalTets    
    %p(1:NpT,tetK ) = invMThat*(b(1:NpT,tetK)./JTet);
    p(1:NpT,tetK ) = f(xT,yT,zT);
end
if ~useLSC % invert full wedge mass matrix
    for ee = 1:length(wedgK)
        e = wedgK(ee);
        p(1:NpW,e) = MW{ee}\b(1:NpW,e);
    end
end
u = zeros(NpMax,K); v = zeros(NpMax,K); w = zeros(NpMax,K); % velocity

% get info of u at all cubature points
[ps us vs ws] = SurfaceInterp(p,u,v,w);

%% error in init cond
uexc = zeros(NcMax,K);
time = 0;
uexc(1:NcH,hexK)  = uexf(xqH, yqH, zqH, time);
uexc(1:NcW,wedgK) = uexf(xqW, yqW, zqW, time);
uexc(1:NcP,pyrK)  = uexf(xqP, yqP, zqP, time);
uexc(1:NcT,tetK)  = uexf(xqT, yqT, zqT, time);

pc = zeros(NcMax,K);
pc(1:NcH,hexK ) = VH*p(1:NpH,hexK );
if useLSC
    pc(1:NcW,wedgK) = (VW*p(1:NpW,wedgK)).*sqJW;
else
    pc(1:NcW,wedgK) = VW*p(1:NpW,wedgK);
end
pc(1:NcP,pyrK ) = VP*p(1:NpP,pyrK );
pc(1:NcT,tetK ) = VT*p(1:NpT,tetK );

err2 = 0;
for e = 1:K       
    diff = uexc(:,e) - pc(:,e);
    err2 = err2 + wJ(:,e)'*(diff.^2);
end
init_err = sqrt(err2);
disp(sprintf('Init cond L2 err = %f\n',init_err))

%% 
resp = zeros(NpMax,K); resu = zeros(NpMax,K); resv = zeros(NpMax,K); resw = zeros(NpMax,K);

% compute time step size using trace inequality of pyramid
[C CK] = estimateCFL(); % CK = local CFLs
dt = 1/C
% dt = 1/max(((N+1)*(N+3))*(2/3)*Fscale(:));

FinalTime = 5*dt;

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
        
        [ps us vs ws] = SurfaceInterp(p,u,v,w);
    end
         
    if (mod(tstep,50)==0)
        up(1:NcH,hexK ) = VH*p(1:NpH,hexK );
        if useLSC
            up(1:NcW,wedgK) = (VW*p(1:NpW,wedgK)).*sqJW;
        else
            up(1:NcW,wedgK) = VW*p(1:NpW,wedgK);
        end        
        up(1:NcP,pyrK ) = VP*p(1:NpP,pyrK );
        up(1:NcT,tetK ) = VT*p(1:NpT,tetK );
                
        clf
        drawMesh
        color_line3(x(ids),y(ids),z(ids),up(ids),'.');
        caxis([-1 1]); view(30,30);  xlabel('x'); ylabel('y'); colorbar
        title(sprintf('t = %f, max val = %f',time,max(abs(up(:)))))        
        axis equal
        drawnow        
        fprintf('Timestep %d/%d\n',tstep,round(FinalTime/dt))
    end
    
    % Increment time
    time = time+dt; tstep = tstep+1;                   

end

% compute error
uexc = zeros(NcMax,K);
uexc(1:NcH,hexK)  = uexf(xqH, yqH, zqH, time);
uexc(1:NcW,wedgK) = uexf(xqW, yqW, zqW, time);
uexc(1:NcP,pyrK)  = uexf(xqP, yqP, zqP, time);
uexc(1:NcT,tetK)  = uexf(xqT, yqT, zqT, time);

pc = zeros(NcMax,K);
pc(1:NcH,hexK ) = VH*p(1:NpH,hexK );
if useLSC
    pc(1:NcW,wedgK) = (VW*p(1:NpW,wedgK)).*sqJW;
else
    pc(1:NcW,wedgK) = VW*p(1:NpW,wedgK);
end
pc(1:NcP,pyrK ) = VP*p(1:NpP,pyrK );
pc(1:NcT,tetK ) = VT*p(1:NpT,tetK );

err2 = 0;
err2K = zeros(K,1);
for e = 1:K
    diff = uexc(:,e) - pc(:,e);    
    err2K(e) = wJ(:,e)'*(diff.^2);
    err2 = err2 + err2K(e);
end
err = sqrt(err2);
disp(sprintf('L2 err = %f\n',err))

% computes RHS for all element types 
function [rhsp rhsu rhsv rhsw] = RHSall(p,u,v,w,ps,us,vs,ws)

hybridgGlobals3D
hybridgGlobalFlags

% hex
if useSkewRHS
    [rhspK rhsuK rhsvK rhswK] = SkewRHS(p,u,v,w,ps,us,vs,ws,VH,VHr,VHs,VHt,VHf,hexK );
else
    [rhspK rhsuK rhsvK rhswK] = RHS(p,u,v,w,ps,us,vs,ws,VH,VHr,VHs,VHt,VHf,hexK );
end
rhsp(1:NpH,hexK) = rhspK; rhsu(1:NpH,hexK) = rhsuK;
rhsv(1:NpH,hexK) = rhsvK; rhsw(1:NpH,hexK) = rhswK;

% wedges
if useLSC
    [rhspK rhsuK rhsvK rhswK] = LSC_RHS(p,u,v,w,ps,us,vs,ws,...
        VW,VWr,VWs,VWt,VWf,wedgK, wW, sqJWf, JWr, JWs, JWt);
else
    [rhspK rhsuK rhsvK rhswK] = Wedge_RHS(p,u,v,w,ps,us,vs,ws,VW,VWr,VWs,VWt,VWf,wedgK);
end
rhsp(1:NpW,wedgK) = rhspK; rhsu(1:NpW,wedgK) = rhsuK;
rhsv(1:NpW,wedgK) = rhsvK; rhsw(1:NpW,wedgK) = rhswK;

% pyrs - always use skew form, same cost but stable if useSEM=1
% if useSkewRHS
    [rhspK rhsuK rhsvK rhswK] = SkewRHS(p,u,v,w,ps,us,vs,ws,VP,VPr,VPs,VPt,VPf,pyrK );
% else
%     [rhspK rhsuK rhsvK rhswK] = RHS(p,u,v,w,ps,us,vs,ws,VP,VPr,VPs,VPt,VPf,pyrK );
% end
rhsp(1:NpP,pyrK ) = rhspK; rhsu(1:NpP,pyrK ) = rhsuK;
rhsv(1:NpP,pyrK ) = rhsvK; rhsw(1:NpP,pyrK ) = rhswK;

% tet - switch w/nodal DG
if useNodalTets
    [rhspK rhsuK rhsvK rhswK] = NodalRHS(p,u,v,w,ps,us,vs,ws,VTr,VTs,VTt,VTf);
else
    [rhspK rhsuK rhsvK rhswK] = RHS(p,u,v,w,ps,us,vs,ws,VT,VTr,VTs,VTt,VTf,tetK );
end
rhsp(1:NpT,tetK ) = rhspK; rhsu(1:NpT,tetK ) = rhsuK;
rhsv(1:NpT,tetK ) = rhsvK; rhsw(1:NpT,tetK ) = rhswK;

% nodal RHS for tets: "collocation" DG
function [rhsp rhsu rhsv rhsw] = NodalRHS(p,u,v,w, ...
    p_surface,u_surface,v_surface,w_surface,...
    Vr,Vs,Vt,Vf)

hybridgGlobals3D;

if nnz(tetK)==0
    rhsp = []; rhsu = []; rhsv = []; rhsw = [];
    return
end

% element-type specifics
Np = NpT; Nfc = size(Vf,1);

% face cubature
wsJK = wsJ(1:Nfc,tetK); 
nxK = nx(1:Nfc,tetK); nyK = ny(1:Nfc,tetK); nzK = nz(1:Nfc,tetK);

% nodal dofs
pK = p(1:Np,tetK); uK = u(1:Np,tetK); vK = v(1:Np,tetK); wK = w(1:Np,tetK);

[p_jump nu_jump] = SurfaceFlux(p_surface,u_surface,v_surface,w_surface,Vf,tetK);
flux_p = .5*(p_jump - nu_jump);
flux_u = .5*(nu_jump - p_jump);

% volume derivative operators
dudx = rxT.*(Vr*uK) + sxT.*(Vs*uK) + txT.*(Vt*uK);
dvdy = ryT.*(Vr*vK) + syT.*(Vs*vK) + tyT.*(Vt*vK);
dwdz = rzT.*(Vr*wK) + szT.*(Vs*wK) + tzT.*(Vt*wK);
divU = dudx + dvdy + dwdz;
rhsp = -divU + invMThat*((Vf'*(wsJK.*flux_p))./JTet);

pr = Vr*pK; ps = Vs*pK; pt = Vt*pK;
dpdx = rxT.*pr + sxT.*ps + txT.*pt;
dpdy = ryT.*pr + syT.*ps + tyT.*pt;
dpdz = rzT.*pr + szT.*ps + tzT.*pt;

flux_uc = wsJK.*flux_u;
rhsu = -dpdx + invMThat*((Vf'*(nxK.*flux_uc))./JTet);
rhsv = -dpdy + invMThat*((Vf'*(nyK.*flux_uc))./JTet);
rhsw = -dpdz + invMThat*((Vf'*(nzK.*flux_uc))./JTet);

return

% strong form RHS with quadrature - OK for hex + pyramids + tets. 
% can do better w/nodal tets
function [rhsp rhsu rhsv rhsw] = RHS(p,u,v,w, ...
    p_surface,u_surface,v_surface,w_surface,...
    V,Vr,Vs,Vt,Vf,typeK)

hybridgGlobals3D;

if nnz(typeK)==0
    rhsp = []; rhsu = []; rhsv = []; rhsw = [];
    return
end

% element-type specifics
Nc = size(V,1); Np = size(V,2); Nfc = size(Vf,1);

invMK = invM(1:Np,typeK); wJK = wJ(1:Nc,typeK);  wsJK = wsJ(1:Nfc,typeK);

rxK = rx(1:Nc,typeK); sxK = sx(1:Nc,typeK); txK = tx(1:Nc,typeK);
ryK = ry(1:Nc,typeK); syK = sy(1:Nc,typeK); tyK = ty(1:Nc,typeK);
rzK = rz(1:Nc,typeK); szK = sz(1:Nc,typeK); tzK = tz(1:Nc,typeK);

nxK = nx(1:Nfc,typeK); nyK = ny(1:Nfc,typeK); nzK = nz(1:Nfc,typeK);
pK = p(1:Np,typeK); uK = u(1:Np,typeK); vK = v(1:Np,typeK); wK = w(1:Np,typeK);

[p_jump nu_jump] = SurfaceFlux(p_surface,u_surface,v_surface,w_surface,Vf,typeK);
flux_p = .5*(p_jump - nu_jump);
flux_u = .5*(nu_jump - p_jump);

% volume derivative operators
dudx = rxK.*(Vr*uK) + sxK.*(Vs*uK) + txK.*(Vt*uK);
dvdy = ryK.*(Vr*vK) + syK.*(Vs*vK) + tyK.*(Vt*vK);
dwdz = rzK.*(Vr*wK) + szK.*(Vs*wK) + tzK.*(Vt*wK);
divU = dudx + dvdy + dwdz;
rhsp = invMK.*(-V'*(wJK.*divU) + Vf'*(wsJK.*flux_p));

pr = Vr*pK; ps = Vs*pK; pt = Vt*pK;
dpdx = rxK.*pr + sxK.*ps + txK.*pt;
dpdy = ryK.*pr + syK.*ps + tyK.*pt;
dpdz = rzK.*pr + szK.*ps + tzK.*pt;

flux_uc = wsJK.*flux_u;
rhsu = invMK.*(-V'*(wJK.*dpdx) + Vf'*(nxK.*flux_uc));
rhsv = invMK.*(-V'*(wJK.*dpdy) + Vf'*(nyK.*flux_uc));
rhsw = invMK.*(-V'*(wJK.*dpdz) + Vf'*(nzK.*flux_uc));

return

function [rhsp rhsu rhsv rhsw] = SkewRHS(p,u,v,w, ...
    p_surface,u_surface,v_surface,w_surface,...
    V,Vr,Vs,Vt,Vf,typeK)

hybridgGlobals3D;

if nnz(typeK)==0
    rhsp = []; rhsu = []; rhsv = []; rhsw = [];
    return
end

% element-type specifics
Nc = size(V,1); Np = size(V,2); Nfc = size(Vf,1);

invMK = invM(1:Np,typeK); wJK = wJ(1:Nc,typeK);  wsJK = wsJ(1:Nfc,typeK);

rxK = rx(1:Nc,typeK); sxK = sx(1:Nc,typeK); txK = tx(1:Nc,typeK);
ryK = ry(1:Nc,typeK); syK = sy(1:Nc,typeK); tyK = ty(1:Nc,typeK);
rzK = rz(1:Nc,typeK); szK = sz(1:Nc,typeK); tzK = tz(1:Nc,typeK);

nxK = nx(1:Nfc,typeK); nyK = ny(1:Nfc,typeK); nzK = nz(1:Nfc,typeK);
pK = p(1:Np,typeK); uK = u(1:Np,typeK); vK = v(1:Np,typeK); wK = w(1:Np,typeK);

[p_jump nu_jump nu_avg] = SurfaceFlux(p_surface,u_surface,v_surface,w_surface,Vf,typeK);
% flux_p = .5*(p_jump - nu_jump);
flux_u = .5*(nu_jump - p_jump);
flux_p = .5*p_jump - nu_avg;

% integrated by parts pressure equation
uc = wJK.*(V*uK); vc = wJK.*(V*vK); wc = wJK.*(V*wK);
Vxu = Vr'*(rxK.*uc) + Vs'*(sxK.*uc) + Vt'*(txK.*uc);
Vyv = Vr'*(ryK.*vc) + Vs'*(syK.*vc) + Vt'*(tyK.*vc);
Vzw = Vr'*(rzK.*wc) + Vs'*(szK.*wc) + Vt'*(tzK.*wc);
gradvU = Vxu + Vyv + Vzw;
rhsp = invMK.*(gradvU + Vf'*(wsJK.*flux_p));

% % volume derivative operators
% dudx = rxK.*(Vr*uK) + sxK.*(Vs*uK) + txK.*(Vt*uK);
% dvdy = ryK.*(Vr*vK) + syK.*(Vs*vK) + tyK.*(Vt*vK);
% dwdz = rzK.*(Vr*wK) + szK.*(Vs*wK) + tzK.*(Vt*wK);
% divU = dudx + dvdy + dwdz;
% rhsp = invMK.*(-V'*(wJK.*divU) + Vf'*(wsJK.*flux_p));

pr = Vr*pK; ps = Vs*pK; pt = Vt*pK;
dpdx = rxK.*pr + sxK.*ps + txK.*pt;
dpdy = ryK.*pr + syK.*ps + tyK.*pt;
dpdz = rzK.*pr + szK.*ps + tzK.*pt;

flux_uc = wsJK.*flux_u;
rhsu = invMK.*(-V'*(wJK.*dpdx) + Vf'*(nxK.*flux_uc));
rhsv = invMK.*(-V'*(wJK.*dpdy) + Vf'*(nyK.*flux_uc));
rhsw = invMK.*(-V'*(wJK.*dpdz) + Vf'*(nzK.*flux_uc));

return

function [rhsp rhsu rhsv rhsw] = LSC_RHS(p,u,v,w, ...
    p_surface,u_surface,v_surface,w_surface,...
    V,Vr,Vs,Vt,Vf,typeK,...
    wqK, sqJKf, JKr, JKs, JKt)

hybridgGlobals3D;

if nnz(typeK)==0
    rhsp = []; rhsu = []; rhsv = []; rhsw = [];
    return
end

% element-type specifics
Nc = size(V,1); Np = size(V,2); Nfc = size(Vf,1);

invMK = invM(1:Np,typeK); % could also just take typeK(1) - same mass matrix
% wJK = wJ(1:Nc,typeK);  
wsJK = wsJ(1:Nfc,typeK);

rxK = rx(1:Nc,typeK); sxK = sx(1:Nc,typeK); txK = tx(1:Nc,typeK);
ryK = ry(1:Nc,typeK); syK = sy(1:Nc,typeK); tyK = ty(1:Nc,typeK);
rzK = rz(1:Nc,typeK); szK = sz(1:Nc,typeK); tzK = tz(1:Nc,typeK);

nxK = nx(1:Nfc,typeK); nyK = ny(1:Nfc,typeK); nzK = nz(1:Nfc,typeK);
pK = p(1:Np,typeK); uK = u(1:Np,typeK); vK = v(1:Np,typeK); wK = w(1:Np,typeK);

[p_jump nu_jump nu_avg] = SurfaceFlux(p_surface,u_surface,v_surface,w_surface,Vf,typeK);
flux_p = .5*p_jump - nu_avg;
flux_u = .5*(nu_jump - p_jump);

% % integrated by parts pressure equation
% LSC approach - chain rule derivatives
uc = diag(wqK)*(V*uK);
Vrxu = Vr'*(rxK.*uc) - V'*(JKr.*rxK.*uc);
Vsxu = Vs'*(sxK.*uc) - V'*(JKs.*sxK.*uc);
Vtxu = Vt'*(txK.*uc) - V'*(JKt.*txK.*uc);

vc = diag(wqK)*(V*vK);
Vryv = Vr'*(ryK.*vc) - V'*(JKr.*ryK.*vc);
Vsyv = Vs'*(syK.*vc) - V'*(JKs.*syK.*vc);
Vtyv = Vt'*(tyK.*vc) - V'*(JKt.*tyK.*vc);

wc = diag(wqK)*(V*wK);
Vrzw = Vr'*(rzK.*wc) - V'*(JKr.*rzK.*wc);
Vszw = Vs'*(szK.*wc) - V'*(JKs.*szK.*wc);
Vtzw = Vt'*(tzK.*wc) - V'*(JKt.*tzK.*wc);

Vxu = Vrxu + Vsxu + Vtxu;
Vyv = Vryv + Vsyv + Vtyv;
Vzw = Vrzw + Vszw + Vtzw;
gradvU = Vxu + Vyv + Vzw;
rhsp = invMK.*(gradvU + Vf'*(sqJKf.*wsJK.*flux_p)); % invM = ones for LSC-DG 

% tmp1 = (sqwJK.*dpdx);
% 
% pr = (Vr*pK - (V*pK).*JKr); 
% ps = (Vs*pK - (V*pK).*JKs); 
% pt = (Vt*pK - (V*pK).*JKt);
% dpdx1 = rxK.*pr + sxK.*ps + txK.*pt;
% norm(tmp1-diag(wW)*dpdx1,'fro')

% volume derivative operators
% pr = Vr*pK; ps = Vs*pK; pt = Vt*pK;

pr = (Vr*pK - (V*pK).*JKr);
ps = (Vs*pK - (V*pK).*JKs);
pt = (Vt*pK - (V*pK).*JKt);
dpdx = rxK.*pr + sxK.*ps + txK.*pt;
dpdy = ryK.*pr + syK.*ps + tyK.*pt;
dpdz = rzK.*pr + szK.*ps + tzK.*pt;

flux_uc = sqJKf.*wsJK.*flux_u;

rhsu = invMK.*(-V'*(diag(wqK)*dpdx) + Vf'*(nxK.*flux_uc)); % invM = ones for LSC-DG
rhsv = invMK.*(-V'*(diag(wqK)*dpdy) + Vf'*(nyK.*flux_uc));
rhsw = invMK.*(-V'*(diag(wqK)*dpdz) + Vf'*(nzK.*flux_uc));


return

% strong form RHS with full mass inversion for wedges
function [rhsp rhsu rhsv rhsw] = Wedge_RHS(p,u,v,w, ...
    p_surface,u_surface,v_surface,w_surface,...
    V,Vr,Vs,Vt,Vf,typeK)

hybridgGlobals3D;

if nnz(typeK)==0
    rhsp = []; rhsu = []; rhsv = []; rhsw = [];
    return
end

% element-type specifics
Nc = size(V,1); Np = size(V,2); Nfc = size(Vf,1);

% invMK = invM(1:Np,typeK); 
wJK = wJ(1:Nc,typeK);  wsJK = wsJ(1:Nfc,typeK);

rxK = rx(1:Nc,typeK); sxK = sx(1:Nc,typeK); txK = tx(1:Nc,typeK);
ryK = ry(1:Nc,typeK); syK = sy(1:Nc,typeK); tyK = ty(1:Nc,typeK);
rzK = rz(1:Nc,typeK); szK = sz(1:Nc,typeK); tzK = tz(1:Nc,typeK);

nxK = nx(1:Nfc,typeK); nyK = ny(1:Nfc,typeK); nzK = nz(1:Nfc,typeK);
pK = p(1:Np,typeK); uK = u(1:Np,typeK); vK = v(1:Np,typeK); wK = w(1:Np,typeK);

[p_jump nu_jump] = SurfaceFlux(p_surface,u_surface,v_surface,w_surface,Vf,typeK);
flux_p = .5*(p_jump - nu_jump);
flux_u = .5*(nu_jump - p_jump);

% volume derivative operators
dudx = rxK.*(Vr*uK) + sxK.*(Vs*uK) + txK.*(Vt*uK);
dvdy = ryK.*(Vr*vK) + syK.*(Vs*vK) + tyK.*(Vt*vK);
dwdz = rzK.*(Vr*wK) + szK.*(Vs*wK) + tzK.*(Vt*wK);
divU = dudx + dvdy + dwdz;
% rhsp = invMK.*(-V'*(wJK.*divU) + Vf'*(wsJK.*flux_p));

pr = Vr*pK; ps = Vs*pK; pt = Vt*pK;
dpdx = rxK.*pr + sxK.*ps + txK.*pt;
dpdy = ryK.*pr + syK.*ps + tyK.*pt;
dpdz = rzK.*pr + szK.*ps + tzK.*pt;

flux_uc = wsJK.*flux_u;
% rhsu = invMK.*(-V'*(wJK.*dpdx) + Vf'*(nxK.*flux_uc));
% rhsv = invMK.*(-V'*(wJK.*dpdy) + Vf'*(nyK.*flux_uc));
% rhsw = invMK.*(-V'*(wJK.*dpdz) + Vf'*(nzK.*flux_uc));

% full mass matrix inversion
rhsp = (-V'*(wJK.*divU) + Vf'*(wsJK.*flux_p));
rhsu = (-V'*(wJK.*dpdx) + Vf'*(nxK.*flux_uc));
rhsv = (-V'*(wJK.*dpdy) + Vf'*(nyK.*flux_uc));
rhsw = (-V'*(wJK.*dpdz) + Vf'*(nzK.*flux_uc));
for e = 1:length(typeK)
    rhse = MW{e}\[rhsp(:,e) rhsu(:,e) rhsv(:,e) rhsw(:,e)];
    rhsp(:,e) = rhse(:,1);    rhsu(:,e) = rhse(:,2);
    rhsv(:,e) = rhse(:,3);    rhsw(:,e) = rhse(:,4);
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
mapBK = mapMK==mapPK;
p_jump(mapBK) = -2*p_surface(mapMK(mapBK));


function [ps us vs ws] = SurfaceInterp(p,u,v,w)

hybridgGlobals3D
hybridgGlobalFlags

ps = zeros(NfcMax,K);
us = zeros(NfcMax,K);
vs = zeros(NfcMax,K);
ws = zeros(NfcMax,K);

ps(1:NfcH,hexK ) = VHf*p(1:NpH,hexK );
us(1:NfcH,hexK ) = VHf*u(1:NpH,hexK );
vs(1:NfcH,hexK ) = VHf*v(1:NpH,hexK );
ws(1:NfcH,hexK ) = VHf*w(1:NpH,hexK );

if useLSC
    % wedges LSC interp
    ps(1:NfcW,wedgK) = (VWf*p(1:NpW,wedgK)).*sqJWf;
    us(1:NfcW,wedgK) = (VWf*u(1:NpW,wedgK)).*sqJWf;
    vs(1:NfcW,wedgK) = (VWf*v(1:NpW,wedgK)).*sqJWf;
    ws(1:NfcW,wedgK) = (VWf*w(1:NpW,wedgK)).*sqJWf;
else
    ps(1:NfcW,wedgK) = VWf*p(1:NpW,wedgK);
    us(1:NfcW,wedgK) = VWf*u(1:NpW,wedgK);
    vs(1:NfcW,wedgK) = VWf*v(1:NpW,wedgK);
    ws(1:NfcW,wedgK) = VWf*w(1:NpW,wedgK);    
end

ps(1:NfcP,pyrK ) = VPf*p(1:NpP,pyrK );
us(1:NfcP,pyrK ) = VPf*u(1:NpP,pyrK );
vs(1:NfcP,pyrK ) = VPf*v(1:NpP,pyrK );
ws(1:NfcP,pyrK ) = VPf*w(1:NpP,pyrK );

ps(1:NfcT,tetK ) = VTf*p(1:NpT,tetK );
us(1:NfcT,tetK ) = VTf*u(1:NpT,tetK );
vs(1:NfcT,tetK ) = VTf*v(1:NpT,tetK );
ws(1:NfcT,tetK ) = VTf*w(1:NpT,tetK );

