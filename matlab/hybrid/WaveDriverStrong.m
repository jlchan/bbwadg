function WaveDriverStrong

% Driver script for solving the 3D advection equations

hybridgGlobals3D
hybridgGlobalFlags
useSEM = 1;

hybrid_mesh
% pyr_mesh
% prism_mesh
% hex_mesh
K = size(EToV,1);
N = 5;

hybridgStartUp

%% set initial conditions

ids = y > -1;

f = @(x,y,z) cos(pi/2*x).*cos(pi/2.*y).*cos(pi/2.*z);
uexf = @(x,y,z,time) cos(pi/2*x).*cos(pi/2.*y).*cos(pi/2.*z)*cos(sqrt(3)*pi/2*time);

b(1:NpH,hexK ) = VH'*(wJ(1:NcH,hexK ).*f(xH,yH,zH));
b(1:NpW,wedgK) = VW'*(wJ(1:NcW,wedgK).*f(xW,yW,zW));
b(1:NpP,pyrK ) = VP'*(wJ(1:NcP,pyrK ).*f(xP,yP,zP)); 
b(1:NpT,tetK ) = VT'*(wJ(1:NcT,tetK ).*f(xT,yT,zT)); 

p = invM.*b;
u = zeros(NpMax,K); v = zeros(NpMax,K); w = zeros(NpMax,K); % velocity

% get info of u at all cubature points
[ps us vs ws] = SurfaceInterp(p,u,v,w);

%%

resp = zeros(NpMax,K); resu = zeros(NpMax,K); resv = zeros(NpMax,K); resw = zeros(NpMax,K);
rhsp = zeros(NpMax,K); rhsu = zeros(NpMax,K); rhsv = zeros(NpMax,K); rhsw = zeros(NpMax,K);

% compute time step size
dt = .5/max(((N+1)*(N+3))*(2/3)*Fscale(:));

FinalTime = .25;

% outer time step loop
figure
time = 0;
tstep = 0;
while time < FinalTime
    
    if(time+dt>FinalTime), 
        dt = FinalTime-time;  
    end
    
    % low storage RK    
    for INTRK = 1:5        
        [rhspK rhsuK rhsvK rhswK] = RHS(p,u,v,w,ps,us,vs,ws,VH,VHr,VHs,VHt,VHf,hexK );
        rhsp(1:NpH,hexK) = rhspK; rhsu(1:NpH,hexK) = rhsuK;
        rhsv(1:NpH,hexK) = rhsvK; rhsw(1:NpH,hexK) = rhswK;
        
        [rhspK rhsuK rhsvK rhswK] = RHS(p,u,v,w,ps,us,vs,ws,VW,VWr,VWs,VWt,VWf,wedgK);
        rhsp(1:NpW,wedgK) = rhspK; rhsu(1:NpW,wedgK) = rhsuK;
        rhsv(1:NpW,wedgK) = rhsvK; rhsw(1:NpW,wedgK) = rhswK;
        
        [rhspK rhsuK rhsvK rhswK] = RHS(p,u,v,w,ps,us,vs,ws,VP,VPr,VPs,VPt,VPf,pyrK );
        rhsp(1:NpP,pyrK ) = rhspK; rhsu(1:NpP,pyrK ) = rhsuK;
        rhsv(1:NpP,pyrK ) = rhsvK; rhsw(1:NpP,pyrK ) = rhswK;
                
        [rhspK rhsuK rhsvK rhswK] = RHS(p,u,v,w,ps,us,vs,ws,VT,VTr,VTs,VTt,VTf,tetK );
        rhsp(1:NpT,tetK ) = rhspK; rhsu(1:NpT,tetK ) = rhsuK;
        rhsv(1:NpT,tetK ) = rhsvK; rhsw(1:NpT,tetK ) = rhswK;        
                                   
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
    
    % Increment time
    time = time+dt; tstep = tstep+1;      
    
    if (mod(tstep,25)==0)        
        up(1:NcH,hexK ) = VH*p(1:NpH,hexK ); 
        up(1:NcW,wedgK) = VW*p(1:NpW,wedgK);
        up(1:NcP,pyrK ) = VP*p(1:NpP,pyrK );
        up(1:NcT,tetK ) = VT*p(1:NpT,tetK );
        
        clf
        drawMesh
        color_line3(x(ids),y(ids),z(ids),up(ids),'.');        
        caxis([-1 1])
        view(30,30);
        xlabel('x'); ylabel('y')
        title(sprintf('time t = %f, max value of u = %i',time,max(abs(up(:)))))
        colorbar
        axis equal
        drawnow
    end           
end

% compute error
uexc = uexf(x,y,z,time);
pc = zeros(NcMax,K);
pc(1:NcH,hexK ) = VH*p(1:NpH,hexK );
pc(1:NcW,wedgK) = VW*p(1:NpW,wedgK);
pc(1:NcP,pyrK ) = VP*p(1:NpP,pyrK );
pc(1:NcT,tetK ) = VT*p(1:NpT,tetK );

err2 = 0;
for e = 1:K
    diff = uexc(:,e) - pc(:,e);
    err2 = err2 + wJ(:,e)'*(diff.^2);
end
err = sqrt(err2);
disp(sprintf('L2 err = %f\n',err))
keyboard


function [rhsp rhsu rhsv rhsw] = RHS(p,u,v,w, ... 
    p_surface,u_surface,v_surface,w_surface,...
    V,Vr,Vs,Vt,Vf,typeK)

hybridgGlobals3D;

Nc = size(V,1); Np = size(V,2); Nfc = size(Vf,1);

% element-type specifics
invMK = invM(1:Np,typeK);
wJK = wJ(1:Nc,typeK);  

rxK = rx(1:Nc,typeK); sxK = sx(1:Nc,typeK); txK = tx(1:Nc,typeK);
ryK = ry(1:Nc,typeK); syK = sy(1:Nc,typeK); tyK = ty(1:Nc,typeK);
rzK = rz(1:Nc,typeK); szK = sz(1:Nc,typeK); tzK = tz(1:Nc,typeK);

pK = p(1:Np,typeK); 
uK = u(1:Np,typeK); vK = v(1:Np,typeK); wK = w(1:Np,typeK);

wsJK = wsJ(1:Nfc,typeK); 
nxK = nx(1:Nfc,typeK); nyK = ny(1:Nfc,typeK); nzK = nz(1:Nfc,typeK);

mapMK = mapM(1:Nfc,typeK); mapPK = mapP(1:Nfc,typeK);

% flux gather
p_jump = p_surface(mapPK) - p_surface(mapMK);
u_jump = u_surface(mapPK) - u_surface(mapMK);
v_jump = v_surface(mapPK) - v_surface(mapMK);
w_jump = w_surface(mapPK) - w_surface(mapMK);
nu_jump = nxK.*u_jump + nyK.*v_jump + nzK.*w_jump;

% apply BCs, todo: try p^+ = - p ^ -, which is energy stable
mapBK = mapMK==mapPK;
p_jump(mapBK) = -2*p_surface(mapMK(mapBK));

flux_p = .5*(p_jump - nu_jump);
flux_u = .5*(nu_jump - p_jump);

% volume derivative operators
dpdx = rxK.*(Vr*pK) + sxK.*(Vs*pK) + txK.*(Vt*pK);
dpdy = ryK.*(Vr*pK) + syK.*(Vs*pK) + tyK.*(Vt*pK);
dpdz = rzK.*(Vr*pK) + szK.*(Vs*pK) + tzK.*(Vt*pK);

dudx = rxK.*(Vr*uK) + sxK.*(Vs*uK) + txK.*(Vt*uK);
dvdy = ryK.*(Vr*vK) + syK.*(Vs*vK) + tyK.*(Vt*vK);
dwdz = rzK.*(Vr*wK) + szK.*(Vs*wK) + tzK.*(Vt*wK);
divU = dudx + dvdy + dwdz;
rhsp = invMK.*(-V'*(wJK.*divU) + Vf'*(wsJK.*flux_p));

flux_uc = wsJK.*flux_u;
rhsu = invMK.*(-V'*(wJK.*dpdx) + Vf'*(nxK.*flux_uc));
rhsv = invMK.*(-V'*(wJK.*dpdy) + Vf'*(nyK.*flux_uc));
rhsw = invMK.*(-V'*(wJK.*dpdz) + Vf'*(nzK.*flux_uc));

return

function [ps us vs ws] = SurfaceInterp(p,u,v,w)

hybridgGlobals3D

ps = zeros(NfcMax,K);
ps(1:NfcH,hexK ) = VHf*p(1:NpH,hexK );
ps(1:NfcW,wedgK) = VWf*p(1:NpW,wedgK);
ps(1:NfcP,pyrK ) = VPf*p(1:NpP,pyrK );
ps(1:NfcT,tetK ) = VTf*p(1:NpT,tetK );

us = zeros(NfcMax,K);
us(1:NfcH,hexK ) = VHf*u(1:NpH,hexK );
us(1:NfcW,wedgK) = VWf*u(1:NpW,wedgK);
us(1:NfcP,pyrK ) = VPf*u(1:NpP,pyrK );
us(1:NfcT,tetK ) = VTf*u(1:NpT,tetK );

vs = zeros(NfcMax,K);
vs(1:NfcH,hexK ) = VHf*v(1:NpH,hexK );
vs(1:NfcW,wedgK) = VWf*v(1:NpW,wedgK);
vs(1:NfcP,pyrK ) = VPf*v(1:NpP,pyrK );
vs(1:NfcT,tetK ) = VTf*v(1:NpT,tetK );

ws = zeros(NfcMax,K);
ws(1:NfcH,hexK ) = VHf*w(1:NpH,hexK );
ws(1:NfcW,wedgK) = VWf*w(1:NpW,wedgK);
ws(1:NfcP,pyrK ) = VPf*w(1:NpP,pyrK );
ws(1:NfcT,tetK ) = VTf*w(1:NpT,tetK );