function AdvecDriver

% Driver script for solving the 3D advection equations

hybridgGlobals3D

% hybrid_mesh
prism_mesh2
K = size(EToV,1);
N = 4;

hybridgStartUp

%%
% make periodic bcs - match node positions
inIds = find(abs(xf+1)<NODETOL); outIds = find(abs(xf-1)<NODETOL);
tmp = ones(1,nnz(inIds));
yMp = yf(inIds)*tmp;  zMp = zf(inIds)*tmp;
yPp = yf(outIds)*tmp; zPp = zf(outIds)*tmp;

D = (yMp-yPp').^2 + (zMp-zPp').^2;
[id1, id2] = find(abs(D)<NODETOL);
inIds = inIds(id1); outIds = outIds(id2);

mapP(inIds) = mapM(outIds); mapP(outIds) = mapM(inIds);

%% set initial conditions

ids = y > 0;

f = @(x,y,z) exp(-5.0*(x.^2 + y.^2 + z.^2));

b(1:NpH,hexK ) = VH'*(wJ(1:NcH,hexK ).*f(xH,yH,zH));
b(1:NpW,wedgK) = VW'*(wJ(1:NcW,wedgK).*f(xW,yW,zW));
b(1:NpP,pyrK ) = VP'*(wJ(1:NcP,pyrK ).*f(xP,yP,zP)); 
b(1:NpT,tetK ) = VT'*(wJ(1:NcT,tetK ).*f(xT,yT,zT)); 

u = invM.*b;

% get info of u at all cubature points
u_surface = zeros(NfcMax,K);
u_surface(1:NfcH,hexK ) = VHf*u(1:NpH,hexK );
u_surface(1:NfcW,wedgK) = VWf*u(1:NpW,wedgK);
u_surface(1:NfcP,pyrK ) = VPf*u(1:NpP,pyrK );
u_surface(1:NfcT,tetK ) = VTf*u(1:NpT,tetK );

% for plotting
up = zeros(NcMax,K);
up(1:NcH,hexK ) = VH*u(1:NpH,hexK ); 
up(1:NcW,wedgK) = VW*u(1:NpW,wedgK);
up(1:NcP,pyrK ) = VP*u(1:NpP,pyrK); 
up(1:NcT,tetK ) = VT*u(1:NpT,tetK);

%%

resu = zeros(NpMax,K);

% compute time step size
dt = .125/max(((N+1)*(N+3))*(2/3)*Fscale(:));

FinalTime = .5;

% outer time step loop
time = 0;
tstep = 0;
while time < FinalTime
    
    if(time+dt>FinalTime), 
        dt = FinalTime-time;  
    end
    
    % low storage RK    
    for INTRK = 1:5        
        rhsu(1:NpH,hexK ) = RHS(u,u_surface,VH,VHr,VHs,VHt,VHf,hexK );
        rhsu(1:NpW,wedgK) = RHS(u,u_surface,VW,VWr,VWs,VWt,VWf,wedgK);
        rhsu(1:NpP,pyrK ) = RHS(u,u_surface,VP,VPr,VPs,VPt,VPf,pyrK );
        rhsu(1:NpT,tetK ) = RHS(u,u_surface,VT,VTr,VTs,VTt,VTf,tetK );
                
        resu = rk4a(INTRK)*resu + dt*rhsu;
        u = u + rk4b(INTRK)*resu;
        
        u_surface = zeros(NfcMax,K);
        u_surface(1:NfcH,hexK ) = VHf*u(1:NpH,hexK );
        u_surface(1:NfcW,wedgK) = VWf*u(1:NpW,wedgK);
        u_surface(1:NfcP,pyrK ) = VPf*u(1:NpP,pyrK );
        u_surface(1:NfcT,tetK ) = VTf*u(1:NpT,tetK );
    end
    
    % Increment time
    time = time+dt; tstep = tstep+1;      
    
    if (mod(tstep,25)==0)        
        up(1:NcH,hexK ) = VH*u(1:NpH,hexK ); 
        up(1:NcW,wedgK) = VW*u(1:NpW,wedgK);
        up(1:NcP,pyrK ) = VP*u(1:NpP,pyrK );
        up(1:NcT,tetK ) = VT*u(1:NpT,tetK );
        
        clf
        color_line3(x(ids),y(ids),z(ids),up(ids),'.');        
        view(30,30);
        xlabel('x'); ylabel('y')
        title(sprintf('time t = %f, max value of u = %i',time,max(abs(up(:)))))
        colorbar
        axis equal
        drawnow
    end   
    
    
end
uexc = exp(-5.0*((x-time).^2 + y.^2 + z.^2));
uex(1:NpH,hexK ) = (VH'*(wJ(1:NcH,hexK).*uexc(1:NcH,hexK)));
uex(1:NpW,wedgK) = (VW'*(wJ(1:NcW,wedgK).*uexc(1:NcW,wedgK)));
uex(1:NpP,pyrK ) = (VP'*(wJ(1:NcP,pyrK).*uexc(1:NcW,pyrK)));
uex(1:NpT,tetK ) = (VT'*(wJ(1:NcT,tetK).*uexc(1:NcT,tetK)));

diff = u-invM.*uex;
M = (1./invM); M(isnan(M) | isinf(M))=0;
err = M.*(diff.^2);
sqrt(sum(err(:)))
keyboard

function rhsu = RHS(u,u_surface,V,Vr,Vs,Vt,Vf,typeK)

hybridgGlobals3D;

Nc = size(V,1); Np = size(V,2); Nfc = size(Vf,1);

% flux gather
alpha = 1;
u_jump = u_surface(mapM(1:Nfc,typeK)) - u_surface(mapP(1:Nfc,typeK));
flux = .5*(nx(1:Nfc,typeK)-alpha*abs(nx(1:Nfc,typeK))).*u_jump; 

% get u on element type
uK = u(1:Np,typeK);

% volume derivative operators
dudx = rx(1:Nc,typeK).*(Vr*uK) + sx(1:Nc,typeK).*(Vs*uK) + tx(1:Nc,typeK).*(Vt*uK);

% integration steps
rhsu = invM(1:Np,typeK).*(-V'*(wJ(1:Nc,typeK).*dudx) + Vf'*(wsJ(1:Nfc,typeK).*flux));

return

