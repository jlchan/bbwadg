
% clear -global *
% clear

Globals2D;

N = 2;
nref = 0;
useJprojection = 1;
FinalTime = 1;

Nq = 2*N+1;
Nfq = N;

K1D = 2;
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);

% This builds the nodal DG stuff
StartUp2D;


    a = .125;
    % warp mesh
    Lx = 1; Ly = 1;
    x0 = 0; y0 = 0;
    x = x + Lx*a*cos(1/2*pi*(x-x0)/Lx).*cos(3/2*pi*(y-y0)/Ly);
    y = y + Ly*a*sin(3/2*pi*(x-x0)/Lx).*cos(1/2*pi*(y-y0)/Ly);


cinfo = BuildCurvedOPS2D_opt(Nq,Nfq);
%     return
% PlotMesh2D;return

% jesse custom curved driver
global Vq wq Vrq Vsq
global Vfq wfq VfqFace
global rxq sxq ryq syq Jq nxq nyq sJq
global Pq Pfq Prq Psq


[rp sp] = EquiNodes2D(25); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

[rq sq wq] = Cubature2D(Nq);
Vq = Vandermonde2D(N,rq,sq)/V;
[Vrq Vsq] = GradVandermonde2D(N,rq,sq);
Vrq = Vrq/V; Vsq = Vsq/V;

% interp nodal to quadrature
[rq1D, wq1D] = JacobiGQ(0, 0, Nfq);
rfq = [rq1D, -rq1D, -ones(size(rq1D))];
sfq = [-ones(size(rq1D)), rq1D, -rq1D];
rfq = rfq(:); sfq = sfq(:);
wfq = repmat(wq1D,Nfaces,1);

Vfq = Vandermonde2D(N,rfq,sfq)/V;
[Vrfq Vsfq] = GradVandermonde2D(N,rfq,sfq);
Vrfq = Vrfq/V;
Vsfq = Vsfq/V;

V1D = Vandermonde1D(N,JacobiGL(0,0,N));
VfqFace = Vandermonde1D(N,rq1D)/V1D;
VfqFace = blkdiag(VfqFace,VfqFace,VfqFace); % repeat for 3 faces

% project down
Pq = (V*V')*(Vq'*diag(wq));
Prq = (V*V')*(Vrq'*diag(wq));
Psq = (V*V')*(Vsq'*diag(wq));
Pfq = (V*V')*(Vfq'*diag(wfq));

Nc = length(wq);
[rxq,sxq,ryq,syq,Jq] = GeometricFactors2D(x,y,Vq*Dr,Vq*Ds);
    
Nfc = length(cinfo(1).gnx(:));
nxq = zeros(Nfc,K); nyq = zeros(Nfc,K);
sJq = zeros(Nfc,K);
for e = 1:K  
    nxq(:,e) = cinfo(e).gnx(:);
    nyq(:,e) = cinfo(e).gny(:);
    sJq(:,e) = cinfo(e).gsJ(:);
end

aa = [nxq(:) nyq(:) sJq(:)],return

global tau;
tau = 0;


%%
if 1
    
    e = zeros(3*Np*K,1);
    A = zeros(3*Np*K);
    for i = 1:3*Np*K
        e(i) = 1;
        ids = 1:Np*K;
        p = reshape(e(ids),Np,K);
        u = reshape(e(ids + Np*K),Np,K);
        v = reshape(e(ids + 2*Np*K),Np,K);
%         [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p);
        [rhsu, rhsv, rhsp] = acousticsRHS2D_skew(u,v,p);
%         [rhsu, rhsv, rhsp] = acousticsRHS2D_JC(u,v,p);
        
        A(:,i) = [rhsp(:);rhsu(:);rhsv(:)];
        e(i) = 0;
        if mod(i,100)==0
            disp(sprintf('on column %d out of %d\n',i,3*Np*K))
        end
    end
    
    %     [rq sq wq] = Cubature2D(Nq);
    %     invMhat = V*V';
    %     for e = 1:K
    %         invMK = invMhat*Vq'*diag(wq./Jq(:,e))*Vq*invMhat;
    %         Mblk{e} = inv(invMK); %kron(eye(3),MK);
    %     end
    %     Mblk = blkdiag(Mblk{:});
    %     Mblk = blkdiag(Mblk,Mblk,Mblk); % for 3 fields
    %     MA = Mblk*A;
    
    lam = eig(A);
    hold on;
    if useJprojection==1
        plot(lam,'o')
    else
        plot(lam,'*')
    end
    title(sprintf('Largest real part = %e',max(real(lam))))
    return
    
end
%%

% Set initial conditions
% First 6 modes of eigenmodes with 6 azimuthal periods
alpha = [2.40482555769577, 5.52007811028631, 8.65372791291101, 11.7915, 14.9309]; % https://math.dartmouth.edu/archive/m23f09/public_html/drum.pdf

% choose radial mode
alpha0 = alpha(2); rad = sqrt(x.^2+y.^2);

p = besselj(0, alpha0*rad);

% d = 25; x0 = -1/3; y0 = 1/3;
% p = exp(-d*sqrt((x-x0).^2 + (y-y0).^2)) + exp(-10*sqrt((x+x0).^2 + (y+y0).^2));
p = exp(-25*(x.^2 + y.^2));
u = zeros(Np, K); v = zeros(Np, K);


% setup
resu = zeros(Np,K); resv = zeros(Np,K); resp = zeros(Np,K);
rLGL = JacobiGQ(0,0,N); rmin = abs(rLGL(1)-rLGL(2));
% dtscale = dtscale2D;
CN = (N+1)^2/2;
h = 1/min(min(sJ./J(Fmask(:),:)));
dtscale = h/CN;

dt = .5*min(dtscale)*rmin;

% outer time step loop
time = 0; tstep = 0;
while (time<FinalTime)
    if(time+dt>FinalTime), dt = FinalTime-time; end
    
    for INTRK = 1:5
        % compute right hand side of TM-mode acoustics's equations        
%         [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p);
        %[rhsu, rhsv, rhsp] = acousticsRHS2D_JC(u,v,p);                
        [rhsu, rhsv, rhsp] = acousticsRHS2D_skew(u,v,p);                
                
        % initiate and increment Runge-Kutta residuals
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;        
        
        % update fields
        u = u+rk4b(INTRK)*resu;
        v = v+rk4b(INTRK)*resv;
        p = p+rk4b(INTRK)*resp;
    end;
        
    
    if 1 && mod(tstep,10)==0
        clf
        pp = Vp*p;
        color_line3(xp,yp,pp,pp,'.');
        view(2)
        %         axis([-1 1 -1 1 -.5 .5])
        title(sprintf('time = %f',time));
        drawnow
    end
    
    % Increment time
    time = time+dt; tstep = tstep+1;
    if mod(tstep,10) ==0
        disp(sprintf('on timestep %d out of %d\n',tstep,round(FinalTime/dt)))
    end
end

xq = Vq*x;
yq = Vq*y;
radq = sqrt(xq.^2 + yq.^2);
pex = besselj(0,alpha0*radq)*cos(alpha0*time);

L2err = 0.0;
for e = 1:K
    errK = sum(wq.*Jq(:,e).*(Vq*p(:,e)-pex(:,e)).^2);
    L2err = L2err + errK;
end
L2err = sqrt(L2err);
href = max(1./Fscale(:));
disp(sprintf('L2 error = %e\n',L2err))





function [rhsu, rhsv, rhsp] = acousticsRHS2D_skew(u,v,p)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

global Vq Vrq Vsq VfqFace
global rxq sxq ryq syq Jq nxq nyq sJq
global Pq Pfq Prq Psq

Globals2D;

% Define field differences at faces
dp = zeros(Nfp*Nfaces,K); dp(:) = p(vmapP)-p(vmapM);
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapP)-u(vmapM);
dv = zeros(Nfp*Nfaces,K); dv(:) = v(vmapP)-v(vmapM);
uavg = zeros(Nfp*Nfaces,K); uavg(:) = .5*(u(vmapP) + u(vmapM));
vavg = zeros(Nfp*Nfaces,K); vavg(:) = .5*(v(vmapP) + v(vmapM));

% Impose reflective boundary conditions (p+ = -p-)
du(mapB) = 0; dv(mapB) = 0; dp(mapB) = -2*p(vmapB);

% can interp numerical fluxes *after* since mult by sJ is higher order
dp = VfqFace*dp;
du = VfqFace*du;
dv = VfqFace*dv;
uavg = VfqFace*uavg;
vavg = VfqFace*vavg;

% evaluate upwind fluxes
ndotdU = nxq.*du + nyq.*dv;
ndotavgU = nxq.*uavg + nyq.*vavg;

global tau
% fluxp =  tau*dp - ndotdU;

fluxp =  tau*dp - 2*ndotavgU;
fluxu =  (tau*ndotdU - dp).*nxq;
fluxv =  (tau*ndotdU - dp).*nyq;

% local derivatives of fields
pr = Vrq*p; ps = Vsq*p;
ur = Vrq*u; us = Vsq*u;
vr = Vrq*v; vs = Vsq*v;

dpdx = pr.*rxq + ps.*sxq;
dpdy = pr.*ryq + ps.*syq;
dudx = ur.*rxq + us.*sxq;
dvdy = vr.*ryq + vs.*syq;
divU = dudx + dvdy;

% compute right hand sides of the PDE's
%rhsp =  Pq*(-divU.*Jq) + Pfq*(fluxp.*sJq/2.0);
uq = Vq*u;
vq = Vq*v;
Ur = rxq.*uq + ryq.*vq;
Us = sxq.*uq + syq.*vq;
rhsp =  Prq*(Ur.*Jq) + Psq*(Us.*Jq) + Pfq*(fluxp.*sJq/2.0);
rhsu =  Pq*(-dpdx.*Jq) + Pfq*(fluxu.*sJq/2.0);
rhsv =  Pq*(-dpdy.*Jq) + Pfq*(fluxv.*sJq/2.0);

% apply inverse mass matrix
rhsp = Pq*((Vq*rhsp)./Jq);
rhsu = Pq*((Vq*rhsu)./Jq);
rhsv = Pq*((Vq*rhsv)./Jq);

end


