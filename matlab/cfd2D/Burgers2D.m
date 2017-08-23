clear
Globals2D

N = 5;
K1D = 8;
FinalTime = 1;
CFL = .5;

[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
StartUp2D;

BuildPeriodicMaps2D(2,2);
% StartUp2D;

% plotting nodes
[rp sp] = EquiNodes2D(50); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

global Vq Pq Lq Lqf Vfqf Vfq Pfqf
global rxJ sxJ ryJ syJ nxJ nyJ nrJ nsJ
global mapPq
Nq = 2*N+1;
[rq sq wq] = Cubature2D(Nq); % integrate u*v*c
Vq = Vandermonde2D(N,rq,sq)/V;
M = Vq'*diag(wq)*Vq;
Pq = M\(Vq'*diag(wq)); % J's cancel out
xq = Vq*x; yq = Vq*y;
Jq = Vq*J;

rxJ = rx.*J; sxJ = sx.*J;
ryJ = ry.*J; syJ = sy.*J;
rxJ = Vq*rxJ; sxJ = Vq*sxJ;
ryJ = Vq*ryJ; syJ = Vq*syJ;

[rq1D wq1D] = JacobiGQ(0,0,N+1);
rfq = [rq1D; -rq1D; -ones(size(rq1D))];
sfq = [-ones(size(rq1D)); rq1D; -rq1D];
wfq = [wq1D; wq1D; wq1D];
Vq1D = Vandermonde1D(N,rq1D)/Vandermonde1D(N,JacobiGL(0,0,N));
% plot(rfq,sfq,'o')
Nfq = length(rq1D);

Vfq = Vandermonde2D(N,rfq,sfq)/V;
Vfqf = kron(eye(3),Vq1D);
Mf = Vfq'*diag(wfq)*Vfq;
Lq = M\(Vfq'*diag(wfq));

Pq1D = (Vq1D'*diag(wq1D)*Vq1D) \ (Vq1D'*diag(wq1D));
Pfqf = kron(eye(3),Pq1D);

nxJ = Vfqf*(nx.*sJ);
nyJ = Vfqf*(ny.*sJ);
% Fscaleq = Vfqf*Fscale;

nrJ = [-zeros(size(rq1D)); ones(size(rq1D)); -ones(size(rq1D)); ];
nsJ = [-ones(size(rq1D)); ones(size(rq1D)); -zeros(size(rq1D)); ];

% % assume affine for now
% (Dx*u,vw) + .5*<uP - P(u),vw> + .5*<uM - P(u),P(vw)>
% (rxJ.*(Dr*u) + sxJ.*(Ds*u),vw)  + .5*(L*(uM - P(u))*nxJ,vw) + .5*<uP - P(u),vw*nxJ>
% Vq*Dr*(rxJ.*u) - .5*Vq*L*P(u) + inv(M)*.5*<uP+uM - P(u),vw>

Drq = Vq*(Dr - .5*Lq*diag(nrJ)*Vfq)*Pq;
diag(wq)*Drq

global fS
fS = @(ux,uy) (ux.^2 + ux.*uy + uy.^2)/6;

if 0
    u = (1:Np)';
    uq = Vq*u;
    uf = Vfq*u;
    
    [ux uy] = meshgrid(uq);
    [ufx ufy] = meshgrid(uf,uq);
    
    FS = fS(ux,uy);
    FSf = fS(ufx,ufy);
    
    Pq*sum((diag(wq)*Vq*Lq).*FSf,2)
    sum((diag(wq)*Vq*Lq).*FSf,1)*(M\Vfq')'
    M\Vfq'*sum((diag(wfq)*Vfq*Pq).*FSf',2)
end

%% make quadrature face maps

if 1
    xf = Vfq*x;
    yf = Vfq*y;
    
    mapMq = reshape(1:nnz(xf),Nfq*Nfaces,K);
    mapPq = mapMq;
    
    for e = 1:K
        for f = 1:Nfaces
            enbr = EToE(e,f);
            if e ~= enbr % if it's a neighbor
                fnbr = EToF(e,f);
                id1 = (1:Nfq) + (f-1)*Nfq;
                id2 = (1:Nfq) + (fnbr-1)*Nfq;
                x1 = xf(id1,e); y1 = yf(id1,e);
                x2 = xf(id2,enbr); y2 = yf(id2,enbr);
                
                [X1 Y1] = meshgrid(x1,y1);
                [X2 Y2] = meshgrid(x2,y2);
                DX = (X1-X2').^2;
                DY = (Y1-Y2').^2;
                D = DX + DY;
                [p,~] = find(D<1e-8);
                
                if length(p) == 0
                    %                     keyboard
                    % assume periodic boundary, find match in x,y
                    [px,~] = find(DX<1e-8);
                    [py,~] = find(DY<1e-8);
                    if length(px)==0
                        p = py;
                    elseif length(py)==0
                        p = px;
                    else
                        keyboard
                    end
                    
                    %                     plot(x1,y1,'o')
                    %                     hold on
                    %                     plot(x2,y2,'x')
                    %                     axis([-1,1,-1,1])
                    %                     keyboard
                    
                end
                mapPq(id1,e) = id2(p) + (enbr-1)*(Nfq*Nfaces);
            end
        end
    end
    
    %     keyboard
    
    %     mapB = find(mapMq==mapPq);
    %     xb = xf(mapB); yb = yf(mapB);
    %
    %     plot(xb,yb,'o')
    %     return
else
    
    mapPq = reshape(mapP,Nfq*Nfaces,K);
    
end

% hold on
% for i = 1:nnz(xf)
%     id = mapPq(i) ;
%     plot(xf(i),yf(i),'o')
%     plot(xf(id),yf(id),'x')
%    pause
% end
% return

%% params setup

x0 = 0; y0 = 0;

% u = Pq*exp(-5^2*((xq-x0).^2 + (yq-y0).^2));
u = Pq*-sin(pi*(xq+yq));

%% testing

global wJq
wJq = diag(wq)*(Vq*J);

% u = randn(Np,K);
% uq = Vq*u;
%
% dudx = rx.*(Dr*u) + sx.*(Ds*u);
% rhs = J.*dudx - .5*Lq*(nxJ.*(Vfq*u));
% % rhsu = (1/3)*Pq*(dfdx + uq.*dudx) + Lq*(nxJ.*((1/3)*Vfq*fproj));
% rhs = rhs./J;
% global wJq
% sum(sum(wJq.*((Vq*rhs).*uq)))
% return

% xf = Vfq*x; yf = Vfq*y;
% norm(xf-xf(mapPq),'fro')
% norm(yf-yf(mapPq),'fro')
% for e = 1:K
%     clf
%     quiver(xf(:,e),yf(:,e),nxJ(:,e),nyJ(:,e))
%     pause
% end
% return

%%

% Runge-Kutta residual storage
resu = zeros(Np,K); resv = zeros(Np,K); resp = zeros(Np,K);

% compute time step size
CN = (N+1)^2/2; % guessing...
CNh = max(CN*max(Fscale(:)));
dt = CFL*2/CNh;

Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

figure(1)
for i = 1:Nsteps
    for INTRK = 1:5
        
        timelocal = dt*i + rk4c(INTRK)*dt;
        
        rhsu  = advecRHS2D(u,timelocal);
        
        if INTRK==5
            uq = Vq*u;
            rhstest(i) = sum(sum(wJq.*(Vq*rhsu).*uq));
        end
        
        resu = rk4a(INTRK)*resu + dt*rhsu;
        u = u+rk4b(INTRK)*resu;
        
    end;
    
    uq = Vq*u;
    energy(i) = sum(sum(wJq.*uq.^2));
    
    if mod(i,10)==0 || i==Nsteps
        clf
        pp = u;
        vv = Vp*pp;
        color_line3(xp,yp,vv,vv,'.');
        axis equal
        axis tight
        colorbar
        title(sprintf('time = %f',dt*i))
        %         view(3)
        drawnow
        
    end
    
end
figure(2)
semilogy(dt*(1:Nsteps),energy,'--')
hold on
semilogy(dt*(1:Nsteps),abs(rhstest),'x')


function [rhsu] = advecRHS2D(u,time)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;

global Vq Pq Lq Lqf Vfqf Vfq Pfqf
global rxJ sxJ ryJ syJ nxJ nyJ nrJ nsJ
global mapPq
global fS

% operators
% Pq
% Vq
% Vfq Lq
Drq = (Vq*Dr*Pq - .5*Vq*Lq*diag(nrJ)*Vfq*Pq);
Dsq = (Vq*Ds*Pq - .5*Vq*Lq*diag(nsJ)*Vfq*Pq);
Lrq = (.5*Vq*Lq*diag(nrJ));
Lsq = (.5*Vq*Lq*diag(nsJ));

uq = Vq*u;
uM = Vfq*u;

% tau = 1;
uP = uM(mapPq);
fSf = fS(uP,uM);

rhsu = zeros(Np,K);
for e = 1:K
    
    [ux uy] = meshgrid(uq(:,e));
    [ufx ufy] = meshgrid(uM(:,e),uq(:,e));
    
    FS = fS(ux,uy);
    FSf = fS(ufx,ufy);
    
    fr = sum(Drq.*FS,2) + sum(Lrq.*(FSf),2); 
    fs = sum(Dsq.*FS,2) + sum(Lsq.*(FSf),2);
    
    % apply J*G^T - dfdx here.
    divF = rxJ(:,e).*fr + sxJ(:,e).*fs + ryJ(:,e).*fr + syJ(:,e).*fs;
    
    fSproj = sum((Vfq*Pq).*FSf',2);
    rhsu(:,e) =  Pq*divF + Lq*(.5*(nxJ(:,e)+nyJ(:,e)).*(fSf(:,e) - fSproj));
end

rhsu = -2*rhsu./J;

end

