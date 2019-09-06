clear
    
N = 6;
K1D = 32;
CFL = .25;
FinalTime = .1;

global Qr Qs wJq Fmask E wq
global x y rxJ sxJ ryJ syJ nx ny sJ mapPq 
global QrFV QsFV

%% low order ops

[r s wq wf QrFV QsFV VfFV Fmask] = build_meshfree_ops(N);

%% SBP ops
load ~/Downloads/QuadratureRules.mat

quadpts = Q_GaussLegendre{N};
[r1D,w1D] = JacobiGQ(0,0,N);
% quadpts = Q_GaussLobatto{N};
% [r1D,w1D] = JacobiGL(0,0,N+1); % 2*(N+1)-1 = 2*N+1

Nfaces = 3;
Nfp = length(r1D);
NODETOL = 1e-8;

rs = quadpts.Points;
wq = quadpts.Weights;

rq = rs(:,1);
sq = rs(:,2);

% vol nodes
Vq = Vandermonde2D(N,rq,sq);
[Vr Vs] = GradVandermonde2D(N,rq,sq);
Dr = Vq\Vr; Ds = Vq\Vs;

M = (Vq'*diag(wq)*Vq);
Pq = M\(Vq'*diag(wq));

% face nodes
ef = ones(size(r1D));
rf = [r1D; r1D; -ef];
sf = [-ef; -r1D; r1D];
wf = [w1D; w1D; w1D];
Vf = Vandermonde2D(N,rf,sf);
Ef = Vf*Pq;

% redefine Vf nodally
Vf = zeros(length(rf),length(rq));
Fmask = [];
for i = 1:length(rf)
    fid = find(abs(rq-rf(i))+abs(sq-sf(i))<1e-10);
    Fmask = [Fmask; fid];
    Vf(i,fid) = 1;
end
Fmask = reshape(Fmask,length(Fmask)/3,3);

e = ones(length(rf)/3,1);
nrJ = [0*e; e; -e];
nsJ = [-e; e; 0*e];

Qr = diag(wq)*Vq*Dr*Pq;
Qs = diag(wq)*Vq*Ds*Pq;

VN = [eye(length(wq));Vf];
Br = diag(nrJ.*wf);
Bs = diag(nsJ.*wf);
QNr = [Qr-.5*Ef'*Br*Ef .5*Ef'*Br;
    -.5*Br*Ef .5*Br];
QNs = [Qs-.5*Ef'*Bs*Ef .5*Ef'*Bs;
    -.5*Bs*Ef .5*Bs];

Qr = VN'*QNr*VN;
Qs = VN'*QNs*VN;

% conversion from meshfree
r = rq;
s = sq;

wfvol = zeros(size(r));
wfvol(Fmask(:)) = wf;
nrvol = zeros(size(r));
nrvol(Fmask(:)) = nrJ;
nsvol = zeros(size(r));
nsvol(Fmask(:)) = nsJ;

if 0
    norm(Qr+Qr' - diag(wfvol.*nrvol),'fro')
    norm(Qs+Qs' - diag(wfvol.*nsvol),'fro')
end

rf = Vf*r;
sf = Vf*s;
Dr = diag(1./wq)*Qr;
Ds = diag(1./wq)*Qs;
L = diag(1./wq)*Vf'*diag(wf);
E = Vf'*diag(wf);

ef = ones(length(Fmask(:))/3,1);
nr = [0*ef; ef; -ef];
ns = [-ef; ef; 0*ef];
sJref = sqrt(nr.^2+ns.^2);
nr = nr./sJref;
ns = ns./sJref;


%%

[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);

iids = find(abs(abs(VX)-1)>1e-8 & abs(abs(VY)-1)>1e-8);

% VX(iids) = VX(iids) + .25/K1D*randn(size(iids));
% VY(iids) = VY(iids) + .25/K1D*randn(size(iids));

va = EToV(:,1)'; vb = EToV(:,2)'; vc = EToV(:,3)';
x = 0.5*(-(r+s)*VX(va)+(1+r)*VX(vb)+(1+s)*VX(vc));
y = 0.5*(-(r+s)*VY(va)+(1+r)*VY(vb)+(1+s)*VY(vc));

% Build connectivity matrix
[EToE, EToF] = tiConnect2D(EToV);

Nfaces = 3;
Nfp = length(rf)/Nfaces;
NODETOL = 1e-10;

% plot(x,y,'.')
% return

% make periodic
xperiod = max(VX)-min(VX);
yperiod = max(VY)-min(VY);
for k1=1:K
    for f1=1:Nfaces
        cx1 = sum(x(Fmask(:,f1), k1))/Nfp;
        cy1 = sum(y(Fmask(:,f1), k1))/Nfp;
        
        k2 = EToE(k1,f1); f2 = EToF(k1,f1);
        if(k2==k1)
            for k2=1:K
                if(k1~=k2)
                    for f2=1:Nfaces
                        if(EToE(k2,f2)==k2)
                            
                            cx2 = sum(x(Fmask(:,f2), k2))/Nfp;
                            cy2 = sum(y(Fmask(:,f2), k2))/Nfp;
                            
                            dx = sqrt( (abs(cx1-cx2)-xperiod)^2 + (cy1-cy2)^2);
                            dy = sqrt( (cx1-cx2)^2 + (abs(cy1-cy2)-yperiod)^2);
                            
                            if(dx<NODETOL | dy<NODETOL)
                                EToE(k1,f1) = k2;  EToE(k2,f2) = k1;
                                EToF(k1,f1) = f2;  EToF(k2,f2) = f1;
                            end
                        end
                    end
                end
            end
        end
    end
end

[rx,sx,ry,sy,J] = GeometricFactors2D(x,y,Dr,Ds);
rxJ = rx.*J; ryJ = ry.*J;
sxJ = sx.*J; syJ = sy.*J;

nrJ = nr.*sJref;
nsJ = ns.*sJref;
nxJ = rxJ(Fmask(:),:).*repmat(nrJ,1,K) + sxJ(Fmask(:),:).*repmat(nsJ,1,K);
nyJ = ryJ(Fmask(:),:).*repmat(nrJ,1,K) + syJ(Fmask(:),:).*repmat(nsJ,1,K);

sJ = sqrt(nxJ.^2 + nyJ.^2);
nx = nxJ./sJ; ny = nyJ./sJ;
Fscale = sJ./J(Fmask(:),:);

wJq = diag(wq)*J;

%%

% make quadrature face maps
xf = Vf*x;
yf = Vf*y;

Nfaces = 3;
mapMq = reshape(1:length(xf(:)),Nfp*Nfaces,K);
mapPq = mapMq;

for e = 1:K
    for f = 1:Nfaces
        enbr = EToE(e,f);
        if e ~= enbr % if it's a neighbor
            fnbr = EToF(e,f);
            id1 = (1:Nfp) + (f-1)*Nfp;
            id2 = (1:Nfp) + (fnbr-1)*Nfp;
            x1 = xf(id1,e); y1 = yf(id1,e);
            x2 = xf(id2,enbr); y2 = yf(id2,enbr);
            
            [X1 Y1] = meshgrid(x1,y1);
            [X2 Y2] = meshgrid(x2,y2);
            DX = (X1-X2').^2;
            DY = (Y1-Y2').^2;
            D = DX + DY;
            [p,~] = find(D<1e-8);
            
            if length(p) == 0 % if no matches                
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
                
            end
            mapPq(id1,e) = id2(p) + (enbr-1)*(Nfp*Nfaces);
        end
    end
end

% for i = 1:length(mapPq(:))
%     plot(xf,yf,'.')
%     hold on
%     plot(xf(i),yf(i),'o','markersize',16)
%     plot(xf(mapPq(i)),yf(mapPq(i)),'x','markersize',16)
%     hold off    
%     pause
% end

%% simulation

if 1
    % modify to make first order skew VFV operators
    % seems to make no difference? 
    Br = diag(sum(QrFV,1));
    Bs = diag(sum(QsFV,1));
    QrFV = .5*(QrFV-QrFV')+.5*Br;
    QsFV = .5*(QsFV-QsFV')+.5*Bs;
end

rk4a = [            0.0 ...
    -567301805773.0/1357537059087.0 ...
    -2404267990393.0/2016746695238.0 ...
    -3550918686646.0/2091501179385.0  ...
    -1275806237668.0/842570457699.0];
rk4b = [ 1432997174477.0/9575080441755.0 ...
    5161836677717.0/13612068292357.0 ...
    1720146321549.0/2090206949498.0  ...
    3134564353537.0/4481467310338.0  ...
    2277821191437.0/14882151754819.0];
rk4c = [             0.0  ...
    1432997174477.0/9575080441755.0 ...
    2526269341429.0/6820363962896.0 ...
    2006345519317.0/3224310063776.0 ...
    2802321613138.0/2924317926251.0];

ssprk = [0 3/4 1/3];

u = exp(-25*(x.^2+y.^2));
u = 0*x;
u(:,abs(mean(x))<.25 & abs(mean(y))<.25) = 1;
% u = 1+sin(pi*x).*sin(pi*y); 

dx = 2/K1D;
CN = (N+1)*(N+2)/2;
dt = CFL*dx / CN;

Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;
res = zeros(size(x));

figure(1)
for i = 1:Nsteps
    
    uold = u;
    for INTRK = 1:3
        
        [rhsH, rhsL] = brhs(u,dt);
        
        % low storage SSP 
        uH = u + dt*rhsH;
        uL = u + dt*rhsL;
        
        r = uH;
        ids = min(uH)<0;
        t = max(0,min(1,min(uL,[],1)./(min(uH,[],1) - min(uL,[],1))));
        t = repmat(t(ids),size(u,1),1);
        
        r(:,ids) = uL(:,ids) + t.*(uH(:,ids)-uL(:,ids));
%         t = 1;
%         u = uL + t.*(uH-uL);
        
        u = ssprk(INTRK)*uold + (1-ssprk(INTRK))*r;
        
    end
            
    if 1
        ufunc(i) = sum(sum((diag(wq)*J).*u.^2));
    end
    
    if i==1 || mod(i,10) == 0 || i==Nsteps
        clf
        vv = u;
        color_line3(x,y,vv,vv,'.');
        title(sprintf('time = %f, step %d / %d\n',dt*i,i));
%         view(3)
        drawnow
    end
end

figure(2)
semilogy(ufunc,'o--')
hold on


function [rhsH, rhsL] = brhs(u,dt)

global Qr Qs wJq Fmask E wq
global QrFV QsFV
global x y rxJ sxJ ryJ syJ nx ny sJ mapPq 

% rhsL = zeros(size(u));
% rhsH = zeros(size(u));

uM = u(Fmask(:),:);
uP = uM(mapPq);
flux = nx.*(uP-uM);

% both low/high order schemes use same E matrix
rhsL = .5*E*(sJ.*flux);
rhsH = rxJ.*(Qr*u) + sxJ.*(Qs*u) + .5*E*(sJ.*flux);

% upwinding 
dflux = abs(nx).*(uP-uM);
rhsH = rhsH - .5*E*(sJ.*dflux);
rhsL = rhsL - .5*E*(sJ.*dflux);

% add local graph viscosity contributions
Lmax = zeros(size(Qr));
for e = 1:size(u,2)        
    
    % wavespeed = [1;0]
    QxFV = rxJ(1,e)*QrFV + sxJ(1,e)*QsFV;
    QyFV = ryJ(1,e)*QrFV + syJ(1,e)*QsFV;
    D2 = sqrt(QxFV.^2 + QyFV.^2);
    ids = D2(:)>1e-10;
    Lmax(ids) = abs(QxFV(ids)./D2(ids)); % | beta dot nij |
    D = max(Lmax.*D2,(Lmax.*D2)');
    D = D - diag(sum(D,2)); % make it so D*1 = 0
    
    FL = (QxFV-D)*u(:,e);
    
    rhsL(:,e) = rhsL(:,e) + FL;
end

% invert lumped mass afterwards
rhsH = -rhsH./wJq;
rhsL = -rhsL./wJq;
end

