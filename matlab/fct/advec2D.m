clear
    
N = 4;
K1D = 16;
CFL = .33;
FinalTime = .25;

opt = 2; % 1 = meshfree, 2 = SBP, 3 = nodal polynomial
useAV = 0;

if opt==1        
    
    load ~/Downloads/QuadratureRules.mat
    
    quadpts = Q_GaussLegendre{N}; [r1D w1D] = JacobiGQ(0,0,N);
    % quadpts = Q_GaussLobatto{N}; [r1D w1D] = JacobiGL(0,0,N);
    xy = quadpts.Points;
    wq  = quadpts.Weights;
    
    r = xy(:,1);
    s = xy(:,2);
    
    Nfaces = 3;
    Nfp = length(r1D);
    NODETOL = 1e-8;
    
    % face nodes
    ef = ones(size(r1D));
    xf = [r1D; r1D; -ef];
    yf = [-ef; -r1D; r1D];
    %wf = [w1D; sqrt(2)*w1D; w1D];
    wf = [w1D; w1D; w1D];
    nrJ = [0*ef; ef; -ef];
    nsJ = [-ef; ef; 0*ef];
    sJref = sqrt(nrJ.^2+nsJ.^2);
    nr = nrJ./sJref;
    ns = nsJ./sJref;
    
    fid = [];
    for i = 1:length(xf)
        fid = [fid; find(abs(r-xf(i))+abs(s-yf(i))<1e-10)];
    end
    
    % extract face-vol data
    wfvol = 0*r;
    nrvol = 0*r;
    nsvol = 0*r;
    wfvol(fid) = wf;
    nrvol(fid) = nrJ;
    nsvol(fid) = nsJ;
    
    ep = 3/(N+1);
    [re se] = rstoxy(r,s);
%     re = r; se = s;
    
    adj = zeros(length(r));
    for i = 1:length(r)        
        d2 = (re(i)-re).^2 + (se(i) - se).^2;
        [d2sort,p] = sort(d2);
        adj(i,d2 <= 1.01*d2sort(3)) = 1; % minimum nbrs for linear recon
        adj(i,i) = 0;        
        rad(i) = sqrt(1.01*d2sort(3));
    end
    adj = (adj + adj' > 0);
    
    if min(sum(adj,2)) < 2 % 2 = min num nbrs for recon
        min(sum(adj,2))
        keyboard
    end

    if 0
        plot(re,se,'o')
        hold on
        text(re+.05,se,num2str((1:length(re))'))
        t = linspace(-pi,pi,50);
        rc = cos(t); sc = sin(t);
        for i = 1:length(re)
            nbrs = find(adj(i,:));
            for j = 1:length(nbrs)
                plot([re(i); re(nbrs(j))],[se(i); se(nbrs(j))],'k--')
            end
        end
        
        for i = 1:length(re)
            plot(re(i)+rad(i)*rc,se(i)+rad(i)*sc)
            pause
        end
        return
    end

    %% build graph laplacians
    
    % basis for p1
    Vfun = @(x,y) [x.^0 x.^1 y(:).^1];
    dxVfun = @(x,y) [0*x x.^0 0*x];
    dyVfun = @(x,y) [0*x 0*x x.^0];
    
    % build graph laplacians
    L = {};
    fx = {};
    for k = 1:3
        L{k} = zeros(length(r));
        fx{k} = zeros(size(r));
        fy{k} = zeros(size(r));
    end
    
    dxphi = dxVfun(r,s);
    dyphi = dyVfun(r,s);
    phif = zeros(length(r),3);
    phif(fid,:) = Vfun(r(fid),s(fid));
    for k = 1:3
        fx{k} = wq.*dxphi(:,k) - (wfvol.*nrvol).*phif(:,k);
        fy{k} = wq.*dyphi(:,k) - (wfvol.*nsvol).*phif(:,k);
    end
    
    for i = 1:length(r)
        inbr = find(adj(i,:));
        for jj = 1:length(inbr)
            phia = Vfun(.5*(r(i) + r(inbr(jj))), .5*(s(i) + s(inbr(jj))));
            for k = 1:3
                L{k}(i,i) = L{k}(i,i) + phia(k);
                L{k}(i,inbr(jj)) = L{k}(i,inbr(jj)) - phia(k);
            end
        end
    end
    
    % solve for potentials
    psix = [pinv(L{1})*fx{1} pinv(L{2})*fx{2} pinv(L{3})*fx{3}]';
    psiy = [pinv(L{1})*fy{1} pinv(L{2})*fy{2} pinv(L{3})*fy{3}]';
    
    %% apply GMLS one node at a time
    
    u = 1 + (r + s);
    
    Qr = zeros(length(u));
    Qs = zeros(length(u));
    u = 0*u;
    for ii = 1:length(u)
        u(ii) = 1;
        
        % apply operator
        duxf = zeros(size(r));
        duyf = zeros(size(r));
        for i = 1:length(u)
            
            inbr = find(adj(i,:));
            
            % compute GMLS recon at point
            idi = [i inbr];
            ci = Vfun(r(idi),s(idi)) \ u(idi);
            
            for jj = 1:length(inbr)
                j = inbr(jj);
                
                % compute GMLS recon at nbr point
                idj = [j find(adj(j,:))];
                cj = Vfun(r(idj),s(idj)) \ u(idj);
                
                % basis/coeffs at virtual face
                cij = .5*(ci + cj);
                phia = Vfun(.5*(r(i) + r(j)),.5*(s(i) + s(j)))';
                
                muxij = phia .* (psix(:,i) - psix(:,j));
                duxf(i) = duxf(i) + muxij'*cij;
                
                muyij = phia .* (psiy(:,i) - psiy(:,j));
                duyf(i) = duyf(i) + muyij'*cij;
            end
        end
        duxf = duxf + wfvol.*nrvol.*u; % add boundary correction back
        duyf = duyf + wfvol.*nsvol.*u; % add boundary correction back
        %     dudx = (duf)./w;
        
        u(ii) = 0;
        
        Qr(:,ii) = duxf;
        Qs(:,ii) = duyf;
    end
    
    Br = diag(wfvol.*nrvol);
    Bs = diag(wfvol.*nsvol);

    % conversion from meshfree
    rq = r;
    sq = s;    
    
    rf = r(fid);
    sf = s(fid);
    
    e = ones(length(rf)/3,1);
    nrJ = [0*e; e; -e];
    nsJ = [-e; e; 0*e];
    
    Fmask = reshape(fid,length(fid)/3,3);
    E = zeros(length(rf),length(r));
    for i = 1:length(rf)
        fid = find(abs(rq-rf(i))+abs(sq-sf(i))<1e-10);        
        E(i,fid) = 1;
    end    
    
    wfvol = zeros(size(r));
    wfvol(Fmask(:)) = wf;
    nrvol = zeros(size(r));
    nrvol(Fmask(:)) = nrJ;
    nsvol = zeros(size(r));
    nsvol(Fmask(:)) = nsJ;

    
    norm(Qr*(1+r+s)./wq-1)
    norm(Qs*(1+r+s)./wq-1)
    
    plot(re,se,'o')
    hold on
    for i = 1:length(re)
        nbrs = find(adj(i,:));
        for j = 1:length(nbrs)
            plot([re(i); re(nbrs(j))],[se(i); se(nbrs(j))],'--')
        end
    end
%     return
%     
%     keyboard                    

elseif opt==2
    % SBP ops
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
    
    E = zeros(length(rf),length(rq));
    Fmask = [];
    for i = 1:length(rf)
        fid = find(abs(rq-rf(i))+abs(sq-sf(i))<1e-10);
        Fmask = [Fmask; fid];
        E(i,fid) = 1;
    end
    Fmask = reshape(Fmask,length(Fmask)/3,3);
    
    e = ones(length(rf)/3,1);
    nrJ = [0*e; e; -e];
    nsJ = [-e; e; 0*e];
    
    Qr = diag(wq)*Vq*Dr*Pq;
    Qs = diag(wq)*Vq*Ds*Pq;
    
    VN = [eye(length(wq));E];
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
    
    norm(Qr+Qr' - diag(wfvol.*nrvol),'fro')
    norm(Qs+Qs' - diag(wfvol.*nsvol),'fro')       
end


Vf = E;
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

% if opt==1
%     % modify to make first order skew VFV operators
%     % seems to make no difference? 
%     Qrsbp = .5*(Qr-Qr')+.5*Br;
%     Qssbp = .5*(Qs-Qs')+.5*Bs;
%     Dr = diag(1./wq)*Qrsbp;
%     Ds = diag(1./wq)*Qssbp;    
% end

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

% u = exp(-100*(x.^2+y.^2));
u = 0*x;
u(:,abs(mean(x))<.25 & abs(mean(y))<.25) = 1;
% u = sin(pi*x).*sin(pi*y); 

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
        
        rhs = arhs(u,Qr,Qs,wJq,Fmask,E,rxJ,sxJ,ryJ,syJ,nx,ny,sJ,mapPq,useAV);
                
        u = ssprk(INTRK)*uold + (1-ssprk(INTRK))*(u + dt*rhs);
        
%         res = rk4a(INTRK)*res + dt*rhs;
%         u   = u  + rk4b(INTRK)*res;
        
    end
            
    if 1
        %             ufunc(i) = sum(sum((diag(wq)*J).*u));
        ufunc(i) = sum(sum((diag(wq)*J).*u.^2));
%         rhs = arhs(u,Qr,Qs,wJq,Fmask,E,rxJ,sxJ,ryJ,syJ,nx,ny,sJ,mapPq,useAV);
%         ufunc(i) = abs(sum(sum((diag(wq)*J).*(u.*rhs))));
    end
    
    if i==1 || mod(i,10) == 0 || i==Nsteps
        clf
        vv = u;
        color_line3(x,y,vv,vv,'.')
        title(sprintf('time = %f, step %d / %d\n',dt*i,i));
%         view(3)
        drawnow
    end
end

figure(2)
semilogy(ufunc,'o--')
hold on


%function rhs = arhs(u,Dr,Ds,Fmask,L,Fscale,rx,sx,ry,sy,nx,ny,mapPq,opt)
function rhs = arhs(u,Qr,Qs,wJq,Fmask,E,rxJ,sxJ,ryJ,syJ,nx,ny,sJ,mapPq,useAV)

uM = u(Fmask(:),:);
uP = uM(mapPq);
flux = nx.*(uP-uM);

dudxJ = rxJ.*(Qr*u) + sxJ.*(Qs*u);
rhs = (dudxJ + .5*E*(sJ.*flux));

if useAV==1
    % add graph viscosity
    Lmax = zeros(size(Qr));
    for e = 1:size(u,2)
        Qx = rxJ(1,e)*Qr + sxJ(1,e)*Qs;
        Qy = ryJ(1,e)*Qr + syJ(1,e)*Qs;
        
        % wavespeed = [1;0]        
        D2 = sqrt(Qx.^2 + Qy.^2);
        ids = D2(:)>1e-10;        
        Lmax(ids) = abs(Qx(ids)./D2(ids)); % | beta dot nij |
        D = max(Lmax.*D2,(Lmax.*D2)'); 
        D = D - diag(sum(D,2)); % make D conservative locally
        rhs(:,e) = rhs(:,e) - D*u(:,e);
    end   
end

% upwinding 
dflux = abs(nx).*(uP-uM);
rhs = rhs - .5*E*(sJ.*dflux);

% invert lumped mass afterwards
rhs = -rhs./wJq;

end