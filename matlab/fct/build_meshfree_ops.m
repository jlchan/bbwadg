% clear
% N = 1;

function [rq sq wq wf Qr Qs Vf Fmask] = build_meshfree_ops(N)

load ~/Downloads/QuadratureRules.mat

quadpts = Q_GaussLegendre{N}; [r1D w1D] = JacobiGQ(0,0,N);
% quadpts = Q_GaussLobatto{N}; [r1D w1D] = JacobiGL(0,0,N);
xy = quadpts.Points;
wq  = quadpts.Weights;

rq = xy(:,1);
sq = xy(:,2);

% face nodes
ef = ones(size(r1D));
rf = [r1D; r1D; -ef];
sf = [-ef; -r1D; r1D];
wf = [w1D; sqrt(2)*w1D; w1D]; % rescale later for compatibility with DG code
nr = [0*ef; ef; -ef];
ns = [-ef; ef; 0*ef];
sJref = sqrt(nr.^2+ns.^2);
nr = nr./sJref;
ns = ns./sJref;

fid = [];
for i = 1:length(rf)
    fid = [fid; find(abs(rq-rf(i))+abs(sq-sf(i))<1e-10)];
end

% extract face-vol data
wfvol = 0*rq;
nrvol = 0*rq;
nsvol = 0*rq;
wfvol(fid) = wf;
nrvol(fid) = nr;
nsvol(fid) = ns;

ep = 3/(N+1);
[re se] = rstoxy(rq,sq);
% xe = x; ye = y;

adj = zeros(length(rq));
for i = 1:length(rq)
    d2 = (re(i)-re).^2 + (se(i) - se).^2;
    [d2sort,p] = sort(d2);
    r2 = 1.01*d2sort(3);
    adj(i,d2 <= r2) = 1; % minimum nbrs for linear recon
    adj(i,i) = 0;
    rad(i) = sqrt(r2);
end
adj = (adj + adj' > 0);

if min(sum(adj,2)) < 2 % 2 = min num nbrs for recon
    min(sum(adj,2))
    keyboard
end

if 0
    tc = linspace(0,2*pi,100)';
    xc = cos(tc); yc = sin(tc);
    
    plot(re,se,'o','markersize',12)
    hold on
    text(re+.05,se,num2str((1:length(re))'),'fontsize',16)
    for i = 1:length(re)
        idi = find(adj(i,:));
        for j = 1:length(idi)
            plot([re(i) re(idi(j))],[se(i) se(idi(j))],'k--')
        end
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
    L{k} = zeros(length(rq));
    fx{k} = zeros(size(rq));
    fy{k} = zeros(size(rq));
end

dxphi = dxVfun(rq,sq);
dyphi = dyVfun(rq,sq);
phif = zeros(length(rq),3);
phif(fid,:) = Vfun(rq(fid),sq(fid));
for k = 1:3
    fx{k} = wq.*dxphi(:,k) - (wfvol.*nrvol).*phif(:,k);
    fy{k} = wq.*dyphi(:,k) - (wfvol.*nsvol).*phif(:,k);
end

for i = 1:length(rq)
    inbr = find(adj(i,:));
    for jj = 1:length(inbr)
        phia = Vfun(.5*(rq(i) + rq(inbr(jj))), .5*(sq(i) + sq(inbr(jj))));
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

u = 1 + (rq + sq);

Qr = zeros(length(u));
Qs = zeros(length(u));
u = 0*u;
for ii = 1:length(u)
    u(ii) = 1;
    
    % apply operator
    duxf = zeros(size(rq));
    duyf = zeros(size(rq));
    for i = 1:length(u)
        
        inbr = find(adj(i,:));
        
        % compute GMLS recon at point
        idi = [i inbr];
        ci = Vfun(rq(idi),sq(idi)) \ u(idi);
        
        for jj = 1:length(inbr)
            j = inbr(jj);
            
            % compute GMLS recon at nbr point
            idj = [j find(adj(j,:))];
            cj = Vfun(rq(idj),sq(idj)) \ u(idj);
            
            % basis/coeffs at virtual face
            cij = .5*(ci + cj);
            phia = Vfun(.5*(rq(i) + rq(j)),.5*(sq(i) + sq(j)))';
            
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

Vf = zeros(length(rf),length(rq));
Fmask = [];
for i = 1:length(rf)
    fid = find(abs(rq-rf(i))+abs(sq-sf(i))<1e-10);
    Fmask = [Fmask; fid];
    Vf(i,fid) = 1;
end
Fmask = reshape(Fmask,length(rf)/3,3); 

% rescale wf for compatibility with DG and sJ
wf = [w1D; w1D; w1D];

% check properties are satisfied
if 0
    norm((Qr*(1+rq+sq))./wq-1)
    norm(Qs*(1+rq+sq)./wq-1)
    norm(sum(Qr,1) - (wfvol.*nrvol)')
    norm(sum(Qs,1) - (wfvol.*nsvol)')
end

