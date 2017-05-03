function Wave2D_spectra

Globals2D
global bx by Vq Pq Vfq Pfq Vrq Vsq Prq Psq


% Polynomial order used for approximation
N = 3;

% Read in Mesh
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Grid/Other/squarereg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Grid/Other/squareireg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Grid/Other/block2.neu');

K1D = 2;
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);

% Initialize solver and construct grid and metric
StartUp2D;

% % rebuild maps for periodic
% BuildPeriodicMaps2D(2,2);

Ne = 25;
[re se] = EquiNodes2D(Ne); [re se] = xytors(re,se);
re = re*(1-2/(Ne+1)); se = se*(1-2/(Ne+1));
Ve = Vandermonde2D(N,re,se)/V;
xe = Ve*x; ye = Ve*y;

[rp sp] = EquiNodes2D(50); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

[rq sq wq] = Cubature2D(3*N);
% [rq sq wq] = QNodes2D(N); [rq sq] = xytors(rq,sq);
Vq = Vandermonde2D(N,rq,sq)/V;
[Vrq Vsq] = GradVandermonde2D(N,rq,sq);
Vrq = Vrq/V; Vsq = Vsq/V;
xq = Vq*x; yq = Vq*y;
M = Vq'*diag(wq)*Vq;
Pq = M\(Vq'*diag(wq));
Prq = M\(Vrq'*diag(wq));
Psq = M\(Vsq'*diag(wq));


[r1D, w1D] = JacobiGQ(0,0,2*N);
% rfq = [r1D; r1D; -ones(size(r1D))];
% sfq = [-ones(size(r1D)); -r1D; r1D];
wfq = [w1D; w1D; w1D];

Vfq1 = Vandermonde1D(N,r1D)/Vandermonde1D(N,r(Fmask(:,1)));
Vfq2 = Vandermonde1D(N,r1D)/Vandermonde1D(N,r(Fmask(:,2)));
Vfq3 = Vandermonde1D(N,r1D)/Vandermonde1D(N,s(Fmask(:,3)));
Vfq = blkdiag(Vfq1,Vfq2,Vfq3);
Pfq = (Vfq'*diag(wfq)*Vfq)\(Vfq'*diag(wfq));

%% build projection to CG

R = getCGRestriction()';
P = R*((R'*kron(diag(J(1,:)),M)*R)\(R'*kron(diag(J(1,:)),M)));
P = blkdiag(P,eye(Np*K*2));

%% eigs

rLGL = JacobiGQ(0,0,N); rmin = abs(rLGL(1)-rLGL(2));
dtscale = dtscale2D; dt = min(dtscale)*rmin*2/3

%% compute matrices

e = zeros(3*Np*K,1);
A = zeros(3*Np*K);

tau = 0;
for i = 1:3*Np*K
    e(i) = 1;
    ids = 1:Np*K;
    p = reshape(e(ids),Np,K);
    u = reshape(e(ids + Np*K),Np,K);
    v = reshape(e(ids + 2*Np*K),Np,K);
    
    % dudt + dudx= 0
    [rhsp rhsu rhsv] = WaveRHS(p,u,v,tau);
    A(:,i) = [rhsp(:); rhsu(:); rhsv(:)];
    e(i) = 0;
end

tau = 1;
A1 = zeros(3*Np*K);
for i = 1:3*Np*K
    e(i) = 1;
    ids = 1:Np*K;
    p = reshape(e(ids),Np,K);
    u = reshape(e(ids + Np*K),Np,K);
    v = reshape(e(ids + 2*Np*K),Np,K);
    
    % dudt + dudx= 0
    [rhsp rhsu rhsv] = WaveRHS(p,u,v,tau);
    A1(:,i) = [rhsp(:); rhsu(:); rhsv(:)];
    e(i) = 0;
end
S = A1-A;
S(abs(S)<1e-8) = 0;

%% plotting eigenmodes

if 1
    
    
    load waveEigTrackerData
    tauvec = [0:1e-1:100];
    [~,id0] = min(abs(tauvec-.1));
    [~,id1] = min(abs(tauvec-1));
    [~,id2] = min(abs(tauvec-100));
    tauids = [id0 id1 id2];
    
    for i = 1:length(tauids)
        tau = tauvec(tauids(i));
        [W D] = eig(A + tau*S);
        d = diag(D);
        [~,p] = sort(imag(d),'descend');
        d = d(p); W = W(:,p);
        
        all_points = vertcat(points{:});
        %track = adjacency_tracks{51};
        track = adjacency_tracks{46};
        track_points = all_points(track, :);
        figure(1)
        plot(track_points(end,1), track_points(end, 2), 'bo','markersize',14,'linewidth',2)
        hold on
        plot(track_points(:,1), track_points(:, 2), 'k-','linewidth',1)
        %         plot(all_points(track(tauids(i)),1), all_points(track(tauids(i)), 2), 'bs','markersize',14,'linewidth',1)
        
        
        len = size(points{1},1);
        i1 = mod(track(tauids(i)),len);
        d(i1)
        plot(real(d(i1)),imag(d(i1)),'r*','markersize',14,'linewidth',2)
        set(gca,'fontsize',15)
        grid on
        
        
        w = real(W(:,i1));
        v = reshape(w(1:Np*K),Np,K);
        [a,p] = max(abs(v(:)));
        v = v*sign(v(p))/a;
        %         vv = Vp*v;
        %         PlotField2D(25,x,y,v);
        %         color_line3(xp,yp,vv,vv,'.');
        % view(3)
        w1 = reshape(w((1:Np*K)+(Np*K)),Np,K);
        w2 = reshape(w((1:Np*K)+2*(Np*K)),Np,K);
        [a,p] = max(abs(w1(:)));        w1 = w1*sign(w1(p))/a;
        [a,p] = max(abs(w2(:)));        w2 = w2*sign(w2(p))/a;
        v1 = Ve*w1; v2 = Ve*w2;
        figure
        PlotMesh2D;hold on;quiver(xe,ye,v1,v2)
        
        set(gca,'fontsize',15)
        grid on
        
        
                print(gcf,'-dpng',sprintf('waveEig2D_u%d.png',i))
        %         keyboard
    end
    
    keyboard
    
    hold on
    plot(real(d),imag(d),'^','markersize',10,'linewidth',2,'DisplayName',sprintf('\\tau = %d',round(tau)))
    legend show
    axis([-400 50 -20 20])
    % ax = axis;
    % axis(1.05*ax);
    % axis equal
    set(gca,'fontsize',15)
    grid on
    %     print(gcf,'-dpng','waveEigs2.png')
    
end
%%

keyboard

%%

tauvec = [0:1e-1:100];
% tauvec = [0:5e-3:4];

if 1
    for ii = 1:length(tauvec)
        tau = tauvec(ii);
        
        if mod(ii,10)==0
            disp(sprintf('computing mode %d out of %d\n',ii,length(tauvec)))
        end
        [W D] = eig(A + tau*S);
        d = diag(D);
        [~,p] = sort(imag(d),'descend');
        d = d(p); W = W(:,p);
        
        points{ii} = [real(d) imag(d)];
        
    end
    
    max_linking_distance = 10;
    max_gap_closing = Inf;
    debug = true;
    disp('running simple tracker')
    [ tracks, adjacency_tracks ] = simpletracker(points,...
        'MaxLinkingDistance', max_linking_distance, ...
        'MaxGapClosing', max_gap_closing, ...
        'Debug', debug);
    
    save waveEigTrackerData points tracks adjacency_tracks
    keyboard
    
    
else
    load waveEigTrackerData
end

%%


plot(points{1}(:,1),points{1}(:,2),'bo','markersize',10,'DisplayName','\tau = 0')
hold on
[~,id] = min(abs(tauvec-1));
plot(points{id}(:,1),points{id}(:,2),'r^','markersize',10,'DisplayName','\tau = 1.0')
plot(points{end}(:,1),points{end}(:,2),'gs','markersize',10,'DisplayName',sprintf('\\tau = %1.1f',tauvec(end)))

legend show

hold on
n_tracks = numel(tracks);
all_points = vertcat(points{:});
for i_track = 1 : n_tracks
    
    track = adjacency_tracks{i_track};
    track_points = all_points(track, :);
    %     if abs(track_points(end,1)) < 5  && track_points(end,2) < 11 && track_points(end,2) > 4
    if (i_track==51)
        plot(track_points(:,1), track_points(:, 2), 'k-','linewidth',1)
    end
    
    %     if abs(track_points(end,1)) < 150 && abs(track_points(end,2))<1
    %         plot(track_points(:,1), track_points(:, 2), 'k-','linewidth',2)
    %     end
end
plot(points{end}(12,1),points{end}(12,2),'m*','markersize',40)
axis equal
set(gca,'fontsize',15)
grid on

axis([-7 5 4 15])

% axis([-40 10 -30 30])
% plot([-7 5],[4 4],'k-','linewidth',2)
% plot([-7 5],[15 15],'k-','linewidth',2)
% plot([-7 -7],[4 15],'k-','linewidth',2)
% plot([5 5],[4 15],'k-','linewidth',2)

plot(points{1}(:,1),points{1}(:,2),'bo','markersize',10)
[~,id] = min(abs(tauvec-1));
plot(points{id}(:,1),points{id}(:,2),'r^','markersize',10)
plot(points{end}(:,1),points{end}(:,2),'gs','markersize',10)


% print(gcf,'-dpng','waveEigs2D.png')
print(gcf,'-dpng','trackedWaveEigs2D.png')

%%
keyboard

%%
FinalTime = 5.0;

time = 0;

% p = exp(-10*(x.^2 + y.^2));
k = 1; % frequency of solution
W = (2*k-1)/2*pi;
p = cos(W*x).*cos(W*y);

u = zeros(Np, K); v = zeros(Np, K);


% Runge-Kutta residual storage
resu = zeros(Np,K); resv = zeros(Np,K); resp = zeros(Np,K);

% compute time step size
CN = (N+1)^2/2; % guessing...
%dt = 1/(CN*max(abs(sJ(:)))*max(abs(1./J(:))));
dt = .1/(CN*max(Fscale(:)));

% outer time step loop
tstep = 0;
% figure
ids = 1:Np*K;
while (time<FinalTime)
    if(time+dt>FinalTime), dt = FinalTime-time; end
    
    for INTRK = 1:5
        
        [rhsp, rhsu, rhsv] = WaveRHS(p,u,v,tau);
        rhs = P*[rhsp(:); rhsu(:); rhsv(:)];
        a = 4/50;
        rhsp = (a)*rhsp + (1-a)*reshape(rhs(ids),Np,K);
        rhsu = (a)*rhsu + (1-a)*reshape(rhs(ids + Np*K),Np,K);
        rhsv = (a)*rhsv + (1-a)*reshape(rhs(ids + 2*Np*K),Np,K);
        
        % initiate and increment Runge-Kutta residuals
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        
        % update fields
        u = u+rk4b(INTRK)*resu;
        v = v+rk4b(INTRK)*resv;
        p = p+rk4b(INTRK)*resp;
        
    end;
    
    if 1 && nargin==0 && mod(tstep,10)==0
        clf
        vv = Vp*p;
        color_line3(xp,yp,vv,vv,'.');
        axis([-1 1 -1 1 -1.25 1.25]);
        %         PlotField2D(N+1, x, y, p); view(2)
        
        drawnow
    end
    
    % Increment time
    time = time+dt; tstep = tstep+1;
end

function rhsu = AdvecRHS(u,tau)

Globals2D
global bx by Vq Pq Vfq Pfq Vrq Vsq Prq Psq

du = zeros(Nfp*Nfaces,K);
bn = bx(vmapM).*nx(:) + by(vmapM).*ny(:);
bnf = reshape((bn - tau*abs(bn(:))),Nfp*Nfaces,K);
ujump = reshape(u(vmapP)-u(vmapM),Nfp*Nfaces,K);

du(:) = .5*Pfq*((Vfq*ujump).*(Vfq*bnf)); % project ujump/bnf
uq = Vq*u;
fx = Pq*((Vq*bx).*uq);
fy = Pq*((Vq*by).*uq);

% fx = bx.*u; fy = by.*u;
% du(:) = 0.5*(ujump).*bnf;

ux = rx.*(Dr*fx) + sx.*(Ds*fx);
uy = ry.*(Dr*fy) + sy.*(Ds*fy);
rhsu = -((ux + uy) + LIFT*(Fscale.*(du)));

function [rhsp, rhsu, rhsv] = WaveRHS(p,u,v,tau)

% function [rhsu, rhsv, rhsp] = acousticsRHS2D(u,v,p)
% Purpose  : Evaluate RHS flux in 2D acoustics TM form

Globals2D;

% Define field differences at faces
dp = zeros(Nfp*Nfaces,K); dp(:) = p(vmapP)-p(vmapM);
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapP)-u(vmapM);
dv = zeros(Nfp*Nfaces,K); dv(:) = v(vmapP)-v(vmapM);

% Impose reflective boundary conditions (p+ = -p-)
du(mapB) = 0; dv(mapB) = 0;
dp(mapB) = -2*p(vmapB);

% evaluate upwind fluxes
ndotdU = nx.*du + ny.*dv;
fluxp =  tau*dp - ndotdU;

fluxu =  (tau*ndotdU - dp).*nx;
fluxv =  (tau*ndotdU - dp).*ny;

pr = Dr*p; ps = Ds*p;
dpdx = rx.*pr + sx.*ps;
dpdy = ry.*pr + sy.*ps;
divU = Dr*(u.*rx + v.*ry) + Ds*(u.*sx + v.*sy);

% compute right hand sides of the PDE's
rhsp =  -divU + LIFT*(Fscale.*fluxp)/2.0;
rhsu =  -dpdx + LIFT*(Fscale.*fluxu)/2.0;
rhsv =  -dpdy + LIFT*(Fscale.*fluxv)/2.0;
