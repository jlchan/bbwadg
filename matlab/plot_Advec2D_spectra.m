function Advec2D_cfl

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

% rebuild maps for periodic
BuildPeriodicMaps2D(2,2);

[re se] = EquiNodes2D(N); [re se] = xytors(re,se);
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

% % Set initial conditions
cx = .1; cy = .1;
D = (x-cx).^2 + (y-cy).^2;
u = exp(-D*5^2).*(1-x.^2).*(1-y.^2);
u = sin(pi*x).*sin(pi*y);

bx = ones(size(x));
by = ones(size(x))*0;
% by = sin(pi*x).*sin(pi*y);

% plot3(x,y,-pi*cos(pi*x).*sin(pi*y),'o');return

% divB = rx.*(Dr*bx) + sx.*(Ds*bx) + ry.*(Dr*by) + sy.*(Ds*by);

FinalTime = 2;
%%

R = getCGRestriction()';
P = R*((R'*kron(diag(J(1,:)),M)*R)\(R'*kron(diag(J(1,:)),M)));

%% eigs

rLGL = JacobiGQ(0,0,N); rmin = abs(rLGL(1)-rLGL(2));
dtscale = dtscale2D; dt = min(dtscale)*rmin*2/3

e = zeros(Np*K,1);
A = zeros(Np*K);
tau = 0;
for i = 1:Np*K
    e(i) = 1;
    ids = 1:Np*K;
    u = reshape(e(ids),Np,K);
    
    % dudt + dudx= 0
    rhsu = AdvecRHS(u,tau);
    A(:,i) = rhsu(:);
    
    e(i) = 0;
    if mod(i,100)==0
        disp(sprintf('on column %d out of %d\n',i,Np*K))
    end
end
M = kron(diag(J(1,:)),inv(V*V'));

S = zeros(Np*K);
tau = 1;
for i = 1:Np*K
    e(i) = 1;
    ids = 1:Np*K;
    u = reshape(e(ids),Np,K);
    
    % dudt + dudx= 0
    rhsu = AdvecRHS(u,tau);
    S(:,i) = rhsu(:);
    
    e(i) = 0;
    if mod(i,100)==0
        disp(sprintf('on column %d out of %d\n',i,Np*K))
    end
end

S = S-A; S(abs(S)<1e-8) = 0;
keyboard
%% plot eigenmodes

if 1
    load eigTrackAdvec
    [~,id0] = min(abs(tauvec-0));
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
        track = adjacency_tracks{17};
        track_points = all_points(track, :);
        
        figure(1)
        plot(track_points(end,1), track_points(end, 2), 'bo','markersize',14,'linewidth',2)
        hold on
        plot(track_points(:,1), track_points(:, 2), 'k-','linewidth',1)
        plot(all_points(track(tauids(i)),1), all_points(track(tauids(i)), 2), 'bs','markersize',14,'linewidth',1)
        pause
        
        if 1
            len = size(points{1},1);
            i1 = mod(track(tauids(i)),len);
            d(i1)
                        
            w = real(W(:,i1));
            v = reshape(w(1:Np*K),Np,K);
            [a,p] = max(abs(v(:)));
            v = v*sign(v(p))/a;
            vv = Vp*v;
            
            figure
            PlotField2D(25,x,y,v);
%             color_line3(xp,yp,vv,vv,'.');
            view(3)
            
            set(gca,'fontsize',15)
            print(gcf,'-dpng',sprintf('advecEig2D%d.png',i))
        end
    end
    keyboard
end


%% plot eigenvalues
if 1
    tmax = 50;
    tauvec = [0:1e-2:tmax];
    tauvec = tmax*tauvec.^2/(tmax)^2;
    
    for ii = 1:length(tauvec)
        tau = tauvec(ii);
        if mod(ii,50)==0
            disp(sprintf('computing mode %d out of %d\n',ii,length(tauvec)))
        end
        
        [W D] = eig(A + tau*S);
        d = diag(D);
        [~,p] = sort(imag(d),'descend');
        d = d(p); W = W(:,p);
        
        points{ii} = [real(d) imag(d)];
        %     modes{ii} = W;
        
    end
    
    max_linking_distance = 10;
    max_gap_closing = Inf;
    debug = true;
    disp('running simple tracker')
    [ tracks, adjacency_tracks ] = simpletracker(points,...
        'MaxLinkingDistance', max_linking_distance, ...
        'MaxGapClosing', max_gap_closing, ...
        'Debug', debug);
    
%     save eigTrackAdvec tracks adjacency_tracks tauvec points
%     keyboard
else
%     load eigTrackAdvec
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
    plot(track_points(:,1), track_points(:, 2), 'k-','linewidth',2)

        
end
axis equal
set(gca,'fontsize',15)
grid on

plot(points{1}(:,1),points{1}(:,2),'bo','markersize',10,'DisplayName','\tau = 0')
[~,id] = min(abs(tauvec-1));
plot(points{id}(:,1),points{id}(:,2),'r^','markersize',10,'DisplayName','\tau = 1.0')
plot(points{end}(:,1),points{end}(:,2),'gs','markersize',10,'DisplayName',sprintf('\\tau = %1.1f',tauvec(end)))


axis([-55 5 -20 20])
print(gcf,'-dpng','trackedAdvecEigs.png')

%%
keyboard


%%
time = 0;

% Runge-Kutta residual storage
resu = zeros(Np,K);

% compute time step size
rLGL = JacobiGQ(0,0,N); rmin = abs(rLGL(1)-rLGL(2));
dtscale = dtscale2D; dt = min(dtscale)*rmin*2/3

% outer time step loop
while (time<FinalTime)
    
    if(time+dt>FinalTime), dt = FinalTime-time; end
    
    for INTRK = 1:5
        
        rhsu = AdvecRHS(u,1);
        
        % initiate and increment Runge-Kutta residuals
        resu = rk4a(INTRK)*resu + dt*rhsu;
        
        % update fields
        u = u+rk4b(INTRK)*resu;
    end
    
    clf;
    up = Vp*u;
    color_line3(xp,yp,up,up,'.');drawnow
    
    % Increment time
    time = time+dt;
end
return

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

