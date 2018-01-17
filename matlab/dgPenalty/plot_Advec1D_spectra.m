function plot_Advec1D_spectra

% function [u] = Advec1D(u, FinalTime)
% Purpose  : Integrate 1D advection until FinalTime starting with
%            initial the condition, u
Globals1D;
N = 3;

% Generate simple mesh
K1D = 8;
[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

% Initialize solver and construct grid and metric
StartUp1D;

vmapP(1) = vmapM(end); % make periodic
vmapP(end) = vmapM(1);


M = inv(V*V');

rp = linspace(-1,1,100)';
Vp = Vandermonde1D(N,rp)/V;
xp =  Vp*x;


%% make CG projection 

R = zeros((N+1)*K,N*K);
offset = 0;
for e = 1:K
    rr = (1:N+1) + (e-1)*(N+1);
    cc = (1:N+1) + offset;
    R(rr,cc) = eye(N+1);
    offset = offset + N;
end
R(:,end) = []; 
R((N+1)*K,1) = 1; % make periodic

P = R*((R'*kron(diag(J(1,:)),M)*R) \ (R'*kron(diag(J(1,:)),M)));
% P = diag(1./sum(R,2))*R*R';
% return

%%
%M = kron(diag(J(1,:)),inv(V*V'));
T = kron(eye(K),V);
% M = blkdiag(M,M);
dhist = [];

dtau = 1e-2;
tmax = 8;
tauvec = 0:dtau:tmax;
% tauvec = tmax*tauvec.^2/(tmax)^2;
% tauvec = [100 ];
for ii = 1:length(tauvec)
    tau = tauvec(ii);
    
    Q = zeros(N+1,K);
    A = zeros((N+1)*K);
    for i = 1:(N+1)*K
        Q(i) = 1;
        rhs = AdvecRHS1D(Q,tau);
        Q(i) = 0;
        A(:,i) = rhs(:);
    end    
%     M = kron(diag(J(1,:)),inv(V*V'));
    [W D] = eig(A); 
    d = diag(D);    
        
%     keyboard
%     figure
%     for j = 1:size(W,2)
%         vv = Vp*reshape(W(:,j),N+1,K);
%         plot(xp,real(vv),'r')
%         hold on
%         plot(xp,imag(vv),'b')        
%     end
%     title(sprintf('\\tau = %f',tau))
            
    [~,p] = sort(imag(d),'descend');
    d = d(p); 
    W = W(:,p);    
    
    points{ii} = [real(d) imag(d)];
    modes{ii} = W;
    
    if mod(ii,50) == 0
        disp(sprintf('on eig %d out of %d\n',ii,length(tauvec)))
    end
end

%%
max_linking_distance = 10;
max_gap_closing = Inf;
debug = true;
disp('running simple tracker')
[ tracks, adjacency_tracks ] = simpletracker(points,...
    'MaxLinkingDistance', max_linking_distance, ...
    'MaxGapClosing', max_gap_closing, ...
    'Debug', debug);

figure
plot(points{1}(:,1),points{1}(:,2),'bo','markersize',10,'linewidth',2,'DisplayName','\tau = 0')
hold on
[~,id] = min(abs(tauvec-1));
plot(points{id}(:,1),points{id}(:,2),'r^','markersize',10,'DisplayName','\tau = 1.0')
[~,id] = min(abs(tauvec-4));
plot(points{id}(:,1),points{id}(:,2),'md','markersize',10,'DisplayName','\tau = 4.0')
plot(points{end}(:,1),points{end}(:,2),'gs','markersize',10,'DisplayName',sprintf('\\tau = %1.1f',tauvec(end)))

% % exact eigs
% lam_ex = pi/2*(1:Np*K); lam_ex = [lam_ex -lam_ex]; 
% plot(zeros(size(lam_ex)),lam_ex,'s','markersize',11,'DisplayName','Exact eigenvalues')

legend show
n_tracks = numel(tracks);
all_points = vertcat(points{:});
for i_track = 1 : n_tracks
   
    % We use the adjacency tracks to retrieve the points coordinates. It
    % saves us a loop.
    
    track = adjacency_tracks{i_track};
    track_points = all_points(track, :);
    
    plot(track_points(:,1), track_points(:, 2), 'k','linewidth',2)    
end
axis equal
set(gca,'fontsize',15)
grid on
%print(gcf,'-dpng','trackedEigs.png')
% keyboard

%% sort eigs by adjacency ids

eignum = 3; 
% eignum = 2*Np*K;
[~,id0] = min(abs(tauvec-0));
[~,id1] = min(abs(tauvec-1));
[~,id2] = min(abs(tauvec-10));
adj_ids = [id0 id1 id2];
sk = 1;
for i = adj_ids %1:2:length(adjacency_tracks{eignum})
    ii = adjacency_tracks{eignum}(i);
    cell_id = floor(ii/(Np*K))+1;
    eig_id = ii - (cell_id-1)*Np*K;        
    
    w{sk} = modes{cell_id}(:,eig_id); % take p mode
    points{cell_id}(eig_id,:)
    sk = sk + 1;
end

%%

c = modes{id1} \ w{1}; 
ids = find(abs(c)<1e-7);
% c(ids) = nan;
[pts, p] = sort(points{id1}(:,1));

% c = c/max(abs(c)); pts = pts/max(abs(pts));
% hold on
semilogy(abs(c(p)),'o','markersize',8,'linewidth',2);
hold on; 
semilogy(abs(pts),'x','markersize',8,'linewidth',2);
legend('c^0_j','-{\rm Re}(\lambda^1_j)','Interpreter','latex');
set(gca,'fontsize',15)
grid on
xlabel('Eigenvalue index','fontsize',15)
ylabel('Magnitude of expansion coefficient','fontsize',15)

ylim([1e-3 5e1])
print(gcf,'-dpng','ceigs0.png')


figure
c = modes{id1} \ w{3}; 
ids = find(abs(c)<1e-7);
% c(ids) = nan;
[pts, p] = sort(points{id1}(:,1));

% c = c/max(abs(c)); pts = pts/max(abs(pts));
% hold on
semilogy(abs(c(p)),'o','markersize',8,'linewidth',2);
hold on; 
semilogy(abs(pts),'x','markersize',8,'linewidth',2);
legend('c^{100}_j','-{\rm Re}(\lambda^1_j)','Interpreter','latex');
set(gca,'fontsize',15)
grid on
xlabel('Eigenvalue index','fontsize',15)
ylabel('Magnitude of expansion coefficient','fontsize',15)
ylim([1e-3 5e1])
print(gcf,'-dpng','ceigs100.png')

% set(h,'Interpreter','LaTeX')

%%
field = 0;
eignum = 3; 
% eignum = 2*Np*K;
[~,id0] = min(abs(tauvec-0));
[~,id1] = min(abs(tauvec-1));
[~,id2] = min(abs(tauvec-10));
adj_ids = [id0 id1 id2];
for i = adj_ids %1:2:length(adjacency_tracks{eignum})
    ii = adjacency_tracks{eignum}(i);
    mod_id = floor(ii/(Np*K))+1;
    eig_id = ii - (mod_id-1)*Np*K;        
    
    w = modes{mod_id}(:,eig_id); % take p mode
    
    w = P*w;
    
    vv = Vp*reshape(w,Np,K);
    [a,idv] = max(abs(vv));
    a = a.*sign(vv(idv));
    vv = vv*diag(1./a);        
        
    subplot(1,2,1)       
    plot(xp,real(vv),'b');
    hold on
    plot(xp,imag(vv),'r');
    hold off
    axis([-1 1 -2 2])
    
    subplot(1,2,2)
    hold on;
    % draw paths
    n_tracks = numel(tracks);
    all_points = vertcat(points{:});
    for i_track = 1 : n_tracks
        track = adjacency_tracks{i_track};
        track_points = all_points(track, :);
        plot(track_points(:,1), track_points(:, 2), 'k','linewidth',1)
    end
    xlim([-100 5]);    
    
    plot(points{end}(:,1),points{end}(:,2),'bx','markersize',10)
    plot(points{id}(:,1),points{id}(:,2),'rs')
    plot(points{1}(:,1),points{1}(:,2),'go')
    plot(all_points(ii,1),all_points(ii,2),'ro','markersize',24)
    
    %     xlim([-25 1])
    

    hold off
    title(sprintf('\\tau = %f',tauvec(i)))
    
    pause
end

%% for paper

field = 0;
eignum = 3; 
% eignum = 2*Np*K;

[~,id0] = min(abs(tauvec-0));
[~,id1] = min(abs(tauvec-1));
[~,id2] = min(abs(tauvec-10));
adj_ids = [id0 id1 id2];
sk = 1;
for i = adj_ids
    
    ii = adjacency_tracks{eignum}(i);
    mod_id = floor(ii/(Np*K))+1;
    eig_id = ii - (mod_id-1)*Np*K;
    w = modes{mod_id}(:,eig_id); % take p mode
    vv = Vp*reshape(w,Np,K);
    [a,idv] = max(abs(vv));
    a = a.*sign(vv(idv));
    vv = vv*diag(1./a);        
        
    figure(1)
    xpa = [xp;nan(1,K)]; xpa = xpa(:);
    vva = [vv;nan(1,K)]; vva = vva(:);
    plot(xpa,real(vva),'b-','linewidth',2,'DisplayName','Real part');
    hold on
    plot(xpa,imag(vva),'r--','linewidth',2,'DisplayName','Imaginary part');    
    hold off
    axis([-1 1 -2 2])
    legend show
    set(gca,'fontsize',15)
    grid on
    print(gcf,'-dpng',sprintf('tauModes%d.png',sk))

    figure(2)
    
%     plot(points{end}(:,1),points{end}(:,2),'x','markersize',10)
      plot(points{i}(:,1),points{i}(:,2),'bo','markersize',10,'linewidth',2)     
%       plot(all_points(ii,1),all_points(ii,2),'o','markersize',10,'linewidth',2)      
    hold on;
    plot(points{i}(eig_id,1),points{i}(eig_id,2),'r*','markersize',15,'linewidth',2)
            
    % draw paths
    n_tracks = numel(tracks);
    colors = hsv(n_tracks);
    all_points = vertcat(points{:});
    for i_track = 1 : n_tracks
        track = adjacency_tracks{i_track};
        track_points = all_points(track, :);        
        plot(track_points(:,1), track_points(:, 2), 'k','linewidth',2)
    end
    hold off
    axis equal
    set(gca,'fontsize',15)
%     xlim([-450 10]) % for diverge
        axis([-24.5 1 5 27])
    %     title(sprintf('\\tau = %f',tauvec(i)))
    grid on
    print(gcf,'-dpng',sprintf('tauEigs%d.png',sk))
    pause
    sk = sk + 1;
end


function [rhsu] = AdvecRHS1D(u,alpha)

% function [rhsu] = AdvecRHS1D(u,time)
% Purpose  : Evaluate RHS flux in 1D advection

Globals1D;

% form field differences at faces
%alpha = 0;
du = zeros(Nfp*Nfaces,K);
du(:) = (u(vmapM)-u(vmapP)).*(nx(:)-alpha*abs(nx(:)))/2;

% compute right hand sides of the semi-discrete PDE
rhsu = -rx.*(Dr*u) + LIFT*(Fscale.*(du));
return

