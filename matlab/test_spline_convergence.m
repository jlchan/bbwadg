clear
Globals1D;

smoothKnots = 20;

NB = 3;
Ksub = 1;
K1D = 1;

uex = @(x) exp(x).*sin(1+2*pi*x);
k = 5;
uex = @(x) sin(.1+k*pi*x); 
% uex = @(x) exp(.5*x);

for i = 2:8
    
    NB = i;
    Ksub = 2*NB;
    
    N = NB+Ksub-1; % number of sub-mesh splines
    
    % Generate simple mesh
    [Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);
    
    ndofs = K*(N+1); % total # dofs
    
    % Initialize solver and construct grid and metric
    StartUp1D;
    
    rp = linspace(-1,1,50)';
    Vp = Vandermonde1D(N,rp)/V;
    xp = Vp*x;
    
    % make splines
    rB = JacobiGL(0,0,N);
    xB = Vandermonde1D(N,rB)/V * x;            
    
    % [rq wq] = JacobiGQ(0,0,N);
    [rgq wgq] = JacobiGQ(0,0,2*N+1);    
    if Ksub> 1
        [~, ~, ~, ~, ~, ~, ~, ~, VX] = bsplineVDM(NB,Ksub,rB,smoothKnots); % VDM for interp, mass, M\S
        h = @(r) repmat(diff(VX(:)'),length(r),1);
        map = @(r) reshape(h(r).*(repmat(r,1,Ksub)+1)/2 + repmat(VX(1:end-1),length(r),1),length(r)*Ksub,1);
        
        rq = map(rgq);
        wq = repmat(wgq,1,Ksub).*h(rgq)/2;
        rq = rq(:); wq = wq(:);
%         keyboard
    else
        rq = rgq(:); wq = wgq(:);
    end
            
    
    xq = Vandermonde1D(N,rq)/V * x;
    
    h = 2/(K1D);
    hB = h/Ksub;
    
    % spline projection
    [Bq] = bsplineVDM(NB,Ksub,rq,smoothKnots); % VDM for interp, mass, M\S
    M = Bq'*diag(wq)*Bq;
    PBq = M\(Bq'*diag(wq));
    
    % polynomial projection
    Vq = Vandermonde1D(NB,rq)/Vandermonde1D(NB,JacobiGL(0,0,NB));
    Pq = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq));
    
%     [Vp] = bsplineVDM(NB,Ksub,rp,smoothKnots);
    
    % set fxn to approximate
    
%     uex = @(x) exp(-(5*sin(1+x)).^2);
%     uex = @(x) tanh(2*5*x);
%     plot(xp,uex(xp));return
    uB = PBq*uex(xq); % initialize w/lsq fit
    u = Pq*uex(xq); % initialize w/lsq fit
    
    err = Bq*uB-uex(xq);
    err = diag(wq)*(err.^2).*J(1);
    L2errB(i) = sqrt(sum(err(:)));
    
    err = Vq*u-uex(xq);
    err = diag(wq)*(err.^2).*J(1);
    
    L2errP(i) = sqrt(sum(err(:)));
    hvec(i) = h;
    hBvec(i) = hB;
    dofs(i) = (NB+1)*K1D;
    dofsB(i) = (N+1)*K1D;
end


% loglog(hvec,L2errP,'o--','markersize',8,'linewidth',2)
% hold on;
% loglog(hBvec,L2errB,'x--','markersize',8,'linewidth',2)

loglog(dofsB,L2errB,'x--','markersize',8,'linewidth',2)
hold on;
% print_pgf_coordinates(dofsB,L2errB)
% loglog(dofs,L2errP,'o--','markersize',8,'linewidth',2)

% legend('Polynomials','Local splines')


