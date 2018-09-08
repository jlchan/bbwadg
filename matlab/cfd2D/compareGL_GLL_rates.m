clear
clear -globals
Globals2D

N = 4;
a = .25; % warping factor

plotMesh = 0;

Kvec = 2.^(1:5);

max_lam = zeros(length(Kvec),1);
min_lam = ones(length(Kvec),1);
sk = 1;
for K1D = Kvec
    [Nv, VX, VY, K, EToV] = QuadMesh2D(K1D,K1D);
    
    StartUp2D;
    
    % vol nodes
    [rq1D wq1D] = JacobiGL(0,0,N);
    [rq sq] = meshgrid(rq1D);
    rq = rq(:); sq = sq(:);
    [wr ws] = meshgrid(wq1D);
    wq = wr(:).*ws(:);
    
    [rq1D wq1D] = JacobiGQ(0,0,N+1);
    [rq2 sq2] = meshgrid(rq1D);
    rq2 = rq2(:); sq2 = sq2(:);
    [wr ws] = meshgrid(wq1D);
    wq2 = wr(:).*ws(:);
    
    Vq = Vandermonde2D(N,rq,sq)/V;    
    Vq2 = Vandermonde2D(N,rq2,sq2)/V;
    
    M = Vq'*diag(wq)*Vq;
    Pq = M\(Vq'*diag(wq)); % J's cancel out
        
    % make curvilinear mesh
    Lx = 1; Ly = 1;
    x0 = 0; y0 = 0;
    dx = cos(1/2*pi*(x-x0)/Lx).*cos(3/2*pi*(y-y0)/Ly);
    x = x + Lx*a*dx;
    dy = cos(3/2*pi*(x-x0)/Lx).*cos(1/2*pi*(y-y0)/Ly);
    y = y + Ly*a*dy;
    
    xq = Vq*x; yq = Vq*y;
    
    [rx sx ry sy J] = GeometricFactors2D(x,y,Vq*Dr,Vq*Ds);
    J = Vq*J;
    J2 = Vq2*J;
    
    for e = 1:K        
        MGLL = Vq'*diag(wq.*J(:,e))*Vq;        
        MGL = Vq2'*diag(wq2.*J2(:,e))*Vq2;        
        
        lamK = eig(MGL,MGLL);
        max_lam(sk) = max(max_lam(sk),max(lamK));
        min_lam(sk) = max(min_lam(sk),min(lamK));        
    end        
    
    sk = sk + 1;
    
    if plotMesh
        rp1D = linspace(-1,1,100)';
        Vp1D = Vandermonde1D(N,rp1D)/Vandermonde1D(N,JacobiGL(0,0,N));
        Vfp = kron(eye(Nfaces),Vp1D);
        xfp = Vfp*x(Fmask(:),:);
        yfp = Vfp*y(Fmask(:),:);
        plot(xfp,yfp,'k.')
        hold on
        %plot(x,y,'b.','markersize',14)
        for e = 1:K
            xe = reshape(x(:,e),N+1,N+1);
            ye = reshape(y(:,e),N+1,N+1);
            for i = 1:N+1
                plot(Vp1D*xe(:,i),Vp1D*ye(:,i),'k-');
                plot(Vp1D*(xe(i,:)'),Vp1D*(ye(i,:)'),'k-');
            end
        end
        axis equal
        L2err = nan;
        axis off
        return
    end
end

plot(min_lam,'o--')
hold on
plot(max_lam,'x--')