clear
Globals2D;

K1D = 1;
useQuadrature = 1;

x0 = .1; y0 = .1;
pex = @(x,y,t) exp(-5^2*((x-x0).^2 + (y-y0).^2));

k = 10;
pex = @(x,y,t) cos(k*pi*x/2).*cos(k*pi*y/2);

% [xp yp] = meshgrid(linspace(-1,1,100));
% vv = pex(xp,yp);
% color_line3(xp,yp,vv,vv,'.');return

% pex = @(x,y,t) sin(k*pi*x).*sin(k*pi*y);
% pex = @(x,y,t) exp(x.*y).*sin(k*pi*x).*sin(k*pi*y);
% pex = @(x,y,t) exp(x.*y);

a = 1/8;
a = .28;

% figure
K1D = 1;
NB = 4;
Kvec = 32;[4 8 16 32];
smoothKnots = 0;

clear L2err ndofs
sk = 1;
for Ksub = Kvec
    Ksub
    N = NB+Ksub-1;
    dofs = (N+1)^2*K1D^2;
    
    % Read in Mesh
    [Nv, VX, VY, K, EToV] = QuadMesh2D(K1D);
    StartUpQuad2D;
    
    % non-affine mappings
    if a > 1e-10
        x = x + a*cos(3/2*pi*y).*cos(pi/2*x);
        y = y + a*sin(3/2*pi*x).*cos(pi/2*y);
        %     x = x + a*(1-y).*(1+y).*(1-x).*(1+x);
        %     y = y + a*(1-x).*(1+x).*(1-y).*(1+y);        
        if 0
            e = 1;
            xK = reshape(x(:,e),N+1,N+1);
            yK = reshape(y(:,e),N+1,N+1);
            rp1D = linspace(-1,1,100)';
            Vp = Vandermonde1D(N,rp1D)/Vandermonde1D(N,JacobiGL(0,0,N));
            hold on
            for i = 1:N+1
                plot(Vp*xK(:,i),Vp*yK(:,i),'k-','linewidth',2);
            end
            for i = 1:N+1
                plot(Vp*xK(i,:)',Vp*yK(i,:)','k-','linewidth',2);
            end
            plot(xK(:),yK(:),'ko','markersize',10,'linewidth',2,'MarkerFaceColor',[.49 1 .63]);hold on
            axis tight
            axis off
            axis equal  
            return
        end
    end
    
    % isoparametric - interpolate geofacs
    [rx,sx,ry,sy,J] = GeometricFactorsQuad2D(x,y,Dr,Ds);
    
    [rq1D wq1D] = JacobiGQ(0,0,N+1);
    if (Ksub>1)
        [~, ~, ~, ~, ~, ~, ~, ~, VX] = bsplineVDM(NB,Ksub,rq1D,smoothKnots);
        [rgq wgq] = JacobiGQ(0,0,NB+1);
        
        h = @(r) repmat(diff(VX(:)'),length(r),1);
        map = @(r) reshape(h(r).*(repmat(r,1,Ksub)+1)/2 + repmat(VX(1:end-1),length(r),1),length(r)*Ksub,1);
        
        rq1D = map(rgq);
        wq1D = repmat(wgq,1,Ksub).*h(rgq)/2;
        rq1D = rq1D(:);
        wq1D = wq1D(:);
    end
    
    [r1D] = JacobiGL(0,0,N);
    
    [rq sq] = meshgrid(rq1D); rq = rq(:); sq = sq(:);
    [wrq wsq] = meshgrid(wq1D); wq = wrq(:).*wsq(:);
    V1D = Vandermonde1D(N,r1D);
    Vq1D = Vandermonde1D(N,rq1D)/V1D;
    Nq = length(wq1D);
    xq = zeros(Nq^2,K);
    yq = zeros(Nq^2,K);
    Jq = zeros(Nq^2,K);
    for e = 1:K
        xq(:,e) = reshape((Vq1D*reshape(x(:,e),N+1,N+1))*Vq1D',(Nq)^2,1);
        yq(:,e) = reshape((Vq1D*reshape(y(:,e),N+1,N+1))*Vq1D',(Nq)^2,1);
        Jq(:,e) = reshape((Vq1D*reshape(J(:,e),N+1,N+1))*Vq1D',(Nq)^2,1);
    end
    wJq = spdiag(wq)*Jq;
%     min(Jq(:))
    if min(Jq(:))<1e-10
        keyboard
    end
%     return

    
    
    [BVDM M1D] = bsplineVDM(NB,Ksub,r1D,smoothKnots); % VDM for interp, mass, M\S
    %             Vp1D = bsplineVDM(NB,Ksub,rp1D,smoothKnots);
    Vq1D = bsplineVDM(NB,Ksub,rq1D,smoothKnots);
    invM1D = inv(M1D);
    Pq1D = invM1D*(Vq1D'*diag(wq1D));
    
    % L2 projection
    Vq = kron(sparse(Vq1D),sparse(Vq1D));
    M = Vq'*spdiag(wJq)*Vq;
    M = sparse(M);
    pproj = reshape(M\(Vq'*(wJq.*pex(xq,yq))),N+1,N+1);
    
    % wadg
    xq = reshape(xq,Nq,Nq);
    yq = reshape(yq,Nq,Nq);
    Jq = reshape(Jq,Nq,Nq);
    wJq = reshape(wJq,Nq,Nq);
    pq = pex(xq,yq);
    pwadg = Pq1D*(((Vq1D*Pq1D)*(pq.*Jq)*(Vq1D*Pq1D)')./Jq)*Pq1D';
    
    L2err_proj(sk) = sqrt(sum(sum((pq-Vq1D*pproj*Vq1D').^2.*wJq)));
    L2err_wadg(sk) = sqrt(sum(sum((pq-Vq1D*pwadg*Vq1D').^2.*wJq)));
    difference(sk) = sqrt(sum(sum((Vq1D*pwadg*Vq1D'-Vq1D*pproj*Vq1D').^2.*wJq)));
    
    hh(sk) = 2/Ksub;
    ndofs(sk) = dofs;
    sk = sk + 1;
end

loglog(hh,L2err_proj,'o--','linewidth',2,'DisplayName',sprintf('L2 projection'))
hold on
loglog(hh,L2err_wadg,'x--','linewidth',2,'DisplayName',sprintf('WADG'))
loglog(hh,difference,'s--','linewidth',2,'DisplayName',sprintf('WADG'))
print_pgf_coordinates(hh,L2err_proj)
print_pgf_coordinates(hh,L2err_wadg)
print_pgf_coordinates(hh,difference)

grid on
legend show