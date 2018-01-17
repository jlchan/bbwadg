clear
Globals1D;

f = @(x) (x < -pi/8) + (x > pi/8).*exp(x).*(sin(1+pi*x));
w = @(x) 2 + sin(pi*x);
Kvec = [4 8 16 32 64 128];
h = 2./Kvec;
sk = 1;
for K1D = Kvec
    N = 1;
    
    r = JacobiGL(0,0,N);    
    
    [rq wq] = JacobiGQ(0,0,N+2);
    [rq2 wq2] = JacobiGQ(0,0,N+4);
            
    [Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);
    
    StartUp1D;
    
    V = Vandermonde1D(N,r);
    Vq = Vandermonde1D(N,rq)/V;
    M = Vq'*diag(wq)*Vq;
    
    Pq = M\(Vq'*diag(wq));
    xq =  Vandermonde1D(N,rq)/Vandermonde1D(N,JacobiGL(0,0,N))*x;
    
    Vq2 = Vandermonde1D(N,rq2)/V;
    wJq = diag(wq2)*(Vq2*J);
    
    for e = 1:K
        fp(:,e) = (Vq'*diag(wq.*w(xq(:,e)))*Vq) \ (Vq'*diag(wq.*w(xq(:,e)))*f(xq(:,e)));
    end
    err = (Vq2*(fp - Pq*((Vq*Pq*(w(xq).*f(xq)))./w(xq)))).^2;   
    L2err(sk) = sqrt(sum(sum(wJq(:,1).*err(:,1))));
    sk = sk + 1;
end

loglog(h,L2err,'x--')
hold on
loglog(h,2*L2err(1)/h(1).^(N+2)*h.^(N+2),'--')

%%
clear
Globals2D

f = @(x,y) (1 + (x+y > sin(pi*x))).*exp(.5*(x+y)).*sin(1+pi*x).*sin(1+pi*y);

plotMesh = 0;
Kvec = [4 8 16 32 64];
h = 2./Kvec;

sk = 1;
for K1D = Kvec
    
    
    N = 5;
    a = 1/4; % curv warping
    
    
    [Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
    StartUp2D;    
    
    Nq = 2*N+1;
    
    [rq sq wq] = Cubature2D(Nq); % integrate u*v*c
    Vq = Vandermonde2D(N,rq,sq)/V;
    M = Vq'*diag(wq)*Vq;
    Pq = M\(Vq'*diag(wq)); % J's cancel out
    xq = Vq*x; yq = Vq*y;        
        
    Nq2 = Nq+4;
    [rq2 sq2 wq2] = Cubature2D(Nq2); % integrate u*v*c
    Vq2 = Vandermonde2D(N,rq2,sq2)/V;
    Pq = (Vq2'*diag(wq2)*Vq2)\(Vq'*diag(wq)); % J's cancel out
    
    % make curvilinear mesh
    
    x0 = 0; y0 = 0;
    x = x + a*cos(1/2*pi*x).*cos(3/2*pi*y);
    y = y + a*cos(3/2*pi*x).*cos(1/2*pi*y);
    
    xq = Vq*x; yq = Vq*y;
    xq2 = Vq2*x; yq2 = Vq2*y;
        
    
    Jq = zeros(length(rq),K);
    Jq2 = zeros(length(rq2),K);
    for e = 1:K
        [rxk,sxk,ryk,syk,Jk] = GeometricFactors2D(x(:,e),y(:,e),Vq*Dr,Vq*Ds);
        Jq(:,e) = Jk;
        
        [rxk,sxk,ryk,syk,Jk] = GeometricFactors2D(x(:,e),y(:,e),Vq2*Dr,Vq2*Ds);
        Jq2(:,e) = Jk;
    end
    
    if plotMesh
        rp1D = linspace(-1,1,100)';
        Vp1D = Vandermonde1D(N,rp1D)/Vandermonde1D(N,JacobiGL(0,0,N));
        Vfp = kron(eye(Nfaces),Vp1D);
        xfp = Vfp*x(Fmask(:),:);
        yfp = Vfp*y(Fmask(:),:);
        plot(xfp,yfp,'k.')
        axis off
        axis equal
        title(sprintf('min J = %f\n',min(Jq(:))));
        return
    end
    
    wJq = diag(wq)*Jq;
    wJq2 = diag(wq2)*Jq2;
    
    fproj = zeros(Np,K);
    for e = 1:K
        fproj(:,e) = (Vq'*diag(wJq(:,e))*Vq)\(Vq'*(f(xq(:,e),yq(:,e)).*wJq(:,e)));
    end
    err1 = Vq2*fproj - f(xq2,yq2);
    L2err_proj(sk) = sqrt(sum(sum(wJq2.*err1.^2)));
    
    fwadg = Pq*((Vq*Pq*(f(xq,yq).*Jq))./Jq);   
    err2 = Vq2*fwadg - f(xq2,yq2);        
    L2err_wadg(sk) = sqrt(sum(sum(wJq2.*err2.^2)));
    
    err3 = Vq2*(fproj-fwadg);
    L2err_diff(sk) = sqrt(sum(sum(wJq2.*err3.^2)));
%     err3 = Vq*(fproj-fwadg);
%     L2err_diff(sk) = sqrt(sum(sum(wJq.*err3.^2)));
    
%     vdiff = Vq*(fproj-fwadg);
%     vdiff = vdiff - Vq*Pq*((Vq*Pq*(vdiff.*Jq))./Jq);
%     ip(sk) = sqrt(abs(sum(sum(f(xq,yq).*wJq.*vdiff))));
    
    sk = sk + 1;
    
end

print_pgf_coordinates(h,L2err_proj)
print_pgf_coordinates(h,L2err_wadg)
print_pgf_coordinates(h,L2err_diff)

%
loglog(h,L2err_proj,'bo--','linewidth',2,'markersize',16,'MarkerFaceColor',[.49 1 .63])
hold on
loglog(h,L2err_wadg,'rx--','linewidth',2,'markersize',16)
loglog(h,L2err_diff,'k^--','linewidth',2,'markersize',16,'MarkerFaceColor',[.49 1 .63])
hold on
loglog(h,h.^(N+2),'b-','linewidth',2)
loglog(h,h.^(2),'b-','linewidth',2)
grid on
% 'MarkerFaceColor',[.49 1 .63]
