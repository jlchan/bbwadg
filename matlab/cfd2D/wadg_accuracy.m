clear
Globals1D;

f = @(x) (x < -pi/8) + (x > pi/8).*exp(x).*(sin(1+pi*x));
f = @(x) 0*(x > .1) + exp(x).*sin(1+pi*x);
w = @(x) 2 + sin(pi*x);
Kvec = [2 4 8 16 32 64 128 256];
h = 2./Kvec;
sk = 1;
for K1D = Kvec
    N = 4;
    
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
    xq2 = Vq2*x;
    wJq = diag(wq2)*(Vq2*J);
    
    for e = 1:K
        fp(:,e) = (Vq'*diag(wq.*w(xq(:,e)))*Vq) \ (Vq'*diag(wq.*w(xq(:,e)))*f(xq(:,e)));
%         fw(:,e) = M\(Vq'*diag(wq./w(xq(:,e)))*Vq)*Pq*(w(xq(:,e)).*f(xq(:,e))
    end
    err = (Vq2*(fp - Pq*((Vq*Pq*(w(xq).*f(xq)))./w(xq)) )).^2;
%     err = (Vq2*(fp - fw )).^2;
    L2err(sk) = sqrt(sum(sum(wJq.*err)));    
    
%     u = Pq*exp(xq);
%     for e = 1:K
%         Tinvwu(:,e) = (Vq'*diag(wq./w(xq(:,e)))*Vq)\(M*u(:,e));
%     end
%     err = (Vq2*u - Vq2*(Tinvwu)./w(xq2)).^2;
%     L2err(sk) = sqrt(sum(sum(wJq.*err)));
    
    sk = sk + 1;
end

loglog(h,L2err,'bx--')
hold on
r = N+2;
loglog(h,5*L2err(1)/h(1).^r*h.^r,'r--')

%%

clear
Globals1D;

f = @(x) (x < -pi/8) + (x > pi/8).*exp(x).*(sin(1+pi*x));
f = @(x) 0*(x > .1) + exp(x).*sin(1+pi*x);
w = @(x) 2 + sin(pi*x);
Kvec = [2 4 8 16 32 64 128 256];
h = 2./Kvec;
sk = 1;
for K1D = Kvec
    N = 4;
    
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
    xq2 = Vq2*x;
    wJq = diag(wq2)*(Vq2*J);
    
    for e = 1:K
        fp(:,e) = (Vq'*diag(wq.*w(xq(:,e)))*Vq) \ (Vq'*diag(wq)*f(xq(:,e)));
%         fw(:,e) = M\(Vq'*diag(wq./w(xq(:,e)))*Vq)*Pq*(w(xq(:,e)).*f(xq(:,e))
    end
    err = (Vq2*(fp - Pq*((Vq*Pq*(f(xq)))./w(xq)) )).^2;
%     err = (Vq2*(fp - fw )).^2;
    L2err(sk) = sqrt(sum(sum(wJq.*err)));    
    
%     u = Pq*exp(xq);
%     for e = 1:K
%         Tinvwu(:,e) = (Vq'*diag(wq./w(xq(:,e)))*Vq)\(M*u(:,e));
%     end
%     err = (Vq2*u - Vq2*(Tinvwu)./w(xq2)).^2;
%     L2err(sk) = sqrt(sum(sum(wJq.*err)));
    
    sk = sk + 1;
end

loglog(h,L2err,'bx--')
hold on
r = N+2;
loglog(h,5*L2err(1)/h(1).^r*h.^r,'r--')

%%
clear
Globals2D

f = @(x,y) exp((x+y)).*sin(pi*x).*sin(pi*y) + 0*(x+y > sin(pi*x));
% f = @(x,y) 1.0*((x+y > sin(pi*x))) + ;
% [xp yp] = meshgrid(linspace(-1,1,100));
% F = f(xp,yp);
% color_line3(xp(:),yp(:),F(:),F(:),'.')
% return

plotMesh = 0;
Kvec = [4 8 16 32];% 64];
h = 2./Kvec;

sk = 1;
for K1D = Kvec
    
    
    N = 4;
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
    rxJ = zeros(length(rq),K);
    Jq2 = zeros(length(rq2),K);
    for e = 1:K
        [rxk,sxk,ryk,syk,Jk] = GeometricFactors2D(x(:,e),y(:,e),Vq*Dr,Vq*Ds);
        rxJ(:,e) = rxk.*Jk;
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
%     wJq = diag(wq)*ones(size(Jq));
%     wJq2 = diag(wq2)*ones(size(Jq2));
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
    
    fq = Vq*Pq*(exp(xq+yq).*sin(pi*xq).*sin(pi*yq));
    fq = exp(xq+yq).*sin(pi*xq).*sin(pi*yq);    
    for e = 1:K
        %TinvJfwadg(:,e) = (Vq'*diag(wq./Jq(:,e))*Vq)\(M*fwadg(:,e));
        TinvJfwadg(:,e) = (Vq'*diag(wq./Jq(:,e))*Vq)\(M*Pq*fq(:,e));
    end    
    %err3 = (Vq*TinvJfwadg - Jq.*(Vq*fwadg));
    err3 = (Vq*TinvJfwadg - Jq.*fq);
    L2err_diff(sk) = sqrt(abs(sum(sum(diag(wq)*err3.^2)))); % N+2 rate
    
    L2err_diff(sk) = .5*wq'*(Jq.*fq-Jq.*(Vq*fwadg)) / (.5*wq'*Jq); % N+2 rate
    %     L2err_diff(sk) = sqrt(abs(sum(sum(wJq.*err3.^2)))); % N+2 rate
    %     L2err_diff(sk) = sqrt(abs(sum(sum(diag(wq)*err3.*(Vq*(fproj-fwadg))))));
    
    %     err3 = Vq*(fproj-fwadg);
    %     L2err_diff(sk) = sqrt(sum(sum(wJq.*err3.^2)));
    
    %     vdiff = Vq*(fproj-fwadg);
    %     vdiff = vdiff - Vq*Pq*((Vq*Pq*(vdiff.*Jq))./Jq);
    %     ip(sk) = sqrt(abs(sum(sum(f(xq,yq).*wJq.*vdiff))));
    
    geo(sk) = sqrt(sum(sum(wJq.*rxJ.^2)));     
    
    sk = sk + 1;
    
end

print_pgf_coordinates(h,L2err_proj)
print_pgf_coordinates(h,L2err_wadg)
print_pgf_coordinates(h,L2err_diff)

% loglog(h,geo,'o--');hold on;loglog(h,h.^2,'--');return

%
loglog(h,L2err_proj,'bo--','linewidth',2,'markersize',16,'MarkerFaceColor',[.49 1 .63])
hold on
loglog(h,L2err_wadg,'rx--','linewidth',2,'markersize',16)
loglog(h,L2err_diff,'k^--','linewidth',2,'markersize',16,'MarkerFaceColor',[.49 1 .63])
hold on
loglog(h,.0025*h.^(2*N+1),'b-','linewidth',2)
% loglog(h,.025*h.^(N+1/2),'b-','linewidth',2)
grid on
% 'MarkerFaceColor',[.49 1 .63]

%% wadg weighted convergence (not curved)

clear
Globals2D

k = 1;
f = @(x,y) exp((x+y)).*sin(k*pi*x).*sin(k*pi*y) + 0*(x+y > sin(pi*x));

Kvec = [4 8 16 32 64];
h = 2./Kvec;

sk = 1;
for K1D = Kvec
    
    
    N = 3;
        
    [Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);
    StartUp2D;    
    
    Nq = 2*N+1;
    
    [rq sq wq] = Cubature2D(Nq); 
    Vq = Vandermonde2D(N,rq,sq)/V;
    M = Vq'*diag(wq)*Vq;
    Pq = M\(Vq'*diag(wq)); 
    xq = Vq*x; yq = Vq*y;        
    
    w = 2 + exp(xq+yq).*sin(pi*xq).*sin(pi*yq);
    
    Jq = Vq*J;
    wJq = diag(wq)*Jq;
    
    
    fproj = zeros(Np,K);
    for e = 1:K
        fproj(:,e) = (Vq'*diag(wJq(:,e).*w(:,e))*Vq) \ (Vq'*(f(xq(:,e),yq(:,e)).*wJq(:,e).*w(:,e)));
    end
    fproj = f(x,y);
    
    err1 = Vq*fproj - f(xq,yq);
    L2err_proj(sk) = sqrt(sum(sum(wJq.*err1.^2)));
    
    fwadg = Pq*((Vq*Pq*(f(xq,yq).*w))./w);   
    fwadg = Pq*((Vq*Pq*((Vq*fproj).*w))./w);   
    err2 = Vq*fwadg - f(xq,yq);        
    L2err_wadg(sk) = sqrt(sum(sum(wJq.*err2.^2)));
    
    err3 = Vq*(fproj-fwadg);
%     err3 = f(xq,yq)-Vq*fwadg;

    L2err_diff(sk) = sqrt(sum(sum(wJq.*err3.^2)));    
    
    sk = sk + 1;
    
end

%
loglog(h,L2err_diff,'k^--','linewidth',2,'markersize',16,'MarkerFaceColor',[.49 1 .63])
hold on
loglog(h,1*h.^(N+2),'b-','linewidth',2)
grid on


%% 3D

clear
Globals3D
N = 3;

f = @(x,y,z) exp(x+y+z).*sin(pi*x).*sin(pi*y).*sin(pi*z);

for kk = 1:4
    filename = ['Grid/cube' num2str(2+kk) '.msh']
    [Nv, VX, VY, VZ, K, EToV] = MeshReaderGmsh3D(filename);
    VX = VX*2; VY = VY*2; VZ = VZ*2; % map -1,1
    StartUp3D;
    
    % make curvilinear
    a = .1;
    x = x + a*cos(.5*pi*x).*cos(.5*pi*y).*cos(.5*pi*z);
    y = y + a*cos(.5*pi*x).*cos(.5*pi*y).*cos(.5*pi*z);
    z = z + a*cos(.5*pi*x).*cos(.5*pi*y).*cos(.5*pi*z);
    
    % a = .05;
    % x = r - a*(r+s).^N;
    % y = s + a/2*(s+t).^N;
    % z = t - a/3*(r+t).^N;
    [rp sp tp] = EquiNodes3D(25);
    Vp = Vandermonde3D(N,rp,sp,tp)/V;
    % plot3(Vp*x,Vp*y,Vp*z,'.')
    
    [rq sq tq wq] = tet_cubature(2*N+1);
    Vq = Vandermonde3D(N,rq,sq,tq)/V;
    M =(Vq'*diag(wq)*Vq);
    Pq = M\(Vq'*diag(wq));
    
    xq = Vq*x;
    yq = Vq*y;
    zq = Vq*z;
    
    Drq = Vq*Dr;
    Dsq = Vq*Ds;
    Dtq = Vq*Dt;
    
    Nfp = (N+1)*(N+2)/2;
    [rqtri sqtri wqtri] = Cubature2D(2*N);
    Vfqf = Vandermonde2D(N,rqtri,sqtri)/Vandermonde2D(N,r(Fmask(:,1)),s(Fmask(:,1)));
    Vfqf = kron(eye(4),Vfqf);
    rfq = Vfqf*r(Fmask(:));
    sfq = Vfqf*s(Fmask(:));
    tfq = Vfqf*t(Fmask(:));
    wfq = repmat(wqtri,4,1);
    Vfq = Vandermonde3D(N,rfq,sfq,tfq)/V;
    Lq = M\(Vfq'*diag(wfq));
    
    Drfq = Vfq*Dr;
    Dsfq = Vfq*Ds;
    Dtfq = Vfq*Dt;
    
    Nfqf = length(rqtri);
    e = ones(Nfqf,1); zz = zeros(Nfqf,1);
    nrJ = [zz;zz;e;-e];
    nsJ = [zz;-e;e;zz];
    ntJ = [-e;zz;e;zz];           
    
    % original cross product form
    xrq = Drq*x; xsq = Dsq*x; xtq = Dtq*x;
    yrq = Drq*y; ysq = Dsq*y; ytq = Dtq*y;
    zrq = Drq*z; zsq = Dsq*z; ztq = Dtq*z;
    Jq = xrq.*(ysq.*ztq-zsq.*ytq) - yrq.*(xsq.*ztq-zsq.*xtq) + zrq.*(xsq.*ytq-ysq.*xtq);

    rxJ =  (ysq.*ztq - zsq.*ytq);
    sxJ = -(yrq.*ztq - zrq.*ytq);
    txJ =  (yrq.*zsq - zrq.*ysq);
    
    ryJ = -(xsq.*ztq - zsq.*xtq);
    syJ =  (xrq.*ztq - zrq.*xtq);
    tyJ = -(xrq.*zsq - zrq.*xsq);
    
    rzJ = (xsq.*ytq - ysq.*xtq);
    szJ = -(xrq.*ytq - yrq.*xtq);
    tzJ = (xrq.*ysq - yrq.*xsq);
    
    % make normals
    xrq = Drfq*x; xsq = Dsfq*x; xtq = Dtfq*x;
    yrq = Drfq*y; ysq = Dsfq*y; ytq = Dtfq*y;
    zrq = Drfq*z; zsq = Dsfq*z; ztq = Dtfq*z;
    
    Jf = xrq.*(ysq.*ztq-zsq.*ytq) - yrq.*(xsq.*ztq-zsq.*xtq) + zrq.*(xsq.*ytq-ysq.*xtq);
    
    rxJf =  (ysq.*ztq - zsq.*ytq);
    sxJf = -(yrq.*ztq - zrq.*ytq);
    txJf =  (yrq.*zsq - zrq.*ysq);
    
    ryJf = -(xsq.*ztq - zsq.*xtq);
    syJf =  (xrq.*ztq - zrq.*xtq);
    tyJf = -(xrq.*zsq - zrq.*xsq);
    
    rzJf = (xsq.*ytq - ysq.*xtq);
    szJf = -(xrq.*ytq - yrq.*xtq);
    tzJf = (xrq.*ysq - yrq.*xsq);        
    
    nxJ = rxJf.*nrJ + sxJf.*nsJ + txJf.*ntJ;
    nyJ = ryJf.*nrJ + syJf.*nsJ + tyJf.*ntJ;
    nzJ = rzJf.*nrJ + szJf.*nsJ + tzJf.*ntJ;
    
    nx = nxJ./Jf;
    ny = nyJ./Jf;
    nz = nzJ./Jf;
    
    sJ = sqrt(nx.^2 + ny.^2 + nz.^2);
    nx = nx./sJ;
    ny = ny./sJ;
    nz = nz./sJ;
    sJ = sJ.*Jf;
        
    wJq = diag(wq)*(Jq);
    for e = 1:K
        MM = Vq'*diag(wq.*Jq(:,e))*Vq;
        uproj(:,e) = MM\(Vq'*(wq.*Jq(:,e).*f(xq(:,e),yq(:,e),zq(:,e))));        
    end
    uwadg = Pq*((Vq*Pq*(f(xq,yq,zq).*Jq))./Jq);
    L2err1(kk) = sqrt(sum(sum(wJq.*(f(xq,yq,zq)-Vq*uproj).^2)));
    L2err2(kk) = sqrt(sum(sum(wJq.*(f(xq,yq,zq)-Vq*uwadg).^2)));
    err = Vq*(uproj-uwadg);
    udiff(kk) = sqrt(sum(sum(wJq.*err.^2)));
    
%     fq = f(xq,yq,zq);
% %     fq = Vq*Pq*fq;    
%     for e = 1:K
%         TinvJfwadg(:,e) = (Vq'*diag(wq./Jq(:,e))*Vq)\(M*Pq*fq(:,e));        
%     end        
%     err3 = (Vq*TinvJfwadg - Jq.*fq);
%     udiff(kk) = sqrt(abs(sum(sum(wJq.*err3.^2)))); % N+2 rate
    
    %     hvec(kk) = max(max(Jf./sJ));
    geo(kk) = sqrt(sum(sum(wJq.*rxJ.^2)));
    JJ(kk) = sqrt(sum(sum(wJq.*Jq.^2)));
end
%%
hvec = .5.^(1:length(L2err1));
% loglog(hvec,J,'o--');hold on;loglog(hvec,hvec.^2);return

loglog(hvec,L2err1,'bo--','markersize',10);
hold on
loglog(hvec,L2err2,'bx--');
loglog(hvec,.5*hvec.^(N+1),'b--');

loglog(hvec,udiff,'ro--','markersize',10);
% loglog(hvec,hvec.^(N+1).*hvec.^3,'k--');
loglog(hvec,1e-2*hvec.^(N+2),'k--');

