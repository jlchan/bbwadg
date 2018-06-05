clear
Globals3D
N = 4;

f = @(x,y) exp(x+y).*sin(pi*x).*sin(pi*y);

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
    
    DNr = [Vq*Dr*Pq - .5*Vq*Lq*diag(nrJ)*Vfq*Pq .5*Vq*Lq*diag(nrJ);
        -.5*diag(nrJ)*Vfq*Pq .5*diag(nrJ)];
    DNs = [Vq*Ds*Pq - .5*Vq*Lq*diag(nsJ)*Vfq*Pq .5*Vq*Lq*diag(nsJ);
        -.5*diag(nsJ)*Vfq*Pq .5*diag(nsJ)];
    DNt = [Vq*Dt*Pq - .5*Vq*Lq*diag(ntJ)*Vfq*Pq .5*Vq*Lq*diag(ntJ);
        -.5*diag(ntJ)*Vfq*Pq .5*diag(ntJ)];
    WN = diag([wq;wfq]);
    
    DNr_skew = [Vq*Dr*Pq - .5*Vq*Lq*diag(nrJ)*Vfq*Pq .5*Vq*Lq*diag(nrJ);
        -.5*diag(nrJ)*Vfq*Pq 0*.5*diag(nrJ)];
    QNr_skew = WN*DNr_skew;
    
    QNr = WN*DNr;
    QNs = WN*DNs;
    QNt = WN*DNt;
    
    % A = rand(size(QNr)); A = A+A';
    
    % WN = diag([wq;wfq]);
    % Qr = [Drq - .5*Vq*Lq*diag(nrJ)];
    
    %%
    
    % original cross product form
    xrq = Drq*x; xsq = Dsq*x; xtq = Dtq*x;
    yrq = Drq*y; ysq = Dsq*y; ytq = Dtq*y;
    zrq = Drq*z; zsq = Dsq*z; ztq = Dtq*z;
    J = xrq.*(ysq.*ztq-zsq.*ytq) - yrq.*(xsq.*ztq-zsq.*xtq) + zrq.*(xsq.*ytq-ysq.*xtq);

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
    
    % Nfqf = length(rqtri);
    % e = ones(Nfqf,1); z = zeros(Nfqf,1);
    % nrJ = [z;z;e;-e];
    % nsJ = [z;-e;e;z];
    % ntJ = [-e;z;e;z];
    
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
    
    %%
    wJq = diag(wq)*(J);
    for e = 1:K
        MM = Vq'*diag(wq.*J(:,e))*Vq;
        uproj(:,e) = MM\(Vq'*(wq.*J(:,e).*f(xq(:,e),yq(:,e))));        
    end
    uwadg = Pq*((Vq*Pq*(f(xq,yq).*J))./J);
    L2err1(kk) = sqrt(sum(sum(wJq.*(f(xq,yq)-Vq*uproj).^2)));
    L2err2(kk) = sqrt(sum(sum(wJq.*(f(xq,yq)-Vq*uwadg).^2)));
    udiff(kk) = sqrt(sum(sum(wJq.*(Vq*(uproj-uwadg)).^2)));
%     hvec(kk) = max(max(Jf./sJ));
end

hvec = .5.^(1:length(L2err1));
loglog(hvec,L2err1,'bo--','markersize',10);
hold on
loglog(hvec,L2err2,'bx--');
loglog(hvec,.5*hvec.^(N+1),'b--');
loglog(hvec,udiff,'r*--','markersize',10);
loglog(hvec,1e-2*hvec.^(N+2),'r--');