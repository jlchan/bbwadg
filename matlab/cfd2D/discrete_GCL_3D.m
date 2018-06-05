
clear
Globals3D
N = 3;

for kk = 1%1:4
    filename = ['Grid/cube' num2str(2+kk) '.msh']
    [Nv, VX, VY, VZ, K, EToV] = MeshReaderGmsh3D(filename);
    VX = VX*2; VY = VY*2; VZ = VZ*2; % map -1,1
    StartUp3D;
    
    % make curvilinear
    a = .125;
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
    
    %% conservative curl form: interpolate to P(N+1) then down to PN
    
    [r2 s2 t2] = Nodes3D(N+1); [r2 s2 t2] = xyztorst(r2,s2,t2);
    V2 = Vandermonde3D(N+1,r2,s2,t2); % interp to degree N+1
    [V2r V2s V2t] = GradVandermonde3D(N+1,r2,s2,t2);
    D2r = V2r/V2; D2s = V2s/V2; D2t = V2t/V2;
    
    xx = (Vandermonde3D(N,r2,s2,t2)/V)*x;
    yy = (Vandermonde3D(N,r2,s2,t2)/V)*y;
    zz = (Vandermonde3D(N,r2,s2,t2)/V)*z;
    
    xr = Dr*x; xs = Ds*x; xt = Dt*x;
    yr = Dr*y; ys = Ds*y; yt = Dt*y;
    zr = Dr*z; zs = Ds*z; zt = Dt*z;

    rxJ1 = Dt*(ys.*z) - Ds*(yt.*z);
    sxJ1 = Dr*(yt.*z) - Dt*(yr.*z);
    txJ1 = Ds*(yr.*z) - Dr*(ys.*z);
    ryJ1 = -(Dt*(xs.*z) - Ds*(xt.*z));
    syJ1 = -(Dr*(xt.*z) - Dt*(xr.*z));
    tyJ1 = -(Ds*(xr.*z) - Dr*(xs.*z));    
    rzJ1 = -(Dt*(ys.*x) - Ds*(yt.*x));
    szJ1 = -(Dr*(yt.*x) - Dt*(yr.*x));
    tzJ1 = -(Ds*(yr.*x) - Dr*(ys.*x));     
    
    xr = D2r*xx; xs = D2s*xx; xt = D2t*xx;
    yr = D2r*yy; ys = D2s*yy; yt = D2t*yy;
    zr = D2r*zz; zs = D2s*zz; zt = D2t*zz;
    rxJnew = D2t*(ys.*zz) - D2s*(yt.*zz);
    sxJnew = D2r*(yt.*zz) - D2t*(yr.*zz);
    txJnew = D2s*(yr.*zz) - D2r*(ys.*zz);
    
    ryJnew = -(D2t*(xs.*zz) - D2s*(xt.*zz));
    syJnew = -(D2r*(xt.*zz) - D2t*(xr.*zz));
    tyJnew = -(D2s*(xr.*zz) - D2r*(xs.*zz));
    
    rzJnew = -(D2t*(ys.*xx) - D2s*(yt.*xx));
    szJnew = -(D2r*(yt.*xx) - D2t*(yr.*xx));
    tzJnew = -(D2s*(yr.*xx) - D2r*(ys.*xx));
    
    rxJnew = (Vandermonde3D(N+1,r,s,t)/V2)*rxJnew; % interp to PN
    sxJnew = (Vandermonde3D(N+1,r,s,t)/V2)*sxJnew;
    txJnew = (Vandermonde3D(N+1,r,s,t)/V2)*txJnew;
    ryJnew = (Vandermonde3D(N+1,r,s,t)/V2)*ryJnew; % interp to PN
    syJnew = (Vandermonde3D(N+1,r,s,t)/V2)*syJnew;
    tyJnew = (Vandermonde3D(N+1,r,s,t)/V2)*tyJnew;
    rzJnew = (Vandermonde3D(N+1,r,s,t)/V2)*rzJnew; % interp to PN
    szJnew = (Vandermonde3D(N+1,r,s,t)/V2)*szJnew;
    tzJnew = (Vandermonde3D(N+1,r,s,t)/V2)*tzJnew;
    
    
    %norm(Dr*rxJnew + Ds*sxJnew + Dt*txJnew,'fro')
    
    wJq = diag(wq)*ones(size(J));
    rxJsize(kk) = sqrt(sum(sum(wJq.*rxJ.^2)));
    
    err2 = (rxJ - Vq*rxJnew).^2+(sxJ - Vq*sxJnew).^2+(txJ - Vq*txJnew).^2;
    err2 = err2 + (ryJ - Vq*ryJnew).^2+(syJ - Vq*syJnew).^2+(tyJ - Vq*tyJnew).^2;
    err2 = err2 + (rzJ - Vq*rzJnew).^2+(szJ - Vq*szJnew).^2+(tzJ - Vq*tzJnew).^2;
    %     wJq = diag(wfq)*ones(length(wfq),K);
    %   nxJnew = nrJ.*(Vfq*rxJnew) + nsJ.*(Vfq*sxJnew) + ntJ.*(Vfq*txJnew);  err2 = (nxJ - nxJnew).^2;
    L2errnew(kk) = sqrt(sum(sum(wJq.*err2)));
    
    err2 = (rxJ - Vq*rxJ1).^2+(sxJ - Vq*sxJ1).^2+(txJ - Vq*txJ1).^2;
    err2 = err2 + (ryJ - Vq*ryJ1).^2+(syJ - Vq*syJ1).^2+(tyJ - Vq*tyJ1).^2;
    err2 = err2 + (rzJ - Vq*rzJ1).^2+(szJ - Vq*szJ1).^2+(tzJ - Vq*tzJ1).^2;
    %     nxJ1 = nrJ.*(Vfq*rxJ1) + nsJ.*(Vfq*sxJ1) + ntJ.*(Vfq*txJ1);  err2 = (nxJ - nxJ1).^2;
    L2err(kk) = sqrt(sum(sum(wJq.*err2)));
    
    xfq = Vfq*x;
    yfq = Vfq*y;
    zfq = Vfq*z;
    [mapM, mapP] = BuildNodeMaps3D(xfq,yfq,zfq,EToE,EToF);
    
    nxJnew = diag(nrJ)*(Vfq*rxJnew) + diag(nsJ)*(Vfq*sxJnew) + diag(ntJ)*(Vfq*txJnew);
    nyJnew = diag(nrJ)*(Vfq*ryJnew) + diag(nsJ)*(Vfq*syJnew) + diag(ntJ)*(Vfq*tyJnew);
    nzJnew = diag(nrJ)*(Vfq*rzJnew) + diag(nsJ)*(Vfq*szJnew) + diag(ntJ)*(Vfq*tzJnew);
    
    idM = mapM(mapM~=mapP);
    idP = mapP(mapM~=mapP);
    err1 = norm(nxJ(idM)+nxJ(idP),'fro');
    fprintf('error between nxJ+/-: regular %g, dGCL: %g\n',err1,norm(nxJnew(idM)+nxJnew(idP),'fro'))
    
    gcl_res = norm(DNr*[Vq;Vfq]*rxJnew + DNs*[Vq;Vfq]*sxJnew + DNt*[Vq;Vfq]*txJnew,'fro');
    gcl_res = gcl_res + norm(DNr*[Vq;Vfq]*ryJnew + DNs*[Vq;Vfq]*syJnew + DNt*[Vq;Vfq]*tyJnew,'fro');
    gcl_res = gcl_res + norm(DNr*[Vq;Vfq]*rzJnew + DNs*[Vq;Vfq]*szJnew + DNt*[Vq;Vfq]*tzJnew,'fro');
    fprintf('dGCL residual: %g\n',gcl_res)
    % rxJ = Dt*(ys.*z) - Ds*(yt.*z);
    % sxJ = Dr*(yt.*z) - Dt*(yr.*z);
    % txJ = Ds*(yr.*z) - Dr*(ys.*z);
    %
    % ryJ = -(Dt*(xs.*z) - Ds*(xt.*z));
    % syJ = -(Dr*(xt.*z) - Dt*(xr.*z));
    % tyJ = -(Ds*(xr.*z) - Dr*(xs.*z));
    %
    % rzJ = -(Dt*(ys.*x) - Ds*(yt.*x));
    % szJ = -(Dr*(yt.*x) - Dt*(yr.*x));
    % tzJ = -(Ds*(yr.*x) - Dr*(ys.*x));    
    
    
    
    hvec(kk) = max(max(Jf./sJ));
end

%%

loglog(hvec,rxJsize,'o--');hold on;loglog(hvec,.65*hvec.^(.5),'k--');return

loglog(hvec,L2errnew,'o--');hold on;
r = (N+2);
scale = 1.*L2errnew(1)/hvec(1).^r;
loglog(hvec,scale*hvec.^r,'k--')
print_pgf_coordinates(hvec,L2errnew)
print_pgf_coordinates(hvec,scale*hvec.^r)


loglog(hvec,L2err,'x--');hold on;
r = (N+1);
scale = 1.5*L2err(1)/hvec(1).^r;
loglog(hvec,scale*hvec.^r,'k--')

print_pgf_coordinates(hvec,L2err)
print_pgf_coordinates(hvec,scale*hvec.^r)
