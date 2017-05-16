function Wave3D_IGA

Globals3D

NB = 2;
Ksub = 2;
smoothKnots = 0;

if Ksub==1
    CN = (NB+1)^2/2;
    dt = 2/CN; % size of domain
else
    CN = (NB+2)*Ksub;
    dt = 2/CN;
end

N = NB+Ksub-1;

Np = (N+1)^3;
Nfp = (N+1)^2;

FinalTime = .5;


%% reference elem stuff

global rxJ sxJ txJ ryJ syJ tyJ rzJ szJ tzJ
global nx ny nz Jq sJ
global DVq1D Vq1D Pq1D Prq1D D1D D1Dt LIFT Fmask
global mapM mapP mapB Nfaces

Nfaces = 6;
if Ksub==1
    [r1D w1D] = JacobiGL(0,0,N);
else
    VX = linspace(-1,1,Ksub+1);
    t = [VX(1)*ones(1,NB) VX VX(end)*ones(1,NB)]; % open knot vec
    for i = 1:N+1
        r1D(i,1) = mean(t((i+1):(i+NB))); % greville
    end
    %     plot(r1D,r1D*0,'o')
end

[r s t] = meshgrid(r1D); 
r = r(:); s = s(:); t = t(:); 

Fmask1 = find(abs(r+1)<1e-8); 
Fmask2 = find(abs(r-1)<1e-8); 
Fmask3 = find(abs(s+1)<1e-8); 
Fmask4 = find(abs(s-1)<1e-8); 
Fmask5 = find(abs(t+1)<1e-8); 
Fmask6 = find(abs(t-1)<1e-8); 
Fmask = [Fmask1(:) Fmask2(:) Fmask3(:) Fmask4(:) Fmask5(:) Fmask6(:)];

rp1D = linspace(-1,1,100);
[rp sp tp] = meshgrid(rp1D);
rp = rp(:); sp = sp(:); tp = tp(:);

% initialize TP operators
[BVDM M1D D1D R rq1D wq1D Vq1D DVq1D] = bsplineVDM(NB,Ksub,r1D,smoothKnots); % VDM for interp, mass, M\S
Pq1D = M1D\(Vq1D'*diag(wq1D));
Prq1D = M1D\(DVq1D'*diag(wq1D));

if Ksub >= NB
    Vq1D = sparse(Vq1D);
    Pq1D = sparse(Pq1D);
    Prq1D = sparse(Prq1D);
end

(N+1)
Nq = size(Vq1D,1)

D1Dt = M1D\(D1D'*M1D); % for skew-symmetric formulation

[rq sq tq] = meshgrid(rq1D); rq = rq(:); sq = sq(:); tq = tq(:);
[wrq wsq wtq] = meshgrid(wq1D); wq = wrq(:).*wsq(:).*wtq(:);

% make face quadrature nodes
rfq = zeros(length(rq1D)^2,Nfaces);
sfq = zeros(length(rq1D)^2,Nfaces);
tfq = zeros(length(rq1D)^2,Nfaces);

[rq2 sq2] = meshgrid(rq1D); rq2 = rq2(:); sq2 = sq2(:); 
[wr2 ws2] = meshgrid(wq1D); w2 = wr2(:).*ws2(:); w2 = repmat(w2,Nfaces,1);
% r = -/+ 1
rfq(:,1) = -1; sfq(:,1) = rq2; tfq(:,1) = sq2; 
rfq(:,2) = 1;  sfq(:,2) = rq2; tfq(:,2) = sq2; 
% s = -/+ 1
rfq(:,3) = rq2; sfq(:,3) = -1; tfq(:,3) = sq2; 
rfq(:,4) = rq2; sfq(:,4) = 1;  tfq(:,4) = sq2; 
% t = -/+ 1
rfq(:,5) = rq2; sfq(:,5) = sq2; tfq(:,5) = -1; 
rfq(:,6) = rq2; sfq(:,6) = sq2; tfq(:,6) =  1; 

rfq = rfq(:); sfq = sfq(:); tfq = tfq(:);

% sum(w2)
% plot3(rfq,sfq,tfq,'o'); return 

% invM1D = inv(M1D);
Mf1D = zeros(N+1,2);
Mf1D(1,1) = 1;
Mf1D(N+1,2) = 1;
LIFT1D = M1D\Mf1D;

% make lift operator
ids = 1:(N+1)^3;
ids = reshape(ids,N+1,N+1,N+1);

LIFT = sparse((N+1)^3,6*(N+1)^2); % initialize as spares matrix

% first/second face (r = +/- 1)
sk = 0;
for i = 1:(N+1)
    iid = squeeze(ids(i,:,:));
    for j = 1:(N+1)
        LIFT(iid(:,j),sk + (j-1)*(N+1) + i) = LIFT1D(:,1); % left-to-right
        LIFT(iid(:,j),sk + (N+1)^2 + (j-1)*(N+1) + i) = LIFT1D(:,2); % right-to-left
    end
end

% third/fourth face (s = +/- 1)
sk = sk + 2*(N+1)^2;
for i = 1:(N+1)
    iid = squeeze(ids(:,:,i));
    for j = 1:(N+1)
        LIFT(iid(:,j),sk + (i-1)*(N+1) + j) = LIFT1D(:,1); % left-to-right
        LIFT(iid(:,j),sk + (N+1)^2 + (i-1)*(N+1) + j) = LIFT1D(:,2); % right-to-left
    end
end

% third/fourth face (s = +/- 1)
sk = sk + 2*(N+1)^2;
for i = 1:(N+1)
    iid = squeeze(ids(:,i,:));
    for j = 1:(N+1)
        LIFT(iid(j,:),sk + (i-1)*(N+1) + j) = LIFT1D(:,1); % left-to-right
        LIFT(iid(j,:),sk + (N+1)^2 + (i-1)*(N+1) + j) = LIFT1D(:,2); % right-to-left
    end
end

%% make mesh and connectivity
if 0 % use preset mesh
    
    addpath IGA_DG   
    addpath IGA_DG/Elbow3D
    
    %% recover vertices and connectivity
    re1D = [-1; 1];
    [re se te] = meshgrid(re1D); re = re(:); se = se(:); te = te(:);
    [vx vy vz] = Elbow3Dgeo(re,se,te);
    [VXYZ idr idc] = uniquetol([vx(:) vy(:) vz(:)],1e-7,'ByRows',true);
    VX = VXYZ(:,1); VY = VXYZ(:,2); VZ = VXYZ(:,3);
    EToV = reshape(idc,8,12)';
    
    [EToE,EToF]= tiConnect3DHex(EToV);    
    
    
    %% get greville nodes for connectivity
    [x y z geo] = Elbow3Dgeo(r,s,t);    
    
    K = size(x,2);
    
    xf = reshape(x(Fmask(:),:),Nfp*Nfaces,K);
    yf = reshape(y(Fmask(:),:),Nfp*Nfaces,K);
    zf = reshape(z(Fmask(:),:),Nfp*Nfaces,K);
    
    mapM = reshape(1:Nfp*Nfaces*K,Nfp*Nfaces,K);
    mapP = reshape(1:Nfp*Nfaces*K,Nfp*Nfaces,K);
    
    NODETOL = 1e-6;
    for e = 1:K
        for f = 1:Nfaces
            enbr = EToE(e,f);
            fnbr = EToF(e,f);
            idM = (1:Nfp) + (f-1)*Nfp;
            idP = (1:Nfp) + (fnbr-1)*Nfp;
            
            if e==enbr
                mapP(idM,e) = mapM(idM,e);
            else
                
                xM = xf(idM,e); yM = yf(idM,e); zM = zf(idM,e);
                xP = xf(idP,enbr); yP = yf(idP,enbr); zP = zf(idP,enbr);
                
                % find find volume node numbers of left and right nodes
                x1 = xf(idM,e); y1 = yf(idM,e); z1 = zf(idM,e);
                x2 = xf(idP,enbr); y2 = yf(idP,enbr); z2 = zf(idP,enbr);
                
                one = ones(1, Nfp);
                x1 = x1*one;  y1 = y1*one;  z1 = z1*one;
                x2 = x2*one;  y2 = y2*one;  z2 = z2*one;
                
                % Compute distance matrix
                D = (x1 -x2').^2 + (y1-y2').^2 + (z1-z2').^2;
                [idM, idP] = find(sqrt(abs(D)) < NODETOL);
                
                mapP(idM + (f-1)*Nfp,e) = idP + (fnbr-1)*Nfp + (enbr-1)*Nfaces*Nfp;
            end
        end
    end
    
    mapB = find(mapM==mapP);
    
    vmapM = reshape(1:Np*K,Np,K);
    vmapM = vmapM(Fmask(:),:);
    vmapP = vmapM(mapP);
%     keyboard
    
    err1 = norm(xf(mapM) - xf(mapP),'fro')+norm(yf(mapM) - yf(mapP),'fro')+norm(zf(mapM) - zf(mapP),'fro');
    err2 = norm(x(vmapM) - x(vmapP),'fro')+norm(y(vmapM) - y(vmapP),'fro')+norm(z(vmapM) - z(vmapP),'fro');
    if err1 + err2 > 1e-10
        keyboard
    end

%     e = 2; plot3(xf(:,e),yf(:,e),zf(:,e),'o');axis equal;return
    
    %% get quadrature nodes
    [xq yq zq geo] = Elbow3Dgeo(rq,sq,tq);
    %     plot3(xq,yq,zq,'o');  axis equal;   return
    rxJ = geo.rxJ; sxJ = geo.sxJ; txJ = geo.txJ;
    ryJ = geo.ryJ; syJ = geo.syJ; tyJ = geo.tyJ;
    rzJ = geo.rzJ; szJ = geo.szJ; tzJ = geo.tzJ;
    Jq = geo.J;
    
    %% surface quadrature nodes and normals
    
    Nfq = length(rfq)/Nfaces;
    [xfq yfq zfq fgeo] = Elbow3Dgeo(rfq,sfq,tfq);
    %     plot3(xfq,yfq,zfq,'o');  axis equal;   return
    rxJf = fgeo.rxJ; sxJf = fgeo.sxJ; txJf = fgeo.txJ;
    ryJf = fgeo.ryJ; syJf = fgeo.syJ; tyJf = fgeo.tyJ;
    rzJf = fgeo.rzJ; szJf = fgeo.szJ; tzJf = fgeo.tzJ;
    Jf = fgeo.J;        
    
    nr = zeros(Nfq*Nfaces,1);
    ns = zeros(Nfq*Nfaces,1);
    nt = zeros(Nfq*Nfaces,1);    
    fids = @(f) (1:Nfq) + (f-1)*Nfq;
    f = 1; nr(fids(f)) = -1;
    f = 2; nr(fids(f)) = 1;
    f = 3; ns(fids(f)) = -1;
    f = 4; ns(fids(f)) = 1;
    f = 5; nt(fids(f)) = -1;
    f = 6; nt(fids(f)) = 1;
    
    nx = zeros(Nfq*Nfaces,K);
    ny = zeros(Nfq*Nfaces,K);
    nz = zeros(Nfq*Nfaces,K);
    sJ = zeros(Nfq*Nfaces,K);
    
    for e = 1:K
        
%         rxf = rxJf(:,e); ryf = ryJf(:,e); rzf = rzJf(:,e)
%         sxf = sxJf(:,e); syf = syJf(:,e); szf = szJf(:,e)
%         txf = txJf(:,e); tyf = tyJf(:,e); tzf = tzJf(:,e)
%         
%         nx(:,e) = rxf.*nr + sxf.*ns + txf.*nt;
%         ny(:,e) = ryf.*nr + syf.*ns + tyf.*nt;
%         nz(:,e) = rzf.*nr + szf.*ns + tzf.*nt;
        
                rxf = rxJf(:,e)./Jf(:,e); ryf = ryJf(:,e)./Jf(:,e); rzf = rzJf(:,e)./Jf(:,e);
                sxf = sxJf(:,e)./Jf(:,e); syf = syJf(:,e)./Jf(:,e); szf = szJf(:,e)./Jf(:,e);
                txf = txJf(:,e)./Jf(:,e); tyf = tyJf(:,e)./Jf(:,e); tzf = tzJf(:,e)./Jf(:,e);
                for f = 1:Nfaces
        
                    fids = (1:Nfq) + (f-1)*Nfq;
                    if f==1
                        nx(fids,e) = -rxf(fids); ny(fids,e) = -ryf(fids); nz(fids,e) = -rzf(fids);
                    elseif f==2
                        nx(fids,e) = rxf(fids); ny(fids,e) = ryf(fids); nz(fids,e) = rzf(fids);
                    elseif f==3
                        nx(fids,e) = -sxf(fids); ny(fids,e) = -syf(fids); nz(fids,e) = -szf(fids);
                    elseif f==4
                        nx(fids,e) = sxf(fids); ny(fids,e) = syf(fids); nz(fids,e) = szf(fids);
                    elseif f==5
                        nx(fids,e) = -txf(fids); ny(fids,e) = -tyf(fids); nz(fids,e) = -tzf(fids);
                    elseif f==6
                        nx(fids,e) = txf(fids); ny(fids,e) = tyf(fids); nz(fids,e) = tzf(fids);
                    end
                end
    end
    
    sJ = sqrt(nx.^2 + ny.^2 + nz.^2);
    nx = nx./sJ; ny = ny./sJ; nz = nz./sJ;
    
    keyboard
    
    if 1
        fids = 1:Nfq*Nfaces;
        for e = 1:K;
            clf
            plot3(xfq(fids,e),yfq(fids,e),zfq(fids,e),'o');
            hold on
            quiver3(xfq(fids,e),yfq(fids,e),zfq(fids,e),nx(fids,e),ny(fids,e),nz(fids,e))
            axis equal
            pause
        end
        keyboard
    end
    if min(Jq(:)) < 1e-8
        keyboard
    end
    
        

    

else % use two-element mesh
    
%     VX = [-1 1 1 -1 -1 1 1 -1];
%     VY = [-1 -1 1 1 -1 -1 1 1];
%     a = 0*.125;
%     VZ = [-1 -1 -1 -1 a -a a -a];
%     EToV = [1 4 2 3 5 8 6 7];
    
    VX = [-1 1 1 -1 -1 1 1 -1 -1 1 1 -1];
    VY = [-1 -1 1 1 -1 -1 1 1 -1 -1 1 1];
    a = .125;
    VZ = [-1 -1 -1 -1 a -a a -a 1 1 1 1];    
    EToV = [1 4 2 3 5 8 6 7; 5 8 6 7 9 12 10 11];
    [EToE,EToF]= tiConnect3DHex(EToV);
    
    K = size(EToV,1);
    
    V11D = Vandermonde1D(1,r1D)/Vandermonde1D(1,[-1 1]);
    x = zeros(Np,K);
    y = zeros(Np,K);
    z = zeros(Np,K);
    for e = 1:K
        x(:,e) = matvec(V11D,(VX(EToV(e,:))'),'all');
        y(:,e) = matvec(V11D,(VY(EToV(e,:))'),'all');
        z(:,e) = matvec(V11D,(VZ(EToV(e,:))'),'all');
    end
    
    % keyboard
    % make connectivity
    xf = reshape(x(Fmask(:),:),Nfp*Nfaces,K);
    yf = reshape(y(Fmask(:),:),Nfp*Nfaces,K);
    zf = reshape(z(Fmask(:),:),Nfp*Nfaces,K);
    
    mapM = reshape(1:Nfp*Nfaces*K,Nfp*Nfaces,K);
    mapP = reshape(1:Nfp*Nfaces*K,Nfp*Nfaces,K);
    
    for e = 1:K
        for f = 1:Nfaces
            enbr = EToE(e,f);
            fnbr = EToF(e,f);
            idM = (1:Nfp) + (f-1)*Nfp;
            idP = (1:Nfp) + (fnbr-1)*Nfp;
            
            if e==enbr
                mapP(idM,e) = mapM(idM,e);
            else
                
                xM = xf(idM,e); yM = yf(idM,e); zM = zf(idM,e);
                xP = xf(idP,enbr); yP = yf(idP,enbr); zP = zf(idP,enbr);
                
                % find find volume node numbers of left and right nodes
                x1 = xf(idM,e); y1 = yf(idM,e); z1 = zf(idM,e);
                x2 = xf(idP,enbr); y2 = yf(idP,enbr); z2 = zf(idP,enbr);
                
                one = ones(1, Nfp);
                x1 = x1*one;  y1 = y1*one;  z1 = z1*one;
                x2 = x2*one;  y2 = y2*one;  z2 = z2*one;
                
                % Compute distance matrix
                D = (x1 -x2').^2 + (y1-y2').^2 + (z1-z2').^2;
                [idM, idP] = find(sqrt(abs(D)) < NODETOL);
                
                mapP(idM + (f-1)*Nfp,e) = idP + (fnbr-1)*Nfp + (enbr-1)*Nfaces*Nfp;
            end
        end
    end
    
    mapB = find(mapM==mapP);
    
    vmapM = reshape(1:Np*K,Np,K);
    vmapM = vmapM(Fmask(:),:);
    vmapP = vmapM(mapP);
    
    for e = 1:K
        xr(:,e) = matvec(D1D,x(:,e),1);
        xs(:,e) = matvec(D1D,x(:,e),2);
        xt(:,e) = matvec(D1D,x(:,e),3);
        yr(:,e) = matvec(D1D,y(:,e),1);
        ys(:,e) = matvec(D1D,y(:,e),2);
        yt(:,e) = matvec(D1D,y(:,e),3);
        zr(:,e) = matvec(D1D,z(:,e),1);
        zs(:,e) = matvec(D1D,z(:,e),2);
        zt(:,e) = matvec(D1D,z(:,e),3);
    end
    
    J = xr.*(ys.*zt-zs.*yt) - yr.*(xs.*zt-zs.*xt) + zr.*(xs.*yt-ys.*xt);
    rxJ =  (ys.*zt - zs.*yt); ryJ = -(xs.*zt - zs.*xt); rzJ = (xs.*yt - ys.*xt);
    sxJ = -(yr.*zt - zr.*yt); syJ =  (xr.*zt - zr.*xt); szJ = -(xr.*yt - yr.*xt);
    txJ =  (yr.*zs - zr.*ys); tyJ = -(xr.*zs - zr.*xs); tzJ = (xr.*ys - yr.*xs);
    
    % normals
    nx = zeros(Nfp*Nfaces,K);
    ny = zeros(Nfp*Nfaces,K);
    nz = zeros(Nfp*Nfaces,K);
    sJ = zeros(Nfp*Nfaces,K);
    for e = 1:K
        
%         fids = @(f) (1:Nfp) + (f-1)*Nfp;
%         nr = zeros(Nfp*Nfaces,1);
%         ns = zeros(Nfp*Nfaces,1);
%         nt = zeros(Nfp*Nfaces,1);
%         f = 1; nr(fids(f)) = -1;
%         f = 2; nr(fids(f)) = 1;
%         f = 3; ns(fids(f)) = -1;
%         f = 4; ns(fids(f)) = 1;
%         f = 5; nt(fids(f)) = -1;
%         f = 6; nt(fids(f)) = 1;
%         
%         nx(:,e) = rxJ(Fmask(:),e).*nr + sxJ(Fmask(:),e).*ns + txJ(Fmask(:),e).*nt;
%         ny(:,e) = ryJ(Fmask(:),e).*nr + syJ(Fmask(:),e).*ns + tyJ(Fmask(:),e).*nt;
%         nz(:,e) = rzJ(Fmask(:),e).*nr + szJ(Fmask(:),e).*ns + tzJ(Fmask(:),e).*nt;        

        rxK = rxJ(:,e)./J(:,e); ryK = ryJ(:,e)./J(:,e); rzK = rzJ(:,e)./J(:,e);
        sxK = sxJ(:,e)./J(:,e); syK = syJ(:,e)./J(:,e); szK = szJ(:,e)./J(:,e);
        txK = txJ(:,e)./J(:,e); tyK = tyJ(:,e)./J(:,e); tzK = tzJ(:,e)./J(:,e);
        
        for f = 1:Nfaces
            rxf = rxK(Fmask(:,f)); ryf = ryK(Fmask(:,f)); rzf = rzK(Fmask(:,f));
            sxf = sxK(Fmask(:,f)); syf = syK(Fmask(:,f)); szf = szK(Fmask(:,f));
            txf = txK(Fmask(:,f)); tyf = tyK(Fmask(:,f)); tzf = tzK(Fmask(:,f));
            
            fids = (1:(N+1)^2) + (f-1)*(N+1)^2;
            if f==1
                nx(fids,e) = -rxf; ny(fids,e) = -ryf; nz(fids,e) = -rzf;
            elseif f==2
                nx(fids,e) = rxf; ny(fids,e) = ryf; nz(fids,e) = rzf;
            elseif f==3
                nx(fids,e) = -sxf; ny(fids,e) = -syf; nz(fids,e) = -szf;
            elseif f==4
                nx(fids,e) = sxf; ny(fids,e) = syf; nz(fids,e) = szf;
            elseif f==5
                nx(fids,e) = -txf; ny(fids,e) = -tyf; nz(fids,e) = -tzf;
            elseif f==6
                nx(fids,e) = txf; ny(fids,e) = tyf; nz(fids,e) = tzf;
            end
        end
    end
    
    sJ = sqrt(nx.^2 + ny.^2 + nz.^2);
    nx = nx./sJ; ny = ny./sJ; nz = nz./sJ;
    
%     e = 1;
%     plot3(xf(:,e),yf(:,e),zf(:,e),'o')
%     hold on
%     quiver3(xf(:,e),yf(:,e),zf(:,e),nx(:,e),ny(:,e),nz(:,e),2)
%     keyboard
    
    % [max(J(:)),min(J(:))]
    err1 = norm(xf(mapM) - xf(mapP),'fro')+norm(yf(mapM) - yf(mapP),'fro')+norm(zf(mapM) - zf(mapP),'fro');
    err2 = norm(x(vmapM) - x(vmapP),'fro')+norm(y(vmapM) - y(vmapP),'fro')+norm(z(vmapM) - z(vmapP),'fro');
    if min(J(:)) < 1e-8
        keyboard
    end
    if err1 + err2 > 1e-10
        keyboard
    end
    
    % interp to qpts
    for e = 1:K
        rxJtmp(:,e) = matvec(Vq1D,rxJ(:,e),'all'); sxJtmp(:,e) = matvec(Vq1D,sxJ(:,e),'all'); txJtmp(:,e) = matvec(Vq1D,txJ(:,e),'all');
        ryJtmp(:,e) = matvec(Vq1D,ryJ(:,e),'all'); syJtmp(:,e) = matvec(Vq1D,syJ(:,e),'all'); tyJtmp(:,e) = matvec(Vq1D,tyJ(:,e),'all');
        rzJtmp(:,e) = matvec(Vq1D,rzJ(:,e),'all'); szJtmp(:,e) = matvec(Vq1D,szJ(:,e),'all'); tzJtmp(:,e) = matvec(Vq1D,tzJ(:,e),'all');
        Jqtmp(:,e) = matvec(Vq1D,J(:,e),'all');
        
        nxK = reshape(nx(:,e),(N+1)^2,6);
        nyK = reshape(ny(:,e),(N+1)^2,6);
        nzK = reshape(nz(:,e),(N+1)^2,6);
        sJK = reshape(sJ(:,e),(N+1)^2,6);
        Nfq = size(Vq1D,1)^2;
        nxKq = zeros(Nfq,6); nyKq = zeros(Nfq,6); nzKq = zeros(Nfq,6);
        sJKq = zeros(Nfq,6);
        for f = 1:6
            tmp = (Vq1D * reshape(nxK(:,f),N+1,N+1))*Vq1D';
            nxKq(:,f) = tmp(:);
            tmp = (Vq1D * reshape(nyK(:,f),N+1,N+1))*Vq1D';
            nyKq(:,f) = tmp(:);
            tmp = (Vq1D * reshape(nzK(:,f),N+1,N+1))*Vq1D';
            nzKq(:,f) = tmp(:);
            tmp = (Vq1D * reshape(sJK(:,f),N+1,N+1))*Vq1D';
            sJKq(:,f) = tmp(:);
        end
        nxtmp(:,e) = nxKq(:);
        nytmp(:,e) = nyKq(:);
        nztmp(:,e) = nzKq(:);
        sJtmp(:,e) = sJKq(:);
    end
    
    rxJ = rxJtmp; ryJ = ryJtmp; rzJ = rzJtmp;
    sxJ = sxJtmp; syJ = syJtmp; szJ = szJtmp;
    txJ = txJtmp; tyJ = tyJtmp; tzJ = tzJtmp;
    Jq = Jqtmp;
    
    nx = nxtmp;
    ny = nytmp;
    nz = nztmp;
    sJ = sJtmp;

end


%%
% plotting
Vp1D_interp = bsplineVDM(NB,Ksub,rp1D,smoothKnots)/bsplineVDM(NB,Ksub,r1D,smoothKnots);
Vq_nodal_1D = bsplineVDM(NB,Ksub,rq1D,smoothKnots)/bsplineVDM(NB,Ksub,r1D,smoothKnots);
Vp1D = bsplineVDM(NB,Ksub,rp1D,smoothKnots);

for e = 1:K
    xp(:,e) = matvec(Vp1D_interp,x(:,e),'all');
    yp(:,e) = matvec(Vp1D_interp,y(:,e),'all');
    zp(:,e) = matvec(Vp1D_interp,z(:,e),'all');
    
    xq(:,e) = matvec(Vq_nodal_1D,x(:,e),'all');
    yq(:,e) = matvec(Vq_nodal_1D,y(:,e),'all');
    zq(:,e) = matvec(Vq_nodal_1D,z(:,e),'all');
end


%% initial cond

k = 1;
pex = @(x,y,z,t) cos(k*pi*x/2).*cos(k*pi*y/2).*cos(k*pi*z/2).*cos(sqrt(3)*.5*k*pi*t);

% % pex = @(x,y,z,t) x.*y.*z;
% x0 = -1.4; y0 = -.5; z0 = 0;
% pex = @(x,y,z,t) exp(-5*((x-x0).^2 + (y-y0).^2 + (z-z0).^2));

vmapB = vmapM(find(vmapM==vmapP));
xb = x(vmapB);
yb = y(vmapB);
zb = z(vmapB);
ids = find(abs(xb-2)<1e-6);
global mapI
mapI = ids;

% plot3(xb,yb,zb,'o')
% hold on
% plot3(xb(mapI),yb(mapI),zb(mapI),'x')
% return

% p = matvec(inv(BVDM),pex(x,y,z,0),'all');
p = zeros(Np,K);
for e = 1:K
    p(:,e) = matvec(Pq1D,pex(xq(:,e),yq(:,e),zq(:,e),0),'all');
end
u = zeros(Np, K);
v = zeros(Np, K);
w = zeros(Np, K);

for e = 1:K
    pq(:,e) = matvec(Vq1D,p(:,e),'all');
end
wJq = spdiag(wq)*Jq;
err = wJq.*(pq - pex(xq,yq,zq,0)).^2;
init_cond_L2err = sqrt(sum(err(:)))

%% estimate max timestep

if 0
    U = randn((N+1)^3*K,4);
    for i = 1:10
        Uprev = U;
        [rhsp, rhsu, rhsv, rhsw] = acousticsRHS3Dq(reshape(U(:,1),Np,K),reshape(U(:,2),Np,K),reshape(U(:,3),Np,K),reshape(U(:,4),Np,K));
        U(:,1) = rhsp(:);
        U(:,2) = rhsu(:);
        U(:,3) = rhsv(:);
        U(:,4) = rhsw(:);
        
        lam = Uprev(:)'*U(:) / norm(Uprev(:));
        U = U/norm(U(:));
    end
    dt = .5/abs(lam)
else
    h = min(min(J)./max(sJ));
    dt = .5*h/(NB*Ksub);
end
% return

% test operators
if 0
    keyboard
    r2 = BVDM\r; s2 = BVDM\s; t2 = BVDM\t;
    Dr = kron(kron(eye(N+1),D1D),eye(N+1));
    norm(Dr*r2 - matvec(D1D,r2,1))
    Ds = kron(kron(eye(N+1),eye(N+1)),D1D);
    norm(Ds*s2 - matvec(D1D,s2,2))
    Dt = kron(kron(D1D,eye(N+1)),eye(N+1));
    norm(Dt*t2 - matvec(D1D,t2,3))
    
    g = BVDM\(r + 2*s + 3*t + r.*s + s.*t + r.*s.*t);
    norm(BVDM*(Dr*g) - (1 + s + s.*t))
    norm(BVDM*(Ds*g) - (2 + r + t + r.*t))
    norm(BVDM*(Dt*g) - (3 + s + r.*s))
    keyboard
end

if 1
    A = sparse((N+1)^3*4*K,(N+1)^3*4*K);
    U = zeros((N+1)^3*K,4);
    for i = 1:(N+1)^3*4*K
        U(i) = 1;
        [rhsp, rhsu, rhsv, rhsw] = acousticsRHS3Dq(reshape(U(:,1),(N+1)^3,K),reshape(U(:,2),(N+1)^3,K),reshape(U(:,3),(N+1)^3,K),reshape(U(:,4),(N+1)^3,K));
        A(:,i) = sparse([rhsp(:);rhsu(:);rhsv(:);rhsw(:)]);
        U(i) = 0;
        if (mod(i,100)==0)
            disp(sprintf('on i = %d out of %d\n',i,(N+1)^3*4*K))
        end
    end
    lam = eig(full(A));
    plot(lam,'o')
    title(sprintf('max real lam = %g\n',max(real(lam))))
    keyboard
end



%%

% Runge-Kutta residual storage
resp = zeros(Np,K);
resu = zeros(Np,K);
resv = zeros(Np,K);
resw = zeros(Np,K);

Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

% outer time step loop
ids = find(abs(yp) < 2e-2);
% ids = find(abs(yp) < .1);
for tstep = 1:Nsteps
    for INTRK = 1:5
        
        timelocal = tstep*dt + rk4c(INTRK)*dt;
        
        [rhsp,rhsu,rhsv,rhsw] = acousticsRHS3Dq(p,u,v,w);
        
        % initiate and increment Runge-Kutta residuals
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        resw = rk4a(INTRK)*resw + dt*rhsw;
        
        % update fields
        u = u+rk4b(INTRK)*resu;
        v = v+rk4b(INTRK)*resv;
        w = w+rk4b(INTRK)*resw;
        p = p+rk4b(INTRK)*resp;
        
    end;
    
    if nargin==0 && (mod(tstep,10)==0 || tstep==Nsteps)
        clf
        for e = 1:K
            vv(:,e) = matvec(Vp1D,p(:,e),'all');
        end
        color_line3(xp(ids),yp(ids),zp(ids),vv(ids),'.');
        axis equal; axis tight
        view(3)
        colorbar
        title(sprintf('time = %f',tstep*dt))
        drawnow
    end
    
    if mod(tstep,25)==0
        disp(sprintf('on tstep %d out of %d\n',tstep,Nsteps))
    end
end

% J = matvec(Vq1D,J,'all');
for e = 1:K
    pq(:,e) = matvec(Vq1D,p(:,e),'all');
end
wJq = spdiag(wq)*Jq;
pterr = pq - pex(xq,yq,zq,FinalTime);
err = wJq.*(pterr).^2;
L2err = sqrt(sum(err(:)));
L2err

return
% keyboard

ids = find(abs(yq) < .5);
figure
color_line3(xq(ids),yq(ids),zq(ids),pterr(ids),'.');
axis equal


% Nsteps
% L2err * Nsteps % error * steps (assume same overall work)

% work = ceil(FinalTime/dt)*(N+1)^4; % total work
% L2err * work


function [rhsp rhsu rhsv rhsw] = acousticsRHS3Dq(p,u,v,w)

global rxJ sxJ txJ ryJ syJ tyJ rzJ szJ tzJ
global nx ny nz Jq sJ
global DVq1D Vq1D Pq1D Prq1D D1D D1Dt LIFT Fmask
global mapM mapP mapB Nfaces

N = size(Vq1D,2)-1;
Nq = size(Vq1D,1);
Nfp = (N+1)^2; 
K = size(p,2);

pf = p(Fmask(:),:);
uf = u(Fmask(:),:);
vf = v(Fmask(:),:);
wf = w(Fmask(:),:);

dp   = pf(mapP)-pf(mapM);
uavg = uf(mapP)+uf(mapM);
vavg = vf(mapP)+vf(mapM);
wavg = wf(mapP)+wf(mapM);
du   = uf(mapP)-uf(mapM);
dv   = vf(mapP)-vf(mapM);
dw   = wf(mapP)-wf(mapM);

if 1  % dirichlet
    dp(mapB) = -2*pf(mapB);
    uavg(mapB) = uf(mapB);
    vavg(mapB) = vf(mapB);
    wavg(mapB) = wf(mapB);
    du(mapB) = 0;
    dv(mapB) = 0;
    dw(mapB) = 0;
else % neumann
    dp(mapB) = 0;
    uavg(mapB) = 0;
    vavg(mapB) = 0;
    wavg(mapB) = 0;
    du(mapB) = -2*uf(mapB);
    dv(mapB) = -2*vf(mapB);
    dw(mapB) = -2*wf(mapB);
end

for e = 1:K       
    dpK = reshape(dp(:,e),Nfp,Nfaces);
    uavgK = reshape(uavg(:,e),Nfp,Nfaces);
    vavgK = reshape(vavg(:,e),Nfp,Nfaces);
    wavgK = reshape(wavg(:,e),Nfp,Nfaces);
    duK = reshape(du(:,e),Nfp,Nfaces);
    dvK = reshape(dv(:,e),Nfp,Nfaces);
    dwK = reshape(dw(:,e),Nfp,Nfaces);
    
    rxJK = rxJ(:,e); ryJK = ryJ(:,e); rzJK = rzJ(:,e);
    sxJK = sxJ(:,e); syJK = syJ(:,e); szJK = szJ(:,e);
    txJK = txJ(:,e); tyJK = tyJ(:,e); tzJK = tzJ(:,e);        
    
    Nfq = Nq^2;
    nxK = reshape(nx(:,e),Nfq,Nfaces);
    nyK = reshape(ny(:,e),Nfq,Nfaces);
    nzK = reshape(nz(:,e),Nfq,Nfaces);
    sJK = reshape(sJ(:,e),Nfq,Nfaces);    
    
    pr = matvec(Vq1D,matvec(D1D,p(:,e),1),'all');
    ps = matvec(Vq1D,matvec(D1D,p(:,e),2),'all');
    pt = matvec(Vq1D,matvec(D1D,p(:,e),3),'all');
    
    uq = matvec(Vq1D,u(:,e),'all');
    vq = matvec(Vq1D,v(:,e),'all');
    wq = matvec(Vq1D,w(:,e),'all');
    
    dpdx = matvec(Pq1D,rxJK.*pr + sxJK.*ps + txJK.*pt,'all');
    dpdy = matvec(Pq1D,ryJK.*pr + syJK.*ps + tyJK.*pt,'all');
    dpdz = matvec(Pq1D,rzJK.*pr + szJK.*ps + tzJK.*pt,'all');
    
    ur = matvec(Pq1D,rxJK.*uq + ryJK.*vq + rzJK.*wq,'all');
    ur = matvec(D1Dt,ur,1);
    us = matvec(Pq1D,sxJK.*uq + syJK.*vq + szJK.*wq,'all');
    us = matvec(D1Dt,us,2);
    ut = matvec(Pq1D,txJK.*uq + tyJK.*vq + tzJK.*wq,'all');
    ut = matvec(D1Dt,ut,3);
    
%     ur = reshape(rxJK.*uq + ryJK.*vq + rzJK.*wq,Nq,Nq,Nq);
%     us = reshape(sxJK.*uq + syJK.*vq + szJK.*wq,Nq,Nq,Nq);
%     ut = reshape(txJK.*uq + tyJK.*vq + tzJK.*wq,Nq,Nq,Nq);
%     br = zeros(Nq,Nq,N+1);
%     bs = zeros(Nq,Nq,N+1);
%     bt = zeros(Nq,Nq,N+1);
%     for i = 1:Nq
%         for j = 1:Nq
%             br(i,j,:) = Pq1D*reshape(ur(i,j,:),Nq,1);
%             bs(i,j,:) = Pq1D*reshape(us(i,j,:),Nq,1);
%             bt(i,j,:) = Prq1D*reshape(ut(i,j,:),Nq,1);
%         end
%     end
%     ur = zeros(N+1,N+1,N+1);
%     us = zeros(N+1,N+1,N+1);
%     ut = zeros(N+1,N+1,N+1);
%     for i = 1:N+1
%         ur(:,:,i) = Pq1D*reshape(br(:,:,i),Nq,Nq)*Prq1D';
%         us(:,:,i) = Prq1D*reshape(bs(:,:,i),Nq,Nq)*Pq1D';
%         ut(:,:,i) = Pq1D*reshape(bt(:,:,i),Nq,Nq)*Pq1D';
%     end
%     ur = ur(:);
%     us = us(:);
%     ut = ut(:);
%     norm(ur-ur0)
%     norm(us-us0)    
%     norm(ut-ut0)
%     keyboard

    divU = -(ur + us + ut);
    
    % Define field differences at faces
    tau = 0;
    fluxp = zeros((N+1)^2,Nfaces);
    fluxu = zeros((N+1)^2,Nfaces);
    fluxv = zeros((N+1)^2,Nfaces);
    fluxw = zeros((N+1)^2,Nfaces);
    for f = 1:Nfaces
        dpfq = (Vq1D*reshape(dpK(:,f),N+1,N+1))*Vq1D';
        ufq  = (Vq1D*reshape(uavgK(:,f),N+1,N+1))*Vq1D';
        vfq  = (Vq1D*reshape(vavgK(:,f),N+1,N+1))*Vq1D';
        wfq  = (Vq1D*reshape(wavgK(:,f),N+1,N+1))*Vq1D';        
        ndotUavg = nxK(:,f).*ufq(:) + nyK(:,f).*vfq(:) + nzK(:,f).*wfq(:);
        
        ufq = (Vq1D*reshape(duK(:,f),N+1,N+1))*Vq1D';
        vfq = (Vq1D*reshape(dvK(:,f),N+1,N+1))*Vq1D';
        wfq = (Vq1D*reshape(dwK(:,f),N+1,N+1))*Vq1D';        
        ndotdU = nxK(:,f).*ufq(:) + nyK(:,f).*vfq(:) + nzK(:,f).*wfq(:);
        
        fp = (tau*dpfq(:) - ndotUavg).*sJK(:,f);
        fU = (tau*ndotdU(:) - dpfq(:)).*sJK(:,f);
        
        fluxp(:,f) = reshape((Pq1D*reshape(fp,Nq,Nq))*Pq1D',(N+1)^2,1);
        fluxu(:,f) = reshape((Pq1D*reshape(fU.*nxK(:,f),Nq,Nq))*Pq1D',(N+1)^2,1);
        fluxv(:,f) = reshape((Pq1D*reshape(fU.*nyK(:,f),Nq,Nq))*Pq1D',(N+1)^2,1);
        fluxw(:,f) = reshape((Pq1D*reshape(fU.*nzK(:,f),Nq,Nq))*Pq1D',(N+1)^2,1);               
    end
    
    rp =  -divU + .5*LIFT*(fluxp(:));
    ru =  -dpdx + .5*LIFT*(fluxu(:));
    rv =  -dpdy + .5*LIFT*(fluxv(:));
    rw =  -dpdz + .5*LIFT*(fluxw(:));

%     rp =  -divU;    ru =  -dpdx;
%     rv =  -dpdy;    rw =  -dpdz;
%     
%     rp = .5*LIFT*(fluxp(:));  ru = .5*LIFT*(fluxu(:));
%     rv = .5*LIFT*(fluxv(:));  rw = .5*LIFT*(fluxw(:));
    
    rhsp(:,e) = matvec(Pq1D,matvec(Vq1D,rp,'all')./Jq(:,e),'all');
    rhsu(:,e) = matvec(Pq1D,matvec(Vq1D,ru,'all')./Jq(:,e),'all');
    rhsv(:,e) = matvec(Pq1D,matvec(Vq1D,rv,'all')./Jq(:,e),'all');
    rhsw(:,e) = matvec(Pq1D,matvec(Vq1D,rw,'all')./Jq(:,e),'all');

end


% apply 1D tensor product matvec
function b = matvec(A,x,dim)

N = round((length(x(:)))^(1/3));
M = size(A,1);
x = reshape(x,N,N,N);

if strcmp(dim,'all')
    % apply in all dimensions
    
    b1 = zeros(M,M,N); % apply xy
    b = zeros(M,M,M); % apply z to xy result
    for i = 1:N
        b1(:,:,i) = (A*reshape(x(:,:,i),N,N))*A';
    end
    for i = 1:M
        b(i,:,:) = reshape(b1(i,:,:),M,N)*A';
    end
    
elseif dim==1  % apply only to one dimension
    b = zeros(M,N,N); % apply z to xy result
    for i = 1:N
        b(:,:,i) = reshape(x(:,:,i),N,N)*A';
    end
elseif dim==2
    b = zeros(N,M,N);
    for i = 1:N
        b(:,:,i) = A*reshape(x(:,:,i),N,N);
    end
elseif dim==3
    b = zeros(N,N,M);
    for i = 1:N
        b(i,:,:) = reshape(x(i,:,:),N,N)*A';
    end
end

b = b(:);

return

