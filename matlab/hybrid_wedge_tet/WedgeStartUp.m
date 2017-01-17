% interp/control nodes
[r s t] = wedge_nodes(N);
Np = length(r);

% plotting nodes
[rp sp tp] = wedge_equi_nodes(25);
Nplot = length(rp);
Vp = wedge_sem_basis(N,rp,sp,tp);

% "exact" cubature 
[rq sq tq wq] = wedge_cub(N);
Nq = length(rq);

Nptri = (N+1)*(N+2)/2;
Nfp = 3*(N+1)^2 + 2*Nptri;  % 3 quad faces, 2 triangle faces

% define face nodes
NODETOL = 1e-5;
Nfaces = 5;
Fmask = cell(Nfaces,1);
Fmask{1} = find(abs(t+1)<NODETOL); % t = -1
Fmask{2} = find(abs(s+1)<NODETOL); % s = -1
Fmask{3} = find(abs(r+s)<NODETOL); % r+s = 0
Fmask{4} = find(abs(r+1)<NODETOL); % r = -1
Fmask{5} = find(abs(t-1)<NODETOL); % t = 1
FmaskAll = [Fmask{1};Fmask{2};Fmask{3};Fmask{4};Fmask{5}];

[rf sf tf wf fids] = wedge_surface_cubature(N);
Nfq = length(rf);

%% reference matrices

[~, Dr, Ds, Dt] = wedge_sem_basis(N,r,s,t);

[r2D s2D] = Nodes2D(N); [r2D s2D] = xytors(r2D,s2D);
VWB = Vandermonde2D(N,r2D,s2D);
[Vrtri Vstri] = GradVandermonde2D(N,r2D,s2D);
Drtri = Vrtri/VWB; Dstri = Vstri/VWB;

[t1D] = JacobiGL(0,0,N); 
VSEM = Vandermonde1D(N,t1D(:));
V1D = Vandermonde1D(N,t1D)/VSEM; % sem in vertical direction
Dt1D = GradVandermonde1D(N,t1D)/VSEM; % sem in vertical direction

% reference quadrature VDMs
Vq = wedge_sem_basis(N,rq,sq,tq);
for f = 1:Nfaces
    rff = rf(fids{f}); sff = sf(fids{f}); tff = tf(fids{f});
    Vfq{f} = wedge_sem_basis(N,rff,sff,tff);
end


%% map nodes

% interp nodes
x = zeros(Np,K); y = zeros(Np,K); z = zeros(Np,K);

% geofacs
J = zeros(Np,K);
rx = zeros(Np,K); ry = zeros(Np,K); rz = zeros(Np,K);
sx = zeros(Np,K); sy = zeros(Np,K); sz = zeros(Np,K);
tx = zeros(Np,K); ty = zeros(Np,K); tz = zeros(Np,K);

% cubature nodes
xq = zeros(Nq,K); yq = zeros(Nq,K); zq = zeros(Nq,K); wJ = zeros(Nq,K);

% face nodes
xf = zeros(Nfq,K); yf = zeros(Nfq,K); zf = zeros(Nfq,K); 

% plotting nodes
xp = zeros(Nplot,K); yp = zeros(Nplot,K); zp = zeros(Nplot,K);

nxf = zeros(length(rf),1); nyf = zeros(length(rf),1); nzf = zeros(length(rf),1);
nx = zeros(length(rf),1); ny = zeros(length(rf),1); nz = zeros(length(rf),1);
Fscale = zeros(length(rf),1);
sJ = zeros(length(rf),1);
wsJ = zeros(length(rf),1);
for e = 1:K
    
    % nodal 
    v = EToV(e,:);
    [xe ye ze Je geo] = wedge_map(r,s,t,VX(v),VY(v),VZ(v));
    x(:,e) = xe; y(:,e) = ye; z(:,e) = ze;
    J(:,e) = Je;
    rx(:,e) = geo.rx; ry(:,e) = geo.ry; rz(:,e) = geo.rz;
    sx(:,e) = geo.sx; sy(:,e) = geo.sy; sz(:,e) = geo.sz;
    tx(:,e) = geo.tx; ty(:,e) = geo.ty; tz(:,e) = geo.tz;    
    
    % face quadrature nodes
    [xfe yfe zfe Jf geof] = wedge_map(rf,sf,tf,VX(v),VY(v),VZ(v));
    xf(:,e) = xfe; yf(:,e) = yfe; zf(:,e) = zfe;    
    
    % normals
    rxf = geof.rx; ryf = geof.ry; rzf = geof.rz;
    sxf = geof.sx; syf = geof.sy; szf = geof.sz;
    txf = geof.tx; tyf = geof.ty; tzf = geof.tz;    
    
    % faces 1-5
    f = 1;
    nxf(fids{f}) = -txf(fids{f}); nyf(fids{f}) = -tyf(fids{f}); nzf(fids{f}) = -tzf(fids{f});
    f = 4; % r = -1
    nxf(fids{f}) = -rxf(fids{f}); nyf(fids{f}) = -ryf(fids{f}); nzf(fids{f}) = -rzf(fids{f});            
    f = 3; % r+s = 0
    nxf(fids{f}) = rxf(fids{f}) + sxf(fids{f}); 
    nyf(fids{f}) = ryf(fids{f}) + syf(fids{f}); 
    nzf(fids{f}) = rzf(fids{f}) + szf(fids{f});        
    f = 2; % s = -1
    nxf(fids{f}) = -sxf(fids{f}); nyf(fids{f}) = -syf(fids{f}); nzf(fids{f}) = -szf(fids{f});
    
    f = 5; 
    nxf(fids{f}) = txf(fids{f});  nyf(fids{f}) = tyf(fids{f});  nzf(fids{f}) = tzf(fids{f});
    
    sJf = sqrt(nxf.^2 + nyf.^2 + nzf.^2);

    nx(:,e) = nxf./sJf; ny(:,e) = nyf./sJf; nz(:,e) = nzf./sJf;
    Fscale(:,e) = sJf; 
    sJ(:,e) = sJf.*Jf; 
    wsJ(:,e) = wf(:).*sJ(:,e);
        
    % cubature nodes
    [xe ye ze Jq] = wedge_map(rq,sq,tq,VX(v),VY(v),VZ(v));
    xq(:,e) = xe; yq(:,e) = ye; zq(:,e) = ze; 
    wJ(:,e) = wq(:).*Jq(:);
    
    % plotting nodes
    [xe ye ze] = wedge_map(rp,sp,tp,VX(v),VY(v),VZ(v));
    xp(:,e) = xe; yp(:,e) = ye; zp(:,e) = ze;
    
end

% h = color_line3(xq,yq,zq,wJ,'.'); set(h,'markersize',32); return
%  h = color_line3(xfq,yfq,zfq,wsJ,'.'); set(h,'markersize',32); return

%% lift operators

FmaskQuad = [Fmask{2}; Fmask{3}; Fmask{4}];
LIFTtri = cell(K,1);
LIFTquad = cell(K,1);
% LIFTtri = zeros(Nptri*K,Nptri);
% LIFTquad = zeros(Nptri*K,(N+1));

% build element lift matrices
for e = 1:K
    
    if 0 % if use face-by-face lift
        for f = 1:Nfaces
            
            M = Vq'*diag(wJ(:,e))*Vq;
            Mf = Vfq{f}'*diag(wsJ(fids{f},e))*Vfq{f};
            LIFTf = M\Mf;
            LIFTf(abs(LIFTf)<1e-8) = 0;
            LIFTf = LIFTf(:,Fmask{f});
            faceLIFT{f,e} = LIFTf;
            Mface{f} = Mf;
            keyboard
            if f==1
                %                 LIFT{f,e} = LIFTf(1:Nptri,1:Nptri);
                %                 LIFTtri{e} = LIFTf(1:Nptri,1:Nptri);
                LIFTtri((1:Nptri) + Nptri*(e-1),1:Nptri) = LIFTf(1:Nptri,1:Nptri);
            elseif f==2 || f==3 || f==4
%                 LIFT{f,e} = LIFTf(1:Nptri,1:(N+1));
                %LIFTquad{e} = [LIFTquad{e} LIFTf(1:Nptri,1:(N+1))];                
            elseif f==5
                LIFT{f,e} = sJ(fids{5}(1))/sJ(fids{1}(1)); % LIFTf(Np-Nptri+(1:Nptri),:);                 
                triScale(e) = sJ(fids{5}(1))/sJ(fids{1}(1));
            end            
            %     xff = xf(fids{f},e); yff = yf(fids{f},e); zff = zf(fids{f},e); sJf = sJ(fids{f},e);
            %     nxf = nx(fids{f},e); nyf = ny(fids{f},e); nzf = nz(fids{f},e);
            %     h = color_line3(x(:,e),y(:,e),z(:,e),J,'.'); set(h,'markersize',32); hold on;
            %     h2 = color_line3(xff,yff,zff,sJf,'.'); set(h2,'markersize',64)
            %     quiver3(xff,yff,zff,nxf,nyf,nzf)
            %     title(sprintf('Face %d',f))
            %     keyboard
        end
        keyboard
    else  % full lift matrices
        
        Vsq = wedge_sem_basis(N,rf,sf,tf);
        Ms = Vsq'*diag(wsJ(:,e))*Vsq; % surface mass
        M = (Vq'*diag(wJ(:,e))*Vq);
        LIFTfull = M\Ms;
        LIFTfull(abs(LIFTfull)<1e-6) = 0;
        LIFT{e} = [LIFTfull(:,Fmask{1}) LIFTfull(:,Fmask{2}) LIFTfull(:,Fmask{3}) LIFTfull(:,Fmask{4}) LIFTfull(:,Fmask{5})];              

% %         LL = 0*LIFTfull;
%         for f = 1:5
%             Lf = M\(Vfq{f}'*diag(wsJ(fids{f},e))*Vfq{f});
%             if f==1 || f==5
%                 LIFTf{f,e} = Lf(1:Nptri,Fmask{f}); % tri face 
%             else
%                 LIFTf{f,e} = Lf(1:Nptri,Fmask{f}(1:N+1)); % quad face = block diag
%             end
% %             LL(:,Fmask{f}) = LL(:,Fmask{f}) + Lf(:,Fmask{f});
%         end
% %         keyboard
%         triScale(e) = sJ(fids{5}(1))/sJ(fids{1}(1));
    end
    
end

%h = color_line3(rf,sf,tf,sJ,'.'); set(h,'markersize',64)

% store only one value per face
nxf = zeros(Nfaces,K); nyf = zeros(Nfaces,K); nzf = zeros(Nfaces,K);
for f = 1:5
   nxf(f,:) = nx(fids{f}(1),:); 
   nyf(f,:) = ny(fids{f}(1),:); 
   nzf(f,:) = nz(fids{f}(1),:);  
end
% store per-node nxyz
nx = zeros(Nfp,K); ny = zeros(Nfp,K); nz = zeros(Nfp,K);
off = 0;
nx(off + (1:Nptri),:)   = repmat(nxf(1,:),Nptri,1); off = off + Nptri;
nx(off + (1:(N+1)^2),:) = repmat(nxf(2,:),(N+1)^2,1); off = off + (N+1)^2;
nx(off + (1:(N+1)^2),:) = repmat(nxf(3,:),(N+1)^2,1); off = off + (N+1)^2;
nx(off + (1:(N+1)^2),:) = repmat(nxf(4,:),(N+1)^2,1); off = off + (N+1)^2;
nx(off + (1:Nptri),:)   = repmat(nxf(5,:),Nptri,1); 

off = 0;
ny(off + (1:Nptri),:)   = repmat(nyf(1,:),Nptri,1); off = off + Nptri;
ny(off + (1:(N+1)^2),:) = repmat(nyf(2,:),(N+1)^2,1); off = off + (N+1)^2;
ny(off + (1:(N+1)^2),:) = repmat(nyf(3,:),(N+1)^2,1); off = off + (N+1)^2;
ny(off + (1:(N+1)^2),:) = repmat(nyf(4,:),(N+1)^2,1); off = off + (N+1)^2;
ny(off + (1:Nptri),:)   = repmat(nyf(5,:),Nptri,1); 

off = 0;
nz(off + (1:Nptri),:)   = repmat(nzf(1,:),Nptri,1); off = off + Nptri;
nz(off + (1:(N+1)^2),:) = repmat(nzf(2,:),(N+1)^2,1); off = off + (N+1)^2;
nz(off + (1:(N+1)^2),:) = repmat(nzf(3,:),(N+1)^2,1); off = off + (N+1)^2;
nz(off + (1:(N+1)^2),:) = repmat(nzf(4,:),(N+1)^2,1); off = off + (N+1)^2;
nz(off + (1:Nptri),:)   = repmat(nzf(5,:),Nptri,1); 


%% matching nodes

vnodeids = reshape(1:K*Np,Np,K);
vmapM = zeros(Nfp,K);
vmapP = zeros(Nfp,K);
foff{1} = 0;
for e = 1:K
    off = 0;
    for f = 1:Nfaces
        vmapM(off + (1:length(Fmask{f})),e) = vnodeids(Fmask{f},e);
        off = off + length(Fmask{f});
        foff{f+1} = off;
    end
end

for e = 1:K
%     foff = 0;
    for f = 1:Nfaces
        enbr = EToE(e,f);
        fnbr = EToF(e,f);
        
        if (enbr==e) % boundary node
            idM = 1:length(Fmask{f});
            vmapP(idM + foff{f},e) = vmapM(idM + foff{f},enbr);            
        else                                    
            vidM = vnodeids(Fmask{f},e); vidM = vidM(:);            
            vidP = vnodeids(Fmask{fnbr},enbr); vidP = vidP(:);
                        
%             xM = x(vidM); yM = y(vidM); zM = z(vidM);
%             xP = x(vidP); yP = y(vidP); zP = z(vidP);                        
%             clf; hold on
%             drawWedgeMesh;
%             ids = [1 2 3 1 4 5 6 4 5 2 3 6];
%             for ee = 1:K
%                 v = EToV(ee,:);
%                 text(VX(v)+.1,VY(v)+.1,VZ(v),num2str((1:length(v))'))
%                 v = v(ids);
%                 plot3(VX(v),VY(v),VZ(v),'ko-','linewidth',2);
%             end
%             plot3(xM,yM,zM,'o');hold on
%             plot3(xP,yP,zP,'*')            
% %             keyboard
            
            tmp = ones(1,length(Fmask{f}));
            xM = x(vidM)*tmp; yM = y(vidM)*tmp; zM = z(vidM)*tmp;
            xP = x(vidP)*tmp; yP = y(vidP)*tmp; zP = z(vidP)*tmp;
                                   
            % Compute distance matrix
            D = (xM -xP').^2 + (yM-yP').^2 + (zM-zP').^2;
                        
            [idM, idP] = find(abs(D)<NODETOL);            
            vmapP(idM + foff{f},e) = vmapM(idP + foff{fnbr},enbr);
        end
%         foff = foff + length(Fmask{f}); % face node number offset
%         foff
    end
end
% vmapM = reshape(vmapM,Nfp,K);
% vmapP = reshape(vmapP,Nfp,K);
vmapM = vmapM(:);
vmapP = vmapP(:);

% Create list of boundary nodes
mapB = find(vmapP(:)==vmapM(:));
vmapB = vmapM(mapB);