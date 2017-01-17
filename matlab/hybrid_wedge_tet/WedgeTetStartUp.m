NODETOL=1e-6;

% # dofs for each elem type
NpH = (N+1)^3;
NpW = (N+1)^2*(N+2)/2;
NpP = (N+1)*(N+2)*(2*N+3)/6;
NpT = (N+1)*(N+2)*(N+3)/6;

%% refernece vol/surface 

Nq = N; % cubature order if we want to overintegrate...

% tet face vertex orderings
fvT{1} = [1 3 2]; fvT{2} = [1 2 4]; fvT{3} = [3 1 4]; fvT{4} = [2 3 4];
r1 = [-1  1 -1 -1 ]'; s1 = [-1 -1  1 -1 ]'; t1 = [-1 -1 -1  1 ]';
[rfT sfT tfT wfT fidsT] = surface_cubature(Nq,r1,s1,t1,fvT);
rstT = [r1 s1 t1];
NfcT = length(wfT);

% wedge face vertex ordering
fvW{1} = [1 3 2]; fvW{2} = [4 5 6];
fvW{3} = [1 2 5 4]; fvW{4} = [3 1 4 6]; fvW{5} = [2 3 6 5];
u = [-1 1 -1 -1 1 -1]; v = [-1 -1 1 -1 -1 1]; w = [-1 -1 -1 1 1 1];
r1 = v(:); s1 = w(:); t1 = u(:); % flipping coordinates for Gmsh
[rfW sfW tfW wfW fidsW] = surface_cubature(Nq,r1,s1,t1,fvW);
rstW = [r1 s1 t1];
NfcW = length(wfW);

% ====================== reference volume stuff =========================

[rqT sqT tqT w] = tet_cubature(2*Nq); % integrate order 2N
wT = w; NcT = length(w(:));
[VT VTr VTs VTt] = tet_basis(N,rqT,sqT,tqT);
VTf = tet_basis(N,rfT,sfT,tfT);

% get nodes
[xn yn zn] = Nodes3D(N); [rT sT tT] = xyztorst(xn,yn,zn);
[VTnodal VTrn VTsn VTtn] = tet_basis(N,rT,sT,tT);
% compute nodal VDMs
VT = VT/VTnodal; VTr = VTrn/VTnodal;  VTs = VTsn/VTnodal;  VTt = VTtn/VTnodal;
VTf = VTf/VTnodal; % quadrature VDM

invMThat = VTnodal*VTnodal'; MThat = inv(invMThat);

% wedge VDMs
% [rqW sqW tqW w] = wedge_cubature(Nq);
[rqW tqW sqW w] = wedge_cub(Nq); % USE REDUCED QUADRATURE FOR WEDGE VOL
wW = w; NcW = length(w(:));
[VW VWr VWs VWt] = wedge_basis(N,rqW,sqW,tqW);
VWf = wedge_basis(N,rfW,sfW,tfW);

% get wedge nodes
[rW sW tW] = wedge_nodes(N); 
% [rW tW sW] = wedge_nodes_perm(N); 
[VWnodal VWrn VWsn VWtn] = wedge_basis(N,rW,sW,tW);

% make nodal quadrature-ops
VW = VW/VWnodal; 
VWr = VWr/VWnodal; 
VWs = VWs/VWnodal; 
VWt = VWt/VWnodal;
VWf = VWf/VWnodal;

% make quadrature-free operators
DWr = VWrn/VWnodal; DWs = VWsn/VWnodal; DWt = VWtn/VWnodal;

if nodalLIFT    
    
    % WARNING: HACK - replace surface quadrature nodes/weights with nodal
    % points and weights set artificially to = 1. re-compute surface ops and
    % surface quad points to compute local wedge LIFT matrices.  
    
    % ======================== tet LIFT matrix ==========================
    
    fmask1   = find( abs(1+tT) < NODETOL)';
    fmask2   = find( abs(1+sT) < NODETOL)';
    fmask3   = find( abs(1+rT) < NODETOL)';
    fmask4   = find( abs(1+rT+sT+tT) < NODETOL)';
    FmaskT  = [fmask1;fmask2;fmask3;fmask4]';
    SmaskT = FmaskT(:);
    
    NfpT = (N+1)*(N+2)/2;
    Emat = zeros(NpT, 4*NfpT);
    for face=1:4
        % process face
        if(face==1); faceR = rT(FmaskT(:,1)); faceS = sT(FmaskT(:,1)); end;
        if(face==2); faceR = rT(FmaskT(:,2)); faceS = tT(FmaskT(:,2)); end;
        if(face==3); faceR = sT(FmaskT(:,3)); faceS = tT(FmaskT(:,3)); end;
        if(face==4); faceR = sT(FmaskT(:,4)); faceS = tT(FmaskT(:,4)); end;
        VFace = Vandermonde2D(N, faceR, faceS);
        massFace = inv(VFace*VFace');
        idr = FmaskT(:,face); idc = (face-1)*NfpT+1:face*NfpT;
        Emat(idr, idc) = Emat(idr, idc)+ massFace;
    end
    % inv(mass matrix)*\I_n (L_i,L_j)_{edge_n}
    LIFTT = invMThat*Emat;   
        
    rfT = rT(FmaskT(:)); sfT = sT(FmaskT(:)); tfT = tT(FmaskT(:));    
    ids = 1:NfpT;
    off = 0;
    for f = 1:4
        fidsT{f} = ids + off; 
        off = off + NfpT;  %FmaskT(:,f);
    end
    NfcT = length(FmaskT(:)); wfT = ones(size(rfT));    
    VTf = tet_basis(N,rfT,sfT,tfT)/VTnodal;
    VTf(abs(VTf)<1e-8) = 0;
    
    % ======================== wedge LIFT matrix =========================
    
    FmaskW{1} = find( abs(1+sW) < NODETOL)';
    FmaskW{2} = find( abs(1-sW) < NODETOL)';
    FmaskW{3} = find( abs(1+rW) < NODETOL)';
    FmaskW{4} = find( abs(1+tW) < NODETOL)';
    FmaskW{5} = find( abs(rW+tW) < NODETOL)';
    SmaskW = [FmaskW{1} FmaskW{2} FmaskW{3} FmaskW{4} FmaskW{5}];              
    
    % save quadrature points for wedges
    rfqW = rfW; sfqW = sfW; tfqW = tfW; wfqW = wfW; 
    fidsqW = fidsW;
    
    rfW = rW(SmaskW); sfW = sW(SmaskW); tfW = tW(SmaskW);
    
    % first tri, then quad faces
    off = 0;
    for f = 1:2
        fidsW{f} = (1:(N+1)*(N+2)/2) + off;
        off = off + (N+1)*(N+2)/2;
    end
    for f = 3:5
        fidsW{f} = (1:(N+1)^2) + off;
        off = off + (N+1)^2;
    end
    NfcW = length(SmaskW(:)); wfW = ones(size(rfW));
%     VWf = wedge_basis(N,rfW,sfW,tfW)/VWnodal;
%     VWf(abs(VWf)<1e-8) = 0;
end


%% original hybrid mesher stuff

% hex face vertex ordering
fvH{1} = [1 4 3 2]; fvH{2} = [1 2 6 5]; fvH{3} = [2 3 7 6];
fvH{4} = [3 4 8 7]; fvH{5} = [4 1 5 8]; fvH{6} = [5 6 7 8];
r1 = [-1  1  1 -1 -1  1  1 -1]'; s1 = [-1 -1  1  1 -1 -1  1  1]';
t1 = [-1 -1 -1 -1  1  1  1  1]';
[rfH sfH tfH wfH fidsH] = surface_cubature(Nq,r1,s1,t1,fvH);
rstH = [r1 s1 t1];
NfcH = length(wfH);

% pyramid face vertex ordering
fvP{1} = [1 2 5]; fvP{2} = [4 1 5];
fvP{3} = [2 3 5]; fvP{4} = [3 4 5];
fvP{5} = [1 4 3 2];
r1 = [ -1   1   1  -1  -1 ]';s1 = [ -1  -1   1   1  -1 ]';
t1 = [ -1  -1  -1  -1   1 ]';
[rfP sfP tfP wfP fidsP] = surface_cubature(Nq,r1,s1,t1,fvP);
rstP = [r1 s1 t1];
NfcP = length(wfP);

[rqP sqP tqP w] = pyr_cubature(Nq);
wP = w; NcP = length(w(:));
[VP VPr VPs VPt] = pyr_basis(N,rqP,sqP,tqP);
VPf = pyr_basis(N,rfP,sfP,tfP);

[rqH sqH tqH w] = hex_cubature(Nq);
wH = w; NcH = length(w(:));
[VH VHr VHs VHt] = hex_basis(N,rqH,sqH,tqH);
VHf = hex_basis(N,rfH,sfH,tfH);

%%

NcMax = max(NcH,NcW);
wJ = zeros(NcMax,K);
x = zeros(NcMax,K); y = zeros(NcMax,K); z = zeros(NcMax,K);
rx = zeros(NcMax,K); ry = zeros(NcMax,K); rz = zeros(NcMax,K);
sx = zeros(NcMax,K); sy = zeros(NcMax,K); sz = zeros(NcMax,K);
tx = zeros(NcMax,K); ty = zeros(NcMax,K); tz = zeros(NcMax,K);

NfcMax = NfcH;
xf = zeros(NfcMax,K); yf = zeros(NfcMax,K); zf = zeros(NfcMax,K);
nx = zeros(NfcMax,K); ny = zeros(NfcMax,K); nz = zeros(NfcMax,K);
Fscale = zeros(NfcMax,K); wsJ = zeros(NfcMax,K);  Jf = zeros(NfcMax,K);
Nv = zeros(K,1);

JsB = zeros(NcMax,K); % for eigenvalue bounds only

% for LSC
NpMax = max([NpH NpW NpP NpT]);
invMref = zeros(NpMax,K);
J = zeros(NcMax,K); Jsurf = zeros(NfcMax,K);
Jr = zeros(NcMax,K); Js = zeros(NcMax,K); Jt = zeros(NcMax,K);

for e = 1:K
    disp(sprintf('Working on elem %d out of %d',e,K))
    v = EToV(e,:); v = v(v > 0); NvK = nnz(v);
    
    Nv(e) = NvK;
    switch NvK
        case 4 % tet
            [xq,yq,zq,rxq,sxq,txq,ryq,syq,tyq,rzq,szq,tzq,Jq] = ...
                tet_geom_factors(VX(v),VY(v),VZ(v),rqT,sqT,tqT);
            
            [xfq,yfq,zfq,rxf,sxf,txf,ryf,syf,tyf,rzf,szf,tzf,Jf] = ...
                tet_geom_factors(VX(v),VY(v),VZ(v),rfT,sfT,tfT);
            
            w = wT; wf = wfT; rst = rstT; fv = fvT; fids = fidsT;
            
        case 5 % pyr            
            [xq,yq,zq,rxq,sxq,txq,ryq,syq,tyq,rzq,szq,tzq,Jq] = ...
                pyr_geom_factors(VX(v),VY(v),VZ(v),rqP,sqP,tqP);
            
            [xfq,yfq,zfq,rxf,sxf,txf,ryf,syf,tyf,rzf,szf,tzf,Jf] = ...
                pyr_geom_factors(VX(v),VY(v),VZ(v),rfP,sfP,tfP);
            
            w = wP; wf = wfP; rst = rstP; fv = fvP; fids = fidsP;
            
        case 6 % wedge
            
            [xq,yq,zq,rxq,sxq,txq,ryq,syq,tyq,rzq,szq,tzq,Jq] = ...
                wedge_geom_factors(VX(v),VY(v),VZ(v),rqW,sqW,tqW);
            
            [xfq,yfq,zfq,rxf,sxf,txf,ryf,syf,tyf,rzf,szf,tzf,Jf] = ...
                wedge_geom_factors(VX(v),VY(v),VZ(v),rfW,sfW,tfW);
            
            w = wW; wf = wfW; rst = rstW; fv = fvW; fids = fidsW;
            
            % save Jr, Js, Jt at quad points for LSC-wedge
            [VW1 VWr1 VWs1 VWt1] = wedge_basis(1,rqW,sqW,tqW);
            VWfN = wedge_basis(N,rfW,sfW,tfW);
            
%             J(1:NcW,e) = Jq;
            
            % local projection of Jq onto higher order space.
            invMWref = 1./diag((VW'*spdiag(w)*VW)); % should be = I for modal
            cJ = invMWref.*(VW'*(wW.*Jq));
            
            Jsurf(1:NfcW,e) = VWfN*cJ;
            Jr(1:NcW,e) = VWr*cJ; Js(1:NcW,e) = VWs*cJ; Jt(1:NcW,e) = VWt*cJ;
            
            
        case 8 % hex
            
            [xq,yq,zq,rxq,sxq,txq,ryq,syq,tyq,rzq,szq,tzq,Jq] = ...
                hex_geom_factors(VX(v),VY(v),VZ(v),rqH,sqH,tqH);
            
            [xfq,yfq,zfq,rxf,sxf,txf,ryf,syf,tyf,rzf,szf,tzf,Jf] = ...
                hex_geom_factors(VX(v),VY(v),VZ(v),rfH,sfH,tfH);
            
            w = wH; wf = wfH; rst = rstH; fv = fvH; fids = fidsH;
            
    end
    
    if min(Jq(:))<0
        disp('Negative Jacobian detected')
        drawMesh(e)
        v = EToV(e,:); v = v(v > 0); NvK = nnz(v);
        hold on;text(VX(v)+.05,VY(v),VZ(v),num2str((1:NvK)'));
        keyboard
    end
    
    NcK = length(xq(:)); NfcK = length(xfq);
    
    x(1:NcK,e) = xq(:); y(1:NcK,e) = yq(:); z(1:NcK,e) = zq(:);
    J(1:NcK,e) = Jq; 
    wJ(1:NcK,e) = w(1:NcK).*Jq;
    
    rx(1:NcK,e) = rxq(:); ry(1:NcK,e) = ryq(:); rz(1:NcK,e) = rzq(:);
    sx(1:NcK,e) = sxq(:); sy(1:NcK,e) = syq(:); sz(1:NcK,e) = szq(:);
    tx(1:NcK,e) = txq(:); ty(1:NcK,e) = tyq(:); tz(1:NcK,e) = tzq(:);
    
    xf(1:NfcK,e) = xfq(:); yf(1:NfcK,e) = yfq(:); zf(1:NfcK,e) = zfq(:);
    
    % get normals
    Nfaces = length(fids);
    ref_face_area = zeros(NfcMax,1);
    for f = 1:Nfaces
        
        % assumes face nodes ordered w.r.t. right hand rule for normal
        v1 = rst(fv{f}(2),:) - rst(fv{f}(1),:);
        v2 = rst(fv{f}(3),:) - rst(fv{f}(2),:);
        nc = cross(v1,v2); % norm of nc = parallelogram area
        
        % scale by area of face
        if length(fv{f}) == 3 % if triangle, halve parallelogram area
            nc = nc/2;
        end
        
        ref_face_area(fids{f}) = norm(nc);
        
        nx(fids{f},e) = nc(1)*rxf(fids{f}) + nc(2)*sxf(fids{f}) + nc(3)*txf(fids{f});
        ny(fids{f},e) = nc(1)*ryf(fids{f}) + nc(2)*syf(fids{f}) + nc(3)*tyf(fids{f});
        nz(fids{f},e) = nc(1)*rzf(fids{f}) + nc(2)*szf(fids{f}) + nc(3)*tzf(fids{f});
        
        % normalize so we can scale by sJ after
        wf(fids{f}) = wf(fids{f})/sum(wf(fids{f}));
        
    end
    
    sJK = sqrt(nx(:,e).^2 + ny(:,e).^2 + nz(:,e).^2);
    nx(:,e) = nx(:,e)./sJK; ny(:,e) = ny(:,e)./sJK; nz(:,e) = nz(:,e)./sJK;
    Fscale(:,e) = sJK; % area of face for planar    
    
    wsJ(1:NfcK,e) = wf(:).*sJK(1:NfcK).*Jf(:); % premultiply jacobian with weights
    JsB(1:NfcK,e) = sJK(1:NfcK).*Jf(:)./ref_face_area(1:NfcK); % Jf = Jacobian on face = used bound
    sJ(1:NfcK,e) = sJK(1:NfcK).*Jf(:);    
end
tet_face_area = ref_face_area;
tetK = find(Nv==4); KT = length(tetK); wJT = wJ(1:NcT,tetK);
pyrK = find(Nv==5); KP = length(pyrK); wJP = wJ(1:NcP,pyrK);
wedgK = find(Nv==6); KW = length(wedgK); wJW = wJ(1:NcW,wedgK);
hexK = find(Nv==8); KH = length(hexK); wJH = wJ(1:NcH,hexK);

% % save to quadrature 
% nxq = nx; nyq = ny; nzq = nz;

% compute lift matrices using quadrature for wedge 
if nodalLIFT    
    KW = length(wedgK);
    LIFTW = cell(KW,1);
    
    w = wW; wf = wfW; rst = rstW; fv = fvW; fids = fidsW;
    
    for ee = 1:KW
        e = wedgK(ee);
        disp(sprintf('Computing wedge lift for elem %d out of %d',e,KW))
        v = EToV(e,:); v = v(v > 0); NvK = nnz(v);
        Nv(e) = NvK;
        
        % compute Dmatrices for each element
        MW{ee} = VW'*spdiag(wJ(1:NcW,e))*VW;        
        MW{ee}(abs(MW{ee})<1e-8) = 0;        
        Vx = diag(rx(1:NcW,e))*VWr + diag(sx(1:NcW,e))*VWs + diag(tx(1:NcW,e))*VWt;        
        Dx{ee} = MW{ee}\(VW'*spdiag(wJ(1:NcW,e))*Vx);
        Vy = diag(ry(1:NcW,e))*VWr + diag(sy(1:NcW,e))*VWs + diag(ty(1:NcW,e))*VWt;        
        Dy{ee} = MW{ee}\(VW'*spdiag(wJ(1:NcW,e))*Vy);
        Vz = diag(rz(1:NcW,e))*VWr + diag(sz(1:NcW,e))*VWs + diag(tz(1:NcW,e))*VWt;        
        Dz{ee} = MW{ee}\(VW'*spdiag(wJ(1:NcW,e))*Vz);
      
%         keyboard
        
        % compute lift        
        [xfq,yfq,zfq,rxf,sxf,txf,ryf,syf,tyf,rzf,szf,tzf,Jf] = ...
            wedge_geom_factors(VX(v),VY(v),VZ(v),rfqW,sfqW,tfqW);
                
        NcK = length(xq(:)); NfcK = length(xfq);
        xfqW(1:NfcK,e) = xfq; yfqW(1:NfcK,e) = yfq; zfqW(1:NfcK,e) = zfq;
        
        % get normals
        Nfaces = length(fidsqW);
        ref_face_area = zeros(NfcMax,1);
        for f = 1:Nfaces
            
            % assumes face nodes ordered w.r.t. right hand rule for normal
            v1 = rst(fv{f}(2),:) - rst(fv{f}(1),:);
            v2 = rst(fv{f}(3),:) - rst(fv{f}(2),:);
            nc = cross(v1,v2); % norm of nc = parallelogram area
            
            % scale by area of face
            if length(fv{f}) == 3 % if triangle, halve parallelogram area
                nc = nc/2;
            end
            
            ref_face_area(fidsqW{f}) = norm(nc);
            
            nxW(fidsqW{f},e) = nc(1)*rxf(fidsqW{f}) + nc(2)*sxf(fidsqW{f}) + nc(3)*txf(fidsqW{f});
            nyW(fidsqW{f},e) = nc(1)*ryf(fidsqW{f}) + nc(2)*syf(fidsqW{f}) + nc(3)*tyf(fidsqW{f});
            nzW(fidsqW{f},e) = nc(1)*rzf(fidsqW{f}) + nc(2)*szf(fidsqW{f}) + nc(3)*tzf(fidsqW{f});
            
            % normalize so we can scale by sJ after
            wf(fidsqW{f}) = wfqW(fidsqW{f})/sum(wfqW(fidsqW{f}));            
        end
                
        sJK = sqrt(nxW(:,e).^2 + nyW(:,e).^2 + nzW(:,e).^2); 
        nxW(:,e) = nxW(:,e)./sJK; 
        nyW(:,e) = nyW(:,e)./sJK; 
        nzW(:,e) = nzW(:,e)./sJK; % normalize        
        sJK = sJK(1:NfcK).*Jf(1:NfcK);
        wsJK = wfqW.*sJK; % includes face area scaling   
        
% keyboard
wsJW(1:size(wsJK,1),ee) = wsJK;         
        
        LIFTW{ee} = [];        
        Ms{ee} = []; Vf = zeros(length(rfqW),NpW);
        for f = 1:5
            VfqW = wedge_basis(N,rfqW(fidsqW{f}),sfqW(fidsqW{f}),tfqW(fidsqW{f}))/VWnodal;
            Mf = VfqW'*diag(wsJK(fidsqW{f}))*VfqW; % surface mass            
            LIFTf = MW{ee}\Mf(:,FmaskW{f}); %MW\Mf;            
            LIFTf(abs(LIFTf)<1e-6) = 0;
            LIFTW{ee} = [LIFTW{ee} LIFTf]; %LIFTf(:,FmaskW{f})];
            Ms{ee} = [Ms{ee} Mf(:,FmaskW{f})];
            %             Vf(fidsqW{f},:) = VfqW;
            
           
%             % testing TP wedge lift
%             if (0 && f==3) % quad face
%                 sJ_tp = reshape(sJK(fidsqW{f}),N+1,N+1); 
%                 [r1D, w1D] = JacobiGQ(0,0,N);                
%                 Vq1D = Vandermonde1D(N,r1D)/Vandermonde1D(N,JacobiGL(0,0,N));
%                 M1D1 = Vq1D'*diag(w1D)*Vq1D; 
%                 M1D2 = Vq1D'*diag(w1D.*sJ_tp(1,:)')*Vq1D; % extract 1D-varying sJ
%                 norm(kron(M1D1,M1D2)/4-Mf(FmaskW{f},FmaskW{f}),'fro') % divide by size of quad
%                 
%                 [rqtri,sqtri,wqtri] = Cubature2D(2*N+1); % match # qpts of vol quadrature
%                 [rtri, stri] = Nodes2D(N); [rtri, stri] = xytors(rtri,stri);
%                 Vqtri = Vandermonde2D(N,rqtri,sqtri)/Vandermonde2D(N,rtri,stri);
%                 Mtri = Vqtri'*diag(wqtri.*Jq(1:length(rqtri)))*Vqtri;
%                 norm(kron(M1D1,Mtri)-MW{1},'fro')
%                 
%                 NpTri = (N+1)*(N+2)/2;
%                 M1Df = zeros(NpTri,N+1);
%                 ids = find(abs(rtri+1)<1e-8); % fmasktri for f==3
%                 M1Df(ids,:) = M1D2;
%                 Lquad = (Mtri\M1Df)/4; % divide by size of quad face for scaling in sJ
%                 norm(Lquad-LIFTf(1:NpTri,1:N+1),'fro')
%                 keyboard                                 
%             end
        end
        
        %         Ms(abs(Ms)<1e-8) = 0;
        %         VfqW = wedge_basis(N,rfqW,sfqW,tfqW)/VWnodal;
%         Mf = VfqW'*diag(wsJK)*VfqW(:,SmaskW); % surface mass
%         
%         Mf(abs(Mf)<1e-8) = 0;
%         L = MW\Mf; L(abs(L)<1e-8) = 0;    
    end
end

%% connectivities

[mapM,mapP,mapB] = BuildHybridMaps();


%%
% cubature points for each element type
xqT = x(1:NcT,Nv==4); yqT = y(1:NcT,Nv==4); zqT = z(1:NcT,Nv==4);
xqP = x(1:NcP,Nv==5); yqP = y(1:NcP,Nv==5); zqP = z(1:NcP,Nv==5);
xqW = x(1:NcW,Nv==6); yqW = y(1:NcW,Nv==6); zqW = z(1:NcW,Nv==6);
xqH = x(1:NcH,Nv==8); yqH = y(1:NcH,Nv==8); zqH = z(1:NcH,Nv==8);

% diff volume terms for nodal tets
% if useNodalTets
    
xT = zeros(NpT,length(tetK)); yT = zeros(NpT,length(tetK)); zT = zeros(NpT,length(tetK));
rxT = zeros(NpT,length(tetK)); ryT = zeros(NpT,length(tetK)); rzT = zeros(NpT,length(tetK));
sxT = zeros(NpT,length(tetK)); syT = zeros(NpT,length(tetK)); szT = zeros(NpT,length(tetK));
txT = zeros(NpT,length(tetK)); tyT = zeros(NpT,length(tetK)); tzT = zeros(NpT,length(tetK));
JTet = zeros(NpT,length(tetK));

for ee = 1:length(tetK)
    e = tetK(ee);
    v = EToV(e,:); v = v(v > 0); NvK = nnz(v);
    
    [xTK,yTK,zTK,rxTK,sxTK,txTK,...
        ryTK,syTK,tyTK,rzTK,szTK,tzTK,JTK] = ...
        tet_geom_factors(VX(v),VY(v),VZ(v),rT,sT,tT);
    
    % nodal data
    xT(1:NpT,ee) = xTK; yT(1:NpT,ee) = yTK; zT(1:NpT,ee) = zTK;
    rxT(1:NpT,ee) = rxTK; ryT(1:NpT,ee) = ryTK; rzT(1:NpT,ee) = rzTK;
    sxT(1:NpT,ee) = sxTK; syT(1:NpT,ee) = syTK; szT(1:NpT,ee) = szTK;
    txT(1:NpT,ee) = txTK; tyT(1:NpT,ee) = tyTK; tzT(1:NpT,ee) = tzTK;
    JTet(1:NpT,ee) = JTK;
end    
% end

% nodal wedges
xW = zeros(NpW,length(wedgK)); yW = zeros(NpW,length(wedgK)); zW = zeros(NpW,length(wedgK));
rxW = zeros(NpW,length(wedgK)); ryW = zeros(NpW,length(wedgK)); rzW = zeros(NpW,length(wedgK));
sxW = zeros(NpW,length(wedgK)); syW = zeros(NpW,length(wedgK)); szW = zeros(NpW,length(wedgK));
txW = zeros(NpW,length(wedgK)); tyW = zeros(NpW,length(wedgK)); tzW = zeros(NpW,length(wedgK));

JW = zeros(NpW,length(wedgK));
for ee = 1:length(wedgK)
    e = wedgK(ee);
    v = EToV(e,:); v = v(v > 0); NvK = nnz(v);
    
    [xWK,yWK,zWK,rxWK,sxWK,txWK,...
        ryWK,syWK,tyWK,rzWK,szWK,tzWK,JWK] = ...
        wedge_geom_factors(VX(v),VY(v),VZ(v),rW,sW,tW);
    
    % nodal data
    xW(1:NpW,ee) = xWK; yW(1:NpW,ee) = yWK; zW(1:NpW,ee) = zWK;
    rxW(1:NpW,ee) = rxWK; ryW(1:NpW,ee) = ryWK; rzW(1:NpW,ee) = rzWK;
    sxW(1:NpW,ee) = sxWK; syW(1:NpW,ee) = syWK; szW(1:NpW,ee) = szWK;
    txW(1:NpW,ee) = txWK; tyW(1:NpW,ee) = tyWK; tzW(1:NpW,ee) = tzWK;    
    JW(1:NpW,ee) = JWK;
end

return
%% re-compute cubature for nodal error

%[rqW sqW tqW wW] = wedge_cubature(N+1);
[rqW tqW sqW wW] = wedge_cub(N+1);
VW = wedge_basis(N,rqW,sqW,tqW)/VWnodal;
xqW = VW*xW; 
yqW = VW*yW;
zqW = VW*zW;
NcW = length(rqW);
wJ(1:NcW,wedgK) = diag(wW)*(VW*JW);

% plot3(rqW,sqW,tqW,'o')
% keyboard

return

%% pre-invert mass matrix

NpMax = NpH; % hex has most dofs
NcMax = NcH;
NfcMax = NfcH;

invM = zeros(NpMax,K);
% if useNodalTets
    invM(1:NpT,tetK) = 1; % mass matrix computed using VDMs instead
% else
%     for e = 1:KT
%         invM(1:NpT,tetK(e)) = 1./diag(VT'*spdiag(wJ(1:NcT,tetK(e)))*VT);
%     end
% end

for e = 1:KP
    invM(1:NpP,pyrK(e)) = 1./diag(VP'*spdiag(wJ(1:NcP,pyrK(e)))*VP);
end

% if useLSC
%     for e = 1:KW
%         invM(1:NpW,wedgK(e)) = invMWref;
%     end
% else
    MW = cell(KW,1);
    for e = 1:KW
        %         invM(1:NpW,wedgK(e)) = 1./diag(VW'*spdiag(wJ(1:NcW,wedgK(e)))*VW);
        Me = VW'*spdiag(wJ(1:NcW,wedgK(e)))*VW;
        Me(abs(Me)<1e-6) = 0; 
        MW{e} = Me; %sparse(Me);
    end
% end

for e = 1:KH
    invM(1:NpH,hexK(e)) = 1./diag(VH'*spdiag(wJ(1:NcH,hexK(e)))*VH);
end

%% Compute LSC quantities for wedge

JW  = J(1:NcW,wedgK);  JWf = Jsurf(1:NfcW,wedgK);
sqJW = 1./sqrt(JW);    sqJWf = 1./sqrt(JWf);
JWr = Jr(1:NcW,wedgK); JWr = .5*JWr./JW; % d(u/sqJ)/dr = (dudr - u * Jr/(2*J))/sqJ
JWs = Js(1:NcW,wedgK); JWs = .5*JWs./JW; % JWr = Jr/(2*J);
JWt = Jt(1:NcW,wedgK); JWt = .5*JWt./JW;

%% Fine-sampled points for plotting
%
% Nplot = 3*N;
% rp1D = linspace(-1,1,Nplot);
%
% % hex
% [rpH spH tpH] = meshgrid(rp1D);
% rpH = rpH(:); spH = spH(:); tpH = tpH(:);
%
% % wedge
% [rpTri spTri] = EquiNodes2D(Nplot); [rpTri spTri] = xytors(rpTri,spTri);
% rpW = []; spW = []; tpW = [];
% for i = 1:Nplot
%     rpW = [rpW(:); rpTri(:)];    tpW = [tpW(:); spTri(:)];
%     spW = [spW(:); ones(size(rpTri))*rp1D(i)];
% end
%
% % pyr
% apP = []; bpP = []; cpP = [];
% c = linspace(-1,1,Nplot);
% for level = 1:Nplot
%     if level < Nplot
%         a1D = linspace(-1,1,Nplot-level);
%     else
%         a1D = 0;
%     end
%     [a b] = meshgrid(a1D);
%     apP = [apP; a(:)];  bpP = [bpP; b(:)];
%     cpP = [cpP; c(level)*ones(size(a(:)))];
% end
% [rpP spP tpP] = pyr_abctorst(apP,bpP,cpP);
%
% % tet
% [rpT spT tpT] = EquiNodes3D(Nplot);
%
% VpH = hex_basis(N,rpH,spH,tpH);  VpW = wedge_basis(N,rpW,spW,tpW);
% VpP = pyr_basis(N,rpP,spP,tpP);  VpT = tet_basis(N,rpT,spT,tpT);
%
% NplotMax = Nplot^3;
% xp = nan(NplotMax,K); yp = nan(NplotMax,K); zp = nan(NplotMax,K);
% sqJp = ones(NplotMax,K);
% for e = 1:K
%     sprintf('Working on plotting for elem %d out of %d',e,K)
%     v = EToV(e,:); v = v(v > 0); NvK = nnz(v);
%     Nv(e) = NvK;
%
%     switch NvK
%         case 4 % tet
%             [xpK,ypK,zpK] = tet_geom_factors(VX(v),VY(v),VZ(v),rpT,spT,tpT);
%             NplotK = length(rpT);
%         case 5 % pyr
%
%             [xpK,ypK,zpK] = pyr_geom_factors(VX(v),VY(v),VZ(v),rpP,spP,tpP);
%             NplotK = length(rpP);
%         case 6 % wedge
%
%             [xpK,ypK,zpK, ~,~,~,~,~,~,~,~,~,JpK] = ...
%                 wedge_geom_factors(VX(v),VY(v),VZ(v),rpW,spW,tpW);
%
%             NplotK = length(rpW);
%             sqJp(1:NplotK,e) = sqrt(JpK);
%
%         case 8 % hex
%             [xpK,ypK,zpK] = hex_geom_factors(VX(v),VY(v),VZ(v),rpH,spH,tpH);
%             NplotK = length(rpH);
%     end
%     xp(1:NplotK,e) = xpK;    yp(1:NplotK,e) = ypK;    zp(1:NplotK,e) = zpK;
%
% end

%% GMSH plotting in terms of monomials (in u,v,w coords)

[rpH spH tpH] = hex_nodes(N); 
upH = rpH; vpH = spH; wpH = tpH;
% upH = tpH; vpH = rpH; wpH = spH;

[rpW spW tpW] = wedge_nodes(N); % GMSH u,v,w coords for wedge 
upW = tpW; vpW = rpW; wpW = spW; upW = (upW+1)/2; vpW = (vpW+1)/2;

[rpP spP tpP] = pyr_nodes(N); 
[upP vpP wpP] = biunit_to_sym(rpP,spP,tpP);

[rpT spT tpT] = tet_nodes(N); 
upT = (rpT + 1)/2; vpT = (spT+1)/2; wpT = (tpT+1)/2;

% hexes
PH = zeros(NpH,3); VpH = zeros(length(rpH),NpH);
ind = 1;
for i = 0:N
    for j = 0:N
        for k = 0:N
            PH(ind,1) = i; PH(ind,2) = j; PH(ind,3) = k;
            VpH(:,ind) = (upH.^i).*(vpH.^j).*(wpH.^k); 
            ind = ind + 1;
        end
    end
end
VDMH = hex_basis(N,rpH,spH,tpH);
FH = inv(VpH); %VpH\VDMH;

% wedges
PW = zeros(NpW,3); VpW = zeros(length(rpW),NpW);
ind = 1;
for i = 0:N
    for j = 0:N
        for k = 0:N-i
            PW(ind,1) = k; PW(ind,2) = i; PW(ind,3) = j;
            VpW(:,ind) = (upW.^k).*(vpW.^i).*(wpW.^j); % for gmsh coords
            ind = ind + 1;
        end
    end
end
% for nodal wedges: eval/sqrt(J) at nodal points then interpolate
VDMW = wedge_basis(N,rpW,spW,tpW);
FW = inv(VpW);

% pyr
PP = zeros(NpP,3); 
% VpP = zeros(length(rpP),NpP);
% ind = 1;
% for i = 0:N
%     for j = 0:N
%         for k = 0:N-max(i,j)
%             PP(ind,1) = i; PP(ind,2) = j; PP(ind,3) = k;
%             VpP(:,ind) = (upP.^i).*(vpP.^j).*(wpP.^k);
%             ind = ind+1;
%         end
%     end
% end
% symmetric ref pyramid - get nodal in terms of GMSH's bergot basis
VpP = bergot_pyr_basis(N,upP,vpP,wpP); 
VDMP = pyr_basis(N,rpP,spP,tpP);
FP = inv(VpP); 
% FP = VpP\VDMP;

% tet
PT = zeros(NpT,3); VpT = zeros(length(rpT),NpT);
ind = 1;
for i = 0:N
    for j = 0:N
        for k = 0:N-(i+j)
            PT(ind,1) = i; PT(ind,2) = j; PT(ind,3) = k;
            VpT(:,ind) = (upT.^i).*(vpT.^j).*(wpT.^k);
            ind = ind+1;
        end
    end
end
VDMT = tet_basis(N,rpT,spT,tpT);
FT = inv(VpT); %VpT\VDMT;

