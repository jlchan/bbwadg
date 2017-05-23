NODETOL=1e-6;

% # dofs for each elem type
NpH = (N+1)^3;
NpW = (N+1)^2*(N+2)/2;
NpP = (N+1)*(N+2)*(2*N+3)/6;
NpT = (N+1)*(N+2)*(N+3)/6;

%% reference vol/surface

Nq = N; % cubature order if we want to overintegrate...

% hex face vertex ordering
fvH{1} = [1 4 3 2]; fvH{2} = [1 2 6 5]; fvH{3} = [2 3 7 6];
fvH{4} = [3 4 8 7]; fvH{5} = [4 1 5 8]; fvH{6} = [5 6 7 8];
r1 = [-1  1  1 -1 -1  1  1 -1]'; s1 = [-1 -1  1  1 -1 -1  1  1]';
t1 = [-1 -1 -1 -1  1  1  1  1]';
[rfH sfH tfH wfH fidsH] = surface_cubature(Nq,r1,s1,t1,fvH);
rstH = [r1 s1 t1];
NfcH = length(wfH);

% wedge face vertex ordering
fvW{1} = [1 3 2]; fvW{2} = [4 5 6];
fvW{3} = [1 2 5 4]; fvW{4} = [3 1 4 6]; fvW{5} = [2 3 6 5];
u = [-1 1 -1 -1 1 -1]; v = [-1 -1 1 -1 -1 1]; w = [-1 -1 -1 1 1 1];
r1 = v(:); s1 = w(:); t1 = u(:); % flipping coordinates for Gmsh
[rfW sfW tfW wfW fidsW] = surface_cubature(Nq,r1,s1,t1,fvW);
rstW = [r1 s1 t1];
NfcW = length(wfW);

% pyramid face vertex ordering
fvP{1} = [1 2 5]; fvP{2} = [4 1 5];
fvP{3} = [2 3 5]; fvP{4} = [3 4 5];
fvP{5} = [1 4 3 2];
r1 = [ -1   1   1  -1  -1 ]';s1 = [ -1  -1   1   1  -1 ]';
t1 = [ -1  -1  -1  -1   1 ]';
[rfP sfP tfP wfP fidsP] = surface_cubature(Nq,r1,s1,t1,fvP);
rstP = [r1 s1 t1];
NfcP = length(wfP);

% tet face vertex orderings
fvT{1} = [1 3 2]; fvT{2} = [1 2 4]; fvT{3} = [3 1 4]; fvT{4} = [2 3 4];
r1 = [-1  1 -1 -1 ]'; s1 = [-1 -1  1 -1 ]'; t1 = [-1 -1 -1  1 ]';
[rfT sfT tfT wfT fidsT] = surface_cubature(Nq,r1,s1,t1,fvT);
rstT = [r1 s1 t1];
NfcT = length(wfT);

% reference volume stuff
[rqT sqT tqT w] = tet_cubature(2*Nq); % integrate order 2N
wT = w; NcT = length(w(:));
[VT VTr VTs VTt] = tet_basis(N,rqT,sqT,tqT);
VTf = tet_basis(N,rfT,sfT,tfT);

if useNodalTets
    [xn yn zn] = Nodes3D(N); [rT sT tT] = xyztorst(xn,yn,zn);
    [VTnodal VTrn VTsn VTtn] = tet_basis(N,rT,sT,tT);
    
    VT = VT/VTnodal; % quadrature VDM
    VTf = VTf/VTnodal;
    
    invMThat = VTnodal*VTnodal';
    MThat = inv(invMThat);
    VTr = VTrn/VTnodal;  VTs = VTsn/VTnodal;  VTt = VTtn/VTnodal;
end

[rqP sqP tqP w] = pyr_cubature(Nq);
wP = w; NcP = length(w(:));
[VP VPr VPs VPt] = pyr_basis(N,rqP,sqP,tqP);
VPf = pyr_basis(N,rfP,sfP,tfP);

[rqW sqW tqW w] = wedge_cubature(Nq);
wW = w; NcW = length(w(:));
[VW VWr VWs VWt] = wedge_basis(N,rqW,sqW,tqW);
VWf = wedge_basis(N,rfW,sfW,tfW);

[rqH sqH tqH w] = hex_cubature(Nq);
wH = w; NcH = length(w(:));
[VH VHr VHs VHt] = hex_basis(N,rqH,sqH,tqH);
VHf = hex_basis(N,rfH,sfH,tfH);

%%

NcMax = NcH;
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

NpMax = max([NpH NpW NpP NpT]);
J = zeros(NcMax,K); Jsurf = zeros(NfcMax,K);

for e = 1:K
    sprintf('Working on elem %d out of %d',e,K)
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
        
        nc
    end
    
    sJ = sqrt(nx(:,e).^2 + ny(:,e).^2 + nz(:,e).^2);
    nx(:,e) = nx(:,e)./sJ; ny(:,e) = ny(:,e)./sJ; nz(:,e) = nz(:,e)./sJ;
    Fscale(:,e) = sJ; % area of face for planar
    
    wsJ(1:NfcK,e) = wf(:).*sJ(1:NfcK).*Jf(:); % premultiply jacobian with weights
    JsB(1:NfcK,e) = sJ(1:NfcK).*Jf(:)./ref_face_area(1:NfcK); % Jf = Jacobian on face = used bound
    %     keyboard
end

%% connectivities

[mapM,mapP,mapB] = BuildHybridMaps();

%%
% cubature points for each element type
xqT = x(1:NcT,Nv==4); yqT = y(1:NcT,Nv==4); zqT = z(1:NcT,Nv==4);
xqP = x(1:NcP,Nv==5); yqP = y(1:NcP,Nv==5); zqP = z(1:NcP,Nv==5);
xqW = x(1:NcW,Nv==6); yqW = y(1:NcW,Nv==6); zqW = z(1:NcW,Nv==6);
xqH = x(1:NcH,Nv==8); yqH = y(1:NcH,Nv==8); zqH = z(1:NcH,Nv==8);

tetK = find(Nv==4); KT = length(tetK); wJT = wJ(1:NcT,tetK);
pyrK = find(Nv==5); KP = length(pyrK); wJP = wJ(1:NcP,pyrK);
wedgK = find(Nv==6); KW = length(wedgK); wJW = wJ(1:NcW,wedgK);
hexK = find(Nv==8); KH = length(hexK); wJH = wJ(1:NcH,hexK);

% diff volume terms for nodal tets
if useNodalTets
    
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
end

%% pre-invert mass matrix

NpMax = NpH; % hex has most dofs
NcMax = NcH;
NfcMax = NfcH;

invM = zeros(NpMax,K);
if useNodalTets
    invM(1:NpT,tetK) = 1; % mass matrix computed using VDMs instead
else
    for e = 1:KT
        invM(1:NpT,tetK(e)) = 1./diag(VT'*spdiag(wJ(1:NcT,tetK(e)))*VT);
    end
end

for e = 1:KP
    invM(1:NpP,pyrK(e)) = 1./diag(VP'*spdiag(wJ(1:NcP,pyrK(e)))*VP);
end


MW = cell(KW,1);
for e = 1:KW
    %         invM(1:NpW,wedgK(e)) = 1./diag(VW'*spdiag(wJ(1:NcW,wedgK(e)))*VW);
    Me = VW'*spdiag(wJ(1:NcW,wedgK(e)))*VW;
    Me(abs(Me)<1e-6) = 0;
    MW{e} = sparse(Me);
end

for e = 1:KH
    invM(1:NpH,hexK(e)) = 1./diag(VH'*spdiag(wJ(1:NcH,hexK(e)))*VH);
end


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

