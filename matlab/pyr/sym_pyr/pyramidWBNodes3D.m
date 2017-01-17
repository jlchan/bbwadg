
% beta = interior blend, maybe implement alphaFace = surface
function [r s t] = pyramidWBNodes3D(N,setFaceNodes,beta)
if nargin == 0
    N = 6;
    setFaceNodes = 1; % default to original WB nodes on faces
    beta = 0;
end
if nargin == 1
    setFaceNodes = 1;
end
if nargin < 3 || isempty(beta)
    beta = 0;
end

% target nodes on surface - default to Nodes3D optimized blend for
% conformity with tet nodes
[mapr maps mapt] = PyramidSurfaceNodes3D(N);

if setFaceNodes == 1    % compute new node positions based on surface distributions
    
%     disp('using interpolation of original WB faces')
    % %% find coefficients mapping from equispaced to warped surface
    [rbc sbc tbc] = PyramidSurfaceEquiNodes3D(N);
    [Vbc v_ids etri_ids equad_ids ftri_ids fquad_ids] = JVandermonde3D(N, rbc, sbc, tbc, beta);
    
    ids = [v_ids etri_ids equad_ids ftri_ids fquad_ids];
    maprst = [mapr maps mapt];
    mapcrst = Vbc(:,ids)\maprst;
    %
    % % evaluate map at equispaced volumed nodes
    [req seq teq] = pyramidEquiNodes(N);
    [Veq v_ids etri_ids equad_ids ftri_ids fquad_ids] = JVandermonde3D(N, req, seq, teq, beta);
    
    ids = [v_ids etri_ids equad_ids ftri_ids fquad_ids];
    rst = Veq(:,ids)*mapcrst; %% final coordinates of all nodes in triangle
    
    r = rst(:,1); s = rst(:,2); t = rst(:,3);
    
elseif setFaceNodes == 2 % replace surface nodes with new WBInterp3D optimized tet nodes    
    
    disp('using interpolation of new WB faces')
     
    % replace surface nodes with new optimized surface nodes
    [mapr maps mapt] = PyramidSurfaceNodesWBInterp3D(N);
    
    % %% find coefficients mapping from equispaced to warped surface
    [rbc sbc tbc] = PyramidSurfaceEquiNodes3D(N);
    [Vbc v_ids etri_ids equad_ids ftri_ids fquad_ids] = JVandermonde3D(N, rbc, sbc, tbc, beta);
    
    ids = [v_ids etri_ids equad_ids ftri_ids fquad_ids];
    maprst = [mapr maps mapt];
    mapcrst = Vbc(:,ids)\maprst;
    %
    % % evaluate map at equispaced volumed nodes
    [req seq teq] = pyramidEquiNodes(N);
    [Veq v_ids etri_ids equad_ids ftri_ids fquad_ids] = JVandermonde3D(N, req, seq, teq, beta);
    
    ids = [v_ids etri_ids equad_ids ftri_ids fquad_ids];
    rst = Veq(:,ids)*mapcrst; %% final coordinates of all nodes in triangle
    
    r = rst(:,1); s = rst(:,2); t = rst(:,3);
    
elseif setFaceNodes==10 % option: deform on edges only, then faces for two separate blends. 
    
    beta1 = 0; % edge blend
    beta2 = beta; % face blend
    
    % =============== STEP 1: deform based on edges ======================
    xieq = linspace(-1,1,N+1);   xigll = JacobiGL(0,0,N);
    xieq = xieq(2:end-1); xigll= xigll(2:end-1); % trim endpoints = vertices
    xieq = (xieq+1)/2; xigll = (xigll+1)/2;
    
    edges = [1 5; 2 5; 3 5 ;4 5; 1 2 ; 2 3; 3 4; 4 1];
    rv = [-1 1 1 -1 0]; sv = [-1 -1 1 1 0]; tv = [0 0 0 0 1];
    req = []; seq = []; teq = [];  rgll = []; sgll = []; tgll = [];
    for e = 1:8
        v1 = edges(e,1); v2 = edges(e,2);        
        req = [req (rv(v2)-rv(v1))*xieq + rv(v1)];
        seq = [seq (sv(v2)-sv(v1))*xieq + sv(v1)];
        teq = [teq (tv(v2)-tv(v1))*xieq + tv(v1)];        
        rgll = [rgll (rv(v2)-rv(v1))*xigll + rv(v1)];
        sgll = [sgll (sv(v2)-sv(v1))*xigll + sv(v1)];
        tgll = [tgll (tv(v2)-tv(v1))*xigll + tv(v1)];
    end
    
    % add in interior square quad base nodes
    xieq = linspace(-1,1,N+1);   xigll = JacobiGL(0,0,N);
    xieq = xieq(2:end-1); xigll= xigll(2:end-1); % trim endpoints = vertices
    
   %     only does edges
    xieq = xieq(:); xigll = xigll(:);
    rqe = [xieq xieq.^0 -xieq -xieq.^0];    sqe = [-xieq.^0 xieq xieq.^0 -xieq];
    tqe = zeros(size(rqe));
    rqgll = [xigll xieq.^0 -xigll -xieq.^0];    sqgll = [-xieq.^0 xigll xieq.^0 -xigll];
    tqgll = zeros(size(rqgll));
      
    % add verts/edges/base
    rbc = [rv(:); req(:); rqe(:)]; sbc = [sv(:); seq(:); sqe(:)];  tbc = [tv(:); teq(:); tqe(:)];
    mapr = [rv(:); rgll(:); rqgll(:)]; maps = [sv(:); sgll(:); sqgll(:)];  mapt = [tv(:); tgll(:); tqgll(:)];
       
    % solve interp problem for map
    [Vbc v_ids etri_ids equad_ids ftri_ids fquad_ids] = JVandermonde3D(N, rbc, sbc, tbc, beta1);
%     ids = [v_ids etri_ids equad_ids fquad_ids]; % use vertices, triangle edges, quad edges, and quad face for interp
    ids = [v_ids etri_ids equad_ids]; % HACK FOR EDGE INTERP ONLY
    mapcrst = Vbc(:,ids)\[mapr maps mapt];
    
    % evaluate map at equispaced volumed nodes
    [req,seq,teq] = pyramidEquiNodes(N);
    [Veq v_ids etri_ids equad_ids ftri_ids fquad_ids] = JVandermonde3D(N, req, seq, teq, beta1);
    
%     ids = [v_ids etri_ids equad_ids fquad_ids]; % use vertices, triangle edges, quad edges, and quad face for interp
    ids = [v_ids etri_ids equad_ids]; % HACK FOR EDGE INTERP ONLY
    rst = Veq(:,ids)*mapcrst; % final coordinates of all nodes in triangle
    r = rst(:,1); s = rst(:,2); t = rst(:,3);
    
    % ============== STEP 2: deform based on face base only ===================
    
    % interior square quad base nodes
    xigll = JacobiGL(0,0,N); 
    xigll= xigll(2:end-1); % trim endpoints = vertices       
    [rqgll sqgll] = meshgrid(xigll); tqgll = zeros(size(rqgll));
        
    tol = 1e-8;
    ids = abs(t) < tol & abs(abs(r)-1) > tol & abs(abs(s)-1) > tol; ids = ids(:); % base interior nodes    
    [Vbc v_ids etri_ids equad_ids ftri_ids fquad_ids] = JVandermonde3D(N, r(ids), s(ids), t(ids), beta2);    
    
    maprst = [rqgll(:) sqgll(:) 0*rqgll(:)] - [r(ids) s(ids) t(ids)];
    
    mapcrst = Vbc(:,fquad_ids)\maprst; % construct map only for quad face
        
    % evaluate map at edge-deformed nodes
    [V2 v_ids etri_ids equad_ids ftri_ids fquad_ids] = JVandermonde3D(N, r, s, t, beta2);
    
    rst = V2(:,fquad_ids)*mapcrst; % final coordinates of all nodes in triangle    
    r = r + rst(:,1); s = s + rst(:,2); t = t+ rst(:,3);
    
else  % just do GLL on edges
    
    disp('enforcing GLL on edges only')
    % pyramid
    %
    %    4------3
    %    |\    /|
    %    | \  / |
    %    |  5   |
    %    | /  \ |
    %    |/    \|
    %    1------2
    
    xieq = linspace(-1,1,N+1);   xigll = JacobiGL(0,0,N);
    xieq = xieq(2:end-1); xigll= xigll(2:end-1); % trim endpoints = vertices
    xieq = (xieq+1)/2; xigll = (xigll+1)/2;
    
    edges = [1 5; 2 5; 3 5 ;4 5; 1 2 ; 2 3; 3 4; 4 1];
    rv = [-1 1 1 -1 0]; sv = [-1 -1 1 1 0]; tv = [0 0 0 0 1];
    req = []; seq = []; teq = [];  rgll = []; sgll = []; tgll = [];
    for e = 1:8
        v1 = edges(e,1); v2 = edges(e,2);
        
        req = [req (rv(v2)-rv(v1))*xieq + rv(v1)];
        seq = [seq (sv(v2)-sv(v1))*xieq + sv(v1)];
        teq = [teq (tv(v2)-tv(v1))*xieq + tv(v1)];
        
        rgll = [rgll (rv(v2)-rv(v1))*xigll + rv(v1)];
        sgll = [sgll (sv(v2)-sv(v1))*xigll + sv(v1)];
        tgll = [tgll (tv(v2)-tv(v1))*xigll + tv(v1)];
    end
    
    % add in interior square quad base nodes
    xieq = linspace(-1,1,N+1);   xigll = JacobiGL(0,0,N);
    xieq = xieq(2:end-1); xigll= xigll(2:end-1); % trim endpoints = vertices
    
%     % HACK: only does edges
%     xieq = xieq(:); xigll = xigll(:);
%     rqe = [xieq xieq.^0 -xieq -xieq.^0];    sqe = [-xieq.^0 xieq xieq.^0 -xieq];
%     tqe = zeros(size(rqe));
%     rqgll = [xigll xieq.^0 -xigll -xieq.^0];    sqgll = [-xieq.^0 xigll xieq.^0 -xigll];
%     tqgll = zeros(size(rqgll));
    
    [rqe sqe] = meshgrid(xieq); tqe = zeros(size(rqe));
    [rqgll sqgll] = meshgrid(xigll); tqgll = zeros(size(rqgll));
    
    % add verts/edges/base
    rbc = [rv(:); req(:); rqe(:)]; sbc = [sv(:); seq(:); sqe(:)];  tbc = [tv(:); teq(:); tqe(:)];
    mapr = [rv(:); rgll(:); rqgll(:)]; maps = [sv(:); sgll(:); sqgll(:)];  mapt = [tv(:); tgll(:); tqgll(:)];
    
%     mapr = mapr-rbc; maps = maps-sbc; mapt = mapt-tbc;        
    
    % solve interp problem for map
    [Vbc v_ids etri_ids equad_ids ftri_ids fquad_ids] = JVandermonde3D(N, rbc, sbc, tbc, beta);
    ids = [v_ids etri_ids equad_ids fquad_ids]; % use vertices, triangle edges, quad edges, and quad face for interp
%     ids = [v_ids etri_ids equad_ids]; % HACK FOR EDGE INTERP ONLY
    mapcrst = Vbc(:,ids)\[mapr maps mapt];
    
    % evaluate map at equispaced volumed nodes
    [req,seq,teq] = pyramidEquiNodes(N);
    [Veq v_ids etri_ids equad_ids ftri_ids fquad_ids] = JVandermonde3D(N, req, seq, teq, beta);
    
    ids = [v_ids etri_ids equad_ids fquad_ids]; % use vertices, triangle edges, quad edges, and quad face for interp
%     ids = [v_ids etri_ids equad_ids]; % HACK FOR EDGE INTERP ONLY
    rst = Veq(:,ids)*mapcrst; % final coordinates of all nodes in triangle
    r = rst(:,1); s = rst(:,2); t = rst(:,3);
    
%     r = req + rst(:,1); s = seq + rst(:,2); t = teq + rst(:,3);
    
end

%%
if nargin ==0    
    plot3(mapr,maps,mapt,'ro');hold on
    plot3(r,s,t,'.');hold on;
    maxleb = PyramidLebesgue3D(N,r,s,t,5000,10);
    title(sprintf('leb const = %d',maxleb))
    view(-15,5)
    %     keyboard
end

% currently uses only face blending
function [V v_ids etri_ids equad_ids ftri_ids fquad_ids] = JVandermonde3D(N, r, s, t, beta)

%% pyramid
%%
%%    4------3
%%    |\    /|
%%    | \  / |
%%    |  5   |
%%    | /  \ |
%%    |/    \|
%%    1------2

% alphastore = [0;0;0;0.1002; 1.1332;1.5608;1.3413;1.2577;1.1603;...
%     1.10153;0.6080;0.4523;0.8856;0.8717;0.9655];
% alphaT = alphastore(N);

tol = 1e-10;
V = [];

% vertex functions
V(:,1) = .25*(1-r-s-t+r.*s./(1-t + tol));
V(:,2) = .25*(1+r-s-t-r.*s./(1-t + tol));
V(:,3) = .25*(1+r+s-t+r.*s./(1-t + tol));
V(:,4) = .25*(1-r+s-t-r.*s./(1-t + tol));
V(:,5) = t;

v_ids = 1:5;

sk = 6;

% edge function for triangular edges
edges = [1 5; 2 5; 3 5; 4 5];
etri_ids = [];
for e = 1:4
    i1 = edges(e,1);
    i2 = edges(e,2);
    
    switch e
        case 1
            blend1 = V(:,2);
            blend2 = V(:,4);
        case 2
            blend1 = V(:,1);
            blend2 = V(:,3);
        case 3
            blend1 = V(:,4);
            blend2 = V(:,2);
        case 4
            blend1 = V(:,1);
            blend2 = V(:,3);            
    end
    blend2 = (1 + beta*(blend1).^2).*(1 + beta*(blend2).^2); % blend towards opposite vertex on faces
    
    for i=0:N-2
        xi = V(:,i1)-V(:,i2);
        V(:,sk) = blend2.*V(:,i1).*V(:,i2).*JacobiP(xi, 1, 1, i);
        etri_ids = [etri_ids sk];
        sk = sk+1;
    end
end

% edge functions for base
edges = [1 2; 2 3; 3 4; 4 1];
equad_ids = [];
for e = 1:4
    i1 = edges(e,1);
    i2 = edges(e,2);
    
%     op_vert = setdiff(1:5,edges(e,:));
    %V(:,op_vert(1)).*V(:,op_vert(2)).*V(:,op_vert(3));
    blend2 = (1 + beta*(V(:,5)).^2); % blend towards top vertex 
    
    for i=0:N-2
        xi = V(:,i1)-V(:,i2);
        V(:,sk) = blend2.*V(:,i1).*V(:,i2).*JacobiP(xi, 1, 1, i);
        equad_ids = [equad_ids sk];
        sk = sk+1;
    end
end

%triangular faces
faces = [1 2 5; 2 3 5; 3 4 5; 4 1 5];
ftri_ids = [];
for f = 1:4
    i1 = faces(f,1);
    i2 = faces(f,2);
    i3 = faces(f,3);
    
    % bubble edge blend
    %         op_vert = setdiff(1:5,faces(f,:));  blend2 = V(:,op_vert(1)).*V(:,op_vert(2));
    
    % plane edge blend
    switch f;
        case 1;
            blend2 = ((1+s)-t)/2;
        case 2;
            blend2 = ((1-r)-t)/2;
        case 3;
            blend2 = ((1-s)-t)/2;
        case 4;
            blend2 = ((1+r)-t)/2;
    end        
    %blend2 = 1-V(:,5); % blend to base   
    
    blend2 = 1 + beta*(blend2).^2;
    %blend2 = beta*(V(:,5)-.5); % blend to base - switches negative at midway pt.      
    
    
    L1 = V(:,i1);    L2 = V(:,i2);    L3 = V(:,i3);
    [x y] = eqbarytoxy(L1,L2,L3);
    [rr ss] = xytors(x,y);
    Vf = Vandermonde2D(N-3,rr,ss);
    for i = 1:size(Vf,2)
        V(:,sk) = blend2.*V(:,i1).*V(:,i2).*V(:,i3).*Vf(:,i);
        ftri_ids = [ftri_ids sk];
        sk = sk + 1;
    end
end

% square face on the bottom
blend2 = 1 + beta*(V(:,5)).^2;
fquad_ids = [];
for i = 0:N-2
    for j = 0:N-2
        V(:,sk) = blend2.*V(:,1).*V(:,2).*V(:,3).*V(:,4).*JacobiP(r,1,1,i).*JacobiP(s,1,1,j);
        fquad_ids = [fquad_ids sk];
        sk = sk + 1;
    end
end


return;



