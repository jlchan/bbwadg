clear -globals
clear

% useQuads = 1;mypath

Globals2D;
N = 4;
K1D = 2;

% Read in Mesh
[Nv, VX, VY, K, EToV] = QuadMesh2D(K1D,K1D);

a = .0/K1D;
ids = find(abs(abs(VX)-1) > 1e-8 & abs(abs(VY)-1) > 1e-8);
VX(ids) = VX(ids) + a*randn(size(ids));
VY(ids) = VY(ids) + a*randn(size(ids));

StartUpQuad2D;

% PlotMesh2D
% for e = 1:K
%     vids = EToV(e,:);
%     text(mean(VX(vids)),mean(VY(vids)),num2str(e))
% end
% return

fmask1   = find( abs(s+1) < NODETOL)'; 
fmask2   = find( abs(r-1) < NODETOL)';
fmask3   = find( abs(s-1) < NODETOL)'; fmask3 = fmask3(end:-1:1);
fmask4   = find( abs(r+1) < NODETOL)'; fmask4 = fmask4(end:-1:1);
Fmask  = [fmask1;fmask2;fmask3;fmask4]';

BuildPeriodicMaps2D(2,2);

[rq1D wq1D] = JacobiGL(0,0,N);
[rq sq] = meshgrid(rq1D);
rq = rq(:); sq = sq(:);
[wr ws] = meshgrid(wq1D);
wq = wr(:).*ws(:);
wfq = repmat(wq1D,4,1);

e = ones(size(rq1D));
rf = [rq1D;e;rq1D;-e];
sf = [-e;rq1D;e;rq1D];
nrJ = [0*e;e;0*e;-e];
nsJ = [-e;0*e;e;0*e];
Vf = Vandermonde2D(N,rf,sf)/V;

M = diag(wq);
Lf = M\(Vf'*diag(wfq));

%% curve elements

a = .0;
x = x + a.*sin(pi*x).*sin(pi*y);
y = y + a.*sin(pi*x).*sin(pi*y);

%% lightweight adaptive mesh data structures

global qt activeK Vsplit Vfsplit

% quadtree
qt = cell(K,1);
for e = 1:K
    qt{e} = struct('parent',e,'children',[],'childId',[]);
end

activeK = ones(K,1);

% VDM for split elements and faces
rs = [-1+.5*(1+r); .5*(1+r); .5*(1+r); -1+.5*(1+r)];
ss = [-1+.5*(1+s); -1+.5*(1+s); .5*(1+s); .5*(1+s)];
Vsplit = Vandermonde2D(N,rs,ss)/V;

rq1Dsplit = [-1+.5*(1+rq1D); .5*(1+rq1D)];
wq1Dsplit = .5*[wq1D; wq1D];
Vfsplit = Vandermonde1D(N,rq1Dsplit)/Vandermonde1D(N,rq1D);
Pqsplit = (Vfsplit'*diag(wq1Dsplit)*Vfsplit)\(Vfsplit'*diag(wq1Dsplit));

%% refine elems

hrefine(round(K/2));
% hrefine(6);
% hrefine(7);
% hrefine(10);
% hrefine(11);
% hrefine(10);


%% build non-con maps

% isolate non-conforming faces, find parent elem           
nonconFaces = zeros(K,Nfaces);
sk = 1;
noncon_face_counter = 1;
for e = 1:K
    if activeK(e)==1
        for f = 1:Nfaces
            children = qt{EToE(e,f)}.children;
            if isempty(children) 
            else % if nbr refined, face is non-conforming                
                nonconFaces(e,f) = noncon_face_counter; %noncon_face_counter;
                noncon_face_counter = noncon_face_counter + 1;
            end
        end        
        sk = sk + 1;
    end
end
[fid eid] = find(nonconFaces'); % transpose for ordering
ncfaces = fid + (eid-1)*Nfaces; % list of faces to be split (including non-active elems)
nceid = eid(:); 
ncfid = fid(:); 
num_conf_faces = Nfaces*K;
num_nonconf_faces = length(ncfaces);
num_faces_total = num_conf_faces+2*num_nonconf_faces;

EToFnc = zeros(K,Nfaces);
for i = 1:num_nonconf_faces
    EToFnc(nceid(i),ncfid(i)) = i;
end

xf = reshape(x(Fmask(:),:),Nfp,Nfaces*K);
yf = reshape(y(Fmask(:),:),Nfp,Nfaces*K);
xfnc = reshape(Vfsplit*xf(:,ncfaces),Nfp,num_nonconf_faces*2);
yfnc = reshape(Vfsplit*yf(:,ncfaces),Nfp,num_nonconf_faces*2);

eK = find(activeK);
text(mean(x(:,eK)),mean(y(:,eK)),num2str(eK))
hold on
plot([xf xfnc],[yf yfnc],'o')
xf(:,ncfaces) = Pqsplit*reshape(xfnc,2*Nfp,num_nonconf_faces);
yf(:,ncfaces) = Pqsplit*reshape(yfnc,2*Nfp,num_nonconf_faces);
plot(xf,yf,'^')
return

% match faces to faces
ee = {[1 2],[2 3],[3 4],[4 1]}; % new refined elements adjacent to each face
en = {[4 3],[1 4],[2 1],[3 2]}; % neighbors adjacent to each new elem
fn = [3 4 1 2]; % faces of neighbors adjacent to each new elem

FToF = zeros(num_faces_total,1); % faces to faces (active faces = nonzero ids)
sk = 1;
for e = 1:K
    if activeK(e) == 1        
        for f = 1:Nfaces
            fid = f + (e-1)*Nfaces;                       
            e_level = get_level(qt,e);
            
            % neighbor data
            enbr = EToE(e,f);
            nbr_children = qt{enbr}.children; 
            enbr_level = get_level(qt,enbr);
                                                
            if ~isempty(nbr_children) % neighbor is refined = split to small face, index into conforming
                              
                
                ncfid = EToFnc(e,f); % get noncon fid
                enbrs = nbr_children(en{f});
                nbrfid = fn(f) + (enbrs-1)*Nfaces;
                FToF(2*ncfid-1 + num_conf_faces) = nbrfid(1);
                FToF(2*ncfid + num_conf_faces) = nbrfid(2);

            elseif e_level==enbr_level+1 % small face to split face: index into noncon
                                              
                fnbr = EToFnc(enbr,fn(f)); % neighboring face                
                childId = qt{e}.childId;
                if ee{f}(2)==childId % if first part of face f (reversed ordering)
                    FToF(fid) = 2*fnbr-1 + num_conf_faces;
                elseif ee{f}(1)==childId % if second part of face f
                    FToF(fid) = 2*fnbr  + num_conf_faces;
                end
                
            elseif e_level==enbr_level % conforming (elem, nbr same level)
                
                FToF(fid) = EToF(e,f) + (enbr-1)*Nfaces; 

            else % if all of the above are false, something's wrong!
                
                keyboard
                
            end
        end
        sk = sk + 1;
    end    
end

xf = [xf xfnc];
yf = [yf yfnc];

if 0
    for e = 1:K
        if activeK(e)==1
            for f = 1:Nfaces
                
                clf
                eK = find(activeK);
                text(mean(x(:,eK)),mean(y(:,eK)),num2str(eK))
                hold on
                plot(xf,yf,'*')
                title(sprintf('elem %d, face %d',e,f))
                if isempty(qt{EToE(e,f)}.children) % if nbr is un-refined, conforming face
                    
                    fid = f + (e-1)*Nfaces;
                    plot(xf(:,fid),yf(:,fid),'o','markersize',16)
                    fP = FToF(fid);
                    plot(xf(:,fP),yf(:,fP),'s','markersize',16)
                    
                else % if non-conforming face
                    
                    ncfid = EToFnc(e,f);
                    fid = 2*ncfid-1 + num_conf_faces;
                    
                    plot(xf(:,fid),yf(:,fid),'o','markersize',16)
                    fP = FToF(fid);
                    plot(xf(:,fP),yf(:,fP),'s','markersize',16)
                    pause
                    
                    fid = 2*ncfid + num_conf_faces;
                    plot(xf(:,fid),yf(:,fid),'o','markersize',16)
                    fP = FToF(fid);
                    plot(xf(:,fP),yf(:,fP),'s','markersize',16)
                end
                
                pause
            end
        end
    end
end

mapMq = reshape(1:Nfp*(num_faces_total),Nfp,num_faces_total);
mapPq = mapMq;

for f = 1:(num_faces_total)
    fP = FToF(f);
    
    if fP~=0
        
        tol = 1e-10;
        [X1 Y1] = meshgrid(xf(:,f),yf(:,f));
        [X2 Y2] = meshgrid(xf(:,fP),yf(:,fP));
        DX = abs(X1-X2');
        DY = abs(Y1-Y2');
        D = DX + DY;
        [p,~] = find(D<tol);
                
        % NOTE - does not work if K1D is too small!!
        if length(p) == 0
            % assume periodic boundary, find match in x,y
            [px,~] = find(DX<tol);
            [py,~] = find(DY<tol);
            if length(px)==0
                p = py;
            elseif length(py)==0
                p = px;
            else
                keyboard
            end
        end
        
        mapPq(:,f) = mapMq(p,fP);                       
    end
end

for f = 1:(num_faces_total)
    fP = FToF(f);   
    if fP~=0 
        for i = 1:Nfp
            id = i + (f-1)*Nfp;
            idM = mapMq(id);
            idP = mapPq(id);
            clf
            plot(xf,yf,'k.','markersize',12)
            hold on
            plot(xf(idM),yf(idM),'o','markersize',16,'linewidth',2)
            plot(xf(idP),yf(idP),'^','markersize',16,'linewidth',2)
            pause(.1)
        end
    end
end


%% other routines: refining 1-irregular meshes

%         nbr 3
%      2----7----3
%      | e4 | e3 |
% nbr4 8----9----6  nbr 2
%      | e1 | e2 |
%      1----5----4
%         nbr 1
function hrefine(e)

Globals2D
global qt activeK Vsplit Vfsplit

if activeK(e)==0
    fprintf('Cannot refine element %d: not active\n',e);
    return
end

% TODO: enforce 1-irregularity

% add to quadtree
qt = [qt;cell(4,1)];
qt{e}.children = K+(1:4);
qt{K+1} = struct('parent',e,'children',[],'childId',1);
qt{K+2} = struct('parent',e,'children',[],'childId',2);
qt{K+3} = struct('parent',e,'children',[],'childId',3);
qt{K+4} = struct('parent',e,'children',[],'childId',4);

% update neighbor data
EToE = [EToE;zeros(4,Nfaces)]; % expand neighbor array
EToF = [EToF;zeros(4,Nfaces)]; % expand neighbor array
newnbr = K + (1:4);

% local neighbors
EToE(K+1,2) = newnbr(2); EToE(K+1,3) = newnbr(4);
EToE(K+2,4) = newnbr(1); EToE(K+2,3) = newnbr(3);
EToE(K+3,1) = newnbr(2); EToE(K+3,4) = newnbr(4);
EToE(K+4,1) = newnbr(1); EToE(K+4,2) = newnbr(3);
%         nbr 3
%      2----7----3
%      | e4 | e3 |
% nbr4 8----9----6  nbr 2
%      | e1 | e2 |
%      1----5----4
%         nbr 1
EToF(K+1,2) = 4; EToF(K+1,3) = 1;
EToF(K+2,4) = 2; EToF(K+2,3) = 1;
EToF(K+3,1) = 3; EToF(K+3,4) = 2;
EToF(K+4,1) = 3; EToF(K+4,2) = 4;

% inherit external neighbors: check nbrs of e
ee = {[1 2],[2 3],[3 4],[4 1]}; % new refined elements adjacent to each face
en = {[4 3],[1 4],[2 1],[3 2]}; % neighbors adjacent to each new elem
fn = [3 4 1 2]; % faces of neighbors adjacent to each new elem
for f = 1:4
    for i = 1:2 % neighbors per face
        new_elem = K+ee{f}(i);
        if isempty(qt{EToE(e,f)}.children) 
            % if neighbor is not refined, connect to parent
            parent = qt{EToE(e,f)}.parent;
            EToE(new_elem,f) = parent;
            EToF(new_elem,f) = EToF(parent,f);
        else % if neighbor refined, pick up neighbors from children
            nbr_elem = qt{EToE(e,f)}.children(en{f}(i));
            EToE(new_elem,f) = nbr_elem;
            EToE(nbr_elem,fn(f)) = new_elem;
            
            EToF(new_elem,f) = EToF(e,f); % inherit face connectivity
        end
    end
end

% could also build by sweeping quadtree but this is easier
activeK(e) = 0;
activeK = [activeK; ones(4,1)];

xs = reshape(Vsplit*x(:,e),Np,4);
ys = reshape(Vsplit*y(:,e),Np,4);
x = [x xs];
y = [y ys];

% update number of elements
K = K + 4; % count refined elements here

end

function level = get_level(qt,e)

next_elem = e;
level = 0;
at_top = 0;
while ~at_top
    parent = qt{next_elem}.parent;
    if parent == next_elem
        at_top = 1;
    else
        next_elem = parent;
        level = level + 1;
    end
end
end
