% check_bern_lift
% function check_bern_lift(inN)
% if nargin==0
%     inN = 3;
% end

% Driver script for solving the 3D IPDG acoustic wave equation
Globals3D;

% N = inN;
% Order of polymomials used for approximation
N = 2;

% % % single element
[VX VY VZ] = Nodes3D(1);
[VX VY VZ] = xyztorst(VX,VY,VZ);
K = 1; EToV = 1:length(VX);

% Initialize solver and construct grid and metric
StartUp3D;

LIFTnodal = LIFT;
%% cubature and plotting
[rq sq tq w] = tet_cubature(2*N+1);
Vq = Vandermonde3D(N, rq, sq, tq)/V;
xq = Vq*x; yq = Vq*y; zq = Vq*z; % phys cubature nodes

% for plotting - build coordinates of all the nodes
[rp sp tp] = EquiNodes3D(25);
Vp = Vandermonde3D(N,rp,sp,tp)*invV;
xp = Vp*x; yp = Vp*y; zp = Vp*z;

%% bernstein conversion

Vp = bern_basis_tet(N,rp,sp,tp);
Vq = bern_basis_tet(N,rq,sq,tq);

[V Vr Vs Vt V1 V2 V3 V4 id] = bern_basis_tet(N,r,s,t);
Dr = V\Vr; Dr(abs(Dr)<1e-8) = 0;
Ds = V\Vs; Ds(abs(Ds)<1e-8) = 0;
Dt = V\Vt; Dt(abs(Dt)<1e-8) = 0;

global D1 D2 D3 D4
D1 = V\V1; D1(abs(D1)<1e-8) = 0;
D2 = V\V2; D2(abs(D2)<1e-8) = 0;
D3 = V\V3; D3(abs(D3)<1e-8) = 0;
D4 = V\V4; D4(abs(D4)<1e-8) = 0;

%% operators

M = Vq'*diag(w)*Vq;
invM = inv(M);

%% lift matrix
Emat = zeros(Np, Nfaces*Nfp);
Mf = zeros(Np, Nfp);

for face=1:Nfaces
    %     % process face
    %     if(face==1); faceR = r(Fmask(:,1)); faceS = s(Fmask(:,1)); end;
    %     if(face==2); faceR = r(Fmask(:,2)); faceS = t(Fmask(:,2)); end;
    %     if(face==3); faceR = s(Fmask(:,3)); faceS = t(Fmask(:,3)); end;
    %     if(face==4); faceR = s(Fmask(:,4)); faceS = t(Fmask(:,4)); end;
    
    %     [rt st] = Nodes2D(N); [rt st] = xytors(rt,st);
    %     [faceR,faceS]-[rt,st]
    %VFace = Vandermonde2D(N, faceR, faceS);
    
    % WARNING: ASSUMES ORDERING OF FACE NODES/FMASK
    [rt st] = Nodes2D(N); [rt st] = xytors(rt,st);
    [rfq sfq wf] = Cubature2D(2*N);
    %     VFace = Vandermonde2D(N, rt, st);
    %     Vfq = Vandermonde2D(N, rfq, sfq)/VFace;
    %     massFace = inv(VFace*VFace');
    
    Vfq = bern_basis_tri(N,rfq,sfq);
    massFace = Vfq'*diag(wf)*Vfq;
    
    idr = Fmask(:,face);
    idc = (face-1)*Nfp+1:face*Nfp;
    
    Emat(idr, idc) = Emat(idr, idc) + massFace;
    if face==1
        Mf(idr,idc)=massFace;
    end
end

% inv(mass matrix)*\I_n (L_i,L_j)_{edge_n}
LIFT = invM*Emat;

% DnLift
global LIFTr LIFTs LIFTt
LIFTr = invM*(Dr'*Emat); LIFTr(abs(LIFTr)<1e-8) = 0;
LIFTs = invM*(Ds'*Emat); LIFTs(abs(LIFTs)<1e-8) = 0;
LIFTt = invM*(Dt'*Emat); LIFTt(abs(LIFTt)<1e-8) = 0;

global nr ns nt
nr = nx.*rx(Fmask(:),:) + ny.*ry(Fmask(:),:) + nz.*rz(Fmask(:),:);
ns = nx.*sx(Fmask(:),:) + ny.*sy(Fmask(:),:) + nz.*sz(Fmask(:),:);
nt = nx.*tx(Fmask(:),:) + ny.*ty(Fmask(:),:) + nz.*tz(Fmask(:),:);

global rr rs rt ss st tt
rr = rx.*rx + ry.*ry + rz.*rz;
rs = rx.*sx + ry.*sy + rz.*sz;
rt = rx.*tx + ry.*ty + rz.*tz;
ss = sx.*sx + sy.*sy + sz.*sz;
st = sx.*tx + sy.*ty + sz.*tz;
tt = tx.*tx + ty.*ty + tz.*tz;

LIFT(abs(LIFT)<1e-8) = 0;

%% decompose lift operators

LIFTf = invM*Mf; LIFTf(abs(LIFTf)<1e-8) = 0;
% get permutations of rows for 2nd-4th cols of LIFT
p = zeros(Np,Nfaces-1);
for f = 1:Nfaces-1
    ids = (1:Nfp) + f*Nfp;
    for i = 1:Np % match rows of first cols with second
        diff = kron(ones(Np,1),LIFTf(i,:)) - LIFT(:,ids);
        [~,j] = min(sum(abs(diff),2));
        p(j,f) =i;
    end
end

[r2D s2D] = Nodes2D(N); [r2D s2D] = xytors(r2D,s2D);
Vtri = bern_basis_tri(N,r2D,s2D);

% many ops
Eblk = {};
for i = 0:N
    for j = 0:i
        EDi = bern_basis_tri(N-j,r2D,s2D)\bern_basis_tri(N-i,r2D,s2D);
        EDi(abs(EDi)<1e-8) = 0;
        Eblk{i+1,j+1} = EDi;
    end
end


% 2d degree ops
ED = {};
for i = 0:N
    EDi = bern_basis_tri(N,r2D,s2D)\bern_basis_tri(N-i,r2D,s2D);
    EDi(abs(EDi)<1e-8) = 0;
    ED{i+1} = EDi;
end
LIFTf = LIFT(:,1:Nfp);
EE = blkdiag(ED{:}); EE = EE';

% get l_i scalings
R = EE*kron(ones(N+1,1),LIFTf(1:Nfp,1:Nfp))./LIFTf;  % ratio
block_starts = []; id = 1;
for i = 0:N
    block_starts = [block_starts id];
    Ni = N-i; Nfpi = (Ni+1)*(Ni+2)/2; % decreasing block sizes
    id = id + Nfpi;
end
R = R(:,1);
l = 1./R(block_starts); %R = R(~isnan(R)); [~, perm] = uniquetol(R); l = 1./R(sort(perm)); % scaling constants
EDL = cell(N+1,1);
for i = 0:N
    EDL{i+1} = ED{i+1}*l(i+1);
end

% built in scalings
EL = [];
for i = 0:N
    EL = [EL;EDL{i+1}'];
end
EEL=([EL EL(p(:,1),:) EL(p(:,2),:) EL(p(:,3),:)]);
% EL1 = EL;
% EL2 = EL(p(:,1),:);
% EL3 = EL(p(:,2),:);
% EL4 = EL(p(:,3),:);

%% check application of LIFT
u = randn(size(LIFT,2),1);
uL = LIFT*u;

% bernstein-specific application of LIFT
L0 = LIFTf(1:Nfp,1:Nfp);

L0u = L0*reshape(u,Nfp,Nfaces); uL2 = EEL*L0u(:);
norm(uL-uL2)

% %% explicit L0 form
% T = inv(Vandermonde2D(N,r2D,s2D))*bern_basis_tri(N,r2D,s2D); % map from bern to modal
% 
% [rq sq wq] = Cubature2D(2*N);
% 
% for i = 0:N
%     Vq2D = bern_basis_tri(N-i,rq,sq);
%     M2D = Vq2D'*diag(wq)*Vq2D;
%     Mblk{i+1} = M2D;
% end
% 
% % get matrix of block scaling coeffs
% tmp = M./(EE*kron(ones(N+1),Mblk{1})*EE');
% 
% C = zeros(N+1);
% r = 0;
% for i = 0:N
%     Ni = (N-i);
%     r = r + (Ni+1)*(Ni+2)/2;
%     c = 0;
%     for j = 0:N
%         Nj = (N-j);
%         c = c + (Nj+1)*(Nj+2)/2;
%         C(i+1,j+1) = tmp(r,c);
%     end
% end
% [L D P] = ldl(C);
% iL = inv(L);
% iD = inv(D);
% 
% 
% LL = L0*0;
% for i = 0:N
%     %     keyboard
%     %iL(i+1,1)^2/D(i+1,i+1)
%     LL = LL + (N-i+1.5)*Eblk{i+1,1}*(Mblk{i+1}\(Eblk{i+1,1}'*Mblk{1}));
% end
% norm(L0-LL,'fro')

u = rand(4*Nfp,1);

b1 = LIFT*u;
L0u = L0*reshape(u,Nfp,Nfaces); b2 = EEL*L0u(:);

norm(b1-b2)

% b3 = V*LIFTnodal*(V\u);
% norm(b1-b3)
% norm(b2-b3)

% iL(:,1).^2./diag(D)
return









%% sparser degree elevation ops
% 2d degree ops
ED = {}; ELblk = {};
for i = 0:N
    Vnmi = bern_basis_tri(N-i,r2D,s2D);
    for j = 0:i
        Vnmj = bern_basis_tri(N-j,r2D,s2D);
        Eij = Vnmj\Vnmi;
        Eij(abs(Eij)<1e-8) = 0;
        ELblk{i+1,j+1} = Eij;
    end
    % make top zeros
    for j = i+1:N
        ELblk{i+1,j+1} = zeros(N-j+1,N-i+1);
    end
    
    EDi = bern_basis_tri(N,r2D,s2D)\bern_basis_tri(N-i,r2D,s2D);
    EDi(abs(EDi)<1e-8) = 0;
    ED{i+1} = EDi;
end
LIFTf = LIFT(:,1:Nfp);
EE = blkdiag(ED{:}); EE = EE';

return












% 1d interpolants
E = {};
for i = 0:N
    Vnmi = bern_basis_1D(N-i,r1D);
    for j = 0:i
        Vnmj = bern_basis_1D(N-j,r1D);
        Eij = Vnmj\Vnmi; Eij(abs(Eij)<1e-8) = 0;
        E{i+1,j+1} = Eij;
    end
    for j = i+1:N
        E{i+1,j+1} = zeros(N-j+1,N-i+1);
    end
end



% get ijkl ids
sk = 1;
for i = 0:N
    for j = 0:N-i
        for k = 0:N-i-j
            l = N-i-j-k;
            idi(sk) = i;
            idj(sk) = j;
            idk(sk) = k;
            idl(sk) = l;
            sk = sk + 1;
        end
    end
end
idi = idi(id);
idj = idj(id);
idk = idk(id);
idl = idl(id);
ids = [idi;idj;idk;idl]
ids = [idi;idj;idk;idl]+1;

% 2*LIFT(:,1:Nfp)

return

[rows,cols,vals] = find(LIFT(:,1:Nfp));
[foo,sids] = sort(vals, 'ascend');
rows = rows(sids);
cols = cols(sids);
vals = vals(sids);

row = ids(:,rows)-1;
col = ids(:,cols)-1;
val = LIFT(rows,cols)';
ch = 65*ones(size(row(1,:),2), 1);
aa = [row(1,:)',row(2,:)',row(3,:)', row(4,:)', ch, col(1,:)',col(2,:)', col(3,:)',col(4,:)', ch, vals(:)]

%     [svals,sids] = sort(vals, 'ascend');

%     for j = 1:length(inds)
%         str = sprintf('row = (%d,%d,%d), col = (%d,%d,%d), value = %2.2f',row(1),row(2),row(3),cols(1,sids(j)),cols(2,sids(j)),cols(3,sids(j)),LIFT(i,inds(sids(j))));
%     end
%     disp(str)
%     str = sprintf('row = (%d,%d,%d),',row(1),row(2),row(3));
%     for j = 1:length(inds);
%         str = strcat(str,sprintf('col %d = (%d,%d,%d), value = %2.2f,\n',j,cols(1,j),cols(2,j),cols(3,j),LIFT(i,inds(j))));
%     end



