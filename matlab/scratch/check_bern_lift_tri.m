% check bern lift triangle
clear
Globals2D

N = 3;

[VX VY] = Nodes2D(1); [VX VY] = xytors(VX,VY);
Nfp = N+1; Np = (N+1)*(N+2)/2; Nfaces=3; NODETOL = 1e-12;
K = 1; EToV = 1:3;

StartUp2D


[rq sq w] = Cubature2D(2*N);
[Vq, ~, ~, ~, ~,~,id]  = bern_basis_tri(N,rq,sq);
% Vq = Vq(:,id);
M = Vq'*diag(w)*Vq;
invM = inv(M);

%%
Emat = zeros(Np, Nfp*Nfaces);

% face 1
[r1D w1D] = JacobiGQ(0,0,N);
V1D = bern_basis_1D(N, r1D);
Mf = V1D'*diag(w1D)*V1D;
Emat(Fmask(:,1),1:Nfp) = Mf;

Emm = zeros(Np,Nfp);
Emm(Fmask(:,1),1:Nfp) = Mf;
LIFTf= invM*Emm;
LIFTf(abs(LIFTf)<1e-8) = 0;

% face 2
faceR = r(Fmask(:,2));
Emat(Fmask(:,2),Nfp+1:2*Nfp) = Mf;

% face 3
faceS = s(Fmask(:,3));
Emat(Fmask(:,3),2*Nfp+1:3*Nfp) = Mf;
%%
Mf = Emat(:,1:Nfp);
LIFT = invM*Emat; LIFT(abs(LIFT)<1e-8) = 0;

E = {};
MM{1} = V1D'*diag(w1D)*V1D;
for i = 0:N
    Vnmii = bern_basis_1D(N-i+1,r1D);
    Vnmi = bern_basis_1D(N-i,r1D);
    E{i+1} = V1D\Vnmi; E{i+1}(abs(E{i+1})<1e-8) = 0;
    MM{i+1} = Vnmi'*diag(w1D)*Vnmi;
end
EE = blkdiag(E{:});

% EE'*kron(ones(N+1,1),LIFT(1:Nfp,1:Nfp))./LIFT(:,1:Nfp)

Mblk = {};
Eblk = {};
for i = 0:N
    Vnmi = bern_basis_1D(N-i,r1D);
    Mblk{i+1} = Vnmi'*diag(w1D)*Vnmi;
    for j = 0:i
        Vnmj = bern_basis_1D(N-j,r1D);        
        Eij = Vnmj\Vnmi; 
        Eij(abs(Eij)<1e-8) = 0;        
        Eblk{i+1,j+1} = Eij';
    end
    for j = i+1:N
        Eblk{i+1,j+1} = zeros(N-i+1,N-j+1);
    end
end

% get matrix of block scaling coeffs
tmp = M./(EE'*kron(ones(N+1),Mblk{1})*EE);
r = 0; 
for i = 0:N
    r = r + (N-i+1);
    c = 0;
    for j = 0:N        
        c = c + (N-j+1);
        C(i+1,j+1) = tmp(r,c);        
    end
end
[L D PL] = ldl(C);
iL = inv(L);
iD = inv(D);

% assemble LDL
Mnew = {}; ELnew = {};
for i = 1:N+1
    Mnew{i} = Mblk{i}*D(i,i);
    for j = 1:N+1
        ELnew{i,j} = Eblk{i,j}*L(i,j);
    end
end
EL = cell2mat(ELnew);
MD = blkdiag(Mnew{:});

% iL = inv(EL);

% get l_i scalings
R = EE'*kron(ones(N+1,1),LIFTf(1:Nfp,1:Nfp))./LIFTf;
R = R(:,1);
off = 0;
for i = 1:N+1    
    l(i) = R(i+off);
    off = off + (N+1-i);
end
l = 1./l(:);

% compute L0 (check)
L0 = zeros(N+1);
P = {};
for i = 0:N;     
    %L0 = L0 + (N+1-i)*Eblk{i+1,1}'*(Mblk{i+1}\(Eblk{i+1,1}*Mblk{1}));    
    L0 = L0 + (N+1-i)*Eblk{i+1,1}'*(Mblk{i+1}\(Eblk{i+1,1}*Mblk{1}));    
    P{i+1} = Eblk{i+1,1}'*(Mblk{i+1}\(Eblk{i+1,1}*Mblk{1}));
end

r1D = JacobiGL(0,0,N);
T1 = inv(Vandermonde1D(N,r1D))*bern_basis_1D(N,r1D); % map from bern to modal
r1D = JacobiGL(0,0,N-1);
T2 = inv(Vandermonde1D(N-1,r1D))*bern_basis_1D(N-1,r1D); % map from bern to modal
iL0 = zeros(N+1);
iL1 = zeros(N+1);
iL2 = zeros(N+1);
iL3 = zeros(N+1);
for j = 1:N+1
    EP{j} = Eblk{j,1}'*Eblk{j,1};
    iL0 = iL0 + C(1,j)*l(j)*Eblk{j,1}'*Eblk{j,1};
    iL1 = iL1 + C(2,j)*l(j)*Eblk{j,1}'*Eblk{j,1};
    iL2 = iL2 + C(3,j)*l(j)*Eblk{j,1}'*Eblk{j,1};
%     iL3 = iL3 + C(4,j)*l(j)*Eblk{j,1}'*Eblk{j,1};
end

return








