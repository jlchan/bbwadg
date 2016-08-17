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
Emat = zeros(Np, Nfp);

% face 1
[r1D w1D] = JacobiGQ(0,0,N);
V1D = bern_basis_1D(N, r1D);
Mf = V1D'*diag(w1D)*V1D;
Emat(Fmask(:,1),1:Nfp) = Mf;

LIFTf= invM*Emat;
LIFTf(abs(LIFTf)<1e-8) = 0;

%% make LDL

EE = {};
for i = 0:N
    Vnmi = bern_basis_1D(N-i,r1D);
    M1D{i+1} = Vnmi'*diag(w1D)*Vnmi;
    for j = 0:i
        Vnmj = bern_basis_1D(N-j,r1D);        
        Eij = Vnmj\Vnmi; 
        Eij(abs(Eij)<1e-8) = 0;        
        EE{i+1,j+1} = Eij;
    end
    for j = i+1:N
        EE{i+1,j+1} = zeros(N-j+1,N-i+1);
    end
end
MD = blkdiag(M1D{:});
L = cell2mat(EE')';
% EE = blkdiag(E{:});
return
% Mblk = kron(Mf,Mf);

%Mblk = kron(ones(size(Mf)),Mf);
C = 1./[3 15/2 30;15/2 15/2 15;30 15 15];
Mblk = kron(C,Mf);

d = [1 -1 1/3]';
LIFTf = EE'*kron(d,LIFTf(1:N+1,1:N+1));
e1 = zeros(N+1,1); e1(1) = 1;
Emat = EE'*kron(e1,Mf);
% norm(Emat - EE'*kron(e1,Mf),'fro')

norm(EE'*kron(C,Mf)*EE*EE'*kron(d,LIFTf(1:N+1,1:N+1)) - EE'*kron(e1,Mf),'fro')

% get ijk ids
sk = 1;
for i = 0:N
    for j = 0:N-i
        k = N-i-j;
        idi(sk) = i; idj(sk) = j; idk(sk) = k;
        sk = sk + 1;
    end
end
idi = idi(id); idj = idj(id); idk = idk(id);
ids = [idi; idj; idk];



return

[rp sp] = EquiNodes2D(75); [rp sp] = xytors(rp,sp);
Vp = bern_basis_tri(N,rp,sp);
% Vp = Vp(:,id);
for i = 1:size(Vp,2)
    vp = Vp(:,i);
    clf;color_line3(rp,sp,vp,vp,'.');pause
end


