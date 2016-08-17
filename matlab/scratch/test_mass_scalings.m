% test bb mass matrix scaling
N = 2;
[rq sq tq w] = tet_cubature(2*N+1);
Vq = bern_basis_tet(N,rq,sq,tq);
M = Vq'*diag(w)*Vq;

% 2D ====================

[rq sq w] = Cubature2D(2*N);
Vq = bern_basis_tri(N,rq,sq);
M2 = Vq'*diag(w)*Vq;

[r s] = Nodes2D(N); [r s] = xytors(r,s);
Nfp = length(r);
E2{1} = eye(Nfp);
for i = 1:N
    E2{i+1} = bern_basis_tri(N,r,s)\bern_basis_tri(N-i,r,s);
end

% 1D =====================

[rq w] = JacobiGQ(0,0,N);
Vq = bern_basis_1D(N,rq);
M1 = Vq'*diag(w)*Vq;

Nfp = length(r);
E1{1} = eye(N+1);
for i = 1:N
    E1{i+1} = bern_basis_1D(N,rq)\bern_basis_1D(N-i,rq);
end

sk = 1;
E12 = {};
for k = 0:N
    for j = 0:N-k
        E12{sk} = E1{j+k+1};
        sk = sk + 1;
    end
end
E = blkdiag(E12{:});

M1 = E'*kron(ones(Nfp),M1)*E;
CM = M./M1;

ids = [];
sk = 1;
off = 0;
for k = 0:N
    for j = 0:N-k
        ids(sk) = off + 1;
        off = off + (N - j - k + 1);
        sk = sk + 1;
    end
end
C = CM(ids,ids);


CM = diag(1./(C(1,:)))*C%*diag(1./sqrt(C(1,:)));
% for i = 2:size(C,1);
%    CM(i,:) =(C(1,:)./C(i,1)).*C(i,:);
% end
