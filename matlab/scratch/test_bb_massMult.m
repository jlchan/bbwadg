% test O(N^{d+1}) mass matrix multiply

N = 2
[r s t] = EquiNodes3D(N);
Np = length(r);
plot3(r,s,t,'o','markersize',16)
text(r+.1,s,t,num2str((1:Np)'))

[rq sq tq w] = tet_cubature(2*N);
Vq = bern_basis_tet(N,rq,sq,tq);

% ==== tet mass matrix
M = Vq'*diag(w)*Vq;
M = M./min(M(:));

% ==== make 1D mats

[rq1D w1D] = JacobiGQ(0,0,N);

E = {};
M1D = {};
for i = 0:N
    Vq1D = bern_basis_1D(N-i,rq1D);
    M1D{i+1} = Vq1D'*diag(w1D)*Vq1D;
    if i==0
        E{i+1} = eye(N+1);
    else
        E{i+1} = bern_basis_1D(N,rq1D)\bern_basis_1D(N-i,rq1D);
    end
end

%% test matmult

M1 = M1D{1};
M1 = M1./min(M1(:));

% get scaling factors
Nfp = (N+1)*(N+2)/2;
c = zeros(Nfp);
sk = 1;
off = 0;
for k = 0:N
    for j = 0:N-k
        ids = off + (1:N-k-j+1);                
        off2 = 0;
        sk2 = 1;
        for kk = 0:N
            for jj = 0:N-kk
                row_ids = off2 + (1:N-kk-jj+1);
                C = M(row_ids,ids)./(E{kk+jj+1}'*M1*E{k+j+1});
                c(sk2,sk) = C(1);
                sk2 = sk2 + 1;
                off2 = off2 + (N-kk-jj+1);
            end
        end        
        
        sk = sk + 1;        
        off = off + (N-k-j+1);
    end
end


% u = (1:Np)';
u = randn(Np,1);
b1 = M*u;

% O(N^{d+1}) application of mass mult
b2 = zeros(Np,1);
sk = 1;
off = 0;
for k = 0:N
    for j = 0:N-k
        
        ids = off + (1:N-k-j+1);
        
        % applying M1 = O(N^2). 
        % apply E as a product of sparse 1D operators (2 nnz per row)
        btmp = M1*E{k+j+1}*u(ids); 
        
        % accumulate result in rhs vector
        off2 = 0;
        sk2 = 1;        
        for kk = 0:N
            for jj = 0:N-kk
                row_ids = off2 + (1:N-kk-jj+1);   
                
                % Will need to apply E as a product of sparse 1D operators.                 
                b2(row_ids) = b2(row_ids)  + c(sk2,sk)*E{kk+jj+1}'*btmp;
                sk2 = sk2 + 1;
                off2 = off2 + (N-kk-jj+1);
            end
        end        
        
        % increment - can replace with map
        off = off + (N-k-j+1);
        sk = sk + 1;
    end
end

norm(b1-b2)