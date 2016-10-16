% 2d diagrams
clear
N = 3;
Np = (N+1)*(N+2)/2;
u = (1:Np)';
v = zeros(Np,1);
v(1) = 1;
% v = u;

N2 = 1;
N2p = (N2+1)*(N2+2)/2;
NM = N + N2;

[r s] = Nodes2D(NM); [r s] = xytors(r,s);
Np2 = length(r);
V = bern_basis_tri(N,r,s);
V2 = bern_basis_tri(N2,r,s);
V2N = bern_basis_tri(NM,r,s);
% spy(VM)
% VM

[re se] = EquiNodes2D(NM); [re se] = xytors(re,se);


for i = 1:N2p
    v = zeros(N2p,1);
    v(i) = 1;
    VM = V2N\(diag(V2*v)*V);
    VM(abs(VM)<1e-8) = 0;
    VM
    
    [idr idc] = find(VM);
    row_ids{i} = idr;
    
    entries = VM(abs(VM)>0);
    clf
    plot(re,se,'o');
    hold on
    plot3(re(idr),se(idr),entries,'x');
    pause
end

%% 3d diagrams
clear
N = 3;
Np = (N+1)*(N+2)*(N+3)/6;
u = (1:Np)';
v = zeros(Np,1);
v(1:2) = 1;
% v = u;

[r s t] = Nodes3D(2*N); [r s t] = xyztorst(r,s,t);
V = bern_basis_tet(N,r,s,t);

V2N = bern_basis_tet(2*N,r,s,t);
VM = V2N\(diag(V*v)*V);
VM(abs(VM)<1e-8) = 0;
spy(VM)
return

[re se te] = EquiNodes3D(2*N);

for i = 1:Np
    v = zeros(Np,1);
    v(i) = 1;
    VM = V2N\(diag(V*v)*V);
    VM(abs(VM)<1e-8) = 0;
    
    [idr idc] = find(VM);
    row_ids{i} = idr;
    %     spy(VM);
    %     clf
    %     plot3(re,se,te,'o');
    %     hold on
    %     plot3(re(idr),se(idr),te(idr),'x');
    %     pause
end

%% test 3d mult

clear
N = 7;

[re se te] = EquiNodes3D(N);
Np = length(re);
C = zeros(Np,1);
sk = 1;
for l = 0:N
    for k = 0:N-l
        for j = 0:N-k-l;
            i = N-j-k-l;
            C(sk)=factorial(N)/(factorial(i)*factorial(j)*factorial(k)*factorial(l));
            sk = sk + 1;
        end
    end
end

[re se te] = EquiNodes3D(2*N);
N2p = length(re);
C2 = zeros(N2p,1);
sk = 1;
for l = 0:2*N
    for k = 0:2*N-l
        for j = 0:2*N-k-l;
            i = 2*N-j-k-l;
            C2(sk)=factorial(2*N)/(factorial(i)*factorial(j)*factorial(k)*factorial(l));
            sk = sk + 1;
        end
    end
end

[r s t] = Nodes3D(2*N); [r s t] = xyztorst(r,s,t);
V = bern_basis_tet(N,r,s,t);
V2N = bern_basis_tet(2*N,r,s,t);
for i = 1:Np
    v = zeros(Np,1);
    v(i) = 1;
%     VM = diag(C2)*(V2N\(diag(V*(v./C))*V*diag(1./C)));
    VM = V2N\(diag(V*v)*V);
    VM(abs(VM)<1e-8) = 0;
    
    [idr idc tmp] = find(VM);
    VMvals{i} = tmp(:);
    
    row_ids{i} = idr;
end

% return

% do mult
u = (1:Np)';
v = (0:Np-1)'; 
VMex = V2N\(diag(V*v)*V);
fg_ex = VMex*u;

E = bern_basis_tet(2*N,r,s,t)\bern_basis_tet(N,r,s,t);
E(abs(E)<1e-8) = 0;
rfg_ex = E'*fg_ex;

u = u.*C; v = v.*C;
fg = zeros(size(re));
for j = 1:Np
    ids = row_ids{j};
    fg(ids) = fg(ids) + v(j)*u(:);
end

norm(fg./C2 - fg_ex)/norm(fg_ex)

% cheap but maybe unstable reduction
C2avg = mean(C2);
C2 = C2/C2avg;
fg = fg./(C2.*C2);
rfg = zeros(Np,1);
for j = 1:Np
    ids = row_ids{j};
    rfg(j) = C(j)*C' * fg(ids);
end
rfg = rfg/C2avg^2;

norm(rfg - rfg_ex)/norm(rfg_ex)


%% 2D projection

N = 3;
[rq sq wq] = Cubature2D(5*N);
VqN = bern_basis_tri(N,rq,sq);
VqM = bern_basis_tri(2*N,rq,sq);
M2 = VqN'*diag(wq)*VqM;
M1 = VqN'*diag(wq)*VqN;
PMN = M1\M2;
% PMN = VqN\VqM;

%% 3d projection
% check projection down
[rq sq tq wq] = tet_cubature(4*N);
VqN = bern_basis_tet(N,rq,sq,tq);
VqM = bern_basis_tet(2*N,rq,sq,tq);
M2 = VqN'*diag(wq)*VqM;
M1 = VqN'*diag(wq)*VqN;
PMN = M1\M2;
% TN = Vandermonde3D(N,rq,sq,tq)\VqN;
% TM = Vandermonde3D(2*N,rq,sq,tq)\VqM;
% imagesc(TN*(PMN/TM))
