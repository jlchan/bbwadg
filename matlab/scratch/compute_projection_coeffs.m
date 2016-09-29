% tests projection matrix
clear

N = 5;
N2 = 2*N;
[r s] = Nodes2D(N2); [r s] = xytors(r,s); [a b] = rstoab(r,s);

for j = 0:N
    E1{j+1} = bern_basis_tri(N2,r,s)\bern_basis_tri(N-j,r,s);
    E1{j+1}(abs(E1{j+1})<1e-8) = 0;
    E2{j+1} = bern_basis_tri(N,r,s)\bern_basis_tri(N-j,r,s);
    E2{j+1}(abs(E2{j+1})<1e-8) = 0;
end

% make projection
[rq sq wq] = Cubature2D(2*N2+1);
Vq1 = bern_basis_tri(N,rq,sq);
M1 = Vq1'*diag(wq)*Vq1;
Vq2 = bern_basis_tri(N2,rq,sq);
M2 = Vq1'*diag(wq)*Vq2;

norm(M1\M2-pinv(E1{1}),'fro')

L = zeros(N+1);
for j = 0:N
    Vq = bern_basis_tri(N-j,rq,sq);
    M = Vq'*diag(wq)*Vq;
    L(1:N-j+1,j+1) = sort(uniquetol(real(eig(M))),'descend');
end
Vq = bern_basis_tri(N2,rq,sq);
M = Vq'*diag(wq)*Vq;
lam = sort(uniquetol(real(eig(M))),'descend');
lam = lam(1:N+1);
c = L\lam;

T1 = bern_basis_tri(N,r,s)\Vandermonde2D(N,r,s);
T2 = bern_basis_tri(N2,r,s)\Vandermonde2D(N2,r,s);
P = 0;
for j = 0:N
    Ej{j+1} = E2{j+1}*E1{j+1}';
    P = P + c(j+1)*Ej{j+1};
end

norm(M1\M2-P,'fro') % difference between projections

return

u = (1:size(P,2))';

% E1{1}'*u
% E1{2}' + E1{3}'
sk = 1;
for j = 0:N-1
    for i = 0:N-1-j
        idi(sk,1) = i;
        idj(sk,1) = j;
        sk =sk + 1;
    end
end
ids = [idi idj];

P2 = 0;
for j = 0:N
    P2 = P2 + c(j+1)*E2{j+1}*E2{j+1}';
end

norm(P-P2*E1{1}','fro')

PP = P2'*(M1\P2); % product appears in WADG method

T = bern_basis_tri(N,r,s)\Vandermonde2D(N,r,s);
% T\(PP*T) % should be diag!

T2 = bern_basis_tri(N2,r,s)\Vandermonde2D(N2,r,s);
PP = M2'*(M1\M2);
subplot(1,2,1)
imagesc(T2\(PP*T2))

P2 = E1{1}*E1{1}';
subplot(1,2,2)
imagesc(T2\(P2*T2))

d1 = diag(T2\(PP*T2)); d2 = diag(T2\(P2*T2));
d1(abs(d1)<1e-8) = 0;
d2(abs(d2)<1e-8) = 0;

return
%% test in-place update

% P*u


T = bern_basis_tri(N,r,s)\Vandermonde2D(N,r,s);

% for N = 4
Np = (N+1)*(N+2)/2;
Np1 = (N-1+1)*(N-1+2)/2;
Np2 = (N-2+1)*(N-2+2)/2;
Np3 = (N-3+1)*(N-3+2)/2;
EE1 = bern_basis_tri(4,r,s)\bern_basis_tri(3,r,s);
EE2 = bern_basis_tri(3,r,s)\bern_basis_tri(2,r,s);
EE3 = bern_basis_tri(2,r,s)\bern_basis_tri(1,r,s);
EE4 = bern_basis_tri(1,r,s)\bern_basis_tri(0,r,s);

% P3 = c(1)*eye(size(P,1)) + c(2)*E21*E21' + c(3)*E21*E10*E10'*E21';
P3 = c(1)*eye(Np) + c(2)*EE1*(eye(Np1) + c(3)/c(2)*EE2*(eye(Np2) + c(4)/c(3)*EE3*(eye(Np3) + c(5)/c(4)*EE4*EE4')*EE3')*EE2')*EE1';

u = randn(size(P3,2),1);
b = c(1)*u;
uF = sum(u); % reduce down 
uF = EE3'*EE2'*EE1'*u + c(N+1)/c(N)*EE4*(uF); % apply small operation at lowest level
uF = EE2'*EE1'*u + c(4)/c(3)*EE3*uF;  % elevate up
uF = EE1'*u + c(3)/c(2)*EE2*uF;  % elevate up
b = b + c(2)*EE1*uF;

norm(P3-P2,'fro')
norm(P2*u-b)

return


%%
% P*u - P2*E1{1}'*u
% P2*E1{1}'
u = randn(size(P,1),1);

b = c(1)*u;
for j = 1:N
    Njp1 = (N-j+1+1)*(N-j+1+2)/2;
    Njp2 = (N-j+1)*(N-j+2)/2;
    E = bern_basis_tri(N-j+1,r,s)\bern_basis_tri(N-j,r,s);
    if j==1
        cc = c(j+1);
    else
        cc = c(j+1)/c(j);
    end
    Etu = E'*u(1:Njp1);
    b(1:Njp2) = b(1:Njp2) + cc*E*Etu;
    u(1:Njp2) = cc*E*Eu;
end

b,P2*u
% d = sqrt(diag(P2));
% imagesc(P2)

%% 3D


