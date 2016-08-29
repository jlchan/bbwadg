N = 3;

for i = 0:N
    [r s] = Nodes2D(N); [r s] = xytors(r,s);
    Ei{i+1} = bern_basis_tri(i+1,r,s)\bern_basis_tri(i,r,s);
    E{i+1} = bern_basis_tri(N,r,s)\bern_basis_tri(N-i,r,s);
end

rule = N;
dim_num = 2;
addpath('scratch/gm_quadrature')
[ w, x ] = gm_rule_set ( rule, dim_num, gm_rule_size ( rule, dim_num ) );
r = x(1,:); r = 2*r-1;
s = x(2,:); s = 2*s-1;
% [r s] = rstoxy(r,s);
% plot(r,s,'o')
% text(r+.05,s,num2str((1:length(r))'));
% for i = 1:length(r)
%     clf
%     plot(r,s,'o')
%     text(r+.05,s,num2str((1:length(r))'));
%     hold on;
%     plot(r(i),s(i),'x')
%     pause
% end
% return

Vq = bern_basis_tri(N,r,s);
% imagesc(Vq)
M1 = Vq'*diag(w)*Vq;
% M1 = M1./min(M1(:));
% Vq = Vq./min(Vq(:));

Pq = M1\(Vq'*diag(w));
imagesc(Pq)
% Vq([4 5 7 8],:) = Vq([7 8 4 5],:)

off = 0;
V = {};
for j = 0:N
    Nptri = (N-j+1)*(N-j+2)/2;    
    V{j+1} = Vq(off + (1:Nptri),:); V{j+1} =  V{j+1}./min(V{j+1}(:));
    off = off + Nptri;
end

%% unif nodes
N = 6;

[r s] = Nodes2D(N); [r s] = xytors(r,s);
plot(r,s,'o')
text(r+.05,s,num2str((1:length(r))'));
V = bern_basis_tri(N,r,s);
V(abs(V)<1e-8) = 0;

%%
N = 3;
[r2 s2 w2] = Cubature2D(2*N);
plot(r2,s2,'o')
text(r2+.05,s2,num2str((1:length(r2))'));
Vq2 = bern_basis_tri(N,r2,s2);
M2 = Vq2'*diag(w2)*Vq2;
M2 = M2./min(M2(:));
% norm(M1-M2,'fro')

