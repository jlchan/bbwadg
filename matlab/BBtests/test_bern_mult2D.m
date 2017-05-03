N = 4;
M = 2;
[rq wq] = JacobiGQ(0,0,2*N);
Vq = bern_basis_1D(N,rq);
VM = bern_basis_1D(M,rq);
Vq2 = bern_basis_1D(N+M,rq);
MM = (Vq2'*diag(wq)*Vq2);
Pq = MM \ (Vq2'*diag(wq));

rp = linspace(-1,1,500);
Vp = bern_basis_1D(M,rp);
for i = 2
    L = Pq*(diag(VM(:,i))*Vq);
    L(abs(L)<1e-8) = 0;
    [i j] = find(L);
    re = linspace(-1,1,N+M+1)';
    
    plot(re,re*0,'bo')
    hold on
    plot(re(i),L(abs(L)>1e-8),'ro--','markersize',12)
    
    plot(re(i),L(abs(L)>1e-8)./((1-re(i))/2),'x--')
    plot(re(i),L(abs(L)>1e-8)./((1+re(i))/2),'x--')     
%     plot(rp,Vp,'--')
end

%% 2D

N = 3;
M = 2;
[re se] = EquiNodes2D(N+M); [re se] = xytors(re,se);
[rq sq wq] = Cubature2D(2*(N+M));
Vq = bern_basis_tri(N,rq,sq);

VA = Vandermonde2D(M,rq,sq);
% VA =
VM = bern_basis_tri(M,rq,sq)*(bern_basis_tri(M,rq,sq)\VA);
Vq2 = bern_basis_tri(N+M,rq,sq);
M = (Vq2'*diag(wq)*Vq2);
Pq = M \ (Vq2'*diag(wq));

for i = 1
    L= Pq*(diag(VM(:,i))*Vq);
    L(abs(L)<1e-8) = 0;
    
    [i j] = find(L);
    
    plot3(re,se,zeros(size(re)),'b*')
    hold on
    plot3(re(i),se(i),L(abs(L)>1e-8),'o','markersize',12)
end

