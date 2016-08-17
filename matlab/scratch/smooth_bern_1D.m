N = 3;

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(-1,1,8);

StartUp1D
map = @(r) ones(length(r),1)*VX(va) + .5*(r(:)+1)*(VX(vb)-VX(va));

V = bern_basis_1D(N,r);

f = @(x) sin(3*pi*x+2);
% f = @(x) -x.^3+x.^2+x; 

[rq w] = JacobiGQ(0,0,N+2); 
xq = map(rq);
Vq = bern_basis_1D(N,rq);
u = (Vq'*diag(w)*Vq)\(Vq'*diag(w)*f(xq));
% u = V\f(x);
% plot(xp,f(xp),'r--')
% u = rand(size(x));

% plotting
Nplot = 100;
rp = linspace(-1,1,Nplot); rp = rp(:);
xp = map(rp);
Vp = bern_basis_1D(N,rp);

% rg = equispaced pts
rg = linspace(-1,1,N+1); rg = rg(:);
xg = map(rg);

vv = Vp*u;
% plot(rp,Vp);return
hold on
plot(xp,vv-f(xp),'-','linewidth',2)
return

% control points for C1 continuity
ids = [(N+1)*(1:K-1)-1; (N+1)*(1:K-1)+2];

plot(xg(ids), u(ids),'-','linewidth',2)
plot(xg,u,'o')
