% function check_TV_1D

Globals1D
N = 5;

% Generate simple mesh
K1D = 16;
[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

% Initialize solver and construct grid and metric
StartUp1D;

[rq wq] = JacobiGQ(0,0,N);
re = linspace(-1,1,N+1)';
rp = linspace(-1,1,100)';
r = JacobiGL(0,0,N);
[V Vr] = bern_basis_1D(N,r);
Dr = V\Vr;
Vq = bern_basis_1D(N,rq);
Vp = bern_basis_1D(N,rp);
Ve = bern_basis_1D(N,re);

xp = Vandermonde1D(N,rp)/Vandermonde1D(N,r)*x;
xe = Vandermonde1D(N,re)/Vandermonde1D(N,r)*x;
xq = Vandermonde1D(N,rq)/Vandermonde1D(N,r)*x;

f = @(r) (r > -2/3).*(r < -1/3) + exp(-25*(r-1/3).^2);
%f = @(r) (r > 1/3).*(r-1/3) + (r < -1/3).*(-(r + 1/3));
% f = @(r) (r > -1/3).*(r < 1/3);
% f = @(r) (r > -2/3) + abs(r-2/3).*(r > 2/3) + exp(-25*(r).^2);
f = @(r) (r > -4/5).*(r < -3/5) + (1-abs(5*(r+.04))).*(abs(r+.04) < 1/5) + exp(-36*(r - 7/10).^2);
% f = @(r) sin(pi*r);
% f = @(x) 1.4 + cos(3+1.3*sin(x));
% f = @(r) (r > 0).*r;
u = V\f(x); 
% u = Dr*u;
% u = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq)*f(xq));

plot(xp(:),f(xp(:)),'b-','linewidth',2)
hold on
up = Vp*u;
plot(xp(:),up(:),'r-.','linewidth',2)


[TVK] = TV1D(u);
TV = sum(TVK(:));
TVK = 1.25*TVK/max(TVK);
plot([x(1,:); x(end,:)],repmat(TVK,2,1),'ko-','linewidth',2,'MarkerFaceColor','g')
%legend('Exact function','Polynomial interpolant','Bernstein coefficients','Location','NorthWest')
legend('Exact function','Polynomial interpolant','Shock detector')%,'Location','NorthWest')
ylim([-.5 1.5])
% plot(VX,f(VX),'x','markersize',8)

% title(sprintf('TV = %f\n',TV))
grid on
set(gca,'fontsize',16)
axis equal
% print(gcf,'-dpng','detector1D.png')

%%

% try limiter

Globals1D
N = 7;

% Generate simple mesh
K1D = 1;
[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

% Initialize solver and construct grid and metric
StartUp1D;

[rq wq] = JacobiGQ(0,0,N);
re = linspace(-1,1,N+1)';
rp = linspace(-1,1,100)';
r = JacobiGL(0,0,N);
[V Vr] = bern_basis_1D(N,r);
Dr = V\Vr;
Vq = bern_basis_1D(N,rq);
Vp = bern_basis_1D(N,rp);
Ve = bern_basis_1D(N,re);

xp = Vandermonde1D(N,rp)/Vandermonde1D(N,r)*x;
xe = Vandermonde1D(N,re)/Vandermonde1D(N,r)*x;
xq = Vandermonde1D(N,rq)/Vandermonde1D(N,r)*x;

f = @(r) r > .25;
% f = @(r) (r > 0).*r;
% f = @(x) sin(1+sin(x));

u = V\f(x); 

% u = Dr*u;
% u = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq)*f(xq));

plot(xp(:),f(xp(:)),'b-','linewidth',2)
hold on
up = Vp*u;
plot(xp(:),up(:),'r-.','linewidth',2)

% ylim([-.5 1.5])
% axis([-1.1 1.1 -.5 1.1])
xlim([-1.1 1.1])

grid on
set(gca,'fontsize',16)
% axis equal

% plot(xe(:),u(:),'o--','linewidth',2,'markersize',10,'MarkerFaceColor','g')
%legend('Exact function','Polynomial interpolant','Bernstein coefficients','Location','Best')
legend('Exact function','Bernstein quasi-interpolant','Location','Best')

% print(gcf,'-dpng','bbcoeffs1.png')
print(gcf,'-dpng','bbquasi.png')

