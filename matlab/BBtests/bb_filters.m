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
T = bern_basis_1D(N,r)\Vandermonde1D(N,r);

xp = Vandermonde1D(N,rp)/Vandermonde1D(N,r)*x;
xe = Vandermonde1D(N,re)/Vandermonde1D(N,r)*x;
xq = Vandermonde1D(N,rq)/Vandermonde1D(N,r)*x;

f = @(r) (r > -2/3).*(r < -1/3) + (1-abs(4*(r-1/3))).*(abs(r-1/3) < 1/4);
u = V\f(x); 

plot(xp(:),f(xp(:)),'b-','linewidth',2)
hold on
up = Vp*u;
% plot(xp(:),up(:),'r-.','linewidth',2)
sum(u,1)
W = Ve;
% W = get_BB_P1smoother(N);
% W = get_BB_smoother(N);

% a = 5; s = 2; nc = 2; n = linspace(0,1,N+1);
% d = exp(-a*((n-n(nc))./(1-n(nc))).^s);  d(1:nc) = 1;
% % clf;plot(n,d);return
% W = T*diag(d)/T;

[TVK] = TV1D(u);
TV = sum(TVK(:));
ids = find(TVK > .1*max(TVK));
% ids = 1:K;
TVK = 1.25*TVK/max(TVK);

% umodal = T\u; ulin = T(:,1:2)*umodal(1:2,:);

% u(:,ids) = f(xe(:,ids));
u(:,ids) = W*u(:,ids); 
% u(:,ids) = W*(u(:,ids)-ulin(:,ids));
% umodal = T\u(:,ids); u(:,ids) = u(:,ids) - T(:,1:2)*umodal(1:2,:) + ulin(:,ids);
sum(u,1)
uBp = Vp*f(xe);
plot(xp(:),uBp(:),'r-.','linewidth',2)
% plot(xe(:),u(:),'o--','linewidth',2,'markersize',10,'MarkerFaceColor','g')

return
% plot([x(1,:); x(end,:)],repmat(TVK,2,1),'ko-','linewidth',2,'MarkerFaceColor','g')
%legend('Exact function','Polynomial interpolant','Bernstein coefficients','Location','NorthWest')
% legend('Exact function','Polynomial interpolant','Shock detector')%,'Location','NorthWest')
% ylim([-.5 1.5])
% plot(VX,f(VX),'x','markersize',8)

% title(sprintf('TV = %f\n',TV))
grid on
set(gca,'fontsize',16)
axis equal
% print(gcf,'-dpng','detector1D.png')


