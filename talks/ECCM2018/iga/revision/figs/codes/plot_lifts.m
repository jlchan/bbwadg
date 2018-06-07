NB = 7;
Ksub = NB;
N = NB+Ksub-1;
[Nv, VX, Ksub, EToV] = MeshGen1D(-1,1,Ksub);
map = @(r) reshape((1/Ksub)*(repmat(r,1,Ksub)+1) + ...
    repmat(VX(1:end-1),length(r),1),length(r)*Ksub,1);

% local quadrature
[rqEx, wqEx] = JacobiGQ(0,0,NB+1); % overkill gauss quadrature
rqEx = map(rqEx);
wqEx = repmat(wqEx,1,Ksub)*(1/Ksub); 
wqEx = wqEx(:); 

t = [VX(1)*ones(1,NB) VX VX(end)*ones(1,NB)]; % open knot vec

% exact quad bspline
VqEx = bspline_basismatrix(NB+1,t,rqEx);

[rq wq] = JacobiGL(0,0,N);

Vq = bspline_basismatrix(NB+1,t,rq);
M = VqEx'*diag(wqEx)*VqEx;
Mh = Vq'*diag(wq)*Vq;
% norm(M-Mh,'fro')

invM = inv(M);

rp = linspace(-1,1,1000)';
Bp = bspline_basismatrix(NB+1,t,rp);

e = zeros(N+1,1); e(1) = 1; %e(N+1) = 1;
Mf = diag(e);
[W D] = eig(Mf,M);
[~,p] = sort(diag(D),'descend'); W = W(:,p);
LB = W(:,1)/W(1,1);


plot(rp,Bp*LB,'linewidth',2);
% plot(rp,Vp*invM(:,:));
hold on

L = Vandermonde1D(N,JacobiGR(0,0,N))\[1;zeros(N,1)];
plot(rp,Vandermonde1D(N,rp)*L,'-.','linewidth',2)
plot(JacobiGR(0,0,N),zeros(N+1,1),'ko','markersize',10)

h = legend('Spline lift','Polynomial lift','Gauss-Radau points');
set(h,'fontsize',16)

set(gca,'fontsize',15)
grid on

figdir = '~/Desktop/bbwadg/docs/IGA-DG/figs/';
print(gcf,'-dpng',[figdir 'lift1.png'])

ylim([-.15 .15])
print(gcf,'-dpng',[figdir 'liftZoom1.png'])

%========================================
clf
plot(rp,Bp*LB,'linewidth',2);
% plot(rp,Vp*invM(:,:));
hold on

N = NB;
L = Vandermonde1D(N,JacobiGR(0,0,N))\[1;zeros(N,1)];
plot(rp,Vandermonde1D(N,rp)*L,'-.','linewidth',2)
plot(JacobiGR(0,0,N),zeros(N+1,1),'ko','markersize',10)

h = legend('Spline lift','Polynomial lift','Gauss-Radau points');
set(h,'fontsize',16)

set(gca,'fontsize',15)
grid on

print(gcf,'-dpng',[figdir 'lift2.png'])

ylim([-.15 .15])
print(gcf,'-dpng',[figdir 'liftZoom2.png'])