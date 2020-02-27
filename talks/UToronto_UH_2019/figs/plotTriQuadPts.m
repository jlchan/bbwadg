N = 3;
[rq sq wq] = Cubature2D(2*N);

[rq1D wq1D] = JacobiGQ(0,0,N);
rf = [rq1D; -rq1D; -ones(size(rq1D))];
sf = [-ones(size(rq1D)); rq1D; -rq1D];

[rq sq] = rstoxy(rq,sq);
[rf sf] = rstoxy(rf,sf);

r1 = [-1 1 -1 -1];
s1 = [-1 -1 1 -1];
[r1 s1] = rstoxy(r1,s1);
%plot([-1 1 -1 -1],[-1 -1 1 -1],'k-','linewidth',2)
plot(r1,s1,'k-','linewidth',2)
hold on
plot(rq,sq,'ko','linewidth',2,'markersize',14,'MarkerFaceColor',[.49 1 .63]);
plot(rf,sf,'rs','linewidth',2,'markersize',14,'MarkerFaceColor',[.49 1 .63]);
axis off
% set(h,'markersize',64)
axis equal
print(gcf,'-dpng','~/Desktop/bbwadg/talks/NAHOM_2019/figs/triQuadrature.png')

return

load ~/Downloads/QuadratureRules.mat

quadpts = Q_GaussLegendre{N};
[r1D,w1D] = JacobiGQ(0,0,N);
% quadpts = Q_GaussLobatto{N};
% [r1D,w1D] = JacobiGL(0,0,N+1); % 2*(N+1)-1 = 2*N+1

rs = quadpts.Points;
wq = quadpts.Weights;

rq = rs(:,1);
sq = rs(:,2);

figure
plot([-1 1 -1 -1],[-1 -1 1 -1],'k-','linewidth',2)
hold on
plot(rq,sq,'ko','linewidth',2,'markersize',14,'MarkerFaceColor',[.49 1 .63]);
axis off