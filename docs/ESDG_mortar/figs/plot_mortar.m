N = 3;
r1D = JacobiGQ(0,0,N);
[rq sq] = meshgrid(r1D);

rq1D = JacobiGQ(0,0,N);
e = ones(size(r1D));
r1D2 = JacobiGQ(0,0,N);
r1D2 = [-1+(1+r1D2)/2; (1+r1D2)/2];
e2 = ones(size(r1D2));

rf = e;
sf = r1D;
rf2 = e2;
sf2 = r1D2;

% rf = [r1D; e; r1D; -e];
% sf = [e; r1D; -e; r1D];
% rf2 = [r1D2; e; r1D; -e];
% sf2 = [e2; r1D; -e; r1D];

a = 1;
b = 2;
plot([-1 1 1 -1 -1],[-1 -1 1 1 -1],'k-','linewidth',2)
hold on
plot(rq,sq,'ko','linewidth',2,'markersize',12,'MarkerFaceColor',[.49 1 .63])
plot(rq,sq,'ko','linewidth',2,'markersize',16,'MarkerFaceColor',[.49 1 .63])
plot(a+rf,sf,'k-','linewidth',2)
plot(a+rf,sf,'rs','linewidth',2,'markersize',12,'MarkerFaceColor',[.49 1 .63])
plot(b+rf2(1:N+1),sf2(1:N+1),'k-','linewidth',2)
plot(b+rf2(N+2:end),sf2(N+2:end),'k-','linewidth',2)
plot(b+rf2,sf2,'bx','linewidth',2,'markersize',11,'MarkerFaceColor',[.49 1 .63])
% plot(4+rq,sq,'ko','linewidth',2,'markersize',16,'MarkerFaceColor',[.49 1 .63])
axis equal
axis off
print(gcf,'-dpng','~/Desktop/bbwadg/docs/ESDG_mortar/figs/mortar.png')