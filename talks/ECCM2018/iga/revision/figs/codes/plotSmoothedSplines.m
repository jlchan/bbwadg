NB = 4;
Ksub = 4;
smoothKnots = 'opt';
rp = linspace(-1,1,2500)';

[Vp M Dr R rBq wBq Bq Brq VX] = bsplineVDM(NB,Ksub,rp,smoothKnots);
plot(rp,Vp)
hold on
plot(VX,VX*0,'ko','markersize',16,'MarkerFaceColor',[.49 1 .63])
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca,'fontsize',14)
set(gca,'xticklabel','')
set(gca,'yticklabel','')
grid on

% print(gcf,'-dpng','../talks/UT2017/figs/unifknots.png')
% print(gcf,'-dpng','../talks/UT2017/figs/optknots.png')


