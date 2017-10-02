clear
load inclusion_aligned_N5K50.mat
% load inclusion_unaligned_N5K50.mat

cMap = gray(256); % invert
%dataMax = 3.25;    centerPoint = 2.5; scalingIntensity = 5;
dataMax = .5;    centerPoint = .65;  scalingIntensity = 5;
dataMin = 0;

xx = 1:length(cMap);
xx = xx - (centerPoint-dataMin)*length(xx)/(dataMax-dataMin);
xx = scalingIntensity * xx/max(abs(xx));
xx = sign(xx).* exp(abs(xx));
xx = xx - min(xx); xx = xx*511/max(xx)+1;
newMap = interp1(xx, cMap, 1:512);
% plot(newMap);return

figure
PlotField2D(N+1,x,y,U{5});
axis equal

colormap(flipud(newMap))
caxis([0,.5])

view(2)
colorbar('TickLabelInterpreter','latex')
set(gca,'fontsize',15)

% axis([-.25 .05 .0 .2])
% print(gcf,'-dpng','-r300','inclusion_aligned_zoom.png')
print(gcf,'-dpng','-r300','inclusion_aligned_shear.png')
% print(gcf,'-dpng','-r300','inclusion_aligned.png')
% print(gcf,'-dpng','-r300','inclusion_unaligned.png')
