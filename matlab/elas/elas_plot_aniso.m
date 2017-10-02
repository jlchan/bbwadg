% load anisoN5K64
load anisoN5K128

p = U{2};

vv = Vp*p;
%         vv = abs(vv);
% color_line3(xp,yp,vv,vv,'.');
% color_line3(x,y,U{2},U{2},'.');
PlotField2D(2*N+1,x,y,U{2});
axis tight
axis equal
axis on
set(gca,'fontsize',15)
colormap(gray)
colorbar
caxis([-.1,.2])
view(2)
% print(gcf,'-dpng','aniso1.png')
print(gcf,'-dpng','aniso2.png')