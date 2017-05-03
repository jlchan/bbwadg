[x y] = meshgrid(linspace(-1,1,500));
r = x.^2 + y.^2;
f = 1 + .5*cos(4*pi*r);
f(r > 1) = nan;
pcolor(x,y,f)
shading interp
hold on
plot(0,.5,'ro','linewidth',2,'markersize',13,'MarkerFaceColor',[.49 1 .63])
axis off
axis equal
print(gcf,'-dpng','-r150','../docs/figs/c2sphere.png')