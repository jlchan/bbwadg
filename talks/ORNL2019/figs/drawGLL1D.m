N = 7;
r = JacobiGL(0,0,N);
rp = linspace(-1,1,150)';
V = Vandermonde1D(N,r);
Vp = Vandermonde1D(N,rp)/V;
u = exp(-sin(pi*r));

plot(rp,Vp*u,'k-','linewidth',4)
hold on
plot(r,u,'ko','linewidth',4,'markersize',16,'MarkerFaceColor',[.49 1 .63])
for i = 1:N+1
    plot(rp,Vp(:,i)*u(i),'k--','linewidth',2)
end
grid on
set(gca,'fontsize',15)
xlim([-1.1 1.1 ])
set(gca,'xticklabel',[],'yticklabel',[])

% print(gcf,'-dpng','~/Desktop/bbwadg/talks/Mech2018/figs/gll1D.png')