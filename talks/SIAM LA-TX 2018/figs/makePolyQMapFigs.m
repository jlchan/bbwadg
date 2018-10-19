Globals2D

N = 5;
[r s] = Nodes2D(N); [r s] = xytors(r,s);
[rq sq] = Cubature2D(2*N+1);
V = Vandermonde2D(N,r,s); 
Vq = Vandermonde2D(N,rq,sq)/V;

[rq1D wq1D] = JacobiGQ(0,0,N+1);
rfq = [rq1D; -rq1D; -ones(size(rq1D))];
sfq = [-ones(size(rq1D)); rq1D; -rq1D];
wfq = [wq1D;wq1D;wq1D];
Vfq = Vandermonde2D(N,rfq,sfq)/V;

u = .75+.5*r.^3-r.*s.^2;

[rp sp] = EquiNodes2D(25); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
% color_line3(rp,sp,Vp*u,Vp*u,'.')

PlotSol2D(rp,sp,Vp*u)
hold on
plot([-1 1 -1 -1],[-1 -1 1 -1],'k-','linewidth',3)
[az el] = view(-45,35);
axis equal
axis off
print(gcf,'-dpng','~/Desktop/bbwadg/talks/Mech2018/figs/poly.png')

%%

clf
plot([-1 1 -1 -1],[-1 -1 1 -1],'k-','linewidth',3)
hold on
h = color_line3(rq,sq,Vq*u,Vq*u,'o');%,'linewidth',3,'markersize',14,'MarkerFaceColor',[.49 1 .63])
set(h,'linewidth',3,'markersize',14,'MarkerFaceColor',[.49 1 .63])
axis off
axis equal
view(az,el)
print(gcf,'-dpng','~/Desktop/bbwadg/talks/Mech2018/figs/polyq.png')

%%
clf
plot([-1 1 -1 -1],[-1 -1 1 -1],'k-','linewidth',3)
hold on
h = color_line3(rfq,sfq,Vfq*u,Vfq*u,'o');
set(h,'marker','s','linewidth',3,'markersize',14,'MarkerFaceColor',[.49 1 .63])
axis off
axis equal
view(az,el)
print(gcf,'-dpng','~/Desktop/bbwadg/talks/Mech2018/figs/polyfq.png')