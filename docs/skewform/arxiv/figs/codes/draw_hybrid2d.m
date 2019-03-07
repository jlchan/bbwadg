N = 3;

r1D = JacobiGL(0,0,N);
[r s] = meshgrid(r1D);

plot([-1 1 1 -1 -1],[-1 -1 1 1 -1],'k-','linewidth',2)
hold on
plot(r,s,'ko','linewidth',2,'markersize',12,'MarkerFaceColor',[.49 1 .63])
axis off

rq1D = JacobiGQ(0,0,N);
rf = [r1D; -r1D;-ones(size(r1D))];
sf = [-r1D;-ones(size(r1D));r1D];

rfq = [rq1D; -rq1D;-ones(size(rq1D))];
sfq = [-rq1D;-ones(size(rq1D));rq1D];

a = 2.5;
plot(a+[-1 1 -1 -1],[-1 -1 1 -1],'k-','linewidth',2)
[rq sq] = Cubature2D(2*N);
plot(a+rq,sq,'ko','linewidth',2,'markersize',11,'MarkerFaceColor',[.49 1 .63])
plot(a+rf,sf,'rs','linewidth',2,'markersize',12,'MarkerFaceColor',[.49 1 .63])
% plot(a+rfq,sfq,'rs','linewidth',2,'markersize',11,'MarkerFaceColor',[.49 1 .63])
axis equal

plot([1.1 a-1.1],[r1D r1D],'k--','linewidth',2)
print(gcf,'-dpng','../hybrid2D.png')

%%
N = 3;

r1D = JacobiGL(0,0,N);
[r s] = meshgrid(r1D);

plot([-1 1 1 -1 -1],[-1 -1 1 1 -1],'k-','linewidth',2)
hold on
plot(r,s,'ko','linewidth',2,'markersize',12,'MarkerFaceColor',[.49 1 .63])
axis off

r1D = JacobiGQ(0,0,N);

rq1D = JacobiGQ(0,0,N);
rf = [r1D; -r1D;-ones(size(r1D))];
sf = [-r1D;-ones(size(r1D));r1D];

rfq = [rq1D; -rq1D;-ones(size(rq1D))];
sfq = [-rq1D;-ones(size(rq1D));rq1D];

a = 2.5;
plot(a+[-1 1 -1 -1],[-1 -1 1 -1],'k-','linewidth',2)
[rq sq] = Cubature2D(2*N);
plot(a+rq,sq,'ko','linewidth',2,'markersize',11,'MarkerFaceColor',[.49 1 .63])
plot(a+rf,sf,'rs','linewidth',2,'markersize',12,'MarkerFaceColor',[.49 1 .63])
% plot(a+rfq,sfq,'rs','linewidth',2,'markersize',11,'MarkerFaceColor',[.49 1 .63])
axis equal

plot([1.1 a-1.1],[r1D r1D],'k--','linewidth',2)
plot([1 ],[r1D ],'rs','linewidth',2,'markersize',12,'MarkerFaceColor',[.49 1 .63])
print(gcf,'-dpng','~/Desktop/proposals/doe/figs/hybrid2D_GQ.png')

%% 

% pyr_mesh;
[VZ VX VY] = pyr_nodes(1); EToV = 1:length(VX); K = 1;
VX(5) = 0;
VY(5) = 1;
VZ(5) = 0;
N = 3;
e = 1;
v = EToV(e,:);

[rtri stri] = Cubature2D(2*N);
e = ones(size(rtri));
[rquad squad ] = meshgrid(JacobiGL(0,0,N));
rquad = rquad(:); squad = squad(:); tquad = -ones(size(squad));
rf = [-e; rtri; rtri;  -rtri; rquad];
sf = [rtri;-e; -stri;  stri;  squad];
tf = [stri; stri; stri; rtri; tquad];

[r1 s1 t1] = pyr_nodes(1);
V1 = pyr_basis(1,rf,sf,tf)/pyr_basis(1,r1,s1,t1);
% rf = -rtri;
% sf = stri;
% tf = rtri;

xf = V1*VX(:);
yf = V1*VY(:);
zf = V1*VZ(:);

ids = [1 2 4 3 1 5 2 5 4 5 3];
plot3(VX(ids),VY(ids),VZ(ids),'k-','linewidth',2)
hold on

idbot = find(abs(tf+1)<1e-8);
idnbot = setdiff(1:length(xf),idbot);
% plot3(xf(idnbot),yf(idnbot),zf(idnbot),'ko','linewidth',2,'markersize',8,'MarkerFaceColor',[.49 1 .63])
plot3(xf(idbot),yf(idbot)*1.0,zf(idbot),'rs','linewidth',2,'markersize',14,'MarkerFaceColor',[.49 1 .63])

r1D = JacobiGL(0,0,N);
[r s t] = meshgrid(r1D);
axis equal
axis off
a = 0;
a2 = 3;
plot3([-1 1 1 -1 -1]-a,[-1 -1 1 1 -1]-a2,[-1 -1 -1 -1 -1],'k-','linewidth',2)
plot3([-1 -1 -1 -1 -1]-a,[-1 1 1 -1 -1]-a2,[-1 -1 1 1 -1],'k-','linewidth',2)
plot3(-[-1 -1 -1 -1 -1]-a,[-1 1 1 -1 -1]-a2,[-1 -1 1 1 -1],'k-','linewidth',2)
plot3([-1 1 1 -1 -1]-a,[-1 -1 -1 -1 -1]-a2,[-1 -1 1 1 -1],'k-','linewidth',2)
plot3([-1 1 1 -1 -1]-a,-[-1 -1 -1 -1 -1]-a2,[-1 -1 1 1 -1],'k-','linewidth',2)

ids = (abs(abs(r)-1)<1e-8 | abs(abs(t)-1)<1e-8) & abs(abs(s-1))>1e-8;

rf = r(ids);
sf = s(ids); 
tf = t(ids);
% plot3(rf,sf-a2,tf,'ko','linewidth',2,'markersize',8,'MarkerFaceColor',[.49 1 .63])

ids = abs(s-1)<1e-8;
rf = r(ids);
sf = s(ids); 
tf = t(ids);
for i = 1:(N+1)^2
    plot3(rf(i),sf(i)-a2,tf(i),'ko','linewidth',2,'markersize',12,'MarkerFaceColor',[.49 1 .63])
   plot3([rf(i) rf(i)],[sf(i)-a2+.99 sf(i)-a2],[tf(i) tf(i)],'k--','linewidth',2)
end
view(-65,30)
% print(gcf,'-dpng','~/Desktop/IGA-DG/docs/multipatch/revision/figs/spectraAdvec.png')
axis tight
print(gcf,'-dpng','../hybrid3D.png')


%% draw GLL nodes

r1D = JacobiGL(0,0,N);


