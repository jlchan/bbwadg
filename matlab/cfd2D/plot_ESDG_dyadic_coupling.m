useQuads = 0;mypath

N = 2;
[rq sq] = Cubature2D(2*N);
[rq1D wq1D] = JacobiGQ(0,0,N);
e = ones(size(rq1D));
rfq = [rq1D; -rq1D; -e];
sfq = [-e; rq1D; rq1D];

plot([-1 1 -1 -1],[-1 -1 1 -1],'k-','linewidth',2)
hold on
plot(rq,sq,'bo','linewidth',2,'markersize',15,'MarkerFaceColor',[.49 1 .63]);
plot(rfq,sfq,'rs','linewidth',2,'markersize',15,'MarkerFaceColor',[.49 1 .63]);

axis off

%%

Globals2D
N = 1;

VX = [-1 0 0 1]; VY = [0 -1/2 1/2 0];
EToV = [1 2 3; 2 3 4]; K = 2;
% [Nv, VX, VY, K, EToV] = unif_tri_mesh(1);
StartUp2D

[rq sq] = Cubature2D(2*N);
[rq1D wq1D] = JacobiGQ(0,0,N);
e = ones(size(rq1D));
rfq = [rq1D; -rq1D; -e];
sfq = [-e; rq1D; rq1D];

Vq = Vandermonde2D(N,rq,sq)/V;
Vfq = Vandermonde2D(N,rfq,sfq)/V;

vids = EToV(1,[1 2 3 1]);
hold on;plot(VX(vids),VY(vids),'k-','linewidth',2)
vids = EToV(2,[1 2 3 1]);
hold on;plot(.25+VX(vids),VY(vids),'k-','linewidth',2)
plot(Vq*x(:,1),Vq*y(:,1),'bo','linewidth',2,'markersize',16,'MarkerFaceColor',[.49 1 .63])
plot(Vq*x(:,2)+.25,Vq*y(:,2),'bo','linewidth',2,'markersize',16,'MarkerFaceColor',[.49 1 .63])

plot(Vfq*x(:,1),Vfq*y(:,1),'rs','linewidth',2,'markersize',16,'MarkerFaceColor',[.49 1 .63])
plot(Vfq*x(:,2)+.25,Vfq*y(:,2),'rs','linewidth',2,'markersize',16,'MarkerFaceColor',[.49 1 .63])
axis off
axis equal