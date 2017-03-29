Globals2D;

N = 4;

Nq = 2*N+1;


% for mesh
[VX VY] = EquiNodes2D(1); VX = VX(:)'; VY = VY(:)'; K = 1; EToV = 1:3;
StartUp2D;

[re se] = EquiNodes2D(N); [re se] = xytors(re,se);
Ve = Vandermonde2D(N,re,se)/V;
[rp sp] = EquiNodes2D(200); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
VB = bern_basis_tri(N,r,s);

N2 = N+1;
[r2 s2] = Nodes2D(N2); [r2 s2] = xytors(r2,s2);
[re2 se2] = EquiNodes2D(N2); [re2 se2] = xytors(re2,se2);
E = bern_basis_tri(N2,r2,s2)\bern_basis_tri(N,r2,s2);
V2 = Vandermonde2D(N,r2,s2)/Vandermonde2D(N,r,s);
xe = Ve*x; ye = Ve*y;
xe2 = E*xe; ye2 = E*ye;

xT = V2*x; yT = V2*y; % original xy coords
xeT = xe2; yeT = ye2; % original xy coords

tri = delaunayFixArea(x',y');
tri2 = delaunayFixArea(xT',yT');
triB = delaunayFixArea(xe',ye');
triB2 = delaunayFixArea(xe2',ye2');

nid = 1:Np;
x0 = x(nid);
y0 = y(nid);
for a = .05%.1:.01:.33
    % a = .1;
%     x(nid) = x0-a;
%     y(nid) = y0+a;
    x = x + a*randn(Np,1);
    y = y + a*randn(Np,1);
    xe = VB\x; ye = VB\y;
    xe2 = E*xe; ye2 = E*ye;
        
    x2 = V2*x; y2 = V2*y;
        
    xp = Vp*x; yp = Vp*y;
    [rx,sx,ry,sy,J] = GeometricFactors2D(x,y,Dr,Ds);
    Jp = Vp*J;
    minJp = min(Jp)    
    min_id = [];
    
    figure
    clf
    PlotField2D(50,x,y,J-1.05*max(J(:)));%(xp,yp,'.','markersize',24)
    % caxis([-.1 .25])
    hold on
    plot(xp(min_id),yp(min_id),'ro','markersize',24)
    triplot(triB2,xe2,ye2,'k','linewidth',2)
    plot3(xe2,ye2,.01*ones(size(xe2)),'ko','linewidth',2,'markersize',12,'MarkerFaceColor',[.49 1 .63])
    view(2)

    figure
    clf
    PlotField2D(50,x,y,J-1.05*max(J(:)));%(xp,yp,'.','markersize',24)
    % caxis([-.1 .25])
    hold on
    plot(xp(min_id),yp(min_id),'ro','markersize',24)
    triplot(triB,xe,ye,'k','linewidth',2)
    plot3(xe,ye,.01*ones(size(xe)),'ko','linewidth',2,'markersize',12,'MarkerFaceColor',[.49 1 .63])
    view(2)
end


%  print(gcf,'-dpng','bb_submesh1.png')
% print(gcf,'-dpng','bb_submesh2.png')

%% map BB points instead

clear

Globals2D;

N = 4;

Nq = 2*N+1;

% for mesh
[VX VY] = EquiNodes2D(1); VX = VX(:)'; VY = VY(:)'; K = 1; EToV = 1:3;
StartUp2D;

[re se] = EquiNodes2D(N); [re se] = xytors(re,se);
Ve = Vandermonde2D(N,re,se)/V;
[rp sp] = EquiNodes2D(200); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
VB = bern_basis_tri(N,r,s);

N2 = 2*N;
[r2 s2] = Nodes2D(N2); [r2 s2] = xytors(r2,s2);
[re2 se2] = EquiNodes2D(N2); [re2 se2] = xytors(re2,se2);
E = bern_basis_tri(N2,r2,s2)\bern_basis_tri(N,r2,s2);
E(abs(E)<1e-8) = 0;
V2 = Vandermonde2D(N,r2,s2)/Vandermonde2D(N,r,s);
xe = Ve*x; ye = Ve*y;
xe2 = E*xe; ye2 = E*ye;

xT = V2*x; yT = V2*y; % original xy coords
xeT = xe2; yeT = ye2; % original xy coords

tri = delaunayFixArea(x',y');
tri2 = delaunayFixArea(xT',yT');
triB = delaunayFixArea(xe2',ye2');

xe = VB\x; ye = VB\y;

nid = 3;
x0 = xe(nid);
y0 = ye(nid);
for a = -1%.1:.01:.33
%     xe(nid) = x0-1.25*a;
    ye(nid) = y0-a;
%     xe(nid) = x0-a;
    
    xe2 = E*xe; ye2 = E*ye;
    x = VB*xe; y = VB*ye;        
    x2 = V2*x; y2 = V2*y;
        
    xp = Vp*x; yp = Vp*y;
    [rx,sx,ry,sy,J] = GeometricFactors2D(x,y,Dr,Ds);
    Jp = Vp*J;
    maxJp = max(Jp)
    minJp = min(Jp)    
    min_id = [];
    
    clf
    PlotField2D(50,x,y,J-1.05*max(J(:)));%(xp,yp,'.','markersize',24)
    % caxis([-.1 .25])
    hold on
    plot(xp(min_id),yp(min_id),'ro','markersize',24)
    triplot(triB,xe2,ye2,'k','linewidth',2)
    plot3(xe2,ye2,.01*ones(size(xe2)),'ko','linewidth',2,'markersize',12,'MarkerFaceColor',[.49 1 .63])
    view(2)
%     axis([-1 1 -1 1])
    drawnow
end
%  print(gcf,'-dpng','bb_submesh1.png')
% print(gcf,'-dpng','bb_submesh2.png')