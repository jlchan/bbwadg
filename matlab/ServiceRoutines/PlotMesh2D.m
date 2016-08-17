function PlotMesh2D()

% function PlotMesh2D()
% Purpose: Show unstructured finite element grid

Globals2D;

axis equal
xmax = max(max(x)); xmin = min(min(x));
ymax = max(max(y)); ymin = min(min(y));

Lx = xmax-xmin;
Ly = ymax-ymin;
xmax = xmax+.1*Lx; xmin = xmin-.1*Lx;
ymax = ymax+.1*Ly; ymin = ymin-.1*Ly;

axis([xmin xmax ymin ymax])
% drawnow; pause(.05);

oFx = reshape(Fx, Nfp, Nfaces*K); oFy = reshape(Fy, Nfp, Nfaces*K);

plot(oFx, oFy, 'k-','linewidth',2)
axis equal
axis([xmin xmax ymin ymax])

% for e = 1:K
%     xm = mean(VX(EToV(e,:)));
%     ym = mean(VY(EToV(e,:)));
%     text(xm,ym,num2str(e));
% end
% drawnow; pause(.05);
axis off
return;
