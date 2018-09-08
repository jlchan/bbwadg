function PlotSol2D(xin, yin, uin)

% function [TRI,xout,yout,uout,interp] = PlotField2D(Nplot, xin, yin, uin)
% Purpose: filled contour plot of solution data

% Globals2D;
K = size(uin,2);
for Nplot = 1:50
    Npout = (Nplot+1)*(Nplot+2)/2;
    if Npout == size(xin,1)
        break
    end    
end

% build equally spaced grid on reference triangle
rout = zeros(Npout,1); sout = zeros(Npout,1);
sk = 1;
for n=1:Nplot+1
    for m=1:Nplot+2-n
        rout(sk) = -1 + 2*(m-1)/Nplot;
        sout(sk) = -1 + 2*(n-1)/Nplot;
        counter(n,m) = sk; sk = sk+1;
    end
end

% build triangulation of equally spaced nodes on reference triangle
tri = [];
for n=1:Nplot+1
    for m=1:Nplot+1-n,
        v1 = counter(n,m);   v2 = counter(n,m+1);
        v3 = counter(n+1,m); v4 = counter(n+1,m+1);
        if(v4)
            tri = [tri;[[v1 v2 v3];[v2 v4 v3]]];
        else
            tri = [tri;[[v1 v2 v3]]];
        end
    end
end

% build triangulation for all equally spaced nodes on all elements
% disp('building triangulation')
TRI = zeros(K*size(tri,1),size(tri,2));
ids = 1:size(tri,1);
for k=1:K    
    TRI(ids + (k-1)*size(tri,1),:) = tri+(k-1)*Npout;
end

% interpolate node coordinates and field to equally spaced nodes

trisurf(TRI, xin(:), yin(:), uin(:));
shading interp

% trimesh(TRI, xout(:), yout(:), uout(:));
% keyboard
% clf;
% hold on
% R = getCGRestriction();
% tricontour([R'*(diag(1./sum(R,2))*R*xout(:)) R'*(diag(1./sum(R,2))*R*yout(:))],TRI,R'*(diag(1./sum(R,2))*R*uout(:)),5)

axis off
axis equal
% material shiny
% lighting gouraud
% camlight headlight
return
