function check_TV_2D

Globals2D

N = 4;

[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Grid/Other/block2.neu');

% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Grid/Other/lshape.neu');
% [Nv VX VY K EToV] = unif_tri_mesh(4);
StartUp2D

Hrefine2D(1:K);
StartUp2D

[r s] = Nodes2D(N); [r s] = xytors(r,s);
[re se] = EquiNodes2D(N); [re se] = xytors(re,se);
[rp sp] = EquiNodes2D(N+4); [rp sp] = xytors(rp,sp);

Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;
Ve = Vandermonde2D(N,re,se)/V;
xe = Ve*x; ye = Ve*y;

[V Vr Vs] = bern_basis_tri(N,r,s);
Dr = V\Vr; Ds = V\Vs;
Ve = bern_basis_tri(N,re,se);

f = @(x,y) (y > .5*sin(pi*x)) + exp(sin(pi*x+y)) + sin(2*pi*(x+y));
% f = @(x,y) y > .1;
% f = @(x,y) sin(25*pi*x+exp(y));
% r = @(x,y) sqrt(x.^2+y.^2);
% f = @(x,y) 1./sqrt(r(x,y));

Vp = bern_basis_tri(N,rp,sp);

u = V\f(x,y);

% up = Vp*u;
% color_line3(xp,yp,up,up,'.'); 
% hold on
% plot3(xe,ye,u,'o');
% return

%% plot TV estimate

clf
[TVK TVx TVy] = compute_TV(u);
TV = repmat(TVK,Np,1);
% TVx = repmat(TVx,Np,1);
% TVy = repmat(TVy,Np,1);
tvp = Vp*TV;
% tvp = log(tvp);

% quiver(x,y,TVx,TVy)
% hold on
color_line3(xp,yp,tvp,tvp,'.')

return

ux = rx.*(Dr*u) + sx.*(Ds*u);
uy = ry.*(Dr*u) + sy.*(Ds*u);
ux = V*ux; uy = V*uy;
d2u = sqrt(ux.^2 + uy.^2); % norm(grad(u))
% d2ux = rx.*(Dr*ux) + sx.*(Ds*ux);
% d2uy = ry.*(Dr*uy) + sy.*(Ds*uy);
% duxy = ry.*(Dr*ux) + sy.*(Ds*ux);
% duyx = rx.*(Dr*uy) + sx.*(Ds*uy);
% d2u = abs(d2ux) + abs(d2uy) + abs(duxy) + abs(duyx);
d2u = repmat(max(abs(d2u),[],1),Np,1);
vp = Vp*d2u;

color_line3(xp,yp,vp,vp,'.')
view(2)
colorbar

return
%%

clf
Klim = find(TVK > Np); % ad-hoc lower bound

% u(:,Klim) = Ve*f(xe(:,Klim),ye(:,Klim)); % true BB approx
% up = Vp*u;
% plot3(xp,yp,up,'.')

% a = 2/(N+1);
a  = .25;
u(:,Klim) = (a*eye(Np) + (1-a)*Ve)*Ve*u(:,Klim); % based on discrete representation
up = Vp*u;
hold on
color_line3(xp,yp,up,up,'.')


function [TV TVx TVy] = compute_TV(u)

Globals2D

off = 0;
TVx = 0;
TVy = 0;
for i = 0:N
    Ni = N-i;
    Npi = Ni+1;
    idsx = (1:Npi) + off;
    TVx = TVx + TV1D(N-i,u(idsx,:));

    offj = 0;
    idsy = [];
    for j = 0:N-i
        idsy(j+1) = i + offj + 1;
        offj = offj + (N-j+1);
    end    
    TVy = TVy + TV1D(N-i,u(idsy,:));
    off = off + Npi;
end
TV = TVx+TVy;


function TV = TV1D(N,u)

TV = 0;
for i = 1:N
    TV = TV + abs(u(i,:) - u(i+1,:));
end
