function check_TV_2D

Globals2D

N = 5;

[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Grid/Other/block2.neu');
% [Nv, VY, VX, K, EToV] = MeshReaderGambit2D('Grid/CFD/kovA02.neu');
% VX = VX - .5; VY = VY - .25;

% [Nv VX VY K EToV] = unif_tri_mesh(4);
StartUp2D

nref = 1;
for ref = 1:nref
    Hrefine2D(1:K);
    StartUp2D
end
[r s] = Nodes2D(N); [r s] = xytors(r,s);
[re se] = EquiNodes2D(N); [re se] = xytors(re,se);
[rp sp] = EquiNodes2D(25); [rp sp] = xytors(rp,sp);

Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;
Ve = Vandermonde2D(N,re,se)/V;
xe = Ve*x; ye = Ve*y;

[V Vr Vs] = bern_basis_tri(N,r,s);
Dr = V\Vr; Ds = V\Vs;
Ve = bern_basis_tri(N,re,se);

rad = @(x,y) sqrt(x.^2+y.^2);
% f = @(x,y) (y > .5*sin(pi*x)) + exp(sin(pi*(.5*x+y))).*sin(4*pi*x).*sin(4*pi*y);
% f = @(x,y) y > .1;
% f = @(x,y) sin(25*pi*x+exp(y));
f = @(x,y) (rad(x+2/3,y) < 1/4) + 4*abs(rad(x,y)-1/4).*(rad(x,y) < 1/4) + exp(-36*rad(x-2/3,y).^2);

Vp = bern_basis_tri(N,rp,sp);

u = V\f(x,y);

up = Vp*u;
%color_line3(xp,yp,up,up,'.'); 
color_line3(xp,yp,f(xp,yp),f(xp,yp),'.'); 
% hold on
% plot3(xe,ye,u,'o');
% return
axis off
axis tight
%print(gcf,'-dpng','u_TVBB.png')
view(-25,25)
print(gcf,'-dpng','u_detector2D.png')

% plot TV estimate
% clf
TVK = compute_TV(u);
TV = repmat(TVK,Np,1);
% TVx = repmat(TVx,Np,1);
% TVy = repmat(TVy,Np,1);
tvp = Vp*TV;
% tvp = log(tvp);

% quiver(x,y,TVx,TVy)
% hold on
figure
% PlotMesh2D
color_line3(xp,yp,tvp,tvp,'.')
colorbar
set(gca,'fontsize',16)
% axis equal
axis off
axis tight

% print(gcf,'-dpng','TVBB.png')
print(gcf,'-dpng','detector2D.png')
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

x = linspace(0,1,N);
% w = 4*x.*(1-x);
w = ones(size(x));
TV = 0;
for i = 1:N
%     TV = TV + abs(u(i,:) - u(i+1,:)).*w(i);
end

tol = 1e-2; %max(abs(u - repmat(mean(u,1),N+1,1)))
for i = 2:N
    flag1 = abs(u(i,:)-u(i-1,:)) > tol;
    flag2 = abs(u(i+1,:)-u(i,:)) > tol;
    TV = TV + abs(sign(u(i,:) - u(i-1,:)).*flag1 - sign(u(i+1,:)- u(i,:)).*flag2); % sign variations at each node    
end
