clear

Globals2D

N = 2;
[VX VY] = Nodes2D(1); 
% VX = VX+rand(size(VX));
% [VX VY] = xytors(VX,VY); 
EToV = 1:length(VX); K = 1;

StartUp2D

u = reshape(rand(Np*K,1),Np,K);

dudx = rx.*(Dr*u) + sx.*(Ds*u);
dudy = ry.*(Dr*u) + sy.*(Ds*u);

xr = Dr*x; xs = Ds*x; yr = Dr*y; ys = Ds*y; 
J = -xs.*yr + xr.*ys;

rx = rx(1); sx = sx(1); ry = ry(1); sy = sy(1);
xr = xr(1); xs = xs(1); yr = yr(1); ys = ys(1);
G2 = [xr xs; yr ys]';    
G = [rx ry; sx sy]';  

nhat = [0 -1; 1/sqrt(2) 1/sqrt(2); -1 0]';
nx = reshape(nx,N+1,Nfaces); ny = reshape(ny,N+1,Nfaces);

invGTn = G*nhat;
for i =1:3
    fJ(i) = norm(invGTn(:,i))
end
n = [nx(1,:); ny(1,:)];
norm(G*nhat*diag(1./fJ) - n,'fro')

n = G*(nhat*diag(1./fJ));
dudx = G(1,1)*Dr*u + G(1,2)*Ds*u;
dudy = G(2,1)*Dr*u + G(2,2)*Ds*u;
norm(dudx - (rx*Dr*u + sx*Ds*u))
norm(dudy - (ry*Dr*u + sy*Ds*u))

G2 = G'*G
f = 3;
dudn = n(1,f)*dudx + n(2,f)*dudy
dudnhat = (nhat(1,f)*(G2(1,1)*Dr*u + G2(1,2)*Ds*u) + nhat(2,f)*(G2(2,1)*Dr*u + G2(2,2)*Ds*u))/fJ(f)



