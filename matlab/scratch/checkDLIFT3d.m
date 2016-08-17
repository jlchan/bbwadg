clear

Globals3D

N = 1;
[VX VY VZ] = Nodes3D(1); 
% VX = VX+rand(size(VX));
% [VX VY VZ] = xyztorst(VX,VY,VZ); 
EToV = 1:length(VX); K = 1;

StartUp3D

rx = rx(1); sx = sx(1); tx = tx(1);
ry = ry(1); sy = sy(1); ty = ty(1);
rz = rz(1); sz = sz(1); tz = tz(1);
G = [rx sx tx; ry sy ty; rz sz tz];
GtG = G'*G;

nx = reshape(nx,Nfp,Nfaces); nx = nx(1,:);
ny = reshape(ny,Nfp,Nfaces); ny = ny(1,:);
nz = reshape(nz,Nfp,Nfaces); nz = nz(1,:);

n = [nx;ny;nz]
nhat = [0 0 -1; 0 -1 0; 1/sqrt(3) 1/sqrt(3) 1/sqrt(3); -1 0 0]';
Gn = G*nhat;
for i = 1:size(Gn,2)
    fJ(i) = norm(Gn(:,i)); 
end
%norm(G*nhat*diag(1./fJ)-n,'fro')


% xr = xr(1); xs = xs(1); yr = yr(1); ys = ys(1);
% G2 = [xr xs; yr ys]';    
% G = [rx ry; sx sy]';  

