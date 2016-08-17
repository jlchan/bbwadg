function [nx,ny,nz,sJ] = tet_sgeofacs(VX,VY,VZ,Dr,Ds,Dt)

[rx,sx,tx,ry,sy,ty,rz,sz,tz,J] = tet_geom_factors(VX,VY,VZ,Dr,Ds,Dt);

Nfq = size(Dr,1)/4; % 4 faces
fids1 = (1:Nfq)';
fids2 = fids1 + Nfq;
fids3 = fids2 + Nfq;
fids4 = fids3 + Nfq;

nx = [-tx(fids1,:); -sx(fids2,:); rx(fids3,:) + sx(fids3,:) + tx(fids3,:); -rx(fids4,:)];
ny = [-ty(fids1,:); -sy(fids2,:); ry(fids3,:) + sy(fids3,:) + ty(fids3,:); -ry(fids4,:)];
nz = [-tz(fids1,:); -sz(fids2,:); rz(fids3,:) + sz(fids3,:) + tz(fids3,:); -rz(fids4,:)];

sJ = sqrt(nx.^2 + ny.^2 + nz.^2);

nx = nx./sJ; ny = ny./sJ;nz = nz./sJ;
sJ = sJ.*J;