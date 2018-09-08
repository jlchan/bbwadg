% builds a simple mesh of NK-by-NK quad elements

function [Nv, VX, VY, K, EToV] = QuadMesh2D(Nx,Ny)

if nargin==1
    Ny = Nx;
end

Nxp = Nx+1;
Nyp = Ny+1;
Nv = Nxp*Nyp;
K = Nx*Ny;

x1D = linspace(-1,1,Nxp);
y1D = linspace(-1,1,Nyp);
[x, y] = meshgrid(x1D,y1D);
[I, J] = meshgrid(1:Nxp,1:Nyp);
inds = (I-1)*Ny + (J+I-1);
EToV = zeros(K,4);
k = 1;
for i = 1:Ny
    for j = 1:Nx
%         EToV(k,:) = [inds(i,j) inds(i+1,j) inds(i+1,j+1) inds(i,j+1)  ];
        EToV(k,:) = [inds(i,j) inds(i,j+1) inds(i+1,j+1) inds(i+1,j)  ];
        k = k+1;
    end
end
% EToV-1
VX = x(:)';
VY = y(:)';


