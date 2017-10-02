function [VX VY VZ EToV] = makeHexWedgeMesh(K1D)
if nargin==0
    K1D = 1;
end
hexVX1D = linspace(-1,1,K1D+1);
[hexVX hexVY hexVZ] = meshgrid(hexVX1D);
hexVX = hexVX(:); hexVY = hexVY(:); hexVZ = hexVZ(:);

% make hex mesh
hexInds = 1 + [0, 1, (K1D+1), 1+(K1D+1),...
    (K1D+1)^2, 1+(K1D+1)^2, (K1D+1)^2+(K1D+1),1+(K1D+1)^2+(K1D+1)]; 
hexEToV = zeros(K1D^3,8);
e = 1;
for i = 0:K1D-1
    for j = 0:K1D-1
        for k = 0:K1D-1
            hexEToV(e,:) = hexInds + i + j*(K1D+1) + k*(K1D+1)^2;
            e = e+1;
        end
    end    
end
hexK = size(hexEToV,1);

%% turn each hex into 2 wedges

% augment at end with middle vertices
VX = hexVX; VY = hexVY; VZ = hexVZ;

p = [1 3 4 2 5 7 8 6]; % permute hex inds for pyramid ordering I used before
p = 1:8;
EToV = zeros(hexK*2,6);
for e = 1:hexK    
    inds = hexEToV(e,p);       

    wEToV = [1 2 3 5 6 7; 1 4 3 5 8 7];
    EToV(2*(e-1)+(1:2),:) = inds(wEToV);
end    
K = size(EToV,1);
VX = VX(:); VY = VY(:); VZ = VZ(:);

% e = 1;
% inds = EToV(1,:);
% plot3(VX(inds),VY(inds),VZ(inds),'.')
% text(VX(inds),VY(inds),VZ(inds),num2str((1:length(inds))'))

