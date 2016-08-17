function [Nv VX VY K EToV] = unif_tri_mesh(K)
% K = 2
[VY VX] = meshgrid(linspace(-1,1,K+1));

% plot(VX,VY,'o')
% text(VX(:)+.1,VY(:),num2str((1:length(VX(:)))'))

sk = 1;
for ey = 1:K
    for ex = 1:K
        id = @(ex,ey) ex + (ey-1)*(K+1);
        id1 = id(ex,ey);
        id2 = id(ex+1,ey);
        id3 = id(ex+1,ey+1);
        id4 = id(ex,ey+1);
        VXe = VX([id1 id2 id3 id4]);
        VYe = VY([id1 id2 id3 id4]);
        EToV(2*sk-1,:) = [id1 id2 id3];
        EToV(2*sk,:) = [id3 id4 id1];
        sk = sk + 1;
%         plot(VX,VY,'o')
%         hold on
%         plot(VXe,VYe,'*')
%         pause
%         clf
    end
end
VX = VX(:)'; VY = VY(:)';
Nv = length(VX(:));
K = size(EToV,1);