function [VX VY VZ EToV EToE EToF] = makeHexPyrMesh(K1D)
if nargin==0
    K1D = 2;
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

%% turn each hex into 6 pyramids

% augment at end with middle vertices
VX = hexVX; VY = hexVY; VZ = hexVZ;
for e = 1:hexK
    VX = [VX; mean(hexVX(hexEToV(e,:)))];
    VY = [VY; mean(hexVY(hexEToV(e,:)))];
    VZ = [VZ; mean(hexVZ(hexEToV(e,:)))];
end

EToV = zeros(hexK*6,5);
for e = 1:hexK
    p = [1 3 4 2 5 7 8 6]; % permute hex inds for pyramid ordering I used before
    pyrInds = [hexEToV(e,p(1:4)) e+length(hexVX) hexEToV(e,p(5:8))];
    pyrEToV = [1 2 3 4 5; 1 4 9 6 5;  1 6 7 2 5;  2 7 8 3 5;  3 8 9 4 5;  7 6 9 8 5];
    EToV(6*(e-1)+(1:6),:) = pyrInds(pyrEToV);
end
K = size(EToV,1);
VX = VX(:); VY = VY(:); VZ = VZ(:);

[EToE,EToF]= tiConnectPyr(EToV);

return

% checking by plotting
for e = 1:K
    clf
    plot3(VX,VY,VZ,'.');hold on
    plot3(VX(EToV(e,:)),VY(EToV(e,:)),VZ(EToV(e,:)),'ro')
    
    fids = [1 2 5 1];
    plot3(VX(EToV(e,fids)),VY(EToV(e,fids)),VZ(EToV(e,fids)),'k-','linewidth',2)
    fids = [2 3 5 2];
    plot3(VX(EToV(e,fids)),VY(EToV(e,fids)),VZ(EToV(e,fids)),'k-','linewidth',2)
    fids = [3 4 5 3];
    plot3(VX(EToV(e,fids)),VY(EToV(e,fids)),VZ(EToV(e,fids)),'k-','linewidth',2)
    fids = [4 1 5 4];
    plot3(VX(EToV(e,fids)),VY(EToV(e,fids)),VZ(EToV(e,fids)),'k-','linewidth',2)    
    
    pause
    view(3)
end

function [EToE,EToF]= tiConnectPyr(EToV)

% function [EToE,EToF]= tiConnect3D(EToV)
% Purpose: tetrahedral face connect algorithm due to Toby Isaac

Nfaces=5;
K = size(EToV,1);
Nnodes = max(EToV(:));

% create list of all faces 1, then 2-5
% fvP{1} = [1 2 5]; fvP{2} = [4 1 5];
% fvP{3} = [2 3 5]; fvP{4} = [3 4 5];
% fvP{5} = [1 4 3 2];

id = [];
for e = 1:K
    tnodes = [sort(EToV(e,[1,2,5]));
        sort(EToV(e,[4,1,5]));
        sort(EToV(e,[2,3,5]));
        sort(EToV(e,[3,4,5]))];
    tfids = tnodes*[Nnodes^2;Nnodes;1];    
    qfid = sort(EToV(e,[1 4 3 2]))*[Nnodes^3;Nnodes^2;Nnodes;1];        
    id = [id;tfids;qfid];
end

% uniquely number each set of faces by their node numbers
EToE = ones(Nfaces,1)*(1:K); EToF = (1:Nfaces)'*ones(1,K);
sorted = sortrows([id, (1:Nfaces*K)' EToE(:) EToF(:)],1);

% find matches in the sorted face list
[indices,~]=find( sorted(1:(end-1),1)==sorted(2:end,1) );

% make links reflexive
matchL = [sorted(indices,:)   ; sorted(indices+1,:)];
matchR = [sorted(indices+1,:) ; sorted(indices,:)];

% insert matches
EToE(matchL(:,2)) = matchR(:,3); EToF(matchL(:,2)) = matchR(:,4);
EToE = EToE'; EToF = EToF'; % to get EToE = K by Nfaces



return;
