function [EToE,EToF]= tiConnect3DHex(EToV)

% function [EToE,EToF]= tiConnect3D(EToV)
% Purpose: tetrahedral face connect algorithm due to Toby Isaac

Nfaces=6;
K = size(EToV,1);
Nnodes = max(max(EToV));

% create list of all faces

% f5 = [1 5 7 3];
% f6 = [2 6 8 4];
% f1 = [1 2 4 3];
% f2 = [5 6 8 7];
% f3 = [1 2 6 5];
% f4 = [3 4 8 7];

f1 = [1     2     5     6];
f2 = [3     4     7     8];
f3 = [1     3     5     7];
f4 = [2     4     6     8];
f5 = [1     2     3     4];
f6 = [5     6     7     8];

% r = -/+1
% f1 = [4,1,8,5];
% f2 = [2,3,6,7];
% 
% % s = -/+1
% f3 = [1,2,5,6];
% f4 = [3,4,7,8];
% 
% % t = -/+1
% f5 = [1,2,4,3];
% f6 = [5,6,8,7];

fnodes = [EToV(:,f1);EToV(:,f2);EToV(:,f3);...
    EToV(:,f4);EToV(:,f5);EToV(:,f6)];
fnodes = sort(fnodes,2)-1;

% set up default element to element and Element to faces connectivity
EToE= (1:K)'*ones(1,Nfaces); EToF= ones(K,1)*(1:Nfaces);

% uniquely number each set of three faces by their node numbers
id = fnodes(:,1)*Nnodes*Nnodes + fnodes(:,2)*Nnodes+fnodes(:,3)+1;
spNodeToNode=[id, (1:Nfaces*K)', EToE(:), EToF(:)];

% Now we sort by global face number.
sorted=sortrows(spNodeToNode,1);

% find matches in the sorted face list
[indices,dummy]=find( sorted(1:(end-1),1)==sorted(2:end,1) );

% make links reflexive
matchL = [sorted(indices,:)   ;sorted(indices+1,:)];
matchR = [sorted(indices+1,:) ;sorted(indices,:)];

% insert matches
EToE(matchL(:,2)) = matchR(:,3); EToF(matchL(:,2)) = matchR(:,4);
return;
