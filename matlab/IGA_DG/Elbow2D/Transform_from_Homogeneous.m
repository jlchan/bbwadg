%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: Transform_from_Homogeneous
%
% Input:  d = number of spatial dimensions
%         P_w = array of homogeneous NURBS control points
%         option = 1 if P_w is single-indexed
%                  2 if P_w is multi-indexed
%
% Output: P = array of NURBS control points
%         w = array of NURBS weights
%
% Purpose: Map control points from homogeneous
%
% Notes: N/A
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P,w] = Transform_from_Homogeneous(d,P_w,option)

if option == 1
    w(:) = P_w(:,d+1);
    for i = 1:d
        P(:,i) = P_w(:,i)./w(:);
    end
elseif option == 2
    w(:,:) = P_w(:,:,d+1);
    for i = 1:d
        P(:,:,i) = P_w(:,:,i)./w(:,:);
    end
end