%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: Transform_to_Homogeneous
%
% Input:  d = number of spatial dimensions
%         P = array of NURBS control points
%         w = array of NURBS weights
%         option = 1 if P_w is single-indexed
%                  2 if P_w is multi-indexed
%
% Output: P_w = array of homogeneous NURBS control points
%
% Purpose: Map control points to homogeneous
%
% Notes: N/A
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P_w] = Transform_to_Homogeneous(d,P,w,option)

if option == 1
    for i = 1:d
        P_w(:,i) = P(:,i).*w(:);
    end

    P_w(:,d+1) = w(:);
elseif option == 2
    for i = 1:d
        P_w(:,:,i) = P(:,:,i).*w(:,:);
    end

    P_w(:,:,d+1) = w(:,:);
end