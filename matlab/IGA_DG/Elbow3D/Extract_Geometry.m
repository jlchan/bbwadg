%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: Extract_Geometry
%
% Input:  d = number of spatial dimensions
%         n_el = number of elements
%         C_operators = array of extraction operators
%         P = array of NURBS control points (single-indexed)
%         w = array of NURBS weights (single-indexed)
%
% Output: P_b = array of Bezier control points
%         w_b = array of Bezier weights
%
% Purpose: Convert the geometry to Bernstein-Bezier
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P_b,w_b] = Extract_Geometry(d,n_el,C_operators,IEN,P,w)

n = size(w);

%%%
% Transform NURBS control points to projective B-spline control points

for i = 1:n
    P_w(i,:) = [w(i)*P(i,:),w(i)];
end

%%%
% Extract on each element and transform back to physical space

for e = 1:n_el
    for A = 1:(d+1)
        P_wb(:,A,e) = transpose(C_operators(:,:,e))*P_w(IEN(:,e),A);
    end
    
    w_b(:,e) = P_wb(:,d+1,e);
    for A = 1:d
        P_b(:,A,e) = P_wb(:,A,e)./w_b(:,e);
    end
end