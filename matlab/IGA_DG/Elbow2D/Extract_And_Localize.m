%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: Extract_And_Localize
%
% Input:  n_1 = number of functions in direction 1
%         n_2 = number of functions in direction 2
%         p_1 = polynomial degree in direction 1
%         p_2 = polynomial degree in direction 2
%         Xi_1 = knot vector in direction 1
%         Xi_2 = knot vector in direction 2
%         P = array of NURBS control points (single-indexed)
%         w = array of NURBS weights (single-indexed)
%
% Output: n_el = number of elements
%         C_operators = array of extraction operators
%         IEN = IEN array
%         P_b = array of Bezier control points
%         w_b = array of Bezier weights
%         P_e = array of local NURBS control points
%         w_e = array of local NURBS weights
%
% Purpose: Extract and localize the 2D B-spline basis
%
% Notes: Algorithm follows directly from notes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [n_el,C_operators,IEN,P_b,w_b,P_e,w_e] = Extract_And_Localize(n_1,n_2,p_1,p_2,Xi_1,Xi_2,P,w)

d = size(P,2);

[n_el,C_operators,IEN] = Extract_Basis(p_1,p_2,n_1,n_2,Xi_1,Xi_2);
[P_b,w_b] = Extract_Geometry(2,n_el,C_operators,IEN,P,w);

for e = 1:n_el
    for A = 1:d
        P_e(:,A,e) = P(IEN(:,e),A);
    end
    w_e(:,e) = w(IEN(:,e));
end
