%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: Extraction3D
%
% Input:  n_1 = number of functions in direction 1
%         n_2 = number of functions in direction 2
%         n_3 = number of functions in direction 3
%         p_1 = polynomial degree in direction 1
%         p_2 = polynomial degree in direction 2
%         p_3 = polynomial degree in direction 3
%         Xi_1 = knot vector in direction 1
%         Xi_2 = knot vector in direction 2
%         Xi_3 = knot vector in direction 3
%
% Output: n_el = number of 3D elements
%         C_e = array of 3D extraction operators
%
% Purpose: Extract the 3D B-spline basis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [n_el, C_e] = Extraction3D(n_1,n_2,n_3,p_1,p_2,p_3,Xi_1,Xi_2,Xi_3)

[n_el_1, C_e_1] = Extraction1D(n_1,p_1,Xi_1);
[n_el_2, C_e_2] = Extraction1D(n_2,p_2,Xi_2);
[n_el_3, C_e_3] = Extraction1D(n_3,p_3,Xi_3);

n_el = n_el_1*n_el_2*n_el_3;

for e1 = 1:1:n_el_1
    for e2 = 1:1:n_el_2
        for e3 = 1:1:n_el_3
            e = (e1-1)*n_el_2*n_el_3 + (e2-1)*n_el_3 + e3;
            C_e(:,:,e) = kron(kron(C_e_1(:,:,e1),C_e_2(:,:,e2)),C_e_3(:,:,e3));
        end
    end
end