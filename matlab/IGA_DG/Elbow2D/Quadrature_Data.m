%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: Quadrature_Data
%
% Input:  n_q = number of quadrature points per direction
%         
% Output: xi_q = array of quadrature points
%         w_q = array of quadrature weights
%
% Purpose: Construct the quadrature point/weight arrays
%
% Notes: N/A
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xi_q, w_q] = Quadrature_Data(n_q)

switch n_q
    case 1
        xi_q = [1/2];
        w_q = [1];
    case 2
        xi_q = [1/2-1/2*sqrt(1/3), 1/2+1/2*sqrt(1/3)];
        w_q = [1/2, 1/2];
    case 3
        xi_q = [1/2-1/2*sqrt(3/5), 1/2, 1/2+1/2*sqrt(3/5)];
        w_q = [5/18, 4/9, 5/18];
    case 4
        xi_q = [1/2-1/2*sqrt(3/7-2/7*sqrt(6/5)), 1/2+1/2*sqrt(3/7-2/7*sqrt(6/5)), 1/2-1/2*sqrt(3/7+2/7*sqrt(6/5)), 1/2+1/2*sqrt(3/7+2/7*sqrt(6/5))];
        w_q = [1/2*(18+sqrt(30))/36, 1/2*(18+sqrt(30))/36, 1/2*(18-sqrt(30))/36, 1/2*(18-sqrt(30))/36];    
end
