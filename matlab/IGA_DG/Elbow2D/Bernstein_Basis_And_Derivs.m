%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: Bernstein_Basis_And_Derivs
%
% Input:  xi_1 = first component of evaluation point
%         xi_2 = second component of evaluation point
%         p_1 = polynomial degree in direction 1
%         p_2 = polynomial degree in direction 2
%
% Output: B = values of 2-D Bernstein polynomials at 
%             (xi_1,xi_2)
%         dBdxi = derivatives of 1-D Bernsteins polynomial
%                 at (xi_1,xi_2)
%
% Purpose: Evaluate the 2-D Bernstein polynomials and their
%          derivatives at a point.
%
% Notes: N/A
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [B,dBdxi] = Bernstein_Basis_And_Derivs(xi_1,xi_2,p_1,p_2)

[B1,dBdxi1] = Bernstein_Basis_And_Derivs_1D(xi_1,p_1);
[B2,dBdxi2] = Bernstein_Basis_And_Derivs_1D(xi_2,p_2);

for a1 = 1:p_1+1
    for a2 = 1:p_2+1
        a = (a1-1)*(p_2+1)+a2;
        
        B(a) = B1(a1)*B2(a2);
        dBdxi(a,1) = dBdxi1(a1)*B2(a2);
        dBdxi(a,2) = B1(a1)*dBdxi2(a2);
    end
end