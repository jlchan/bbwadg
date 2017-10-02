%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: Bernstein_Basis_And_Derivs_1D
%
% Input:  xi = evaluation point
%         p = polynomial degree
%
% Output: B = values of 1-D Bernstein polynomials at xi
%         dBdx = derivatives of 1-D Bernsteins polynomial
%                at xi
%
% Purpose: Evaluate the 1-D Bernstein polynomials and their
%          derivatives at a point.
%
% Notes: N/A
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [B,dBdx] = Bernstein_Basis_And_Derivs_1D(xi,p)

basis_p = bernstein_poly(p,xi);
basis_pm1 = bernstein_poly(p-1,xi);

B = basis_p;

dBdx(1) = -p*basis_pm1(1);
for a = 2:p
    dBdx(a) = p*(basis_pm1(a-1)-basis_pm1(a));
end
dBdx(p+1) = p*basis_pm1(p);