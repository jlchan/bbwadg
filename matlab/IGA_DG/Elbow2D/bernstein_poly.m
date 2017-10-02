%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: bernstein_poly
%
% Input:  p = polynomial degree
%         x = evaluation point
%
% Output: bern = value of Bernstein polynomials at x
%
% Purpose: Evaluate the Bernstein polynomials at a point.
%
% Notes: N/A
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bern = bernstein_poly(p,x)

if (p == 0)    
    bern(1) = 1.0;
    
elseif (p > 0)
    
    bern(1) = 1.0 - x;
    bern(2) = x;
    
    for i = 2 : p
        bern(i+1) = x * bern(i);
        for j = i-1 : -1 : 1
            bern(j+1) = x * bern(j) + ( 1.0 - x ) * bern(j+1);
        end
        bern(1) = ( 1.0 - x ) * bern(1);
    end
    
end