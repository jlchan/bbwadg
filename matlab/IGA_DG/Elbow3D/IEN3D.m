%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: IEN3D
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
% Output: IEN3D = 3D IEN array
%
% Purpose: Compute the 2D IEN array
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [IEN] = IEN3D(n_1,n_2,n_3,p_1,p_2,p_3,Xi_1,Xi_2,Xi_3)

IEN1 = IEN1D(n_1,p_1,Xi_1);
IEN2 = IEN1D(n_2,p_2,Xi_2);
IEN3 = IEN1D(n_3,p_2,Xi_3);

n_el_1 = size(IEN1,2);
n_el_2 = size(IEN2,2);
n_el_3 = size(IEN3,2);

for e1 = 1:1:n_el_1
    for a1 = 1:1:(p_1+1)
        i1 = IEN1(a1,e1);
        
        for e2 = 1:1:n_el_2
            for a2 = 1:1:(p_2+1)
                i2 = IEN1(a2,e2);
                
                for e3 = 1:1:n_el_3
                    for a3 = 1:1:(p_3+1)
                        i3 = IEN1(a3,e3);
                
                        e = (e1-1)*n_el_2*n_el_3 + (e2-1)*n_el_3 + e3;
                        a = (a1-1)*(p_2+1)*(p_3+1) + (a2-1)*(p_3+1) + a3;
                        i = (i1-1)*n_2*n_3 + (i2-1)*n_3 + i3;
                        
                        IEN(a,e) = i;
                    end
                end
                
            end
        end
        
    end
end