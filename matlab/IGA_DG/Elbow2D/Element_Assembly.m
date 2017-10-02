%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: Element_Assembly
%
% Input:  e = element number
%         k_e = element stiffness matrix
%         f_e = element forcing vector
%         IEN = IEN array
%         ID = ID array
%         BC = Dirichlet BC array
%         g = vector of Dirichlet control variables
%         K = global stiffness matrix
%         F = global forcing vector
%
% Output: K = global stiffness matrix (updated)
%         F = global forcing vector (updated)
%
% Purpose: Assemble the element stiffness and forcing
%
% Notes: Algorithm follows directly from notes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [K,F] = Element_Assembly(e,k_e,f_e,IEN,ID,BC,g,K,F)

n_loc = size(IEN,1);
d = 2;

for A = 1:d
    for a = 1:n_loc
        p = d*(a-1)+A;
        i = IEN(a,e);
        P = ID(A,i);
        
        if BC(A,i) == 0
            for B = 1:d
                for b = 1:n_loc
                    q = d*(b-1)+B;
                    j = IEN(b,e);
                    Q = ID(B,j);
                    
                    if BC(B,j) == 0
                        K(P,Q) = K(P,Q) + k_e(p,q);
                    else
                        F(P) = F(P) - k_e(p,q)*g(B,j);
                    end
                    
                end                                
            end
            
            F(P) = F(P) + f_e(p);
        end
    end
end