%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: Boundary_Arrays
%
% Input:  problem = problem description index
%         n_1 = number of functions in direction 1
%         n_2 = number of functions in direction 2
%         Xi_1 = knot vector in direction 1
%         Xi_2 = knot vector in direction 2
%
% Output: BC = Dirichlet BC array
%         g = vector of Dirichlet control variables
%         Neumann = Neumann BC array
%
% Purpose: Construct the BC data and arrays
%
% Notes: N/A
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [BC,g,Neumann] = Boundary_Arrays(problem,n_1,n_2,Xi_1,Xi_2)

n_el_1 = length(unique(Xi_1))-1;
n_el_2 = length(unique(Xi_2))-1;

n_el = n_el_1*n_el_2;

d = 2;

switch problem
    %%%
    % Array Declarations for Problem 3
    
    case 1
        BC = zeros(d,n_1*n_2);
        g = zeros(d,n_1*n_2);
        
        for i2 = 1:n_2
            i_left = i2;
            i_right = (n_1-1)*n_2 + i2;
            
            BC(2,i_left) = 1;
            BC(1,i_right) = 1;
        end        
        
        Neumann = zeros(d,4,n_el);
        
        for e2 = 1:n_el_2
            e_left = (1-1)*n_el_2 + e2;
            e_right = (n_el_1-1)*n_el_2 + e2;
            
            Neumann(1,4,e_left) = 1;
            Neumann(2,2,e_right) = 1;
        end
        
        for e1 = 1:n_el_1
            e_bottom = (e1-1)*n_el_2 + 1;
            e_top = (e1-1)*n_el_2 + n_el_2;
            
            Neumann(1,1,e_bottom) = 1;
            Neumann(2,1,e_bottom) = 1;
            
            Neumann(1,3,e_top) = 1;
            Neumann(2,3,e_top) = 1;
        end

    %%%
    % Array Declarations for Problem 4 Part 1
    
    case 2
        BC = zeros(d,n_1*n_2);
        g = zeros(d,n_1*n_2);
        
        for i2 = 1:n_2
            i_left = i2;
            i_right = (n_1-1)*n_2 + i2;
            
            BC(1,i_left) = 1;
            BC(2,i_left) = 1;
            
            BC(1,i_right) = 1;
        end
        
        half_n_1 = ceil(n_1/2);
        
        for i1 = 1:half_n_1
            i_top = (i1-1)*n_2 + n_2;
            
            BC(1,i_top) = 1;
            BC(2,i_top) = 1;
        end
        
        Neumann = zeros(d,4,n_el);
        
        for e2 = 1:n_el_2
            e_right = (n_el_1-1)*n_el_2 + e2;
                        
            Neumann(2,2,e_right) = 1;
        end
        
        for e1 = 1:n_el_1
            e_bottom = (e1-1)*n_el_2 + 1;            
            
            Neumann(1,1,e_bottom) = 1;
            Neumann(2,1,e_bottom) = 1;            
        end
        
        half_n_el_1 = n_el_1/2;
        
        for e1 = half_n_el_1+1:n_el_1
            e_top = (e1-1)*n_el_2 + n_el_2;
            
            Neumann(1,3,e_top) = 1;
            Neumann(2,3,e_top) = 1;
        end
        
    %%%
    % Array Declarations for Problem 4 Part 2
    
    case 3
        BC = zeros(d,n_1*n_2);
        g = zeros(d,n_1*n_2);
        
        for i2 = 1:n_2
            i_left = i2;
            i_right = (n_1-1)*n_2 + i2;
            
            BC(2,i_left) = 1;            
            BC(1,i_right) = 1;
        end
        
        half_n_1 = ceil(n_1/2);
        
        for i1 = 1:half_n_1
            i_top = (i1-1)*n_2 + n_2;
            
            BC(1,i_top) = 1;
        end
        
        Neumann = zeros(d,4,n_el);
        
        for e2 = 1:n_el_2
            e_left = (n_el_1-1)*n_el_2 + e2;
            e_right = (n_el_1-1)*n_el_2 + e2;
            
            Neumann(1,4,e_left) = 1;
            Neumann(2,2,e_right) = 1;
        end
        
        for e1 = 1:n_el_1
            e_bottom = (e1-1)*n_el_2 + 1;
            
            Neumann(1,1,e_bottom) = 1;
            Neumann(2,1,e_bottom) = 1;            
        end
        
        half_n_el_1 = n_el_1/2;
        
        for e1 = 1:half_n_el_1
            e_top = (e1-1)*n_el_2 + n_el_2;
            
            Neumann(2,3,e_top) = 1;
        end
        
        for e1 = half_n_el_1+1:n_el_1
            e_top = (e1-1)*n_el_2 + n_el_2;
            
            Neumann(1,3,e_top) = 1;
            Neumann(2,3,e_top) = 1;
        end
end