%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: L2_Error
%
% Input:  p_1 = polynomial degree in direction 1
%         p_2 = polynomial degree in direction 2
%         n_1 = number of functions in direction 1
%         n_2 = number of functions in direction 2
%         Xi_1 = knot vector in direction 1
%         Xi_2 = knot vector in direction 2
%         P = array of NURBS control points (single-indexed)
%         w = array of NURBS weights (single-indexed)
%         n_q = number of quadrature points per direction
%         d = vector of NURBS control variables
%         r_i = innter radius
%         r_o = outer radius
%         E = Young's modulus
%         nu = Poisson ratio
%         p_i = internal pressure loading
%
% Output: L2x = L2 error of discrete x-displacement field
%         L2y = L2 error of discrete y-displacement field
%
% Purpose: Evaluate the L2 error for Problem 3
%
% Notes: Evaluates the L2 error using quadrature
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [L2x,L2y] = L2_Error(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,n_q,d,r_i,r_o,E,nu,p_i)

n_loc = (p_1+1)*(p_2+1);

%%%
% Extract, Localize, and Declare Quadrature Data

[n_el,C_operators,IEN,P_b,w_b,~,w_e] = Extract_And_Localize(n_1,n_2,p_1,p_2,Xi_1,Xi_2,P,w);
[xi_q,w_q] = Quadrature_Data(n_q);

%%%
% Construct the ID array
n = n_1*n_2;
dim = 2;
ndof = dim*n;

[ID] = Construct_ID(dim,n);

%%%
% Initialize L2 error

L2x = 0;
L2y = 0;

%%%
% Loop over elements and quadrature points

for e = 1:n_el
    for q_1 = 1:n_q
        for q_2 = 1:n_q
            [R,~,x,J] = Shape_Function(xi_q(q_1),xi_q(q_2),p_1,p_2,C_operators(:,:,e),P_b(:,:,e),w_b(:,e),w_e(:,e));
            
            j = det(J);                        
            
            %%%
            % Exact solution
            
            r = sqrt(x(1)^2 + x(2)^2);
            
            rat1 = r_o/r_i;
            rat2 = r_o/r;
            
            u_r = p_i*(1+nu)/E*r/(rat1^2-1)*((1-2*nu)+rat2^2);
            
            u_x = u_r*x(1)/r;
            u_y = u_r*x(2)/r;
            
            %%%
            % Numerical solution
            
            u_xh = 0;
            u_yh = 0;
            for a = 1:n_loc
                u_xh = u_xh + R(a)*d(ID(1,IEN(a,e)));
                u_yh = u_yh + R(a)*d(ID(2,IEN(a,e)));
            end
            
            %%%
            % Add quadrature point contribution to integral
            
            L2x = L2x + (u_x-u_xh)^2*w_q(q_1)*w_q(q_2)*j;
            L2y = L2y + (u_y-u_yh)^2*w_q(q_1)*w_q(q_2)*j;
            
        end
    end
end

%%%
% Square root of integral

L2x = sqrt(L2x);
L2y = sqrt(L2y);