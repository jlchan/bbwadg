%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: Element_Formation
%
% Input:  p_1 = polynomial degree in direction 1
%         p_2 = polynomial degree in direction 2
%         C_e = element extraction matrix
%         P_b = array of element Bezier control points
%         w_b = array of element Bezier weights
%         w_e = array of element NURBS weights
%         n_q = number of quadrature points per direction
%         xi_q = array of quadrature points
%         w_q = array of quadrature weights
%         Neumann = Neumann BC array
%         D = stiffness tensor
%         f = body force
%         h = applied traction
%
% Output: k_e = element stiffness matrix
%         f_e = element forcing vector
%
% Purpose: Construct the element stiffness and forcing
%
% Notes: Algorithm follows directly from notes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [k_e,f_e] = Element_Formation(p_1,p_2,C_e,P_b,w_b,w_e,n_q,xi_q,w_q,Neumann,D,f,h)

n_loc = (p_1+1)*(p_2+1);
d = 2;
ndof_loc = d*n_loc;

%%%
% Initialize element stiffness and forcing

k_e = zeros(ndof_loc,ndof_loc);
f_e = zeros(ndof_loc,1);

%%%
% Loop over element interiors

for q_1 = 1:n_q
    for q_2 = 1:n_q        
        [R,dRdx,x,J] = Shape_Function(xi_q(q_1),xi_q(q_2),p_1,p_2,C_e,P_b,w_b,w_e);
        
        j = det(J);
        D_loc = D(x(1),x(2));
        f_loc = f(x(1),x(2));
        
        for a = 1:n_loc
            B_ae = [dRdx(a,1) 0; 0 dRdx(a,2); dRdx(a,2) dRdx(a,1)];
            
            for b = 1:n_loc
                B_be = [dRdx(b,1) 0; 0 dRdx(b,2); dRdx(b,2) dRdx(b,1)];
                
                k_loc = B_ae'*D_loc*B_be*w_q(q_1)*w_q(q_2)*j;
                
                for A = 1:d
                    for B = 1:d
                        p = d*(a-1)+A;
                        q = d*(b-1)+B;
                        
                        k_e(p,q) = k_e(p,q) + k_loc(A,B);
                    end
                end
                
            end
            
            for A = 1:d
                p = d*(a-1)+A;
                
                f_e(p) = f_e(p) + R(a)*f_loc(A)*w_q(q_1)*w_q(q_2)*j;
            end            
        end
        
    end
end

%%%
% Loop over element boundaries

for A = 1:d
    for side = 1:4
        
        if Neumann(A,side) == 1
            for q = 1:n_q
                switch side
                    case 1
                        xi_1 = xi_q(q);
                        xi_2 = 0;
                        tangent = [1;0];
                    case 2
                        xi_1 = 1;
                        xi_2 = xi_q(q);
                        tangent = [0;1];
                    case 3
                        xi_1 = xi_q(q);
                        xi_2 = 1;
                        tangent = [1;0];
                    case 4
                        xi_1 = 0;
                        xi_2 = xi_q(q);
                        tangent = [0;1];
                end
                
                [R,~,x,J] = Shape_Function(xi_1,xi_2,p_1,p_2,C_e,P_b,w_b,w_e);
                
                j_boundary = norm(J*tangent);
                
                h_loc = h(x(1),x(2));
                
                for a = 1:n_loc
                    p = d*(a-1)+A;
                    
                    f_e(p) = f_e(p) + R(a)*h_loc(A)*w_q(q)*j_boundary;
                end
                
            end
        end
    end
end