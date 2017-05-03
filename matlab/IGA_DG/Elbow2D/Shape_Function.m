%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: Shape_Function
%
% Input:  xi_1 = first component of evaluation point
%         xi_2 = second component of evaluation point
%         p_1 = polynomial degree in direction 1
%         p_2 = polynomial degree in direction 2
%         C_e = element extraction matrix
%         P_b = array of element Bezier control points
%         w_b = array of element Bezier weights
%         w_e = array of element NURBS weights
%
% Output: R = values of NURBS at (xi_1,xi_2)
%         dRdx = derivatives of NURBS at (xi_1,xi_2)
%         x = position at (xi_1,xi_2)
%         J = Jacobian matrix at (xi_1,xi_2)
%
% Purpose: Evaluate NURBS and derivatives at a point.
%
% Notes: Algorithm follows directly from notes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R,dRdx,x,J] = Shape_Function(xi_1,xi_2,p_1,p_2,C_e,P_b,w_b,w_e)

n_loc = (p_1+1)*(p_2+1);

%%%
% Initialize basis functions, derivatives, position, and Jacobian

R = zeros(n_loc,1);
dRdxi = zeros(n_loc,2);

x = zeros(1,2);
J = zeros(2,2);

wb = 0;
dwbdxi = zeros(1,2);

%%%
% Compute Bernstein basis functions and their parametric derivatives

[B,dBdxi] = Bernstein_Basis_And_Derivs(xi_1,xi_2,p_1,p_2);

%%%
% Compute weighting function and its parametric derivatives

for a = 1:n_loc
    wb = wb + w_b(a)*B(a);
    dwbdxi(1) = dwbdxi(1) + w_b(a)*dBdxi(a,1);
    dwbdxi(2) = dwbdxi(2) + w_b(a)*dBdxi(a,2);
end

%%%
% Compute NURBS basis functions and their parametric derivatives

for a = 1:n_loc
    for b = 1:n_loc
        R(a) = R(a) + w_e(a)*C_e(a,b)*B(b)/wb;
        
        dRdxi(a,1) = dRdxi(a,1) + w_e(a)*C_e(a,b)*(dBdxi(b,1)/wb - dwbdxi(1)*B(b)/wb^2);
        dRdxi(a,2) = dRdxi(a,2) + w_e(a)*C_e(a,b)*(dBdxi(b,2)/wb - dwbdxi(2)*B(b)/wb^2);
    end
end

%%%
% Compute physical space quantities

for a = 1:n_loc
    x(1) = x(1) + w_b(a)*P_b(a,1)*B(a)/wb;
    x(2) = x(2) + w_b(a)*P_b(a,2)*B(a)/wb;
    
    J(1,1) = J(1,1) + w_b(a)*P_b(a,1)*(dBdxi(a,1)/wb - dwbdxi(1)*B(a)/wb^2);
    J(1,2) = J(1,2) + w_b(a)*P_b(a,1)*(dBdxi(a,2)/wb - dwbdxi(2)*B(a)/wb^2);
    J(2,1) = J(2,1) + w_b(a)*P_b(a,2)*(dBdxi(a,1)/wb - dwbdxi(1)*B(a)/wb^2);
    J(2,2) = J(2,2) + w_b(a)*P_b(a,2)*(dBdxi(a,2)/wb - dwbdxi(2)*B(a)/wb^2);
end

dxidx = inv(J);

%%%
% Compute physical space derivatives of NURBS basis functions

for a = 1:n_loc
    dRdx(a,1) = dRdxi(a,1)*dxidx(1,1) + dRdxi(a,2)*dxidx(2,1);
    dRdx(a,2) = dRdxi(a,1)*dxidx(1,2) + dRdxi(a,2)*dxidx(2,2);
end