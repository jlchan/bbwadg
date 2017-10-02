%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: Elasticity_Parameters
%
% Input:  problem = problem description index
%
% Output: D = stiffness matrix
%         f = body force
%         h = applied traction
%
% Purpose: Set the 2-D linear elasticity problem parameters
%
% Notes: N/A
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [D, f, h] = Elasticity_Parameters(problem)

%%%
% Geometry parameters for Problem 3

r_i = 75/1000;

%%%
% Geometry parameters for Problem 4

r_t = 5;
h_c = 1.5;
height = r_t + h_c;

%%%
% Tolerance for geometry checks

eps = 1e-10;

switch problem
    
    %%%
    % Problem parameters for Problem 3
    
    case 1
        E = 200*10^9;
        nu = 0.3;
        p_i = 30*10^6;
        
        mu = E/(2*(1+nu));
        lambda = (nu*E)/((1+nu)*(1-2*nu));
        
        D_mat = [2*mu+lambda lambda 0; lambda 2*mu+lambda 0; 0 0 mu];
        
        D = @(x1,x2) D_mat;
        f = @(x1,x2) [0;0];
        h = @(x1,x2) [p_i*x1/sqrt(x1^2+x2^2);p_i*x2/sqrt(x1^2+x2^2)]*heaviside(r_i-sqrt(x1^2+x2^2)+eps);
        
    %%%
    % Problem parameters for Problem 4 Part 1
    
    case 2
        E = 30*10^9;
        nu = 0.2;
        rho = 2400;
        g = 9.81;
        p = 300*10^3;
        
        mu = E/(2*(1+nu));
        lambda = (nu*E)/((1+nu)*(1-2*nu));
        
        D_mat = [2*mu+lambda lambda 0; lambda 2*mu+lambda 0; 0 0 mu];
        
        D = @(x1,x2) D_mat;
        f = @(x1,x2) [0;-rho*g];
        h = @(x1,x2) [0;-p]*heaviside(x2-height+eps);
        
    %%%
    % Problem parameters for Problem 4 Part 2
    
    case 3
        E = 30*10^9;
        nu = 0.2;
        rho = 2400;
        g = 9.81;
        p = 300*10^3;
        
        mu = E/(2*(1+nu));
        lambda = (nu*E)/((1+nu)*(1-2*nu));
        
        D_mat = [2*mu+lambda lambda 0; lambda 2*mu+lambda 0; 0 0 mu];
        
        D = @(x1,x2) D_mat;
        f = @(x1,x2) [0;-rho*g];
        h = @(x1,x2) [0;-p*heaviside(x2-height+eps)]; 
end
        