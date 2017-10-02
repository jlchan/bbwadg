%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: Linear_Elasticity
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
%         problem = problem description index
%
% Output: d = vector of NURBS control variables
%
% Purpose: Solve the 2-D linear elasticity problem
%
% Notes: Algorithm follows directly from notes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [d] = Linear_Elasticity(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,n_q,problem)

n = n_1*n_2;
dim = 2;
ndof = dim*n;

%%%
% Step 1: Initialize the Problem

[n_el,C_operators,IEN,P_b,w_b,~,w_e] = Extract_And_Localize(n_1,n_2,p_1,p_2,Xi_1,Xi_2,P,w);
[ID] = Construct_ID(dim,n);

[D,f,h] = Elasticity_Parameters(problem);
[xi_q,w_q] = Quadrature_Data(n_q);
[BC,g,Neumann] = Boundary_Arrays(problem,n_1,n_2,Xi_1,Xi_2);

%%%
% Step 2: Construct the Matrix System

K = zeros(ndof,ndof); % nzmax = dim*(2*p_1+1)*(2*p_2+1)*ndof; K = spalloc(ndof,ndof,nzmax); % Utilize commented portion for higher levels of refinement
F = zeros(ndof,1);

for e = 1:n_el    
    [k_e,f_e] = Element_Formation(p_1,p_2,C_operators(:,:,e),P_b(:,:,e),w_b(:,e),w_e(:,e),n_q,xi_q,w_q,Neumann(:,:,e),D,f,h);
    [K,F] = Element_Assembly(e,k_e,f_e,IEN,ID,BC,g,K,F);
end

%%%
% Step 2.1: Account for Dirichlet BCs

for i = 1:n
    for A = 1:dim
        if BC(A,i) == 1
            P = ID(A,i);
            K(P,P) = 1;
            F(P) = g(A,i);
        end
    end
end

%%%
% Step 3: Solve the Matrix System

d = K\F;