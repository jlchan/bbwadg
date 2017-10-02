%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: Max_Elasticity
%
% Input:  p_1 = polynomial degree in direction 1
%         p_2 = polynomial degree in direction 2
%         n_el = number of elements
%         P = array of control points
%         w = array of weights
%         E = Young's modulus
%         nu = Poisson ratio
%         state = integer for plane stress/strain
%         field = integer for viz field
%         amp = displacement amplification factor
%
% Output: max_d = maximum displacement magnitude
%         max_vm = maximum von Mises stress
%         max_d_arc = maximum displacement magnitude in arc
%         max_vm_arc = maximum von Mises stress in arc
%
% Purpose: Evaluate max displacement and von Mises stress
%
% Notes: N/A
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [max_d,max_vm,max_d_arc,max_vm_arc] = Max_Elasticity(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,d,E,nu,state)

%%%
% Define Lame parameters
mu = E/(2*(1+nu));
lambda = (nu*E)/((1+nu)*(1-2*nu));

%%%
% Define inner radius
r_i = 5;
tol = 1e-10;

%%%
% Construct stiffness matrix
if state == 1   
    D = [2*mu+lambda lambda 0; lambda 2*mu+lambda 0; 0 0 mu];
elseif state == 2
    D = E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
end

%%%
% Extract the basis, geometry, and temperature field
[~,C_operators,IEN,P_b,w_b,~,w_e] = Extract_And_Localize(n_1,n_2,p_1,p_2,Xi_1,Xi_2,P,w);

n_el_1 = length(unique(Xi_1))-1;
n_el_2 = length(unique(Xi_2))-1;

n_loc = (p_1+1)*(p_2+1);

%%%
% Construct the ID array
n = n_1*n_2;
dim = 2;
ndof = dim*n;

[ID] = Construct_ID(dim,n);

%%%
% Determine the evaluation points on the parent element

x_eval = 0:0.1:1;
y_eval = 0:0.1:1;

%%%
% Initialize max values

max_d = 0;
max_vm = 0;
max_d_arc = 0;
max_vm_arc = 0;

%%%
% Loop over points, evaluate quantities, and find max values

for e1 = 1:n_el_1    
  
    for e2 = 1:n_el_2
        e = (e1-1)*n_el_2 + e2;
        
        for i1 = 1:11
            xi_1 = x_eval(i1);
            
            for i2 = 1:11
                xi_2 = y_eval(i2);
                
                [R,dRdx,x,~] = Shape_Function(xi_1,xi_2,p_1,p_2,C_operators(:,:,e),P_b(:,:,e),w_b(:,e),w_e(:,e));                                
                
                r = sqrt(x(1)^2+x(2)^2);
                
                u = 0;
                v = 0;
                
                Exx = 0;
                Eyy = 0;
                TwoExy = 0;
                
                for a = 1:n_loc
                    u = u + R(a)*d(ID(1,IEN(a,e)));
                    v = v + R(a)*d(ID(2,IEN(a,e)));
                    
                    Exx = Exx + dRdx(a,1)*d(ID(1,IEN(a,e)));
                    Eyy = Eyy + dRdx(a,2)*d(ID(2,IEN(a,e)));
                    TwoExy = TwoExy + dRdx(a,2)*d(ID(1,IEN(a,e))) + dRdx(a,1)*d(ID(2,IEN(a,e)));                                        
                end                
                
                Evec = [Exx;Eyy;TwoExy];
                Sigvec = D*Evec;
                
                Sigmaxx = Sigvec(1);
                Sigmayy = Sigvec(2);
                Sigmaxy = Sigvec(3);
                
                if state == 1
                    Sigmazz = lambda*(Exx+Eyy);
                else
                    Sigmazz = 0;
                end                                
                
                magDisp = sqrt(u^2+v^2);
                vonMises = sqrt(0.5*((Sigmaxx-Sigmayy)^2+(Sigmayy-Sigmazz)^2+(Sigmazz-Sigmaxx)^2+6*Sigmaxy^2));
                
                if magDisp > max_d
                    max_d = magDisp;
                end
                
                if vonMises > max_vm
                    max_vm = vonMises;
                end
                
                if abs(r-r_i) < tol                
                    if magDisp > max_d_arc
                        max_d_arc = magDisp;
                    end
                    
                    if vonMises > max_vm_arc
                        max_vm_arc = vonMises;
                    end
                end
            end
        end
        
    end
    
end