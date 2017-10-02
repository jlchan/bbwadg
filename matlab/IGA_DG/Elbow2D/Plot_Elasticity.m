%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: Plot_Elasticity
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
% Output: N/A
%
% Purpose: Plot the elasticity solution
%
% Notes: N/A
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Plot_Elasticity(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,d,E,nu,state,field,amp)

%%%
% Define Lame parameters
mu = E/(2*(1+nu));
lambda = (nu*E)/((1+nu)*(1-2*nu));

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
% Determine the plot points on the parent element

x_plot = 0:0.1:1;
y_plot = 0:0.1:1;

%%%
% Determine the X, Y, and T matrices for plotting

X_mat = [];
Y_mat = [];
Z_mat = [];

for e1 = 1:n_el_1
    X_row = [];
    Y_row = [];
    Z_row = [];
  
    for e2 = 1:n_el_2
        e = (e1-1)*n_el_2 + e2;        
        
        X = zeros(11,11);
        Y = zeros(11,11);
        Z = zeros(11,11);
        
        for i1 = 1:11
            xi_1 = x_plot(i1);
            
            for i2 = 1:11
                xi_2 = y_plot(i2);
                
                [R,dRdx,x,~] = Shape_Function(xi_1,xi_2,p_1,p_2,C_operators(:,:,e),P_b(:,:,e),w_b(:,e),w_e(:,e));                                
                
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
                
                X(i1,i2) = x(1) + amp*u;
                Y(i1,i2) = x(2) + amp*v;
                
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
                
                Sigmat = [Sigmaxx Sigmaxy; Sigmaxy Sigmayy];                
                Principal = eig(Sigmat);                
                
                vonMises = sqrt(0.5*((Sigmaxx-Sigmayy)^2+(Sigmayy-Sigmazz)^2+(Sigmazz-Sigmaxx)^2+6*Sigmaxy^2));
                
                switch field
                    case 1
                        Z(i1,i2) = u;
                    case 2
                        Z(i1,i2) = v;
                    case 3
                        Z(i1,i2) = sqrt(u^2+v^2);
                    case 4
                        Z(i1,i2) = Sigmaxx;
                    case 5
                        Z(i1,i2) = Sigmaxy;
                    case 6
                        Z(i1,i2) = Sigmayy;
                    case 7
                        Z(i1,i2) = Sigmazz;
                    case 8
                        Z(i1,i2) = max(Principal);
                    case 9
                        Z(i1,i2) = min(Principal);
                    case 10
                        Z(i1,i2) = vonMises;                        
                end
                        
            end            
        end        
        
        X_row = [X_row X];
        Y_row = [Y_row Y];
        Z_row = [Z_row Z];        
    end
    
    X_mat = [X_mat; X_row];
    Y_mat = [Y_mat; Y_row];
    Z_mat = [Z_mat; Z_row];
end

%%%
% Remove repeated rows/columns in the X,Y, and T matrices

index_1 = 1:11;

for e1 = 2:n_el_1
    first = (e1-1)*11+2;
    index_1 = [index_1 first:first+9];
end

index_2 = 1:11;

for e2 = 2:n_el_2
    first = (e2-1)*11+2;
    index_2 = [index_2 first:first+9];
end

X_mat = X_mat(index_1,index_2);
Y_mat = Y_mat(index_1,index_2);
Z_mat = Z_mat(index_1,index_2);

Z_min = min(min(Z_mat));
Z_max = max(max(Z_mat));

%%%
% Contour plot of the desired field with colorbar

figure
hold on
axis equal

contourf(X_mat,Y_mat,Z_mat,20,'EdgeColor','none');
caxis([Z_min Z_max])
colorbar

%%%
% Plot element boundaries (isoparametric lines)

for e1 = 1:n_el_1
  
    for e2 = 1:n_el_2
        e = (e1-1)*n_el_2 + e2;
        
        for i1 = 1:11
            xi_1 = x_plot(i1);
            
            for i2 = 1:11
                xi_2 = y_plot(i2);
                
                [R,~,x,~] = Shape_Function(xi_1,xi_2,p_1,p_2,C_operators(:,:,e),P_b(:,:,e),w_b(:,e),w_e(:,e));                                
                
                u = 0;
                v = 0;
                
                for a = 1:n_loc
                    u = u + R(a)*d(ID(1,IEN(a,e)));
                    v = v + R(a)*d(ID(2,IEN(a,e)));                                                         
                end
                
                X(i1,i2) = x(1) + amp*u;
                Y(i1,i2) = x(2) + amp*v;               
            end
            
        end
        
        for i = 1:10
            line([X(i,1),X(i+1,1)],[Y(i,1),Y(i+1,1)],'Color','k');
            line([X(i,11),X(i+1,11)],[Y(i,11),Y(i+1,11)],'Color','k');
        end
        
        for j = 1:10
            line([X(1,j),X(1,j+1)],[Y(1,j),Y(1,j+1)],'Color','k');
            line([X(11,j),X(11,j+1)],[Y(11,j),Y(11,j+1)],'Color','k');
        end
    end
    
end

% Label axes and establish a title

xlabel('{\it{x}} (in meters)','FontName','Times New Roman','FontSize',20)
ylabel('{\it{y}} (in meters)','FontName','Times New Roman','FontSize',20)

switch field
    case 1
        Desired_Field = 'x-displacement Field (in meters)';
    case 2
        Desired_Field = 'y-displacement Field (in meters)';
    case 3
        Desired_Field = 'Displacement Magnitude (in meters)';
    case 4
        Desired_Field = 'Stress Field Normal to x (in Pascals)';
    case 5
        Desired_Field = 'In-Plane Shear Stress Field (in Pascals)';
    case 6
        Desired_Field = 'Stress Field Normal to y (in Pascals)';
    case 7
        Desired_Field = 'Stress Field Normal to z (in Pascals)';
    case 8
        Desired_Field = 'Maximum Principal Stress (in Pascals)';
    case 9
        Desired_Field = 'Minimum Principal Stress (in Pascals)';
    case 10
        Desired_Field = 'Von Mises Stress (in Pascals)';
end

s = sprintf(Desired_Field, char(176));
title(s,'FontName','Times New Roman','FontSize',20)

set(gca,'FontName','Times New Roman','FontSize',16)

hold off