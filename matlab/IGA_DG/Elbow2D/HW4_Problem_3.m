%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Solution to Homework 4, Problem 3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Problem and Geometry Parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

problem = 1;
dim = 2;

r_i = 75/1000;
t = 15/1000;
r_o = r_i+t;
r_m = 0.5*(r_i+r_o);

E = 200*10^9;
nu = 0.3;
p_i = 30*10^6;

state = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Original Bezier Control Points and Weights
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P_original(:,1) = [-r_i;-r_m;-r_o;-r_i;-r_m;-r_o;0;0;0];
P_original(:,2) = [0;0;0;r_i;r_m;r_o;r_i;r_m;r_o];

w_original = [1;1;1;cos(pi/4);cos(pi/4);cos(pi/4);1;1;1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Convergence Study for p = 2 (6 Levels of Refinement)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
% Geometric parameters for first level

P = P_original;
w = w_original;

p_1 = 2;
p_2 = 2;
n_1 = 3;
n_2 = 3;
Xi_1 = [0,0,0,1,1,1];
Xi_2 = [0,0,0,1,1,1];

%%%
% Number of quadrature points

n_q = 3;

%%%
% Solve linear elasticity problem

d = Linear_Elasticity(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,n_q,problem);

%%%
% Evaluate L2-norm of error

[L2x,L2y] = L2_Error(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,n_q,d,r_i,r_o,E,nu,p_i);
L2x_p2(1) = L2x;
L2y_p2(1) = L2y;
n_p2(1) = dim*n_1*n_2;

%%%
% Iterate to higher levels of refinement

for level = 2:6
    
    %%%
    % Uniformly refine the knot vectors
    
    uXi_1 = unique(Xi_1);
    uXi_2 = unique(Xi_2);
    
    num_add_1 = size(uXi_1,2) - 1;
    num_add_2 = size(uXi_2,2) - 1;
    
    add_Xi_1 = zeros(1,num_add_1);
    for i = 1:num_add_1
        add_Xi_1(i) = 0.5*(uXi_1(i) + uXi_1(i+1));
    end
    
    add_Xi_2 = zeros(1,num_add_2);
    for i = 1:num_add_2
        add_Xi_2(i) = 0.5*(uXi_2(i) + uXi_2(i+1));
    end
    
    [P,w] = Transform_Single_to_Multi(dim,P,w,n_1,n_2);
    [n_1,n_2,Xi_1,Xi_2,P,w] = NURBS_Surface_Refine(dim,add_Xi_1,add_Xi_2,p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w);
    [P,w] = Transform_Multi_to_Single(dim,P,w,n_1,n_2);
    
    %%%
    % Solve linear elasticity problem
    
    d = Linear_Elasticity(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,n_q,problem);
    
    %%%
    % Evaluate L2-norm of error
    
    [L2x,L2y] = L2_Error(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,n_q,d,r_i,r_o,E,nu,p_i);
    L2x_p2(level) = L2x;
    L2y_p2(level) = L2y;
    n_p2(level) = dim*n_1*n_2;
    
    %%%
    % Plot elasticity fields for third level
    
    if level == 3
        Plot_Elasticity(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,d,E,nu,state,1,50)
        Plot_Elasticity(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,d,E,nu,state,2,50)
        Plot_Elasticity(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,d,E,nu,state,8,50)
        Plot_Elasticity(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,d,E,nu,state,9,50)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Convergence Study for p = 3 (6 Levels of Refinement)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
% Geometric parameters for first level - obtained via degree elevation

[p_1,p_2,P,w] = Bezier_Surface_Elevate(dim,p_1,p_2,P_original,w_original);

n_1 = 4;
n_2 = 4;
Xi_1 = [0,0,0,0,1,1,1,1];
Xi_2 = [0,0,0,0,1,1,1,1];

%%%
% Number of quadrature points

n_q = 4;

%%%
% Solve linear elasticity problem

d = Linear_Elasticity(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,n_q,problem);

%%%
% Evaluate L2-norm of error

[L2x,L2y] = L2_Error(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,n_q,d,r_i,r_o,E,nu,p_i);
L2x_p3(1) = L2x;
L2y_p3(1) = L2y;
n_p3(1) = dim*n_1*n_2;

%%%
% Iterate to higher levels of refinement

for level = 2:6
    
    %%%
    % Uniformly refine the knot vectors
    
    uXi_1 = unique(Xi_1);
    uXi_2 = unique(Xi_2);
    
    num_add_1 = size(uXi_1,2) - 1;
    num_add_2 = size(uXi_2,2) - 1;
    
    add_Xi_1 = zeros(1,num_add_1);
    for i = 1:num_add_1
        add_Xi_1(i) = 0.5*(uXi_1(i) + uXi_1(i+1));
    end
    
    add_Xi_2 = zeros(1,num_add_2);
    for i = 1:num_add_2
        add_Xi_2(i) = 0.5*(uXi_2(i) + uXi_2(i+1));
    end
    
    [P,w] = Transform_Single_to_Multi(dim,P,w,n_1,n_2);
    [n_1,n_2,Xi_1,Xi_2,P,w] = NURBS_Surface_Refine(dim,add_Xi_1,add_Xi_2,p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w);
    [P,w] = Transform_Multi_to_Single(dim,P,w,n_1,n_2);
    
    %%%
    % Solve linear elasticity problem
    
    d = Linear_Elasticity(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,n_q,problem);
    
    %%%
    % Evaluate L2-norm of error
    
    [L2x,L2y] = L2_Error(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,n_q,d,r_i,r_o,E,nu,p_i);
    L2x_p3(level) = L2x;
    L2y_p3(level) = L2y;
    n_p3(level) = dim*n_1*n_2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot Convergence Study Results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
% First plot the x-Displacement error
%%%

%%%
% Plot L2-norm of error vs NDOF for p = 2 and p = 3

figure
loglog(n_p2,L2x_p2,n_p3,L2x_p3,'LineWidth',2)
hold on

%%%
% Draw lines illustrating slope of convergence for p = 2 and p = 3

pert = 3;

line([n_p2(5),n_p2(6)],[L2x_p2(5)*pert,L2x_p2(6)*pert],'LineWidth',2,'Color','k')
line([n_p2(5),n_p2(6)],[L2x_p2(5)*pert,L2x_p2(5)*pert],'LineWidth',2,'Color','k')
line([n_p2(6),n_p2(6)],[L2x_p2(5)*pert,L2x_p2(6)*pert],'LineWidth',2,'Color','k')
text((n_p2(5)+n_p2(6))/2,L2x_p2(5)*pert*2,'1','FontName','Times New Roman','FontSize',20)
text(n_p2(6)*1.1,pert*(L2x_p2(5)+L2x_p2(6))/2,'1.50','FontName','Times New Roman','FontSize',20)

line([n_p3(5),n_p3(6)],[L2x_p3(5)*pert,L2x_p3(6)*pert],'LineWidth',2,'Color','k')
line([n_p3(5),n_p3(6)],[L2x_p3(5)*pert,L2x_p3(5)*pert],'LineWidth',2,'Color','k')
line([n_p3(6),n_p3(6)],[L2x_p3(5)*pert,L2x_p3(6)*pert],'LineWidth',2,'Color','k')
text((n_p3(5)+n_p3(6))/2,L2x_p3(5)*pert*2,'1','FontName','Times New Roman','FontSize',20)
text(n_p3(6)*1.1,pert*(L2x_p3(5)+L2x_p3(6))/2,'1.94','FontName','Times New Roman','FontSize',20)

hold off

%%%
% Label axes and apply legend

xlabel('Number of Degrees of Freedom','FontName','Times New Roman','FontSize',20)
ylabel('{\it{L}}^2-norm of the x-Displacement Error','FontName','Times New Roman','FontSize',20)

h = legend('{\it{p}} = 2','{\it{p}} = 3','Location','northeast');

set(h,'FontName','Times New Roman','FontSize',20)
set(gca,'FontName','Times New Roman','FontSize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
% Next plot the y-Displacement error
%%%

%%%
% Plot L2-norm of error vs NDOF for p = 2 and p = 3

figure
loglog(n_p2,L2y_p2,n_p3,L2y_p3,'LineWidth',2)
hold on

%%%
% Draw lines illustrating slope of convergence for p = 2 and p = 3

pert = 3;

line([n_p2(5),n_p2(6)],[L2y_p2(5)*pert,L2y_p2(6)*pert],'LineWidth',2,'Color','k')
line([n_p2(5),n_p2(6)],[L2y_p2(5)*pert,L2y_p2(5)*pert],'LineWidth',2,'Color','k')
line([n_p2(6),n_p2(6)],[L2y_p2(5)*pert,L2y_p2(6)*pert],'LineWidth',2,'Color','k')
text((n_p2(5)+n_p2(6))/2,L2y_p2(5)*pert*2,'1','FontName','Times New Roman','FontSize',20)
text(n_p2(6)*1.1,pert*(L2y_p2(5)+L2y_p2(6))/2,'1.50','FontName','Times New Roman','FontSize',20)

line([n_p3(5),n_p3(6)],[L2y_p3(5)*pert,L2y_p3(6)*pert],'LineWidth',2,'Color','k')
line([n_p3(5),n_p3(6)],[L2y_p3(5)*pert,L2y_p3(5)*pert],'LineWidth',2,'Color','k')
line([n_p3(6),n_p3(6)],[L2y_p3(5)*pert,L2y_p3(6)*pert],'LineWidth',2,'Color','k')
text((n_p3(5)+n_p3(6))/2,L2y_p3(5)*pert*2,'1','FontName','Times New Roman','FontSize',20)
text(n_p3(6)*1.1,pert*(L2y_p3(5)+L2y_p3(6))/2,'1.94','FontName','Times New Roman','FontSize',20)

hold off

%%%
% Label axes and apply legend

xlabel('Number of Degrees of Freedom','FontName','Times New Roman','FontSize',20)
ylabel('{\it{L}}^2-norm of the y-Displacement Error','FontName','Times New Roman','FontSize',20)

h = legend('{\it{p}} = 2','{\it{p}} = 3','Location','northeast');

set(h,'FontName','Times New Roman','FontSize',20)
set(gca,'FontName','Times New Roman','FontSize',16)