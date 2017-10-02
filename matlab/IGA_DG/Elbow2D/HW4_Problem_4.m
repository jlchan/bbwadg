%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Solution to Homework 4, Problem 4
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Problem and Geometry Parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dim = 2;

r_t = 5;
l_a = 3;
h_c = 1.5;
l = r_t+l_a;
h = r_t+h_c;

E = 30*10^9;
nu = 0.2;

state = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Original Composite Bezier Geometry Information
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P_original(:,1) = [-r_t,-0.5*(r_t+l),-l,-r_t,-0.5*(r_t+l),-l,-r_t*cos(pi/4),-0.5*(r_t*cos(pi/4)+l),-l,-r_t*(sqrt(2)-1),-0.5*(r_t*(sqrt(2)-1)+0.5*l),-0.5*l,0,0,0];

P_original(:,2) = [0,0,0,r_t*(sqrt(2)-1),0.5*(r_t*(sqrt(2)-1)+0.5*h),0.5*h,r_t*cos(pi/4),0.5*(r_t*cos(pi/4)+h),h,r_t,0.5*(r_t+h),h,r_t,0.5*(r_t+h),h];
               
w_original = [1,1,1,cos(pi/8),1,1,1,1,1,cos(pi/8),1,1,1,1,1];

p_1 = 2;
p_2 = 2;
n_1_original = 5;
n_2_original = 3;
Xi_1_original = [0,0,0,1,1,2,2,2];
Xi_2_original = [0,0,0,1,1,1];

%%%
% Number of quadrature points

n_q = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% **Part 1** Pin-Supported Sides
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
% Problem type

problem = 2;

%%%
% Geometric parameters for first level

P = P_original;
w = w_original;

n_1 = n_1_original;
n_2 = n_2_original;
Xi_1 = Xi_1_original;
Xi_2 = Xi_2_original;

%%%
% Solve linear elasticity problem

d = Linear_Elasticity(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,n_q,problem);

%%%
% Evaluate max displacement magnitude and von Mises stress

[max_d,max_vm,max_d_arc,max_vm_arc] = Max_Elasticity(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,d,E,nu,state);

max_d_part1(1) = max_d;
max_vm_part1(1) = max_vm;
max_d_arc_part1(1) = max_d_arc;
max_vm_arc_part1(1) = max_vm_arc;
n_part1(1) = dim*n_1*n_2;

%%%
% Iterate to higher levels of refinement

for level = 2:5 % Iterate to 8 for "converged" results
    
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
    % Evaluate max displacement magnitude and von Mises stress
    
    [max_d,max_vm,max_d_arc,max_vm_arc] = Max_Elasticity(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,d,E,nu,state);
    
    max_d_part1(level) = max_d;
    max_vm_part1(level) = max_vm;
    max_d_arc_part1(level) = max_d_arc;
    max_vm_arc_part1(level) = max_vm_arc;
    n_part1(level) = dim*n_1*n_2;
    
    %%%
    % Plot elasticity fields for fifth level
    
    if level == 5
        Plot_Elasticity(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,d,E,nu,state,1,3000)
        Plot_Elasticity(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,d,E,nu,state,2,3000)
        Plot_Elasticity(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,d,E,nu,state,10,3000)
    end
end

%%%
% Plot max displacement vs NDOF

figure
semilogx(n_part1,max_d_part1,'LineWidth',2)
xlabel('Number of Degrees of Freedom','FontName','Times New Roman','FontSize',20)
ylabel('Computed Max Displacement','FontName','Times New Roman','FontSize',20)
title('Max Displacement for Problem 4 Part 1','FontName','Times New Roman','FontSize',20)
set(gca,'FontName','Times New Roman','FontSize',16)

%%%
% Plot max von Mises stress vs NDOF

figure
semilogx(n_part1,max_vm_part1,'LineWidth',2)
xlabel('Number of Degrees of Freedom','FontName','Times New Roman','FontSize',20)
ylabel('Computed Max von Mises Stress','FontName','Times New Roman','FontSize',20)
title('Max von Mises Stress for Problem 4 Part 1','FontName','Times New Roman','FontSize',20)
set(gca,'FontName','Times New Roman','FontSize',16)

%%%
% Plot arc max displacement vs NDOF

figure
semilogx(n_part1,max_d_arc_part1,'LineWidth',2)
xlabel('Number of Degrees of Freedom','FontName','Times New Roman','FontSize',20)
ylabel('Computed Max Displacement','FontName','Times New Roman','FontSize',20)
title('Arc Max Displacement for Problem 4 Part 1','FontName','Times New Roman','FontSize',20)
set(gca,'FontName','Times New Roman','FontSize',16)

%%%
% Plot arc max von Mises stress vs NDOF

figure
semilogx(n_part1,max_vm_arc_part1,'LineWidth',2)
xlabel('Number of Degrees of Freedom','FontName','Times New Roman','FontSize',20)
ylabel('Computed Max von Mises Stress','FontName','Times New Roman','FontSize',20)
title('Arc Max Stress for Problem 4 Part 1','FontName','Times New Roman','FontSize',20)
set(gca,'FontName','Times New Roman','FontSize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% **Part 2** Roller-Supported Sides
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
% Problem type

problem = 3;

%%%
% Geometric parameters for first level

P = P_original;
w = w_original;

n_1 = n_1_original;
n_2 = n_2_original;
Xi_1 = Xi_1_original;
Xi_2 = Xi_2_original;

%%%
% Solve linear elasticity problem

d = Linear_Elasticity(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,n_q,problem);

%%%
% Evaluate max displacement magnitude and von Mises stress

[max_d,max_vm,max_d_arc,max_vm_arc] = Max_Elasticity(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,d,E,nu,state);

max_d_part2(1) = max_d;
max_vm_part2(1) = max_vm;
max_d_arc_part2(1) = max_d_arc;
max_vm_arc_part2(1) = max_vm_arc;
n_part2(1) = dim*n_1*n_2;

%%%
% Iterate to higher levels of refinement

for level = 2:5 % Iterate to 8 for "converged" results
    
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
    % Evaluate max displacement magnitude and von Mises stress
    
    [max_d,max_vm,max_d_arc,max_vm_arc] = Max_Elasticity(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,d,E,nu,state);
    
    max_d_part2(level) = max_d;
    max_vm_part2(level) = max_vm;
    max_d_arc_part2(level) = max_d_arc;
    max_vm_arc_part2(level) = max_vm_arc;
    n_part2(level) = dim*n_1*n_2;
    
    %%%
    % Plot elasticity fields for fifth level
    
    if level == 5
        Plot_Elasticity(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,d,E,nu,state,1,3000)
        Plot_Elasticity(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,d,E,nu,state,2,3000)
        Plot_Elasticity(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,d,E,nu,state,10,3000)
    end
end

%%%
% Plot max displacement vs NDOF

figure
semilogx(n_part2,max_d_part2,'LineWidth',2)
xlabel('Number of Degrees of Freedom','FontName','Times New Roman','FontSize',20)
ylabel('Computed Max Displacement','FontName','Times New Roman','FontSize',20)
title('Max Displacement for Problem 4 Part 2','FontName','Times New Roman','FontSize',20)
set(gca,'FontName','Times New Roman','FontSize',16)

%%%
% Plot max von Mises stress vs NDOF

figure
semilogx(n_part2,max_vm_part2,'LineWidth',2)
xlabel('Number of Degrees of Freedom','FontName','Times New Roman','FontSize',20)
ylabel('Computed Max von Mises Stress','FontName','Times New Roman','FontSize',20)
title('Max von Mises Stress for Problem 4 Part 2','FontName','Times New Roman','FontSize',20)
set(gca,'FontName','Times New Roman','FontSize',16)

%%%
% Plot arc max displacement vs NDOF

figure
semilogx(n_part2,max_d_arc_part2,'LineWidth',2)
xlabel('Number of Degrees of Freedom','FontName','Times New Roman','FontSize',20)
ylabel('Computed Max Displacement','FontName','Times New Roman','FontSize',20)
title('Arc Max Displacement for Problem 4 Part 2','FontName','Times New Roman','FontSize',20)
set(gca,'FontName','Times New Roman','FontSize',16)

%%%
% Plot arc max von Mises stress vs NDOF

figure
semilogx(n_part2,max_vm_arc_part2,'LineWidth',2)
xlabel('Number of Degrees of Freedom','FontName','Times New Roman','FontSize',20)
ylabel('Computed Max von Mises Stress','FontName','Times New Roman','FontSize',20)
title('Arc Max Stress for Problem 4 Part 2','FontName','Times New Roman','FontSize',20)
set(gca,'FontName','Times New Roman','FontSize',16)