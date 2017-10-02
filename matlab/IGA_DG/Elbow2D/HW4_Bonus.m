%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Solution to Homework 4, Bonus
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
problem = 3;
n_q = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Loop Over Vaiable Weights and Find Max Stress for Each
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w_vec = [0.05,0.1,0.2,0.3,0.4,0.5,0.55,0.6,0.65,sqrt(2)/2,0.8,0.9,1,1.2,1.4,1.8,2,2.5,3,4,6];

for itr = 1:size(w_vec,2)
    itr

    w_var = w_vec(itr);
    
    %%%
    % Geometric parameters for first level
    
    [P,w] = HW4_Bonus_Geometry(w_var);
    
    p_1 = 2; p_2 = 2; n_1 = 5; n_2 = 3;
    Xi_1 = [0,0,0,1,1,2,2,2]; Xi_2 = [0,0,0,1,1,1];
    
    %%%
    % Iterate to higher levels of refinement
    
    for level = 2:5
        
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
        
    end
    
    %%%
    % Solve linear elasticity problem for final level
    
    d = Linear_Elasticity(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,n_q,problem);
    
    %%%
    % Evaluate von Mises stress
    
    [~,max_vm,~,~] = Max_Elasticity(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,d,E,nu,state);
        
    max_vonMises(itr) = max_vm
end

[M,I] = min(max_vonMises);

figure
plot(w_vec,max_vonMises,'--k','LineWidth',1.5)
xlabel('Weight Value','FontName','Times New Roman','FontSize',20)
ylabel('Max von Mises Stress','FontName','Times New Roman','FontSize',20)
title('Max von Mises Stress for Bonus','FontName','Times New Roman','FontSize',20)
set(gca,'FontName','Times New Roman','FontSize',16)
hold on

scatter(w_vec(I),max_vonMises(I),'o','filled','SizeData',100)
scatter(w_vec(10),max_vonMises(10),'o','filled','SizeData',100)
h = legend('Interpolant','Optimum','Original','Location','SouthEast');
set(h,'FontName','Times New Roman','FontSize',16)
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot the Solution for the Optimum Weight
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w_var = w_vec(I);

%%%
% Geometric parameters for first level

[P,w] = HW4_Bonus_Geometry(w_var);

p_1 = 2; p_2 = 2; n_1 = 5; n_2 = 3;
Xi_1 = [0,0,0,1,1,2,2,2]; Xi_2 = [0,0,0,1,1,1];

%%%
% Iterate to higher levels of refinement

for level = 2:5
    
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
    
end

%%%
% Solve linear elasticity problem for final level

d = Linear_Elasticity(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,n_q,problem);

%%%
% Plot von Mises stress for final level

Plot_Elasticity(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,d,E,nu,state,10,0)
Plot_Elasticity(p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,w,d,E,nu,state,10,3000)