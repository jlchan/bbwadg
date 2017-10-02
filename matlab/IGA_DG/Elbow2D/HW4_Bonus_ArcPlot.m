%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot Various Arc Configurations Considered for the Bonus
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;

d = 2;

r_t = 5;
l_a = 3;
h_c = 1.5;
l = r_t+l_a;
h = r_t+h_c;

p = 2;
n = 3;
Xi = [0 0 0 1 1 1];

P(:,1) = [-r_t,-r_t,0];
P(:,2) = [0,r_t,r_t];

x_plot = 0:0.01:1;

for i = 1:size(x_plot,2)
    x = x_plot(i);
    B(:,i) = bernstein_poly(p,x);
end

w_vec = [0.1,0.25,0.5,sqrt(2)/2,1,2,4];

x_plot = [];
y_plot = [];

for itr = 1:size(w_vec,2)

    w_var = w_vec(itr);
    w = [1;w_var;1];

    wx_temp = B'*(P(:,1).*w);
    wy_temp = B'*(P(:,2).*w);
    w_temp = B'*w;

    x_plot = [x_plot wx_temp./w_temp];
    y_plot = [y_plot wy_temp./w_temp];
end

figure
plot(x_plot,y_plot,'LineWidth',1.5)
xlabel('{\it{x}} (in meters)','FontName','Times New Roman','FontSize',20)
ylabel('{\it{y}} (in meters)','FontName','Times New Roman','FontSize',20)
title('Variable Arc Designs','FontName','Times New Roman','FontSize',20)
set(gca,'FontName','Times New Roman','FontSize',16)
h = legend('w = 0.1','w = 0.25','w = 0.5','w = 0.7071','w = 1','w = 2','w = 4','Location','SouthEast');
set(h,'FontName','Times New Roman','FontSize',16)