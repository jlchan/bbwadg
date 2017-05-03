%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: HW4_Bonus_Geometry
%
% Input:  w_var = variable weight
%
% Output: P = array of control points for composite Bezier
%         w = array of weights for composite Bezier
%
% Purpose: Determine the bonus geometry control net.
%
% Notes: N/A
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P,w] = HW4_Bonus_Geometry(w_var)

d = 2;

r_t = 5;
l_a = 3;
h_c = 1.5;
l = r_t+l_a;
h = r_t+h_c;

p_1D = 2;
n_1D = 3;
Xi_1D = [0 0 0 1 1 1];

P_1D(:,1) = [-r_t,-r_t,0];
P_1D(:,2) = [0,r_t,r_t];

w_1D = [1;w_var;1];

add_Xi = [0.5 0.5];

[~,~,P_1D,w_1D] = NURBS_Curve_Refine(d,add_Xi,p_1D,n_1D,Xi_1D,P_1D,w_1D);

P(:,1) = [P_1D(1,1),0.5*(P_1D(1,1)-l),-l,...
          P_1D(2,1),0.5*(P_1D(2,1)-l),-l,...
          P_1D(3,1),0.5*(P_1D(3,1)-l),-l,...
          P_1D(4,1),0.5*(P_1D(4,1)-0.5*l),-0.5*l,...
          P_1D(5,1),0.5*(P_1D(5,1)),0];

P(:,2) = [P_1D(1,2),0.5*(P_1D(1,2)),0,...
          P_1D(2,2),0.5*(P_1D(2,2)+0.5*h),0.5*h,...
          P_1D(3,2),0.5*(P_1D(3,2)+h),h,...
          P_1D(4,2),0.5*(P_1D(4,2)+h),h,...
          P_1D(5,2),0.5*(P_1D(5,2)+h),h];

w = [w_1D(1),1,1,w_1D(2),1,1,w_1D(3),1,1,w_1D(4),1,1,w_1D(5),1,1];