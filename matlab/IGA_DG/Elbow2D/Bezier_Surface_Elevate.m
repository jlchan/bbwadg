%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: Bezier_Surface_Elevate
%
% Input:  d = number of spatial dimensions
%         p_1 = polynomial degree in direction 1
%         p_2 = polynomial degree in direction 2
%         P = array of control points (single-indexed)
%         w = array of weights (single-indexed)
%
% Output: p_1_new = new polynomial degree in direction 1
%         p_2_new = new polynomial degree in direction 1
%         P_new = new array of control points (single-indexed)
%         w_new = new array of weights (single-indexed)
%
% Purpose: Elevate the degree of a given Bezier surface.
%
% Notes: Utilizes the formula from the last pages of the
%        notes on Bernstein Polynomials and Bezier curves.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p_1_new,p_2_new,P_new,w_new] = Bezier_Surface_Elevate(d,p_1,p_2,P,w)

%%%
% Transform to Homogeneous Multi-Index

[P_multi,w_multi] = Transform_Single_to_Multi(d,P,w,p_1+1,p_2+1);
[P_w] = Transform_to_Homogeneous(d,P_multi,w_multi,2);

%%%
% Degree elevate in direction 1

p_1_new = p_1 + 1;

Q_w = zeros(p_1_new+1,p_2+1,d+1);

Q_w(1,:,:) = P_w(1,:,:);
for i = 1:p_1
    Q_w(i+1,:,:) = i/(p_1+1)*P_w(i,:,:) + (1-i/(p_1+1))*P_w(i+1,:,:);
end
Q_w(p_1_new+1,:,:) = P_w(p_1+1,:,:);

%%%
% Degree elevate in direction 2

p_2_new = p_2 + 1;

R_w = zeros(p_1_new+1,p_2_new+1,d+1);

R_w(:,1,:) = Q_w(:,1,:);
for i = 1:p_2
    R_w(:,i+1,:) = i/(p_2+1)*Q_w(:,i,:) + (1-i/(p_2+1))*Q_w(:,i+1,:);
end
R_w(:,p_2_new+1,:) = Q_w(:,p_2+1,:);

%%%
% Transform Back to Physical Single-Index

[P_multi_new,w_multi_new] = Transform_from_Homogeneous(d,R_w,2);
[P_new,w_new] = Transform_Multi_to_Single(d,P_multi_new,w_multi_new,p_1_new+1,p_2_new+1);