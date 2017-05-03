%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: NURBS_Pipe_Elbow_Surface
%
% Input:  h = height
%         w = width
%         r = pipe radius
%         R = elbow turning radius
%
% Output: p_1 = polynomial degree in direction 1
%         p_2 = polynomial degree in direction 2
%         n_1 = number of functions in direction 1
%         n_2 = number of functions in direction 2
%         Xi_1 = knot vector in direction 1
%         Xi_2 = knot vector in direction 2
%         P = array of NURBS control points (multi-indexed)
%         w = array of NURBS weights (multi-indexed)
%
% Purpose: Compute the control net for a pipe elbow.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p_1,p_2,n_1,n_2,Xi_1,Xi_2,P,W] = NURBS_Pipe_Elbow_Surface(h,w,r,R)

Xi_1 = [0 0 0 .25 .25 .5 .5 .75 .75 1 1 1];
Xi_2 = [0 0 0 1/3 1/3 2/3 2/3 1 1 1];

p_1 = 2;
p_2 = 2;

n_1 = length(Xi_1)-p_1-1;
n_2 = length(Xi_2)-p_2-1;


P(:,1,:) = [[r,0;r,r;0,r;-r,r;-r,0;-r,-r;0,-r;r,-r;r,0],zeros(n_1,1)]+...
    [ones(n_1,1)*(-R), zeros(n_1,1), ones(n_1,1)*h];
W(:,1) = [1, sqrt(2)/2 1 sqrt(2)/2, 1 sqrt(2)/2, 1 sqrt(2)/2 1]';
P(:,2,:) = [[r,0;r,r;0,r;-r,r;-r,0;-r,-r;0,-r;r,-r;r,0],zeros(n_1,1)]+...
    [ones(n_1,1)*-R,zeros(n_1,1),ones(n_1,1)*h/2];
W(:,2) = [1, sqrt(2)/2 1 sqrt(2)/2, 1 sqrt(2)/2, 1 sqrt(2)/2 1]';
P(:,3,:) = [[r,0;r,r;0,r;-r,r;-r,0;-r,-r;0,-r;r,-r;r,0],zeros(n_1,1)]+...
    [ones(n_1,1)*-R,zeros(n_1,1),ones(n_1,1)*0];
W(:,3) = [1, sqrt(2)/2 1 sqrt(2)/2, 1 sqrt(2)/2, 1 sqrt(2)/2 1]';


P(:,4,:) = [-R+r, 0, -R+r;...
    -R+r, r, -R+r;...
    -(R), r, -(R);...
    -R-r, r, -R-r;...
    -R-r, 0, -R-r;...
    -R-r, -r, -R-r;...
    -(R), -r, -(R);...
    -R+r, -r, -R+r;...
    -R+r, 0, -R+r];


W(:,4) = [sqrt(2)/2 1/2 sqrt(2)/2 1/2 sqrt(2)/2 1/2 sqrt(2)/2 1/2 sqrt(2)/2]';

P(:,5,:) = [zeros(n_1,1), [0,r;r,r;r,0;r,-r;0,-r;-r,-r;-r,0;-r,r;0,r]]+...
    [ones(n_1,1)*0,zeros(n_1,1),ones(n_1,1)*-R];

W(:,5) = [1, sqrt(2)/2 1 sqrt(2)/2, 1 sqrt(2)/2, 1 sqrt(2)/2 1]';

P(:,6,:) = [zeros(n_1,1), [0,r;r,r;r,0;r,-r;0,-r;-r,-r;-r,0;-r,r;0,r]]+...
    [ones(n_1,1)*w/2,zeros(n_1,1),ones(n_1,1)*-R];

W(:,6) = [1, sqrt(2)/2 1 sqrt(2)/2, 1 sqrt(2)/2, 1 sqrt(2)/2 1]';

P(:,7,:) = [zeros(n_1,1), [0,r;r,r;r,0;r,-r;0,-r;-r,-r;-r,0;-r,r;0,r]]+...
    [ones(n_1,1)*w,zeros(n_1,1),ones(n_1,1)*-R];

W(:,7) = [1, sqrt(2)/2 1 sqrt(2)/2, 1 sqrt(2)/2, 1 sqrt(2)/2 1]';