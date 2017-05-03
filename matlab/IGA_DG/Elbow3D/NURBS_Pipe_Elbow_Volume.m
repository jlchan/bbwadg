%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: NURBS_Pipe_Elbow_Volume
%
% Input:  h = height
%         w = width
%         r = pipe radius
%         R = elbow turning radius
%         t = thickness
%
% Output: p_1 = polynomial degree in direction 1
%         p_2 = polynomial degree in direction 2
%         p_3 = polynomial degree in direction 3
%         n_1 = number of functions in direction 1
%         n_2 = number of functions in direction 2
%         n_3 = number of functions in direction 3
%         Xi_1 = knot vector in direction 1
%         Xi_2 = knot vector in direction 2
%         Xi_3 = knot vector in direction 3
%         P = array of NURBS control points (single-indexed)
%         w = array of NURBS weights (single-indexed)
%
% Purpose: Compute the control lattice for a pipe elbow.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p_1,p_2,p_3,n_1,n_2,n_3,Xi_1,Xi_2,Xi_3,P,W] = NURBS_Pipe_Elbow_Volume(h,w,r,R,t)

[p_1,p_2,n_1,n_2,Xi_1,Xi_2,P_1,W_1] = NURBS_Pipe_Elbow_Surface(h,w,r,R);
[~,~,~,~,~,~,P_2,W_2] = NURBS_Pipe_Elbow_Surface(h,w,r+t,R);
[~,~,~,~,~,~,P_3,W_3] = NURBS_Pipe_Elbow_Surface(h,w,r-t,R);

P = zeros(9,7,3,3);

P(:,:,1,1) = P_1(:,:,1);
P(:,:,2,1) = P_2(:,:,1);
P(:,:,3,1) = P_3(:,:,1);

P(:,:,1,2) = P_1(:,:,2);
P(:,:,2,2) = P_2(:,:,2);
P(:,:,3,2) = P_3(:,:,2);

P(:,:,1,3) = P_1(:,:,3);
P(:,:,2,3) = P_2(:,:,3);
P(:,:,3,3) = P_3(:,:,3);

W(:,:,1) = W_1(:,:);
W(:,:,2) = W_2(:,:);
W(:,:,3) = W_3(:,:);

P1 = reshape(permute(P(:,:,:,1),[3,2,1]),[189,1]);
P2 = reshape(permute(P(:,:,:,2),[3,2,1]),[189,1]);
P3 = reshape(permute(P(:,:,:,3),[3,2,1]),[189,1]);

P = [P1 P2 P3];
W = reshape(permute(W,[3,2,1]),[189,1]);

Xi_3 = [0 0 0 1 1 1];
p_3 = 2;
n_3 = 3;