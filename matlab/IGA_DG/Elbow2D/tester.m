% function [

r_i = .5;
t = .5;
r_o = r_i+t;
r_m = 0.5*(r_i+r_o);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Original Bezier Control Points and Weights
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P_original(:,1) = [-r_i;-r_m;-r_o;-r_i;-r_m;-r_o;0;0;0];
P_original(:,2) = [0;0;0;r_i;r_m;r_o;r_i;r_m;r_o];
w_original = [1;1;1;cos(pi/4);cos(pi/4);cos(pi/4);1;1;1];

%%%
% Geometric parameters for first level

P = P_original;
w = w_original;
% P(:,1) = [-1 0 1 -1 0 1 -1 0 1]';
% P(:,2) = [0 0 0 .5 .5 .5 1 1 1]';
% w = ones(size(w));

p_1 = 2;
p_2 = 2;
n_1 = 3;
n_2 = 3;
Xi_1 = [0,0,0,1,1,1];
Xi_2 = [0,0,0,1,1,1];

%%%
% Number of quadrature points

n_q = 3;

if 1
    n = n_1*n_2;
    dim = 2;
    ndof = dim*n;
    [n_el,C_operators,IEN,P_b,w_b,~,w_e] = Extract_And_Localize(n_1,n_2,p_1,p_2,Xi_1,Xi_2,P,w);
    [ID] = Construct_ID(dim,n);
    
    [xi_q,w_q] = Quadrature_Data(n_q);
    
    Nplot = 25;
    xp = linspace(0,1,Nplot);
    
    x = zeros(n_q^2,n_el);
    y = zeros(n_q^2,n_el);
    for e = 1:n_el
        C_e = C_operators(:,:,e);
        sk = 1;
        for q_1 = 1:Nplot
            for q_2 = 1:Nplot
                %[R,dRdx,xij,J] = Shape_Function(xi_q(q_1),xi_q(q_2),p_1,p_2,C_e,P_b,w_b,w_e);
                [R,dRdx,xij,G,wb] = Shape_Function(xp(q_1),xp(q_2),p_1,p_2,C_e,P_b(:,:,e),w_b(:,e),w_e(:,e));                
                
                J(sk,e) = det(G);
                x(sk,e) = xij(1);
                y(sk,e) = xij(2);
                ww(sk,e) = wb;
                sk = sk + 1;
            end
        end
    end
    vv = ww.*y;
    h = color_line3(x,y,vv,vv,'.'); set(h,'markersize',32);
    axis equal
    hold on
    return
end

