function [x y drdx dsdx drdy dsdy J] = Elbow2D(x_in,y_in)

addpath('/Users/jchan985/Desktop/bbwadg/matlab/IGA_DG/Elbow2D')

if nargin==0
    x_in = linspace(-1,1,8);
    y_in = linspace(-1,1,8);
end

% convert to [0,1]
if nargin==1
    y_in = x_in;
end
x_in = (1+x_in)/2;
y_in = (1+y_in)/2;

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

% % testing
% P(:,1) = [1 0 -1 1 0 -1 1 0 -1]';
% a = 2; P(:,2) = [-a -a -a 0 0 0 a a a]' - P(:,1);
% a = 2; P(:,2) = [-a -a -a 0 0 0 a a a]';
% w = ones(size(w));

p_1 = 2; p_2 = 2;
n_1 = 3; n_2 = 3;
Xi_1 = [0,0,0,1,1,1];
Xi_2 = [0,0,0,1,1,1];

%%%
% Number of quadrature points
n_q = 3;

n = n_1*n_2;
dim = 2;
[n_el,C_operators,IEN,P_b,w_b,~,w_e] = Extract_And_Localize(n_1,n_2,p_1,p_2,Xi_1,Xi_2,P,w);
[ID] = Construct_ID(dim,n);

[xi_q,w_q] = Quadrature_Data(n_q);

Np1 = length(x_in(:));
Np2 = length(y_in(:));

x = zeros(Np1*Np2,n_el);
y = zeros(Np1*Np2,n_el);
drdx = zeros(Np1*Np2,n_el);
drdy = zeros(Np1*Np2,n_el);
dsdx = zeros(Np1*Np2,n_el);
dsdy = zeros(Np1*Np2,n_el);

J = zeros(Np1*Np2,n_el);

for e = 1:n_el
    C_e = C_operators(:,:,e);
    sk = 1;
    for q_1 = 1:Np1
        for q_2 = 1:Np2
            [R,dRdx,xij,G] = Shape_Function(x_in(q_1),y_in(q_2),p_1,p_2,C_e,P_b(:,:,e),w_b(:,e),w_e(:,e));
            
            G = fliplr(inv(G')*2); % correction from John - map to [-1,1] and ordering i guess?
            drdx(sk,e) = G(1,1); dsdx(sk,e) = G(1,2);
            drdy(sk,e) = G(2,1); dsdy(sk,e) = G(2,2);

            J(sk,e) = 1/abs(det(G));
            x(sk,e) = xij(1);
            y(sk,e) = xij(2);
            sk = sk + 1;
        end
    end
    if mod(e,10)==0
        disp(sprintf('on element %d out of %d\n',e,n_el))
    end
end
%     plot(x,y,'o')
%     axis equal
return


