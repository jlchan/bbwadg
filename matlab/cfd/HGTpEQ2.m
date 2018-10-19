function [H,D,tL,tR,x]=HGTpEQ2(N)
%==========================================================================
% Purpose: This function computes the operator for
%          p = 2 of degree = 2 for the first derivative with N nodes on
%          the domain [-1,1] for the hybrid Gauss trapezoidal nodal
%          distribution with a=j=2
% Inputs:
%          N = number of nodes, the minimum is 10
%	   construction of the operators
%
% Outputs:
%           H = SBP norm
%           D = first derivative operator
%           tL = projection operator for the left boundary
%           tR = projection operator for the right boundary
%           x = nodal distribution for the operator
%
% Author: David C. Del Rey Fernandez January 2016
%
%==========================================================================
x = transpose(linspace(1,N,N));
x(1,1) = 761/470-(17/470)*sqrt(119);
x(2,1) = 761/470+(17/470)*sqrt(119);

x(N,1) = N-291/470+(17/470)*sqrt(119);
x(N-1,1) = N-291/470-(17/470)*sqrt(119);

g = (N-1)/2;
theta = 1+g;
x = x/g-theta/g;
dx = 2/(N-1);

tL = zeros(N,1);tR = zeros(N,1);

tL(1,1) = 0.738402e6 / 0.2132393e7 + 0.481613142e9 / 0.4313831039e10 * sqrt(0.119e3);
tL(2,1) = -0.481613142e9 / 0.4313831039e10 * sqrt(0.119e3) + 0.738402e6 / 0.2132393e7;
tL(3,1) =  0.321e3 / 0.823e3;
tL(4,1) = -0.214e3 / 0.2591e4;

tR(N,1) = tL(1,1);
tR(N-1,1) = tL(2,1);
tR(N-2,1) = tL(3,1);
tR(N-3,1) = tL(4,1);

%=== non-unity weights of the norm matrix
h11 = 0.19025345e8 / 0.25588716e8 - 0.71518740e8 / 0.4313831039e10 * sqrt(0.119e3);
h22 = 0.19025345e8 / 0.25588716e8 + 0.71518740e8 / 0.4313831039e10 * sqrt(0.119e3);
h33 = 0.2516e4 / 0.2469e4;
h44 = 0.7726e4 / 0.7773e4;

H = eye(N,N);
H(1,1) = h11;H(2,2) = h22;H(3,3) = h33;H(4,4) = h44;
H(N,N) = h11;H(N-1,N-1) = h22;H(N-2,N-2) = h33;H(N-3,N-3) = h44;
H = dx*H;

%== construct the derivative operator Note: for practical implementation
% one would not first construct S and then construct D1 = H^-1*Q
% instead, this multiplication is performed once and for all
% Here, I am seperating them so that those who wish to derive these 
% operators of their own accord can see the values of the various matrices

q11 = -(0.738402e6 / 0.2132393e7 + 0.481613142e9 / 0.4313831039e10 * sqrt(0.119e3)) ^ 2 / 0.2e1;
q12 = 0.2582155325e10 / 0.25882986234e11 * sqrt(0.119e3) - (-0.481613142e9 / 0.4313831039e10 * sqrt(0.119e3) + 0.738402e6 / 0.2132393e7) * (0.738402e6 / 0.2132393e7 + 0.481613142e9 / 0.4313831039e10 * sqrt(0.119e3)) / 0.2e1; 
q13 = 0.483617271e9 / 0.3509918878e10 - 0.62859015829e11 / 0.835360692964e12 * sqrt(0.119e3); 
q14 = -0.197542311e9 / 0.11050060526e11 + 0.111600501617e12 / 0.7889743215564e13 * sqrt(0.119e3);  

q21 = -0.2582155325e10 / 0.25882986234e11 * sqrt(0.119e3) - (-0.481613142e9 / 0.4313831039e10 * sqrt(0.119e3) + 0.738402e6 / 0.2132393e7) * (0.738402e6 / 0.2132393e7 + 0.481613142e9 / 0.4313831039e10 * sqrt(0.119e3)) / 0.2e1;
q22 = -(-0.481613142e9 / 0.4313831039e10 * sqrt(0.119e3) + 0.738402e6 / 0.2132393e7) ^ 2 / 0.2e1;
q23 = 0.62859015829e11 / 0.835360692964e12 * sqrt(0.119e3) + 0.483617271e9 / 0.3509918878e10; 
q24 = -0.111600501617e12 / 0.7889743215564e13 * sqrt(0.119e3) - 0.197542311e9 / 0.11050060526e11; 
 
q31 = -0.957671355e9 / 0.3509918878e10 + 0.450211994765e12 / 0.14201131780388e14 * sqrt(0.119e3);  
q32 = -0.450211994765e12 / 0.14201131780388e14 * sqrt(0.119e3) - 0.957671355e9 / 0.3509918878e10;
q33 = -0.103041e6 / 0.1354658e7;
q34 = 0.18042395e8 / 0.25588716e8; 
q35 = -0.1e1 / 0.12e2; 

q41 = 0.513578367e9 / 0.11050060526e11 - 0.660425978833e12 / 0.134125634664588e15 * sqrt(0.119e3); 
q42 = 0.660425978833e12 / 0.134125634664588e15 * sqrt(0.119e3) + 0.513578367e9 / 0.11050060526e11;
q43  = -0.17218067e8 / 0.25588716e8;
q44 = -0.22898e5 / 0.6713281e7;
q45 = 0.2e1 / 0.3e1; 
q46 = -0.1e1 / 0.12e2; 

Q = zeros(N,N);
Q(1,1) = q11; Q(1,2) = q12;Q(1,3) = q13;Q(1,4) = q14;
Q(2,1) = q21;Q(2,2) = q22;Q(2,3) = q23;Q(2,4) = q24;
Q(3,1) = q31;Q(3,2) = q32;Q(3,3) = q33;Q(3,4) = q34;Q(3,5) = q35;
Q(4,1) = q41;Q(4,2) = q42;Q(4,3) = q43;Q(4,4) = q44;Q(4,5) = q45;Q(4,6) = q46;

%-- populate the interior point operator
for j = 5:N-4
    Q(j,j-2:j+2) = [1/12,-2/3,0,2/3,-1/12];
end

%-- construct the lower portion of Q
for j = 1:4
    for k = 1:6
        Q(N-(j-1),N-(k-1)) = -Q(j,k);
    end
end

D = H\Q;
