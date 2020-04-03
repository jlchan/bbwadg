function [H,D,tL,tR,x]=CSBPpEQ2(N)
%==========================================================================
% Purpose: This function computes the classical FD SBP operator for p = 2
%          of degree = 2 for the first derivative with N nodes on
%          the domain [-1,1]
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

x = transpose(linspace(-1,1,N));
dx = 2/(N-1);

%==== norm matrix
H = dx*eye(N,N);
H(1:4,1:4)=dx*[0.17e2 / 0.48e2 0 0 0; 0 0.59e2 / 0.48e2 0 0; 0 0 0.43e2 / 0.48e2 0; 0 0 0 0.49e2 / 0.48e2;];
% lower poriton
for i=1:4
    for j=1:4
        H(N-(i-1),N-(j-1))=H(i,j);
    end
end
%==== first derivative
D = zeros(N,N);
D(1:4,1:6) = ...
1/dx*[-0.24e2 / 0.17e2 0.59e2 / 0.34e2 -0.4e1 / 0.17e2 -0.3e1 / 0.34e2 0 0; -0.1e1 / 0.2e1 0 0.1e1 / 0.2e1 0 0 0; 0.4e1 / 0.43e2 -0.59e2 / 0.86e2 0 0.59e2 / 0.86e2 -0.4e1 / 0.43e2 0; 0.3e1 / 0.98e2 0 -0.59e2 / 0.98e2 0 0.32e2 / 0.49e2 -0.4e1 / 0.49e2;];

% interior of D1
for i=5:N-4
    D(i,i-2:i+2)=1/dx*[1/12,-2/3,0,2/3,-1/12];
end
% the lower portion of D1
for i=1:6
    for j=1:6
        D(N-(i-1),N-(j-1))=-D(i,j);
    end
end
%===== matrices necessary for the SATs
tL = zeros(N,1);
tL(1,1)=1;
tR = zeros(N,1);
tR(N,1)=1;
