function [H,D,tL,tR,x]=CSBPpEQ1(N)
%==========================================================================
% Purpose: This function computes the classical FD SBP operator for p = 1
%          of degree = 1 for the first derivative with N nodes on
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
H(1,1) = dx/2;
H(N,N) = dx/2;
%==== first derivative
D = zeros(N,N);
D(1,1:2) = 1/dx*[-1,1];

D(N,N-1:N) = 1/dx*[-1,1];

% interior of D
for i=2:N-1
    D(i,i-1:i+1)=1/dx*[-1/2,0,1/2];
end

%===== matricies necessary for the SATs
tL = zeros(N,1);
tL(1,1)=1;
tR = zeros(N,1);
tR(N,1)=1;