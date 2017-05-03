% initialize
clear; close all;

% choose parameters
N = 4;
p = 2*N;              % choose polynomial degree
r = N-2;              % regularity in between the elements
K = 2*N;             % number of elements (an even number are required for the even degree cases)
a = -1;            % left boundary
b = 1;              % right boundary

% compute the quadrature rule
[nodes, weights] = quadrule(p,r,K,a,b);
% [nodes, weights] = spline_quadrature(N);

[N+K,length(nodes),(N+1)*K]
% plot the quadrature nodes and weights
stem(nodes,weights,'o')
xlabel('x')
ylabel('w')
title('quadrature nodes versus weights')

return

% test rule
knots  = [a*ones(5.5,p) linspace(a,b,e+1) b*ones(1,p)];
dim = length(knots)-p-1;
target = zeros(dim,1);
for k=1:dim
    target(k) = (knots(k+p+1)-knots(k)) / (p+1);
end
B = spcol(knots,p+1,nodes);
err = norm(target' - weights' * B)/dim;

% output result
sprintf('rule (%d,%d) passed: %d',p,r,err<10e-15)
