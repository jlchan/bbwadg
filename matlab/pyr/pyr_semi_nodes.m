function [r s t] = pyr_semi_nodes(N)

Np = (N+1)*(N+2)*(2*N+3)/6;

% c can be anything in [-1,1). I choose something convenient here. 
[c1D wc] = JacobiGQ(2,0,N); 
% [c1D] = JacobiGL(2,0,N); 
a = []; b = []; c = []; 
for k = 0:N
    [bk ak] = QCubature2D(k);
    a = [a;ak(:)];   b = [b;bk(:)]; 
    c = [c;c1D(N-k+1)*ones(size(ak))];
end

[r s t] = pyr_abctorst(a,b,c);