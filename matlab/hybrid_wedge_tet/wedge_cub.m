function [r s t w] = wedge_cub(N)

hybridgGlobalFlags

% triangle nodes
[rt st w2D] = Cubature2D(2*N+1); 

% tensor product with 1D SEM nodes
if useSEM
    [t1D] = JacobiGL(0,0,N); V = Vandermonde1D(N,t1D); w1D = sum(inv(V*V'),2);
else
    [t1D w1D] = JacobiGQ(0,0,N);
end
[te re] = meshgrid(t1D,rt); [~,se] = meshgrid(t1D,st);
r = re(:); s = se(:); t = te(:);

[w2D w1D] = meshgrid(w1D,w2D); 
% [w1D w2D] = meshgrid(w1D,w2D); 
w = w2D(:).*w1D(:);
w = 4*w/sum(w);

% [r s wrt] = Cubature2D(2*N+1); % for integrating O(2*N + J = 2*N+1?)
% 
% [t1D ws] = JacobiGQ(0,0,N);
% 
% % % get GLL weights by lumping
% [s r] = meshgrid(t1D,r);
% [~, t] = meshgrid(t1D,s);
% r = r(:); s = s(:); t = t(:);
% [wrt ws] = meshgrid(ws,wrt);
% w = wrt(:).*ws(:);
% 
% w = 4*w/sum(w);
