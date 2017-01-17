function [r s t w] = wedge_cubature(N)

% cheaper for N=1
if N==1
    [r s t w] = wedge_cubature_N3(N);
    return
end

[r t wrt] = Cubature2D(2*N+1); % for integrating O(2*N + J = 2*N+1?)

[s1D ws] = JacobiGQ(0,0,N);

% % get GLL weights by lumping
[s1D ws] = JacobiGQ(0,0,N);


[s r] = meshgrid(s1D,r);
[~, t] = meshgrid(s1D,t);
r = r(:); s = s(:); t = t(:);
[wrt ws] = meshgrid(ws,wrt);
w = wrt(:).*ws(:);

w = 4*w/sum(w);
