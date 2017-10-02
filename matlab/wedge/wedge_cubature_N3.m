function [r s t w] = wedge_cubature_N3(N)

[rgq w] = JacobiGQ(0,0,N);
[rgq1 w1] = JacobiGQ(1,0,N);

% [a b c] = meshgrid(rgq, rgq, rgq1);
% [wa wb wc] = meshgrid(w,w,w1);

% j first for TP convenience
[a c b] = meshgrid(rgq, rgq1, rgq);
[wa wc wb] = meshgrid(w,w1,w);
a = a(:); b = b(:); c = c(:);
wa = wa(:); wb = wb(:); wc = wc(:);
w = wa.*wb.*wc;

[r s t] = wedge_abctorst(a,b,c);

% [r t wrt] = Cubature2D(2*N); % for integrating O(2*N + J = 2*N+1?)
% [s1D ws] = JacobiGQ(0,0,N-1);
% 
% [r s] = meshgrid(r,s1D);
% [t ~] = meshgrid(t,s1D);
% r = r(:); s = s(:); t = t(:);
% [wrt ws] = meshgrid(wrt,ws);
% w = wrt(:).*ws(:);

w = 4*w/sum(w);
