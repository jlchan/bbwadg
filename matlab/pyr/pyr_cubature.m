function [r s t w] = biunit_pyr_cubature(N,alpha)

if nargin==1    
    alpha = 2;
end
[a1D w1D] = JacobiGQ(0,0,N);
[c1D wc1D] = JacobiGQ(alpha,0,N);
% [c1D wc1D] = jagsrd(N+1,alpha,0);

[a c b] = meshgrid(a1D, c1D, a1D);
[wa wc wb] = meshgrid(w1D,wc1D,w1D);

% [a b c] = meshgrid(a1D,a1D,c1D);
% [wa wb wc] = meshgrid(w1D,w1D,wc1D);

a = a(:); b = b(:); c = c(:);
wa = wa(:); wb = wb(:); wc = wc(:);
w = wa.*wb.*wc;
w = (8/3)*w./sum(w); % scale by b*h/3 volume of pyr.  b = 4, h = 2

% convert to abc
[r s t] = pyr_abctorst(a,b,c);