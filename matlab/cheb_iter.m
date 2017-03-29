% chebyshev iter
function [x rvec rhist] = cheb_iter(A,b,x,Lmax,Lmin,tol,maxit)

if nargin < 6
    maxit = 100;
    tol = 1e-7;
end
p = zeros(size(b));

d = (Lmax+Lmin)/2;
c = (Lmax-Lmin)/2;
beta = 0; alpha = 1.0; 

rhist = b-A*x;
for i = 1:maxit
    r = b-A*x;
    p = r + beta*p;    
    
    alpha = 1/(d-beta/alpha);
    x = x + alpha*p;            
    beta = (c*alpha/2)^2;
    
    rvec(i) = norm(r);    
  
    if nargout == 3
        rhist = [rhist x];
    end
end