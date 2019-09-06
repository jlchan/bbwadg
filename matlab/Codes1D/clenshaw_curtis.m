function [r w] = clenshaw_curtis(N,kind)

if nargin==1
    kind=1;
end

j = (0:N)';
if kind==1
    r = cos((2*j+1)*pi/(2*(N+1))); % no boundary points
elseif kind == 2
    r = cos(pi*j/N); % boundary points    
else
    r = cos((2*j)*pi/(2*(N+1))); % nonsymmetric version
end
r = sort(r);

% compute weights in a dumb way
V = Vandermonde1D(N,r);
[rq wq] = JacobiGQ(0,0,N);
Vq = Vandermonde1D(N,rq)/V;
w = Vq'*wq;