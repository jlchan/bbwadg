function An = corner(A,n)

if nargin==1
    n = 5;
end

An = A(1:n,1:n);
    