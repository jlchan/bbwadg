clear
N = 6;
rb = JacobiGL(0,0,N+1); 
rb = rb(2:end-1);
e = ones(size(rb));
x = [-1;  1; -1; rb;  -e; -rb];
y = [-1; -1; 1; -e; rb; rb];
% plot(x,y,'o')
% return
[x y] = rstoxy(x,y);
a = .95;
r = []; 
s = [];
for i = 1:N-1
    [ri si] = xytors(a*x,a*y);
    r = [r; ri];
    s = [s; si];
    a = a * a;
end

% rb = JacobiGL(0,0,N); 
% rb = rb(2:N);
% e = ones(size(rb));
% x = [-1;  1; -1; rb; rb; -e];
% y = [-1; -1; 1; -e; rb; rb];
% [x y] = rstoxy(x,y);
% a = .95;
% for i = 1:N
%     [ri si] = xytors(a*x,a*y);
%     r = [r; ri];
%     s = [s; si];
%     a = a*.75;
% end

plot(r,s,'o')

N = 2*N;
V = Vandermonde2D(N,r,s);
[size(V,2) rank(V),cond(V) ]