function [Vq x xq] = SincBasis(N,Nq)

dx = 2*pi/N;
tol = 1e-12;
sinc = @(x) (sin(pi*x/dx)+tol)./(2*pi/dx*tan(x/2)+tol);

x = linspace(-pi,pi,N+1)';
x = x(1:N);
xq = linspace(-pi,pi,Nq+1); %ceil(3/2*N))';
xq = xq(1:end-1);
xp = linspace(-pi,pi,1000)';
xp = xp(1:end-1);

Vp = zeros(length(xp),N);
Vq = zeros(length(xq),N);
V = zeros(N);

for i = 1:N
    Vp(:,i) = sinc(xp-x(i));
    Vq(:,i) = sinc(xq-x(i));
    V(:,i) = sinc(x-x(i));
end
% norm(Vq*sin(x)-sin(xq)) % low frequency mode 
% norm(Vq*sin((N/2-1)*x)-sin((N/2-1)*xq)) % highest mode representable
% 
% plot(xp,Vp)