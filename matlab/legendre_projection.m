clear

N = 5;
a = 1e4;
f = @(x) 1./(1+a*x.^2);

% projection using Legendre polynomials: use legendreP in Matlab
b = zeros(N+1,1);
for i = 0:N
    b(i+1) = integral(@(x) legendreP(i,x) .* f(x),-1,1);
    M(i+1,i+1) = integral(@(x) legendreP(i,x).*legendreP(i,x),-1,1);
end
cm = b./diag(M);

% plot result
x = linspace(-1,1,100);
uL = 0;
for i = 0:N
    uL = uL + cm(i+1)*legendreP(i,x);
end
plot(x,f(x))
hold on
plot(x,uL,'--')