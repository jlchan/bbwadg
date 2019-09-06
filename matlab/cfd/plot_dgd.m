clear
N = 2;

dx = 1/(N+1);
r = linspace(-1+dx,1-dx,N+1)';
rdual = linspace(-1,1,N+2)';
rp = linspace(-1,1,1000)';

V = Vandermonde1D(N,r);
Vp = Vandermonde1D(N,rp)/V;
rint = linspace(-dx,dx,25)';
Vint = Vandermonde1D(N,rint)/V;

N0 = ceil((N+1)/2);
e = zeros(N+1,1);
hold on
plot(r,r*0,'o')
plot(rdual,rdual*0,'x')
% ylim([-1.5,1.5])
% plot(rint,Vint*e,'o')

iH = @(r,id) r >= rdual(id) - 1e-8;
flag = 0;
for id = 1:N+1
    flag = flag + iH(rp,id);
end

Dr = V\GradVandermonde1D(N,r);
v = zeros(size(rp));
e = zeros(N+1,1);
for i = 1:length(rp)
    id = flag(i);
    e(id) = 1;        
    ri = r(id) - rp(i);
    Vpi = Vandermonde1D(N,ri)*inv(V);
    v(i) = Vpi * e;
    e(id) = 0;
end

plot(rp,v,'.')
% plot(rp,flag,'.')


