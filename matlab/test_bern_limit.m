clear

Nvec = 2:15;

rp = linspace(-1,1,150)';

for i = 1:length(Nvec)
    N = Nvec(i);
    
    r = JacobiGL(0,0,N);
    
    V = bern_basis_1D(N,r);
    re = linspace(-1,1,N+1)';
    Ve = bern_basis_1D(N,re);
    
    Vp = bern_basis_1D(N,rp);
    
    a = 1e1;
    f = @(x) 1./(1+exp(-a*(x-.25)));
%     f = @(x) sin(x);
%     clf
%     plot(rp,f(rp))
%     hold on
%     plot(rp,Vp*f(re),'--')
%     pause
    err1(i) = max(abs(f(re) - Ve*f(re))) / max(abs(f(re)));    
    err2(i) = max(abs(f(re) - Ve*(V\f(r))))/ max(abs(f(re)));    
end

semilogy(Nvec,err1,'o--')
hold on
% semilogy(Nvec,err2,'x--')
% semilogy(Nvec,5.^(-Nvec),'.-')
semilogy(Nvec,5e-2./Nvec,'--')
semilogy(Nvec,.5./sqrt(Nvec),'--')

%% test choice of limiter switch

N = 7;

r = JacobiGL(0,0,N);

V = bern_basis_1D(N,r);
re = linspace(-1,1,N+1)';
Ve = bern_basis_1D(N,re);

rp = linspace(-1,1,150)';
Vp = bern_basis_1D(N,rp);

f = @(x) 1./(1+exp(-10*(x-.1)));
% f = @(x) sin(pi*x);

u = f(re);
err = max(abs(f(re) - Ve*u)) / max(abs(u));

a = min(1,(sqrt(N)*err).^1); % if f smooth, err = O(1/N). if rough, err = O(1/sqrt(N))

% mix interpolant/projection w/BB approximant
u = (1-a)*(V\f(r)) + f(r)*a;
clf
plot(rp,f(rp),'-')
hold on
plot(re,Ve*u,'o')
plot(rp,Vp*u,'--')
drawnow


%% test iterated bernstein approximant - ech, it converges to the equispaced interpolant

N = 10;

r = JacobiGL(0,0,N);

V = bern_basis_1D(N,r);
re = linspace(-1,1,N+1)';
Ve = bern_basis_1D(N,re);

rp = linspace(-1,1,150)';
Vp = bern_basis_1D(N,rp);

f = @(x) 1./(1+exp(-100*(x-.1)));
% f = @(x) sin(pi*x);

err = max(abs(f(r) - Ve*f(re)));

u = f(r);
for i = 1:0
    err = Ve*u - f(re);
    u = u - err;
end
plot(rp,f(rp),'-')
% plot(rp,Vp*(Ve\f(re)),'.-')
hold on
plot(re,Ve*u,'o')
plot(rp,Vp*u,'--')

    
