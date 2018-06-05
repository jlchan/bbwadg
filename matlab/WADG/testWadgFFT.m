% function testWadgFFT
% Fourier series

clear

global N iD
N = 10000;
iD = 1i*pi*[0:N/2 (-N/2+1):-1];

x = linspace(-1+2/N,1,N);

k = 4;
% a = 2+sin(pi*x); 
a = 1 + 10*(abs(x+.25)<.25) + 100*(abs(x-.25)<.5);
a = 0 + exp(sin(k*pi*x)*10);
% plot(x,a);return
x = x(:);

f = sin(pi*x); %1.0*(abs(x)<.5);
ff = fft(f); ff(1) = 0; f = ifft(ff);

u = applyWadg(f,a);
f2 = applyWK(u,a);

maxit = 5;
[u,FLAG,relres,iter,resvec] = gmres(@(u) applyWK(u,a), f, maxit, 1e-11, maxit, @(f) applyWadg(f,a));

% plot(x,f)
% hold on
% plot(x,f2,'x')

figure(1)
plot(x,u,'.-')
hold on
title(sprintf('iter = %d, %d, residual = %g',iter(1),iter(2),relres))

function f = applyWK(u,a)

global N iD

u = u(:)';
uf = fft(u);

% differentiate
duf = iD.*uf;

% transform back to nodal, scale by a
dua = fft(ifft(duf).*a);

% differentiate again
d2ua = conj(iD).*dua;
d2ua(1) = 0;

f = ifft(d2ua);
f = f(:);

end

function u = applyWadg(f,a)
global N iD

f = f(:)';
f_f = fft(f);

% invert laplacian and differentiate
uf = iD.*f_f./(conj(iD).*iD); uf(1) = 0;

% transform back to nodal
ua = fft(ifft(uf)./a);

% differentiate and invert laplacian
uf = conj(iD).*ua./(conj(iD).*iD);
uf(1) = 0;

u = ifft(uf);
u = u(:);
end

