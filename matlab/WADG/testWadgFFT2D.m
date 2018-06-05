% function testWadgFFT2D
global N Dx Dy
N = 1000;

x = linspace(-1 + 2/N,1,N);
[x y] = meshgrid(x,x);
% f = 1.0*(max(abs(x),abs(y))<.5); 
f = sin(x*pi).*sin(y*pi);
% mesh(x,y,f);return

k = 4;
a = 1 + exp(sin(k*x*pi).*sin(k*y*pi)*5);

% c = 100;
% a = c*(1 + 1/c + sin(k*pi*x).*sin(k*pi*y));

a = 2*(2 + sin(k*pi*x).*sin(k*pi*y));
a = 1 + 10.0*(max(abs(x+.25),abs(y+.25))<.5) + 20.0*(max(abs(x-.25),abs(y-.25))<.5);
% mesh(x,y,a);return

ff = fft2(f); ff(1) = 0; f = ifft2(ff);
iD = 1i*pi*[0:N/2 (-N/2+1):-1];
[Dx Dy] = meshgrid(iD);

u0 = applyWadg(f,a);
maxit = 50;
% [u,FLAG,relres,iter,resvec] = pcg(@(u) applyWK(u,a), f(:), 1e-10, maxit);
% [u,FLAG,relres,iter,resvec] = gmres(@(u) applyWK(u,a), f(:), 5, 1e-10, maxit, @(f) applyWadg(f,a));
[u,FLAG,relres,iter,resvec] = pcg(@(u) applyWK(u,a), f(:), 1e-8, maxit, @(f) applyWadg(f,a),[],u0);
[u2,FLAG,relres2,iter2,resvec2] = pcg(@(u) applyWK(u,a), f(:), 1e-8, maxit, @(f) applyInvK(f,a),[],u0);

figure(1)
clf
mesh(x,y,reshape(real(u),N,N))
title(sprintf('iter = %d, residual = %g',iter,relres))

figure(2)
semilogy(resvec,'o--')
hold on
semilogy(resvec2,'x--')
% U = zeros(N);
% A = zeros(N^2);
% P = zeros(N^2);
% for i = 1:N^2
%     U(i) = 1;
%     A(:,i) = applyWK(U,a);
%     P(:,i) = applyWadg(U,a);
%     U(i) = 0;
% end
% lam = eig(P*A);
% figure
% plot(lam,'o')


function f = applyWK(u,a)

global N Dx Dy

u = reshape(u(:),N,N);
uf = fft2(u);
uf(1) = 0;

% differentiate in x,y
dux = Dx.*uf;
duy = Dy.*uf;

% transform back to nodal, scale by a
uax = fft2(ifft2(dux).*a);
uay = fft2(ifft2(duy).*a);

% divergence of result
d2ua = conj(Dx).*uax + conj(Dy).*uay;
% d2ua(1) = 0;

f = ifft2(d2ua);
f = f(:);
end

function u = applyWadg(f,a)

global N Dx Dy

f = reshape(f(:),N,N);
ff = fft2(f);

% invert laplacian and differentiate
uf = ff./(conj(Dx).*Dx+conj(Dy).*Dy); uf(1) = 0; 
ufx = Dx.*uf; 
ufy = Dy.*uf;

% transform back to nodal, scale
uax = fft2(ifft2(ufx)./a);
uay = fft2(ifft2(ufy)./a);

% take divergence and invert laplacian
d2ua = (conj(Dx).*uax + conj(Dy).*uay)./(conj(Dx).*Dx+conj(Dy).*Dy);
d2ua(1) = 0;

u = ifft2(d2ua);
u = u(:);
end

function u = applyInvK(f,a)

global N Dx Dy

f = reshape(f(:),N,N);
ff = fft2(f);

% invert laplacian and differentiate
uf = ff./(conj(Dx).*Dx+conj(Dy).*Dy); uf(1) = 0; 

u = ifft2(uf);
u = u(:);
end