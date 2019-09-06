clear
N = 5;
xq = linspace(-1,1,2*N+2)'; % interpolatory
% xq = linspace(-pi,pi,4*N+2)'; % interpolatory
xp = linspace(-1,1,1000)';
xq = xq(1:end-1);
dx = xq(2)-xq(1);

sk = 1;
D = zeros(2*N+1);
for k = -N:N
    Vq(:,sk) = exp(1i*pi*k*xq)/sqrt(2);
    Vp(:,sk) = exp(1i*pi*k*xp)/sqrt(2);    
    D(sk,sk) = 1i*pi*k;
    sk = sk + 1;
end
M = dx*(Vq'*Vq);
% plot(x,Vq)

Pq = M\(Vq'*dx);

% f = @(x) exp(sin(pi*x));
% df = @(x) exp(sin(x)).*cos(x);
% a = .5;
% f = @(x) (1-x).*(x > a) + -(x+1).*(x < a);
% df = @(x) -1.*(x>a) + -1.*(x<a);
% plot(xp,df(xp))
% hold on
% plot(xp,Vp*D*Pq*f(xq),'r-')
% plot(xq,Vq*D*Pq*f(xq),'o')

Dq = real(Vq*D*Pq); % imag part = 0
% figure
% imagesc(Dq)
% Dq = Dq*sqrt(2*N+1)
