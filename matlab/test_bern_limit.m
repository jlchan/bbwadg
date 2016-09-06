clear

Nvec = 3:10;

rp = linspace(-1,1,150)';

for i = 1:length(Nvec)
    N = Nvec(i);
    
    r = JacobiGL(0,0,N);
    [rq wq] = JacobiGQ(0,0,N);
    Vq = bern_basis_1D(N,rq);
    Pq = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq));
    
    V = bern_basis_1D(N,r);
    re = linspace(-1,1,N+1)';
    Ve = bern_basis_1D(N,re);
    
    Vp = bern_basis_1D(N,rp);
    
    a = 1;
    f = @(x) 1./(1+exp(-a*(x-.25)));
    %     f = @(x) sin(x);
    %     clf
    %     plot(rp,f(rp))
    %     hold on
    %     plot(rp,Vp*f(re),'--')
    %     pause
    err1(i) = max(abs(f(re) - Ve*f(re))) ;    % BB approx
    err2(i) = max(abs(f(re) - Ve*(Pq*f(rq))));
end

semilogy(Nvec,err1,'o--')
hold on
% semilogy(Nvec,err2,'x--')
% semilogy(Nvec,err2,'x--')
% semilogy(Nvec,5.^(-Nvec),'.-')
C1 = 2.5e-2;
C2 = 1.3e-2;
semilogy(Nvec,C1./Nvec,'--')
semilogy(Nvec,C2./sqrt(Nvec),'--')

%% test choice of limiter switch

N = 7;

r = JacobiGL(0,0,N);
[rq wq] = JacobiGQ(0,0,N+4);
Vq = bern_basis_1D(N,rq);
V = bern_basis_1D(N,r);
re = linspace(-1,1,N+1)';
Ve = bern_basis_1D(N,re);

rp = linspace(-1,1,150)';
Vp = bern_basis_1D(N,rp);

errs = [];
sk = 1;
avec = 50:-.1:.5;
avec = 25;
for a = avec
    f = @(x) 1./(1+exp(-a*(x-.1)));
    % f = @(x) exp(-a*(x-.25).^2);
    %     f = @(x) sin(pi/2*x);
    
    %     uB = f(re);
    uB = Ve*(V\f(r)); % BB approx of interp - with Gibbs
    
    % N2 = N+1;
    % re2 = linspace(-1,1,N2+1)';
    % Ve2 = bern_basis_1D(N,re2);
    % E = bern_basis_1D(N2,JacobiGL(0,0,N2))\bern_basis_1D(N,JacobiGL(0,0,N2));
    % uB = E\(Ve2*(V\f(r))); % BB approx on aux space, interp back
    
    % spatial scaling
    %     b = .25*(1-re).*(1+re);
    b = 1;
    err = max(abs(uB - Ve*uB).*b) / max(abs(uB));  % pointwise diff
    alpha = (sqrt(N)*err)^(1); % if f smooth, err = O(1/N). if rougher, err = O(1/sqrt(N)).
    
    % mix interpolant/projection w/BB approximant
    u = (1-alpha)*(V\f(r)) + uB*alpha;
    
    e1 = sqrt(wq'*(f(rq)-Vq*(V\f(r))).^2);
    e2 = sqrt(wq'*(f(rq)-Vq*u).^2);
    e3 = sqrt(wq'*(f(rq)-Vq*f(re)).^2);
    
    errs(sk,:) = [e1 e2 e3];
    sk = sk + 1;
    
    if length(avec)==1
        clf
        plot(rp,f(rp),'-')
        hold on
        plot(rp,Vp*(V\f(r)),'-') % interp
        plot(re,Ve*u,'o') % mixed BB + interp
        plot(re,Ve*uB,'x') % pure bb limiter
        plot(rp,Vp*u,'--')
        plot(rp,Vp*uB,'--') % pure bb limiter
        title(sprintf('L2 error: interp = %f, limited = %f, BB = %f\n',e1,e2,e3));
        legend('exact','interp','limit','BB')
    end
end
if length(avec)>1
    semilogy(errs)
    legend('interp','limit','BB')
end

%%

N = 9;

r = JacobiGL(0,0,N);
[rq wq] = JacobiGQ(0,0,N+4);
Vq = bern_basis_1D(N,rq);
V = bern_basis_1D(N,r);
re = linspace(-1,1,N+1)';
Ve = bern_basis_1D(N,re);

rp = linspace(-1,1,150)';
Vp = bern_basis_1D(N,rp);

a = 100;
f = @(x) 1./(1+exp(-a*(x-.1)));
u = V\f(r); % interp coeffs
uB = Ve*(V\f(r)); % BB approx of interp with Gibbs

TV = 0;
for i = 1:N
    TV = TV + abs(u(i,:) - u(i+1,:));
end
TV = TV/(N*max(abs(u)));

% spatial scaling
%     b = .25*(1-re).*(1+re);
b = 1;
err = max(abs(uB - Ve*uB).*b);  % pointwise diff
alpha = err/sqrt(N); % if f smooth, err = O(1/N). if rougher, err = O(1/sqrt(N)).
alpha = TV;

% mix interpolant/projection w/BB approximant
u = (1-alpha).*(V\f(r)) + uB.*alpha;

plot(rp,f(rp),'-')
hold on
plot(rp,Vp*(V\f(r)),'-') % interp

plot(re,uB,'s') % control points
plot(re,Ve*u,'o') % mixed BB + interp
plot(rp,Vp*u,'--')
% plot(re,Ve*uB,'x') % pure bb limiter
% plot(rp,Vp*uB,'--') % pure bb limiter
title(sprintf('TV = %f\n',TV))
% title(sprintf('L2 error: interp = %f, limited = %f, BB = %f\n',e1,e2,e3));
legend('exact','control pts','interp','limit')

%% pics!

for N = 9;
    
    r = JacobiGL(0,0,N);
    [rq wq] = JacobiGQ(0,0,N+4);
    Vq = bern_basis_1D(N,rq);
    V = bern_basis_1D(N,r);
    re = linspace(-1,1,N+1)';
    Ve = bern_basis_1D(N,re);
    
    rp = linspace(-1,1,150)';
    Vp = bern_basis_1D(N,rp);
    
    a = 1000;
    f = @(x) 1./(1+exp(-a*x));
    uB = f(re);
    plot(rp,f(rp),'--')
    hold on
    plot(re,uB,'o')
    plot(rp,Vp*uB,'-')
%     pause
end

%% 2D

N = 15;

[r s] = Nodes2D(N); [r s] = xytors(r,s);

[rp sp] = EquiNodes2D(75); [rp sp] = xytors(rp,sp);
[re se] = EquiNodes2D(N); [re se] = xytors(re,se);


V = bern_basis_tri(N,r,s);
Vp = bern_basis_tri(N,rp,sp);

uex = @(x,y) x > -.25 + .5*y;
u = uex(re,se);
% u = V\uex(r,s);



vv = Vp*u;
h = color_line3(rp,sp,vv,vv,'.');
set(h,'markersize',24)

hold on
plot3(re,se,u,'o');
