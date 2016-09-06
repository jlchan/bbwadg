% function test_conservation_limit

% function [u] = Advec1D(u, FinalTime)
% Purpose  : Integrate 1D advection until FinalTime starting with
%            initial the condition, u

Globals1D;

% Order of polymomials used for approximation
N = 5;

% Generate simple mesh
K1D = 1;
[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

% Initialize solver and construct grid and metric
StartUp1D;

vmapP(1) = vmapM(end); % make periodic
vmapP(end) = vmapM(1);

rp = linspace(-1,1,100)';
Vp = bern_basis_1D(N,rp);
xp = Vandermonde1D(N,rp)/V*x;

global VB DB VBe xe TB
re = linspace(-1,1,N+1)';
[VB VBr] = bern_basis_1D(N,r);
DB = VB\VBr;
DB(abs(DB)<1e-8) = 0;
VBe = bern_basis_1D(N,re);
xe = Vandermonde1D(N,re)/V * x;

TB = VB\Vandermonde1D(N,r);


% Set initial conditions
d = 1000;
%     uex = @(x) -1./(1 + exp(-d*(x-1/3))) + 1./(1 + exp(-d*(x+1/3)));%x > -1/3 & x < 1/3;%exp(-25*x.^2);
%     uex = @(x) (x > -3/4 & x < -1/4) + 0*abs(x-.4).*(x > 0);
uex = @(x) x > .1;
%     uex = @(x) sin(pi*x);

u = uex(x);


uB = VB\u;
uavg = repmat(sum(uB,1)/(N+1),N+1,1);

u1 = TB\uB; % modal
u1(3:end,:) = 0;
u1 = TB*u1; % convert back to BB
u2 = uB - u1;
uB2 = uex(xe) - u1;

% % plot(xp,(Vp*(u1)),'-');
% hold on
% usmooth = uex(xe);%.2*(uB-u1) + u1;
%
% plot(xp,uex(xp),'--')
% % uB = uex(xe);
% plot(xp,Vp*uB,'-')
% plot(xp,Vp*usmooth,'-')
% plot(xe,uB,'o--')
% plot(xe,usmooth,'s--')

[rq wq] = JacobiGQ(0,0,N+10);
Vq = bern_basis_1D(N,rq);
xq = Vandermonde1D(N,rq)/V * x;
invM = inv(Vq'*diag(wq)*Vq);
L = invM(:,1);

% ||dudx||_L1
TV1 = sum(sum(diag(wq)*abs(Vq*DB*uB)*diag(J(1,:))));

% discrete TV estimate
TV = 0;
for i = 1:N
    TV = TV + abs(uB(i+1,:) - uB(i,:));
end
TV2 = J(1)*sum(TV);

% semilogy(1:N,TV1,'o--')
% hold on;
% semilogy(1:N,TV2,'x-')
% semilogy(1:N,exp(.62*(1:N)),'s--')

% plot(xp,abs(uex(xp)-Vp*VB*uex(xe)),'.-');
% hold on
% plot(xp,abs(Vp*VB*VBe*uex(xe) - 0*uex(xp)),'--');

% axis([-1 1 -2 2])
% axis([-.55 -.2 -2 2])

v = uB(1:3);
m = size(v,1); 
mfunc = zeros(1,size(v,2));
s = sum(sign(v),1)/m;
ids = find(abs(s)==1);
if(~isempty(ids))
    mfunc(ids) = s(ids).*min(abs(v(:,ids)),[],1);
end
