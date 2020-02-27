% clear

N = 3;
K = 64;

tau = 1;
CFL = .5;
FinalTime = 1;

% plotting points
Npts = 500; a = -K/2; b = -a;
rp = linspace(a,b,Npts)';

[Vp VX] = GDVDM(N,K,rp);
VX = VX';
i = K/2+1; 

% VX = 2*VX/K;
h = VX(2)-VX(1);

% piecewise quadrature
rq = []; wq = [];
for e = 1:K    
    [rqe wqe] = JacobiGQ(0,0,N);
    D1D = GradVandermonde1D(N,rqe)/Vandermonde1D(N,rqe);
    rq = [rq; h*(rqe+1)/2 + VX(e)];
    wq = [wq; h/2*wqe];    
end

Vq = GDVDM(N,K,rq);
M = Vq'*diag(wq)*Vq;
VN = GDVDM(N,K,rq);
DN = kron(eye(K),D1D);
rx = 2;
Q = VN'*diag(wq)*(rx*DN*VN); % Galerkin first deriv mat
Q(abs(Q)<1e-8) = 0;
D = M\Q;

% % possibly reduced quadrature versions
% [H, ~, X] = CSBPp4(K+1);
% L = max(VX)-min(VX);
% wq = diag(H)*L/2; rq = X*L/2;
% Vq = GDVDM(N,K,rq);
% M = (Vq'*diag(wq)*Vq);
% Q = M*D;

h = 2/K;
VX = VX*h;
rp = rp*h;

%%

% scaling
M = h*M;

B = zeros(K+1,2);
B(1,1) = 1;
B(K+1,2) = 1;

D = M\Q;
L = M\B;
nx = [-1;1];

M1 = M;
m1 = M1(:,1);
M1(1,:) = 0;
M1(:,1) = 0;
M1(1,1) = 1;

[nm,~] = size(M1);
% M2 = M;
% pskip = N+2;
% inds = pskip:nm-(pskip)+1;
inds = 1:nm-(N+1); % only inflow
% inds = N+2:nm; % only outflow
MDiag = diag(sum(M1,2));
% M1(inds,:) = MDiag(inds,:);

Q1 = Q;
q1 = Q(:,1);
Q1(1,:) = 0;
Q1(:,1) = 0;
Q1(1,1) = 1;
D1 = M1\Q1;

dt = CFL*h;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

u0 = @(x) sin(pi*x); %exp(-5^2*x.^2);
u = u0(VX);

% % plot(rp,Vp*u,'o--'); return
% for t = 0:.01:2
%     plot(rp,sin(pi*(rp-t)));
%     ylim([-1,1])
%     drawnow
% end

%% rk coeffs

rk4a = [            0.0 ...
        -567301805773.0/1357537059087.0 ...
        -2404267990393.0/2016746695238.0 ...
        -3550918686646.0/2091501179385.0  ...
        -1275806237668.0/842570457699.0];
rk4b = [ 1432997174477.0/9575080441755.0 ...
         5161836677717.0/13612068292357.0 ...
         1720146321549.0/2090206949498.0  ...
         3134564353537.0/4481467310338.0  ...
         2277821191437.0/14882151754819.0];
rk4c = [             0.0  ...
         1432997174477.0/9575080441755.0 ...
         2526269341429.0/6820363962896.0 ...
         2006345519317.0/3224310063776.0 ...
         2802321613138.0/2924317926251.0];     

res = 0*u;
for i = 1:Nsteps
    
    for INTRK = 1:5               
                
        t = (i-1)*dt + rk4c(INTRK)*dt;
        
        u1 = sin(pi*(VX(1)-t));
%         du1 = -pi*cos(pi*(VX(1)-t));
% 
%         u(1) = u1;
%         rhs = D1*u + M1\(m1*du1 + u1*q1);
%         rhs(1) = du1;
        
        uM = u([1;K+1]);        
        uP = u([K+1;1]);        
        uP(1) = u1;
        uP(end) = uM(end);
        f = .5*(uP-uM).*nx - .5*tau*(uP-uM);
        rhs = (D*u) + .5*L*f;        
        
        rhs = -rhs;           
        res = rk4a(INTRK)*res + dt*rhs;
        u   = u + rk4b(INTRK)*res;
%         u(1) = u1;
        
    end
    
    if mod(i,10)==0 || i==Nsteps
        plot(rp,Vp*u,'o--')
%         hold on
%         plot(-1,sin(pi*(VX(1)-t)),'x')
%         hold off
        title(sprintf('time = %f\n',i*dt))
        drawnow
    end        
end

title(sprintf('err = %g',norm(u - sin(pi*(VX-FinalTime)),'inf')))
