function burgersSystem1D

% function [u] = Advec1D(u, FinalTime)
% Purpose  : Integrate 1D advection until FinalTime starting with
%            initial the condition, u

Globals1D;

% Order of polymomials used for approximation
N = 5;
K1D = 16;
FinalTime = .1;

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

% Initialize solver and construct grid and metric
StartUp1D;

vmapP(1) = vmapM(end); % make periodic
vmapP(end) = vmapM(1);

rp = linspace(-1,1,25)';
Vp = Vandermonde1D(N,rp)/V;
xp=  Vp*x;

% Set initial conditions
global  Vq Pq
 
Nq = 2*N+2;
Nq = N;
[rq wq] = JacobiGQ(0,0,Nq);
Vq = Vandermonde1D(N,rq)/V;
Pq = (Vq'*diag(wq)*Vq)\(Vq'*diag(wq));
xq = Vq*x;
wJq = diag(wq)*(Vq*J);

u1ex = @(x) cos(pi/2*x);
u2ex = @(x) sin(pi*x);
% a = .3;
% u1ex = @(x) 1+exp(-10*(x-a).^2);
% u2ex = @(x) 1+exp(-10*(x+a).^2);
u1 = Pq*u1ex(xq);
u2 = Pq*u2ex(xq);

time = 0;

% Runge-Kutta residual storage
res1 = zeros(Np,K);
res2 = zeros(Np,K);

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));
dt   = .1*xmin;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps;

M = inv(V*V');
% LIFT = diag(1./sum(M,2))*M*LIFT; % SEM


% outer time step loop

figure(1)
for tstep=1:Nsteps
    
    % low storage RK
    for INTRK = 1:5
        [rhs1, rhs2] = RHS1D(u1, u2);
        
        res1 = rk4a(INTRK)*res1 + dt*rhs1;
        res2 = rk4a(INTRK)*res2 + dt*rhs2;
        
        u1 = u1 + rk4b(INTRK)*res1;
        u2 = u2 + rk4b(INTRK)*res2;        
    end;
    
    % Increment time
    time = time+dt;
    
    if (1 && mod(tstep,5)==0 ) || (tstep==Nsteps)
        
        clf
        hold on;
        plot(xp,Vp*u1,'-');     plot(x,u1,'o')
        plot(xp,Vp*u2,'-');     plot(x,u2,'^')
        hold off
        
        axis([-1 1 -3 3])
        title(sprintf('Time = %f\n',time))
        drawnow
    end
    
    u1q = Vq*u1;
    u2q = Vq*u2;
    enorm(tstep) = sum(sum(wJq.*(u1q.^2 + u2q.^2)));
    
    if mod(tstep,1000)==0
        disp(sprintf('tstep = %d out of %d\n',tstep,Nsteps))
    end
end;

figure(2);
plot(dt*(1:Nsteps),enorm)
hold on


function [rhs1, rhs2] = RHS1D(u1,u2)

Globals1D;
global Vq Pq

jump = @(u) reshape(u(vmapP)-u(vmapM),Nfp*Nfaces,K);
Dh = @(u) rx.*(Dr*u) + .5*LIFT*(Fscale.*jump(u).*nx);

opt=3;
if opt==1
    u1q = Vq*u1;
    u2q = Vq*u2;
    f1 = Pq*(u1q.^2 + u2q.^2);
    f2 = Pq*(2*u1q.*u2q);
    rhs1 = -Dh(f1);
    rhs2 = -Dh(f2);
    
    c = max(u1 + u2, u1 - u2);
    c = abs(c);
    cf = reshape(max(c(vmapM),c(vmapP)),Nfp*Nfaces,K);
    rhs1 = rhs1 + .5*LIFT*(Fscale.*cf.*jump(u1));
    rhs2 = rhs2 + .5*LIFT*(Fscale.*cf.*jump(u2));
    
elseif opt==2
    
    du1 = Dh(u1);
    du2 = Dh(u2);
    u1q = Vq*u1;
    u2q = Vq*u2;
    f1 = Pq*(u1q.^2 + u2q.^2);
    f2 = Pq*(2*u1q.*u2q);
    uDu1 = Pq*(u1q.*(Vq*du1) + u2q.*(Vq*du2));
    uDu2 = Pq*(u2q.*(Vq*du1) + u1q.*(Vq*du2));
    rhs1 = (2/3)*(Dh(f1) + uDu1);
    rhs2 = (2/3)*(Dh(f2) + uDu2);
    
    rhs1 = -rhs1;
    rhs2 = -rhs2;
    
    c = max(u1 + u2, u1 - u2);
    c = abs(c);
    cf = reshape(max(c(vmapM),c(vmapP)),Nfp*Nfaces,K);
    rhs1 = rhs1 + .5*LIFT*(Fscale.*cf.*jump(u1));
    rhs2 = rhs2 + .5*LIFT*(Fscale.*cf.*jump(u2));
    
%     rhs1 = Pq*((Vq*exp(-u1.^2)).*(Vq*rhs1));
%     rhs2 = Pq*((Vq*exp(-u2.^2)).*(Vq*rhs2));

elseif opt==3 % integrated splitting
    
    dfdU11 = @(u1,u2) 2*u1; 
    dfdU12 = @(u1,u2) 2*u2; 
    dfdU22 = @(u1,u2) 2*u1;    
    
    u1q = Vq*u1;
    u2q = Vq*u2;    
    
    Nq = 2;
    [tq wtq] = JacobiGQ(0,0,Nq); 
    tq = (1+tq)/2; wtq = wtq/2;
    
    F11 = 0; F12 = 0; F22 = 0;
    for i = 1:length(tq)
        F11 = F11 + dfdU11(u1q*tq(i),u2q*tq(i))*tq(i)*wtq(i);
        F12 = F12 + dfdU12(u1q*tq(i),u2q*tq(i))*tq(i)*wtq(i);
        F22 = F22 + dfdU22(u1q*tq(i),u2q*tq(i))*tq(i)*wtq(i);
    end
    
    A1 = Dh(Pq*(F11.*u1q + F12.*u2q));
    A2 = Dh(Pq*(F12.*u1q + F22.*u2q));
    
    dU1 = Vq*Dh(u1);
    dU2 = Vq*Dh(u2);
    B1 = Pq*(F11.*dU1 + F12.*dU2);
    B2 = Pq*(F12.*dU1 + F22.*dU2);
    
    % split forms
    rhs1 = -(A1 + B1);
    rhs2 = -(A2 + B2);
    
    c = max(u1 + u2, u1 - u2);
    c = abs(c);
    cf = reshape(max(c(vmapM),c(vmapP)),Nfp*Nfaces,K);
    rhs1 = rhs1 + .5*LIFT*(Fscale.*cf.*jump(u1));
    rhs2 = rhs2 + .5*LIFT*(Fscale.*cf.*jump(u2));
    
end

return


