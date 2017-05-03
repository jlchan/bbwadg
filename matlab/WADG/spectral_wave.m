% single-domain spectral method - used to check
% function [p invVquad] = spectral_wave(Nin,cfun,FinalTime)
% Nin = order
% cfun = @(x) ... (wavespeed)
% FinalTime = time to run until
% p = pressure solution
% invVquad = inverse VDM used to transform the solution for plotting

function [p invVquad] = spectral_wave(Nin,cfun,FinalTime)

Globals2D

if nargin==0
    N = 50;
    %cfun = @(x,y) ones(size(x));
    cfun = @(x,y) 1 + .5*sin(pi*x).*sin(pi*y);
    FinalTime = 1.0;
else
    N = Nin;    
end

r = JacobiGL(0,0,N);
V = Vandermonde1D(N,r);
Dr = GradVandermonde1D(N,r)/V;

if nargin==0
    rp = linspace(-1,1,150);
    [xp yp] = meshgrid(rp);   
    Vp1D = Vandermonde1D(N,rp)/V;
    
end
[x y] = meshgrid(r);
Np = length(x); K = 1;

global c vmapB Dr
vmapB = find(abs(abs(x)-1)<1e-8 | abs(abs(y)-1)<1e-8);
% plot(x,y,'o'); hold on;plot(x(vmapB),y(vmapB),'*');return

c = cfun(x,y);

%% eigs
Np = (N+1)^2;
if 0 && nargin==0 
    e = zeros(3*Np,1);
    A = zeros(3*Np);
    for i = 1:3*Np
        e(i) = 1;
        ids = 1:Np;        
        p = reshape(e(ids),N+1,N+1);
        u = reshape(e(ids + Np),N+1,N+1);
        v = reshape(e(ids + 2*Np),N+1,N+1);
        [rhsp, rhsu, rhsv] = spectral_rhs(p,u,v);
        A(:,i) = [rhsp(:);rhsu(:);rhsv(:)];
        e(i) = 0;
        if mod(i,100)==0
            disp(sprintf('on column %d out of %d\n',i,3*Np))
        end
    end
    lam = eig(A);
    hold on;
    plot(lam,'o')
    title(sprintf('Largest real part = %e',max(real(lam))))    
    keyboard
    
%     load AN3K4; norm(A-blkdiag(blkdiag(Pc{:}),eye(2*Np*K))*Aconst,'fro') % should be 0
    return
end

%%

% initial condition
% p = exp(-100*(x.^2 + (y-.25).^2));
k = 1; W = (2*k-1)/2*pi; p = cos(W*x).*cos(W*y); % frequency of solution
u = zeros(N+1,N+1); v = zeros(N+1);

time = 0;
resp = zeros(N+1); resu = zeros(N+1); resv = zeros(N+1);  % Runge-Kutta residual storage

dt = .1/max(c(:)*(N+1)^2); % compute time step size

% outer time step loop
tstep = 0;

while (time<FinalTime)
    if(time+dt>FinalTime), dt = FinalTime-time; end
    
    for INTRK = 1:5
        
        timelocal = time + rk4c(INTRK)*dt;
        
        [rhsp, rhsu, rhsv] = spectral_rhs(p,u,v,timelocal);
        
        % initiate and increment Runge-Kutta residuals
        resp = rk4a(INTRK)*resp + dt*rhsp;
        resu = rk4a(INTRK)*resu + dt*rhsu;
        resv = rk4a(INTRK)*resv + dt*rhsv;
        
        % update fields
        u = u+rk4b(INTRK)*resu;
        v = v+rk4b(INTRK)*resv;
        p = p+rk4b(INTRK)*resp;
        
    end;
    
    if 0 && mod(tstep,100)==0
        disp(sprintf('on time %f\n',time))
        pcolor(x,y,p); shading interp;
        axis off
        axis equal
        drawnow;       
    end
    
    % Increment time
    time = time+dt; tstep = tstep+1;
end

invVquad = kron(inv(V),inv(V));
% if nargin==0
%     Vp = Vp*invVquad;
%     Vp = kron(Vp1D/,Vp1D/V);
%     pp = Vp*p(:);
%     Vp = Vandermonde2DQuad(N,xp,yp);
%     
%     plot3(xp(:),yp(:),pp(:),'.')
% end

% surf(xp,yp,pp)
% shading interp

function [rhsp rhsu rhsv] = spectral_rhs(p,u,v,time)

global c vmapB Dr

rhsp = -c.*(Dr*u + v*Dr'); % div U
rhsu = -Dr*p;
rhsv = -p*Dr';

rhsp(vmapB) = 0;



