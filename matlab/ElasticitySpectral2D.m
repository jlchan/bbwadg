% single-domain spectral method - used to check
% function [p invVquad] = spectral_wave(Nin,cfun,FinalTime)
% Nin = order
% cfun = @(x) ... (wavespeed)
% FinalTime = time to run until
% p = pressure solution
% invVquad = inverse VDM used to transform the solution for plotting

function [U invVquad] = ElasticitySpectral2D(Nin,FinalTime)

% Low storage Runge-Kutta coefficients
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
     
global N nx ny mu lambda vmapX vmapY D Fmask LIFT M invM Nfld E
Nfld = 5;
     
if nargin==0
    N = 41;
    FinalTime = .5;
else
    N = Nin;
end

r = JacobiGL(0,0,N);
V = Vandermonde1D(N,r);
D = GradVandermonde1D(N,r)/V;
if nargin==0
    rp = linspace(-1,1,150);
    [xp yp] = meshgrid(rp);   
    Vp1D = Vandermonde1D(N,rp)/V;    
end
[y x] = meshgrid(r);
    
Fmask = [find(abs(y(:)+1) < 1e-8), find(abs(x(:)-1)<1e-8), find(abs(y(:)-1)<1e-8), find(abs(x(:)+1)<1e-8)];

nx = zeros(N+1,4);
nx(:,1) = 0;  nx(:,2) = 1;
nx(:,3) = 0;  nx(:,4) = -1;

ny = zeros(N+1,4);
ny(:,1) = -1; ny(:,2) = 0;
ny(:,3) = 1;  ny(:,4) = 0;

nx = nx(:); ny = ny(:);

E = zeros(4*(N+1),(N+1)^2);
for i = 1:4*(N+1)
    E(i,Fmask(i)) = 1;
end
% M = diag(sum(inv(V*V'),2));
% invM = inv(M);
M = inv(V*V');
invM = V*V';

%%

if 0
    x = x(:); y = y(:);
    
    global Nfp Np Nfaces rx ry sx sy K vmapP vmapM Dr Ds mapB vmapB Fscale
    
    mapB = 1:4*(N+1);
    vmapB = Fmask(:);
    Dr = kron(eye(N+1),D);
    Ds = kron(D,eye(N+1));
    
    Fscale = ones(size(Fmask(:)));
    vmapP = Fmask(:); vmapM = Fmask(:);
    K = 1;
    Nfp = N+1; Nfaces = 4; Np = (N+1)^2;
    rx = ones(size(x));
    sy = ones(size(x));
    ry = zeros(size(x));
    sx = zeros(size(x));
    
    LIFT = kron(invM,invM) * E'* blkdiag(M,M,M,M);
    flux = rand(size(LIFT,2),1);
    val1 = LIFT*flux;
    tmp = M*reshape(flux,N+1,4);
    tmp = reshape(E'*tmp(:),N+1,N+1);
    val2 = invM*tmp*invM;
end



%%

a = .25; b = 1;
mu = 1 + a*cos(b*pi*x).*cos(b*pi*y);
lambda = 1 + a*sin(b*pi*x).*sin(b*pi*y);

mu0 = 1;
w = sqrt(2*mu0);
u1 = @(x,y,t) cos(w*pi*t).*cos(pi*x).*sin(pi*y);
u2 = @(x,y,t) -cos(w*pi*t).*sin(pi*x).*cos(pi*y);

u = zeros(N+1);
U{1} = u1(x,y,0);
U{2} = u2(x,y,0);
U{3} = u;
U{4} = u;
U{5} = u;


%%

time = 0;
for fld = 1:Nfld
    res{fld} = zeros(N+1);
%     res{fld} = zeros(Np,K);
end

dt = 1/(N+1)^2; % compute time step size

% outer time step loop
tstep = 0;

while (time<FinalTime)
    if(time+dt>FinalTime), dt = FinalTime-time; end
    
    for INTRK = 1:5
        
        timelocal = time + rk4c(INTRK)*dt;        
        [rhs] = elas_rhs(U,timelocal);
        
        for fld = 1:Nfld 
            res{fld} = rk4a(INTRK)*res{fld} + dt*rhs{fld};
            U{fld} = U{fld}+rk4b(INTRK)*res{fld};
        end        
        
    end;
    
    if nargin==0 && mod(tstep,100)==0
%         keyboard
        disp(sprintf('on time %f\n',time))
        clf
        vv = U{3};
        color_line3(x,y,vv,vv,'.'); % shading interp;
        view(3)
        axis equal
        xlabel('x')
        ylabel('y')
        axis([-1 1 -1 1 -2 2])
        %         axis off
        drawnow
    end
    
    if mod(tstep,1000)==0
        disp(sprintf('on timestep %d out of %d\n',tstep,ceil(FinalTime/dt)))
    end
    
    % Increment time
    time = time+dt; tstep = tstep+1;
end

invVquad = kron(inv(V),inv(V));

if nargin==0
    clf
    vv = U{3};
    color_line3(x,y,vv,vv,'.'); % shading interp;
    view(3)
    axis tight
    axis([-1 1 -1 1 -2 2])
end

% surf(xp,yp,pp)
% shading interp

function [rhs] = elas_rhs(U,time)

global N nx ny mu lambda vmapX vmapY D Fmask LIFT M invM Nfld
global Nfp Np Nfaces rx ry sx sy K vmapP vmapM Dr Ds mapB vmapB Fscale  E

% % Define field differences at faces
% for fld = 1:Nfld
%     u = U{fld};
%     
%     % compute jumps
%     dU{fld} = zeros(Nfp*Nfaces,K);
%     dU{fld}(:) = u(vmapP)-u(vmapM);
%     
% %     ur = Dr*u;
% %     us = Ds*u;
% %     Ux{fld} = rx.*ur + sx.*us;
% %     Uy{fld} = ry.*ur + sy.*us;
% end
% 
% divSx = zeros(Np,K); 
% divSy = zeros(Np,K);
% du1dx = zeros(Np,K);
% du2dy = zeros(Np,K);
% du12dxy = zeros(Np,K);
% 
% divSx(:) = (D*reshape(U{3},N+1,N+1) + reshape(U{5},N+1,N+1)*D');
% divSy(:) = (D*reshape(U{5},N+1,N+1) + reshape(U{4},N+1,N+1)*D');
% du1dx(:) = D*reshape(U{1},N+1,N+1);
% du2dy(:) = reshape(U{2},N+1,N+1)*D';
% du12dxy(:) = (D*reshape(U{2},N+1,N+1) + reshape(U{1},N+1,N+1)*D');
% 
% % divSx = Ux{3} + Uy{5}; % d(Sxx)dx + d(Sxy)dy
% % divSy = Ux{5} + Uy{4}; % d(Sxy)dx + d(Syy)dy
% % du1dx = Ux{1}; % du1dx
% % du2dy = Uy{2}; % du2dy
% % du12dxy = Ux{2} + Uy{1}; % du2dx + du1dy
% 
% for fld = 1:Nfld
%     dU{fld} = zeros(size(nx));
% end
% nSx = nx.*dU{3} + ny.*dU{5};
% nSy = nx.*dU{5} + ny.*dU{4};
%  
% opt=3;
% if opt==1 % traction BCs    
%     nSx(mapB) = -2*(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
%     nSy(mapB) = -2*(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));
%     %     end    
% elseif opt==2 % basic ABCs
%     nSx(mapB) = -(nx(mapB).*U{3}(vmapB) + ny(mapB).*U{5}(vmapB));
%     nSy(mapB) = -(nx(mapB).*U{5}(vmapB) + ny(mapB).*U{4}(vmapB));
%     dU{1}(mapB) = -U{1}(vmapB);
%     dU{2}(mapB) = -U{2}(vmapB);    
% elseif opt==3 % zero velocity
%     dU{1}(mapB) = -2*U{1}(vmapB);
%     dU{2}(mapB) = -2*U{2}(vmapB);    
%     
% end
% 
% % stress fluxes
% nUx = dU{1}.*nx;
% nUy = dU{2}.*ny;
% nUxy = dU{2}.*nx + dU{1}.*ny;
% 
% % evaluate central fluxes
% fc{1} = nSx;
% fc{2} = nSy;
% fc{3} = nUx;
% fc{4} = nUy;
% fc{5} = nUxy;
% 
% % penalization terms - reapply An
% fp{1} = nx.*fc{3} + ny.*fc{5};
% fp{2} = nx.*fc{5} + ny.*fc{4};
% fp{3} = fc{1}.*nx;
% fp{4} = fc{2}.*ny;
% fp{5} = fc{2}.*nx + fc{1}.*ny;
% 
% flux = cell(5,1);
% for fld = 1:Nfld   
%     flux{fld} = zeros(Nfp*Nfaces,K);
%     %flux{fld}(:) = fc{fld}(:) + tau{fld}(vmapM).*fp{fld}(:);
%     flux{fld}(:) = fc{fld}(:) + fp{fld}(:);
% end
% 
% % compute right hand sides of the PDE's
% rr{1} =  divSx   +  LIFT*(flux{1})/2.0;
% rr{2} =  divSy   +  LIFT*(flux{2})/2.0;
% rr{3} =  du1dx   +  LIFT*(flux{3})/2.0;
% rr{4} =  du2dy   +  LIFT*(flux{4})/2.0;
% rr{5} =  du12dxy +  LIFT*(flux{5})/2.0;
% 
% rhs{1} = rr{1};
% rhs{2} = rr{2};
% rhs{3} = (2*mu+lambda).*rr{3} + lambda.*rr{4};
% rhs{4} = lambda.*rr{3} + (2*mu+lambda).*rr{4};
% rhs{5} = (mu) .* rr{5};
% 
% 
% return;
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% return

rhs{1} = (D*U{3} + U{5}*D');
rhs{2} = (D*U{5} + U{4}*D');
U1x = D*U{1};
U2y = U{2}*D';
u12xy = (D*U{2} + U{1}*D');

rhs{3} = U1x;
rhs{4} = U2y;
rhs{5} = u12xy;

% rhs{3} = (2*mu+lambda).*U1x + lambda.*U2y;
% rhs{4} = lambda.*U1x + (2*mu+lambda).*U2y;
% rhs{5} = mu.*u12xy;

% for fld = 1:5
%     rhs{fld}(Fmask(:,3)) = rhs{fld}(Fmask(:,1));
%     rhs{fld}(Fmask(:,4)) = rhs{fld}(Fmask(:,2));
% end
% return

opt=2;
if opt==1
    % need to enforce weak velocity BCs
    dU{1} = -2*U{1}(Fmask(:));
    dU{2} = -2*U{2}(Fmask(:));
    nSx = zeros((N+1)*4,1);
    nSy = zeros((N+1)*4,1);
elseif opt==2
    nSx = -2*(nx.*U{3}(Fmask(:)) + ny.*U{5}(Fmask(:)));
    nSy = -2*(nx.*U{5}(Fmask(:)) + ny.*U{4}(Fmask(:)));
    dU{1} = zeros((N+1)*4,1);
    dU{2} = zeros((N+1)*4,1);
 
end
nUx = dU{1}.*nx;
nUy = dU{2}.*ny;
nUxy = dU{2}.*nx + dU{1}.*ny;

fc{1} = nSx;
fc{2} = nSy;
fc{3} = nUx;
fc{4} = nUy;
fc{5} = nUxy;

% penalization terms - reapply An
fp{1} = nx.*fc{3} + ny.*fc{5};
fp{2} = nx.*fc{5} + ny.*fc{4};
fp{3} = fc{1}.*nx;
fp{4} = fc{2}.*ny;
fp{5} = fc{2}.*nx + fc{1}.*ny;

flux = cell(5,1);
for fld = 1:5      
    flux{fld} = reshape(fc{fld} + fp{fld},N+1,4);
end

ftmp = zeros(N+1);
for fld = 1:5    
    ftmp = M*flux{fld};
    ftmp = reshape(E'*ftmp(:),N+1,N+1);
    rhs{fld} = rhs{fld} + invM*ftmp*invM'; 
end

U1x = rhs{3};
U2y = rhs{4};
u12xy = rhs{5};

rhs{3} = (2*mu+lambda).*U1x + lambda.*U2y;
rhs{4} = lambda.*U1x + (2*mu+lambda).*U2y;
rhs{5} = mu.*u12xy;
