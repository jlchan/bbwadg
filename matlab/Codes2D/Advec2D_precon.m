function Advec2D

Globals2D

% Polynomial order used for approximation 
N = 7;

% Read in Mesh
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Grid/Other/squarereg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Grid/Other/squareireg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Grid/Other/block2.neu');
K1D = 4;
[Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);

% Initialize solver and construct grid and metric
StartUp2D;

% rebuild maps for periodic
BuildPeriodicMaps2D(2,2);

[rp sp] = EquiNodes2D(15); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;


% % Set initial conditions
cx = .1; cy = .1;
D = (x-cx).^2 + (y-cy).^2;
u = exp(-D*5^2).*(1-x.^2).*(1-y.^2);

global bx by
bx = ones(size(x));
by = .5*sin(pi*x); %.*sin(pi*y);

FinalTime = 2;

%% eigs

rLGL = JacobiGQ(0,0,N); rmin = abs(rLGL(1)-rLGL(2));
dtscale = dtscale2D; dt = min(dtscale)*rmin*2/3

if 1 && nargin==0
    e = zeros(Np*K,1);
    A = zeros(Np*K);
    Aavg = zeros(Np*K);
    for i = 1:Np*K
        e(i) = 1;
        ids = 1:Np*K;
        u = reshape(e(ids),Np,K);
        
        % dudt + dudx= 0
        rhsu = AdvecRHS(u,0);      
        A(:,i) = rhsu(:);
        
        rhsuavg = AdvecRHS(u,1);      
        Aavg(:,i) = rhsuavg(:);
        
        e(i) = 0;
        if mod(i,100)==0
            disp(sprintf('on column %d out of %d\n',i,Np*K))
        end
    end
%     lam = eig(A);
%     hold on;
%     plot(lam,'o')
%     title(sprintf('Largest real part = %e',max(real(lam))))
    
    VB = bern_basis_tri(N,r,s);
    if 1
        T1 = kron(eye(K),inv(VB));
        T2 = kron(eye(K),VB);
    else
        T1 = eye(Np*K);
        T2 = T1;
    end
    
    %     imagesc(A)
    h = sqrt(mean(J(:)));
        
    A = T1*A*T2; % transform to bernstein
    A(abs(A)<1e-8) = 0;    
    Aavg = T1*Aavg*T2; % transform to bernstein
    Aavg(abs(Aavg)<1e-8) = 0;    
    [nnz(Aavg),nnz(T2*Aavg*T1)]
    
    A = eye(size(A)) - 1e-1*A;
    Aavg = eye(size(A)) - 1e-1*Aavg;

    P = {};
    for e = 1:K
        ids = (1:Np) + (e-1)*Np;
        %P{e} = Aavg(ids,ids); % block Jacobi
        P{e} = A(ids,ids); % block Jacobi
    end
    P = blkdiag(P{:});    
    P(abs(P)<1e-8) = 0;
    
    cond(T2*(Aavg\A)*T1)        
    cond(T2*(P\A)*T1)        
    
    keyboard
end

%%
time = 0;

% Runge-Kutta residual storage  
resu = zeros(Np,K); 

% compute time step size
rLGL = JacobiGQ(0,0,N); rmin = abs(rLGL(1)-rLGL(2));
dtscale = dtscale2D; dt = min(dtscale)*rmin*2/3

% outer time step loop 
while (time<FinalTime)
  
  if(time+dt>FinalTime), dt = FinalTime-time; end

   for INTRK = 1:5    
      
       rhsu = AdvecRHS(u,1);

      % initiate and increment Runge-Kutta residuals
      resu = rk4a(INTRK)*resu + dt*rhsu;  
        
      % update fields
      u = u+rk4b(INTRK)*resu; 
   end
   
   clf;
   up = Vp*u;
   color_line3(xp,yp,up,up,'.');drawnow
   
   % Increment time
   time = time+dt;
end
return

function rhsu = AdvecRHS(u,take_avg)

Globals2D
global bx by

if nargin==1
    take_avg = 0;
end

if take_avg
    bxavg = (V\bx) / sqrt(2);  bxavg = repmat(bxavg(1,:),Np,1);
    byavg = (V\by) / sqrt(2);  byavg = repmat(byavg(1,:),Np,1);    
    
    alpha = 1; du = zeros(Nfp*Nfaces,K);
    bn = bxavg(vmapM).*nx(:) + byavg(vmapM).*ny(:);
    du(:) = 0.5*(u(vmapM)-u(vmapP)).*(bn - alpha*abs(bn(:)));
    fx = bxavg.*u; fy = byavg.*u;
else
    alpha = 1; du = zeros(Nfp*Nfaces,K);
    bn = bx(vmapM).*nx(:) + by(vmapM).*ny(:);
    du(:) = 0.5*(u(vmapM)-u(vmapP)).*(bn - alpha*abs(bn(:)));
    fx = bx.*u; fy = by.*u;
end

ux = rx.*(Dr*fx)+sx.*(Ds*fx);
uy = ry.*(Dr*fy)+sy.*(Ds*fy);
rhsu = -(ux + uy) + LIFT*(Fscale.*(du));

