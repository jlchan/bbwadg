function Advec2D_cfl

Globals2D
global bx by Vq Pq Vfq Pfq Vrq Vsq Prq Psq


% Polynomial order used for approximation
N = 9;

% Read in Mesh
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Grid/Other/squarereg.neu');
% [Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Grid/Other/squareireg.neu');
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D('Grid/Other/block2.neu');

% K1D = 8;
% [Nv, VX, VY, K, EToV] = unif_tri_mesh(K1D);

% Initialize solver and construct grid and metric
StartUp2D;

% Hrefine2D(1:K);
% StartUp2D


% rebuild maps for periodic
BuildPeriodicMaps2D(2,2);

[re se] = EquiNodes2D(N); [re se] = xytors(re,se);
Ve = Vandermonde2D(N,re,se)/V;
xe = Ve*x; ye = Ve*y;

[rp sp] = EquiNodes2D(50); [rp sp] = xytors(rp,sp);
Vp = Vandermonde2D(N,rp,sp)/V;
xp = Vp*x; yp = Vp*y;

[rq sq wq] = Cubature2D(3*N);
% [rq sq wq] = QNodes2D(N); [rq sq] = xytors(rq,sq);
Vq = Vandermonde2D(N,rq,sq)/V;
[Vrq Vsq] = GradVandermonde2D(N,rq,sq);
Vrq = Vrq/V; Vsq = Vsq/V;
xq = Vq*x; yq = Vq*y;
M = Vq'*diag(wq)*Vq;
Pq = M\(Vq'*diag(wq));
Prq = M\(Vrq'*diag(wq));
Psq = M\(Vsq'*diag(wq));


[r1D, w1D] = JacobiGQ(0,0,2*N);
% rfq = [r1D; r1D; -ones(size(r1D))];
% sfq = [-ones(size(r1D)); -r1D; r1D];
wfq = [w1D; w1D; w1D];

Vfq1 = Vandermonde1D(N,r1D)/Vandermonde1D(N,r(Fmask(:,1)));
Vfq2 = Vandermonde1D(N,r1D)/Vandermonde1D(N,r(Fmask(:,2)));
Vfq3 = Vandermonde1D(N,r1D)/Vandermonde1D(N,s(Fmask(:,3)));
Vfq = blkdiag(Vfq1,Vfq2,Vfq3);
Pfq = (Vfq'*diag(wfq)*Vfq)\(Vfq'*diag(wfq));

[re se] = EquiNodes2D(N); [re se] = xytors(re,se);
VBe = bern_basis_tri(N,re,se);
VB = bern_basis_tri(N,r,s);


% % Set initial conditions
cx = .1; cy = .1;
D = (x-cx).^2 + (y-cy).^2;
u = exp(-D*5^2).*(1-x.^2).*(1-y.^2);
u = sin(pi*x).*sin(pi*y);
rad = @(x,y) sqrt(x.^2 + y.^2);
u = rad(x,y) < sin(pi*x);

TV = TV2D(VB\u);
ids = TV > Np/2; 
vp = Vp*u;
% vp = Vp*repmat(TV,Np,1);
color_line3(xp,yp,vp,vp,'.')
return

bx = ones(size(x));
by = ones(size(x));
% by = sin(pi*x).*sin(pi*y); 

% plot3(x,y,-pi*cos(pi*x).*sin(pi*y),'o');return

% divB = rx.*(Dr*bx) + sx.*(Ds*bx) + ry.*(Dr*by) + sy.*(Ds*by);

FinalTime = .1;

%% eigs

if 0
    e = zeros(Np*K,1);
    A = zeros(Np*K);
    for i = 1:Np*K
        e(i) = 1;
        ids = 1:Np*K;
        u = reshape(e(ids),Np,K);
        
        % dudt + dudx= 0
        rhsu = AdvecRHS(u);
        A(:,i) = rhsu(:);
        
        e(i) = 0;
        if mod(i,100)==0
            disp(sprintf('on column %d out of %d\n',i,Np*K))
        end
    end
    [W D] = eig(A); 
    lam = diag(D);
    
%     ids = find(real(lam) > -5e-1);    
%     lam = lam(ids); W = W(:,ids); 
    
    [~,p] = sort(abs(real(lam)),'descend');
    W = W(:,p); lam = lam(p);
    
    ids = abs(lam) > 2*(N+1)*max(1./Fscale(:));
    hold on
    plot(real(lam),imag(lam),'o');
    plot(real(lam(ids)),imag(lam(ids)),'x')
    %     hold on; plot(pi*[1 1],[min(imag(lam)) max(imag(lam))],'--')    
    %     return
    
    skip = nnz(ids);
    
    Q = W(:,skip+1:end);
    W = eye(size(Q,1)) - Q*Q'; % directly project out high modes    
    
    keyboard
    for i = 1:size(W,2)
        w = real(W(:,i));
        w = reshape(w,Np,K);
        clf
        subplot(1,2,1)
        
        vv = Vp*w;
        color_line3(xp,yp,vv,vv,'.')
        view(3)
        
        subplot(1,2,2)
        plot(lam,'o')
        hold on
        plot(real(lam(i)),imag(lam(i)),'rx')
        hold off
        
        pause
    end
    
    plot(real(lam),imag(lam),'o')
    hold on
    mu = eig(Q'*A*Q);
    plot(real(mu),imag(mu),'.')
    
    return
end

%%
time = 0;

% Runge-Kutta residual storage
resu = zeros(Np,K);

% compute time step size
rLGL = JacobiGQ(0,0,N); rmin = abs(rLGL(1)-rLGL(2));
% dtscale = dtscale2D; dt = min(dtscale)*rmin*2/3
dt = .75*rmin./max(Fscale(:));

% outer time step loop
while (time<FinalTime)
    
    if(time+dt>FinalTime), dt = FinalTime-time; end
    
    for INTRK = 1:5
        
        rhsu = AdvecRHS(u);
        
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

figure

uB = VB\u;
TV = TV2D(uB);
ids = TV > 2*(N+1); 
TVp = repmat(TV,length(rp),1);
color_line3(xp,yp,TVp,TVp,'.')

figure
TVp(:,ids) = nan;
color_line3(xp,yp,TVp,TVp,'.')

figure
u(:,ids) = VB*VBe*(uB(:,ids));
% u = VB*VBe*VBe*(VB\u);
up = Vp*u;
color_line3(xp,yp,up,up,'.');
return

function rhsu = AdvecRHS(u)

Globals2D
global bx by Vq Pq Vfq Pfq Vrq Vsq Prq Psq

alpha = 1; du = zeros(Nfp*Nfaces,K);
bn = bx(vmapM).*nx(:) + by(vmapM).*ny(:);
bnf = reshape((bn - alpha*abs(bn(:))),Nfp*Nfaces,K);
ujump = reshape(u(vmapP)-u(vmapM),Nfp*Nfaces,K);

du(:) = .5*Pfq*((Vfq*ujump).*(Vfq*bnf)); % project ujump/bnf
uq = Vq*u; 
fx = Pq*((Vq*bx).*uq); 
fy = Pq*((Vq*by).*uq); 

% fx = bx.*u; fy = by.*u;
% du(:) = 0.5*(ujump).*bnf;

ux = rx.*(Dr*fx) + sx.*(Ds*fx);
uy = ry.*(Dr*fy) + sy.*(Ds*fy);
rhsu = -((ux + uy) + LIFT*(Fscale.*(du)));

