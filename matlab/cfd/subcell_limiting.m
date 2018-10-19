clear
N = 5;
[rq wq] = JacobiGQ(0,0,N+4);
% [rq wq] = JacobiGL(0,0,N);

opt = 1;
if opt==1
    [r w] = JacobiGQ(0,0,N);
    p = @(x) JacobiP(x,0,0,N+1);
elseif opt==2
    [r w] = JacobiGL(0,0,N);
    p = @(x) (1+x).*(1-x).*GradJacobiP(x,0,0,N);
else % Newton Cotes
    r = linspace(-1,1,N+1)';
    jj = 1:N+1; r = sort(cos((2*jj-1)*pi/(2*(N+1)))');
    %     p = @(x) cos((N+1)*acos(x));
    VNq = Vandermonde1D(N,rq)/Vandermonde1D(N,r); w = VNq'*wq;
%     w = 2/(N+1)*ones(N+1,1);
end
% r_nodal = JacobiGL(0,0,N);
r_nodal = r;

rf = [-1;cumsum(w)-1];

rsub = zeros(length(rq),N+1);
for i = 1:N+1
    a = rf(i);
    b = rf(i+1);
    rsub(:,i) = a + (b-a)*(1+rq)/2;
    wsub(:,i) = (b-a)*wq/2;
end

V = Vandermonde1D(N,r_nodal);
D = GradVandermonde1D(N,r_nodal)/V;
Vq = Vandermonde1D(N,rq)/V;
Vfv = Vandermonde1D(N,rf)/V;

M = Vq'*diag(wq)*Vq;
for i = 1:N+1
    Vqi{i} = Vandermonde1D(N,rsub(:,i))/V;
    B(:,i) = Vqi{i}'*wsub(:,i); % subcell integrals of basis
end
Wfv = M\B; % test functions which extract cell avgs

Pq = M\(Vq'*diag(wq));
Vf = Vandermonde1D(N,[-1;1])/V;

rp = linspace(-1,1,1e3)';
Vp = Vandermonde1D(N,rp)/V;
Vpsub = Vandermonde1D(N,rsub(:))/V;

% maps nodal -> sub-cell averages
% AA = M\(Wfv'*M);
AA = diag(1./w)*(Wfv'*M);

% plot fluxes
Q = M*D;
e = ones(size(r));
DF = diag([-e;0]) + diag(e,1);
DF = DF(1:end-1,:);

% correction functions
rNp1 = JacobiGL(0,0,N+1);
VNp1 = Vandermonde1D(N+1,rNp1);
% VpNp1 = Vandermonde1D(N+1,rp)/VNp1;
VqNp1 = Vandermonde1D(N+1,rq)/VNp1;
DNp1 = GradVandermonde1D(N+1,r)/VNp1;
% VNm1 = Vandermonde1D(N+1,r)/VNp1;
VfvNp1 = Vandermonde1D(N+1,rf)/VNp1;
iD = [DNp1; wq'*VqNp1]\[eye(N+1); 0*wq'*Vq]; % set zero avg
LIFT = M\(Vf'); % lift
C = iD*LIFT; % correction function
QC = M*DNp1*C;
QC2 = (Wfv')\(DF*VfvNp1*C);
VfvC = VfvNp1*C;

% L/R correction fluxes
QCR = (Wfv')\(VfvC(2:N+2,:));
QCL = (Wfv')\(VfvC(1:N+1,:));
QCavg = .5*(QCR + QCL);
QCR = QCR-QCavg;
QCL = QCL-QCavg;

% L/R derivative fluxes
QfvR = (Wfv')\Vfv(2:end,:) ;
QfvL = (Wfv')\Vfv(1:end-1,:);
Qavg = .5*(QfvR + QfvL);
QfvR = QfvR - Qavg + AA; % preserve avg
QfvL = QfvL - Qavg + AA;

norm(QCR - QCL - M*LIFT,'fro')
norm(QCR - QCL - QC,'fro')
norm(QfvR - QfvL - Q,'fro')
norm(AA-(QfvR+QfvL)/2) % should be zero by construction

% testing
uex = @(r) exp(-10^2*(r+.5).^2);
uex = @(r) sin(pi*r);
%uex = @(r) 1e-11 + exp(-10^2*(r+.5).^2) + (abs(r-.5) < .2);
% uex = @(r) (abs(r) < .2);

%% make mesh, ops

K1D = 40;
[Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);

h = VX(2)-VX(1); % assume unif
rxJ = 1;
J = h/2;

for e = 1:K
    x(:,e) = VX(e) + .5*(1+r)*h;
end
% plot(x,x*0,'o')
% return

% make periodic
Nfp = 1; Nfaces = 2;
mapM = reshape(1:Nfp*Nfaces*K,Nfp*Nfaces,K);
mapP = mapM;
for e = 2:K-1
    mapP(1,e) = mapM(2,e-1);
    mapP(2,e) = mapM(1,e+1);
end
if K==1
    mapP(:,1) = [mapM(2,1);mapM(1,1)];
else
    mapP(:,1) = [mapM(2,K);mapM(1,2)];
    mapP(:,K) = [mapM(2,K-1);mapM(1,1)];
end

xf = Vf*x;
xq = Vq*x;
xp = Vp*x;
xfv = Vfv*x;

u = Pq*uex(xq);
umin = min(u(:))-.5;
umax = max(u(:))+.5;

%%

% runge kutta coefficients
rk3a = [0 .75 1/3];
rk3b = [1 .25 2/3];
rk3c = [1 .25 2/3];

FinalTime = 2.0;
CFL = 2;
dt = CFL*h/(N+1)^2; % set time-step
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;

nx = [-1;1];

figure(1)
resu = zeros(size(u));
for i = 1:Nsteps
    uold = u;
    for INTRK = 1:3 
        uf = Vf*u;
        du = uf(mapP)-uf(mapM);
        
        flux = .5*nx.*du - .5*du;
        fR = rxJ*QfvR*u + QCR*flux;
        fL = rxJ*QfvL*u + QCL*flux;

        %         uR = QfvR*u;
        %         uL = QfvL*u;        
        du = [zeros(1,K); diff(AA*u); zeros(1,K)];
        uL = du(1:end-1,:);
        uR = du(2:end,:);
        tau = 0*h*w; 
        rhs = -((fR-fL)./w - tau.*(uR-uL))./J; 
         
        u = rk3a(INTRK)*uold + rk3b(INTRK)*u + rk3c(INTRK)*dt*rhs;
    end
        
    if mod(i,25)==0
        clf
        ua = AA*u;
        plot(xp,Vp*u,'-','linewidth',2,'markersize',16)
        hold on
%         plot(r,u,'o','linewidth',2,'markersize',16)       
        for i = 1:N+1
%             plot(rf(i:i+1),[fL(i,:) fR(i,:)],'--','linewidth',2)
            plot(xfv(i:i+1,:),[1; 1]*ua(i,:),'--','linewidth',2)
        end
        axis([-1,1,umin,umax])
        drawnow
    end
end
figure
plot(xq,Vq*u - uex(xq))
sqrt(sum(sum(wq.*(Vq*u - uex(xq - mod(FinalTime,2))).^2)))
%sqrt(sum(sum(wq.*(Vq*u - uex(rq - FinalTime)).^2)))