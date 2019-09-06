clear

N = 7;
% [rq wq] = JacobiGQ(0,0,N+4); 

opt = 2;
if opt==1
    [rq wq] = JacobiGQ(0,0,N);
    [r w] = JacobiGQ(0,0,N);
    p = @(x) JacobiP(x,0,0,N+1);
elseif opt==2
    [rq wq] = JacobiGL(0,0,N);
    [r w] = JacobiGL(0,0,N);
    p = @(x) (1+x).*(1-x).*GradJacobiP(x,0,0,N);
else % Newton Cotes
    [rq wq] = JacobiGQ(0,0,N+4); 
    r = linspace(-1,1,N+1)';
%     jj = 1:N+1; r = sort(cos((2*jj-1)*pi/(2*(N+1)))');
%     p = @(x) cos((N+1)*acos(x));
    VNq = Vandermonde1D(N,rq)/Vandermonde1D(N,r); w = VNq'*wq;
    w = 2/(N+1)*ones(N+1,1);
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
Lf = M\(Vf');

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
VpNp1 = Vandermonde1D(N+1,rp)/VNp1;
VqNp1 = Vandermonde1D(N+1,rq)/VNp1;
DNp1 = GradVandermonde1D(N+1,r)/VNp1;
VNm1 = Vandermonde1D(N+1,r)/VNp1;
VfvNp1 = Vandermonde1D(N+1,rf)/VNp1;
iD = [DNp1; wq'*VqNp1]\[eye(N+1); 0*wq'*Vq]; % set zero avg
L = M\(Vf'); % lift
C = iD*L; % correction function
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

norm(QCR - QCL - M*L,'fro')
norm(QCR - QCL - QC,'fro')
norm(QfvR - QfvL - Q,'fro')
norm(AA-(QfvR+QfvL)/2) % should be zero by construction

% testing
uex = @(r) exp(r);
uex = @(r) exp(sin(pi*r));
% uex = @(r) 1 + r;
% uex = @(r) (tanh(1e3*(r-.2)));
% uex = @(r) JacobiP(r,0,0,N+1);
u = Pq*uex(rq); % + VNm1*C*(uex([-1;1])-Vf*Pq*uex(rq));

% % plot sub-cell averages
% ua = AA*u;
% clf
% plot(rp,Vp*u,'--','linewidth',2)
% hold on
% plot(r,u,'o','linewidth',2,'markersize',16)
% plot(rf,Vfv*u,'x','linewidth',2,'markersize',16)
% for i = 1:N+1
%     plot(rf(i:i+1),ua(i)*[1 1],'-','linewidth',2)
%     plot(rf(i:i+1),u(i)*[1 1],'--','linewidth',2)
% end

% figure
% clf
% plot(rp,Vp*u,'-','linewidth',2)
% hold on

ua = AA*u; % cell avgs
du = .5*(uex([-1;1]) - Vf*u); % weak deriv flux
% du = 0*du;
Qflux = [QfvL*u + QCL*du, QfvR*u + QCR*du];  % Wfv'*

% 
% % uwfv = Wfv*u;
% % for i = 1:N+1
% %     plot(rf(i:i+1),uwfv(i)*[1 1],'-','linewidth',2)
% % end
% 
% for i = 1:N+1
%     plot(rf(i:i+1),Qflux(i,1:2),'o--','linewidth',2,'markersize',16)
%     
%     % nodal values, cell avgs
%     plot(r(i),u(i),'.','linewidth',2,'markersize',50)
%     plot(rf(i:i+1),ua(i)*[1 1],'-','linewidth',2,'markersize',16)
%     
%     % compare to subcell P1 projection
%     V1q = Vandermonde1D(1,rq)/Vandermonde1D(1,[-1,1]);
%     P1 = (V1q'*diag(wq)*V1q)\(V1q'*diag(wq));
%     P1flux(i,:) = P1*Vqi{i}*u;
% %     plot(rf(i:i+1),P1flux(i,1:2),'x--','linewidth',2,'markersize',16)
% end
% 
% title(sprintf('diff bw recon and P1 proj = %g\n',norm(Qflux-P1flux,'fro')/(N+2)))

figure
hold on
plot(rp,Vp*u)

xc = sum([rf(1:end-1) rf(2:end)],2)/2;

slopes = Qflux(:,2)-Qflux(:,1);
for i = 2:N
    plot(rf(i:i+1),Qflux(i,1:2),'o--','linewidth',2,'markersize',16)
    
    % nodal values, cell avgs
    plot(r(i),u(i),'.','linewidth',2,'markersize',50)
    plot(rf(i:i+1),ua(i)*[1 1],'-','linewidth',2,'markersize',16)
    
    x0 = .5*sum(rf(i:i+1));
    du = (Qflux(i,2)-Qflux(i,1))/w(i);
    %     plot(rf(i:i+1),ua(i) + du*(rf(i:i+1) - x0),'--')    
    %         duL = (Qflux(i-1,2)-Qflux(i-1,1))/w(i-1);
    %         duR = (Qflux(i+1,2)-Qflux(i+1,1))/w(i+1);
    
    % MUSCL minmod limiter
    xL = .5*sum(rf(i-1:i));
    xR = .5*sum(rf(i+1:i+2));
    hL = (x0-xL)/2;
    hR = (xR-x0)/2;
    hL = w(i)/2;
    hR = w(i)/2;
    duL = (ua(i)-ua(i-1))/hL;
    duR = (ua(i+1)-ua(i))/hR;
    du = minmod([duL;du;duR]);
    plot(rf(i:i+1),ua(i) + du*(rf(i:i+1) - x0),'^-','linewidth',2,'markersize',16)
end
return

%% solve for polynomial sub-cell
Q = M*D;
Iwt = .5*abs(DF)';
Iq1 = [DF; (Iwt*w)']\[Q; w'];
Bfv = [-1 0;zeros(N-1,2);0 1];
Iq2 =[Vf(1,:); DF(:,2:end-1)\(Q - Bfv*Vf); Vf(2,:)];

% DCDR version of sub-cell
Mtp = tril(ones(N+2,N+1),-1); Mt = -triu(ones(N+2,N+1)); 
% I = eye(N+1); tL = I(:,1); tR = I(:,end);
tL = Vf(1,:)';
tR = Vf(2,:)';
Iq = .5*((Mtp+Mt)*Q + ones(N+2,1)*(tL'+tR'));

QQ = zeros(N+2,N+1);
QQ(1,:) = QfvL(1,:);
for i = 2:N+1
    QQ(i,:) = .5*(QfvR(i-1,:) + QfvL(i,:));
end
QQ(end,:) = QfvR(end,:);


%%

% look at jumps
figure
semilogy(rp,abs(uex(rp)-Vp*u),'--','linewidth',2,'markersize',16)
hold on
for i = 1:N+1
    % plot Lf type shock indicator
    if i==1
        plot(r(i),abs(Qflux(i,2)-Qflux(i+1,1)),'x','linewidth',2,'markersize',16)
    elseif i==N+1
        plot(r(i),abs(Qflux(i,1)-Qflux(i-1,2)),'x','linewidth',2,'markersize',16)
    else
        absLf = abs(Qflux(i,2)-Qflux(i+1,1))+abs(Qflux(i,1)-Qflux(i-1,2));
        plot(rf(i:i+1),absLf*[1 1],'-','linewidth',2,'markersize',16)
    end
    
    % plot local estimator || u - u_1 ||^2
    V1q = Vandermonde1D(1,rq)/Vandermonde1D(1,[-1,1]);
    plot(rf(i:i+1),sqrt(sum(w(i)/sum(wq)*wq.*(Vqi{i}*u - V1q*Qflux(i,1:2)').^2))*[1 1],'--','linewidth',2,'markersize',16)
end

plot([-1 1],max(w).^2*[1 1],'--','linewidth',2)


