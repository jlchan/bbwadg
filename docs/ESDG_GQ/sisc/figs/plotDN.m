N = 5;
[rq wq] = JacobiGQ(0,0,N);
rp = linspace(-1,1,50)';
Vq = Vandermonde1D(N,rq);
D = GradVandermonde1D(N,rq)/Vq;
Vf = Vandermonde1D(N,[-1;1])/Vq;
Vp = Vandermonde1D(N,rp)/Vq;

B = [-1 0;0 1];
QN = [D - .5*Vf'*B*Vf .5*Vf'*B;
    -.5*B*Vf .5*B];
PN = [eye(N+1) diag(1./wq)*Vf'];

% gex = @(x) exp(x).*sin(pi*x);%exp(-25*x.^2); %
% gex = @(x) 1./(1+25*(x-0).^2);
a = 4;
gex = @(x) exp(-a*x.^2);
dgex = @(x) -2*a*x.*exp(-a*x.^2);
fex = @(x) 1 + 0*sin(pi*x);
f = fex(rq);
g = gex(rq);
fdgex = @(x) fex(x).*dgex(x);

fN = fex([rq;-1;1]);
gN = gex([rq;-1;1]);
u1 = f.*(D*g);
u2 = PN*(fN.*(QN*gN));
hold on
p1 = plot(rp,Vp*u1,'b-','linewidth',2);
p2 = plot(rp,Vp*u2,'r-','linewidth',2);
plot(rq,u1,'bo','linewidth',2,'markersize',15,'MarkerFaceColor',[.49 1 .63])
plot(rq,u2,'rs','linewidth',2,'markersize',15,'MarkerFaceColor',[.49 1 .63])
plot(rp,fdgex(rp),'k--','linewidth',2)
% legend('GSBP','Decoupled SBP')
set(gca,'fontsize',16)
grid on


% print_pgf_coordinates(rq,u1)
% print_pgf_coordinates(rp,Vp*u1)
% print_pgf_coordinates(rq,u2)
% print_pgf_coordinates(rp,Vp*u2)
% print_pgf_coordinates(rp,fdgex(rp))

%% convergence

for N = 1:15
    [rq wq] = JacobiGQ(0,0,N);
    rp = linspace(-1,1,50)';
    Vq = Vandermonde1D(N,rq);
    D = GradVandermonde1D(N,rq)/Vq;
    Vf = Vandermonde1D(N,[-1;1])/Vq;
    Vp = Vandermonde1D(N,rp)/Vq;
    
    [rq2 wq2] = JacobiGQ(0,0,N+4);
    Vq2 = Vandermonde1D(N,rq2)/Vq;
    
    B = [-1 0;0 1];
    QN = [D - .5*Vf'*B*Vf .5*Vf'*B;
        -.5*B*Vf .5*B];
    PN = [eye(N+1) diag(1./wq)*Vf'];
    
    % gex = @(x) exp(x).*sin(pi*x);%exp(-25*x.^2); %
    % gex = @(x) 1./(1+25*(x-0).^2);
    a = 4;
    gex = @(x) exp(-a*x.^2);
    dgex = @(x) -2*a*x.*exp(-a*x.^2);
    fex = @(x) 1 + 0*sin(pi*x);
    f = fex(rq);
    g = gex(rq);
    fdgex = @(x) fex(x).*dgex(x);
    
    fN = fex([rq;-1;1]);
    gN = gex([rq;-1;1]);
    u1 = f.*(D*g);
    u2 = PN*(fN.*(QN*gN));
    
    err1(N) = sqrt(sum(wq2.*(fdgex(rq2)-Vq2*u1).^2));
    err2(N) = sqrt(sum(wq2.*(fdgex(rq2)-Vq2*u2).^2));
end
semilogy(1:N,err1,'o--')
hold on
semilogy(1:N,err2,'x--')
print_pgf_coordinates(1:N,err1)
print_pgf_coordinates(1:N,err2)

%% convergence in h

clear
Globals1D

% Order of polymomials used for approximation
Kvec = [2 4 8 16 32 64 128];
N = 1;
for sk = 1:length(Kvec)
    K1D = Kvec(sk);
    [Nv, VX, K, EToV] = MeshGen1D(-1,1,K1D);
    StartUp1D;
    
    [rq wq] = JacobiGQ(0,0,N);
    rp = linspace(-1,1,50)';
    Vq = Vandermonde1D(N,rq);
    D = GradVandermonde1D(N,rq)/Vq;
    Vf = Vandermonde1D(N,[-1;1])/Vq;
    Vp = Vandermonde1D(N,rp)/Vq;
    
    [rq2 wq2] = JacobiGQ(0,0,N+4);
    Vq2 = Vandermonde1D(N,rq2)/Vq;
    
    B = [-1 0;0 1];
    QN = [D - .5*Vf'*B*Vf .5*Vf'*B;
        -.5*B*Vf .5*B];
    PN = [eye(N+1) diag(1./wq)*Vf'];
    
    xq = (Vq/V)*x;
    xq2 = Vq2*xq;
    
    a = 4;
    gex = @(x) exp(-a*x.^2);
    dgex = @(x) -2*a*x.*exp(-a*x.^2);
    fex = @(x) 1 + 0*sin(pi*x);
    f = fex(xq);
    g = gex(xq);
    fdgex = @(x) fex(x).*dgex(x);
        
    xN = [xq;Vf*xq];
    fN = fex(xN);
    gN = gex(xN);
    u1 = rx.*(f.*(D*g));
    u2 = rx.*(PN*(fN.*(QN*gN)));
    
%     xp = (Vandermonde1D(N,rp)/V)*x;
%     plot(xp,Vp*u1)
%     hold on
%     plot(xp,Vp*u2)
%     return

    wJq2 = diag(wq2)*((Vandermonde1D(N,rq2)/V)*J);
%     err1(sk) = sqrt(sum(sum(wJq2.*((gex(xq2)-Vq2*g).^2))));
    err1(sk) = sqrt(sum(sum(wJq2.*((fdgex(xq2)-Vq2*u1).^2))));
    err2(sk) = sqrt(sum(sum(wJq2.*((fdgex(xq2)-Vq2*u2).^2))));
end
h = 1./Kvec;
loglog(h,err1,'o--')
hold on
loglog(h,err2,'x--')
loglog(h,h.^(N),'k--')