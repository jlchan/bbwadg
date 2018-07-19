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


print_pgf_coordinates(rq,u1)
print_pgf_coordinates(rp,Vp*u1)
print_pgf_coordinates(rq,u2)
print_pgf_coordinates(rp,Vp*u2)
print_pgf_coordinates(rp,fdgex(rp))

