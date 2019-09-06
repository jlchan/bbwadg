clear
K = 10;
%FinalTime = .35;
FinalTime = 2.5;
tau = 1e-3;

xv = linspace(-1,1,K+1)';
x = .5*(xv(1:end-1)+xv(2:end));
xp = linspace(-1,1,1000)';
dx = (max(xv(:))-min(xv(:)))/K;

uex = @(x) exp(sin(1+pi*x)); %x+.5*tanh(5*(x-.25));
u = uex(x);
plot(xp,uex(xp),'r-','linewidth',2)
hold on
for i = 1:K
    plot([xv(i) xv(i+1)],u(i)*[1 1],'b-','linewidth',2)
    plot(x(i),u(i),'bo','linewidth',2,'markersize',12)
end

i = 5;
text(x(i)-.5/K, u(i)+.25,'${\bf u}_{i-1}$','fontsize',36,'Interpreter','latex')
i = 6;
text(x(i), u(i)+.25,'${\bf u}_i$','fontsize',36,'Interpreter','latex')
i = 7;
text(x(i)-.5/K, u(i)+.25,'${\bf u}_{i+1}$','fontsize',36,'Interpreter','latex')

axis([-1,1,0, 3.25])
grid on
xlabel('')
ylabel('')
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
% axis equal