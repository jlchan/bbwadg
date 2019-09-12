% time

N = 1e7;
x = 2*rand(N,1);
y = rand(N,1);
logx = log(x);
logy = log(y);
logmean(x,y,logx,logy);
logmean(x,y);
return 

d = 5;
du = logspace(-d,d,1000);

uL = 10^d*1.1; %rand;
uR = uL + du;

% uL = single(uL);
% uR = single(uR);

uLog = logmean(uL,uR);
uLogNaive = (uL-uR)./(log(uL)-log(uR));
% uavg = .5*(sqrt(uL.*uR)+(uL+uR)/2);
uavg = .5*(uL+uR);
uMax = max(abs(uL),abs(uR));
loglog(du./uMax, abs(uLogNaive - logmean(uL,uR))./uMax,'x','linewidth',2,'markersize',10,'DisplayName','uLog-uLogNaive')
hold on
loglog(du./uMax, abs(uLogNaive - uavg)./uMax,'o','linewidth',2,'markersize',10,'DisplayName','uLogNaive-uavg')
loglog(du./uMax, abs(uLog - uavg)./uMax,'.','linewidth',2,'markersize',10,'DisplayName','uLog-uavg')

h = legend('show');
set(h,'fontsize',16)
xlabel('Relative perturbation','fontsize',16)
set(gca,'fontsize',16)

grid on
axis tight

