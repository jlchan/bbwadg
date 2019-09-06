du = logspace(-8,0,1000);

uL = 1e0*1.1; %rand;
uR = uL + du;

% uL = single(uL);
% uR = single(uR);

uLog = logmean(uL,uR);
uLogNaive = (uL-uR)./(log(uL)-log(uR));
uavg = (uL+uR)/2;
uMax = max(abs(uL),abs(uR));
loglog(du, abs(uLogNaive - logmean(uL,uR))./uMax,'x','linewidth',2,'markersize',10,'DisplayName','uLog-uLogNaive')
hold on
loglog(du, abs(uLogNaive - uavg)./uMax,'o','linewidth',2,'markersize',10,'DisplayName','uLogNaive-uavg')
loglog(du, abs(uLog - uavg)./uMax,'.','linewidth',2,'markersize',10,'DisplayName','uLog-uavg')

h = legend('show');
set(h,'fontsize',16)
xlabel('Relative perturbation','fontsize',16)
set(gca,'fontsize',16)

loglog(du, 1e-11*ones(size(du)),'--','linewidth',2)

grid on