% load romerr25
% load romerr75
load rom_errors

sk = 50;
semilogy(dt*(1:Nsteps),err25,'b-','linewidth',2);
hold on;
semilogy(dt*(1:sk:Nsteps),err25fullquad(1:sk:Nsteps),'bx','linewidth',2,'markersize',10);

semilogy(dt*(1:Nsteps),err75,'r-','linewidth',2);
semilogy(dt*(1:sk:Nsteps),err75fullquad(1:sk:Nsteps),'rx','linewidth',2,'markersize',10);

semilogy(dt*(1:Nsteps),err125,'k-','linewidth',2);
semilogy(dt*(1:sk:Nsteps),err125fullquad(1:sk:Nsteps),'kx','linewidth',2,'markersize',10);

grid on
fontsz = 16;
xlabel('Time','fontsize',fontsz)
h = legend('With hyper-reduction','Without hyper-reduction','Location','Best');
set(h,'fontsize',fontsz)
set(gca,'fontsize',fontsz)
axis tight
text(dt*50,.05,'25 modes','fontsize',24)
text(dt*1200,.018,'75 modes','fontsize',24)
text(dt*1500,.003,'125 modes','fontsize',24)

%%
load rom_entropy 
load fom_entropy 
semilogy(dt*(1:Nsteps),fom_entropy,'b-','linewidth',2);
hold on;
semilogy(dt*(1:Nsteps),rom_entropy25,'k:','linewidth',2);
semilogy(dt*(1:Nsteps),rom_entropy75,'r--','linewidth',2)
grid on
fontsz = 16;
xlabel('Time','fontsize',fontsz)
% ylabel({'Discrete entropy dissipation'},'fontsize',fontsz);
% legend show
% legend location southwest
%legend('FOM','Option 1','Option 2','Option 3','Location','Best')
h = legend('FOM','25 modes','75 modes','Location','Best');
set(h,'fontsize',fontsz)
set(gca,'fontsize',fontsz)

%%

load dtest_rom_25modes
load dtest_rom_75modes
load dtest_fom

Nsteps = length(dtest);
dt = .75/Nsteps;
skip = 20;
plot(dt*(1:skip:Nsteps),dtest(1:skip:Nsteps),'--','linewidth',2,'markersize',10)
hold on
plot(dt*(1:skip:Nsteps),dtest1(1:skip:Nsteps),':','linewidth',2,'markersize',10)
semilogy(dt*(1:skip:Nsteps),dtest2(1:skip:Nsteps),'-.','linewidth',2,'markersize',10)
semilogy(dt*(1:skip:Nsteps),dtest3(1:skip:Nsteps),'.-','linewidth',2,'markersize',10)
grid on
fontsz = 16;
xlabel('Time','fontsize',fontsz)
% ylabel({'Discrete entropy dissipation'},'fontsize',fontsz);
% legend show
% legend location southwest
%legend('FOM','Option 1','Option 2','Option 3','Location','Best')
h = legend('FOM','(5.5)','(5.7)','(5.8)','Location','Best');
set(h,'fontsize',fontsz)
axis tight
set(gca,'fontsize',fontsz)
ylim([.75,4.25])

%%
skip = 5;
% plot(dt*(1:skip:Nsteps),dtest(1:skip:Nsteps),'--','linewidth',2,'markersize',10)
semilogy(dt*(1:skip:Nsteps),abs(dtest1(1:skip:Nsteps)-dtest1(1:skip:Nsteps)),'o','linewidth',2,'markersize',10)
hold on
semilogy(dt*(1:skip:Nsteps),abs(dtest1(1:skip:Nsteps)-dtest2(1:skip:Nsteps)),'x--','linewidth',2,'markersize',10)
semilogy(dt*(1:skip:Nsteps),abs(dtest1(1:skip:Nsteps)-dtest3(1:skip:Nsteps)),'.-','linewidth',2,'markersize',10)
grid on
fontsz = 16;
xlabel('Time','fontsize',fontsz)
% ylabel({'Discrete entropy dissipation'},'fontsize',fontsz);
% legend show
% legend location southwest
set(gca,'fontsize',fontsz)
h = legend('Difference (5.5)-(5.7)','Difference (5.5)-(5.8)','Location','Best');
set(h,'fontsize',fontsz)
