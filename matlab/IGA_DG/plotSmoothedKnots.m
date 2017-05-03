%% try diff approaches: iterate to t = sum(re_j*BNj(t))
clear
NB = 3;
Ksub = 8;
N = NB+Ksub-1;
VX = linspace(-1,1,Ksub+1);
t = [VX(1)*ones(1,NB) VX VX(end)*ones(1,NB)]; % open knot vec
t0 = t;
for i = 1:N+1
    r(i,1) = mean(t((i+1):(i+NB))); % greville
end

re = linspace(-1,1,N+1)';    
R = [];
V = [];
for ii = 1:10    
    R = [R;r(:)'];    
    V = [V;t(:)'];

    Ve = bspline_basismatrix(NB+1, t, t0);
    t = Ve*re; 
    
    % recompute greville
    for i = 1:N+1
        r(i,1) = mean(t((i+1):(i+NB))); 
    end    
end

% N = 3, K = 8
VXopt = [   -1.0000   -0.6661   -0.4443   -0.2222    0.0000    0.2222    0.4443    0.6661    1.0000]; 
topt = [VXopt(1)*ones(1,NB) VXopt VXopt(end)*ones(1,NB)]'; % open knot vec
for i = 1:N+1
    ropt(i,1) = mean(topt((i+1):(i+NB)));
end

figure
hold on
for i = 1:size(V,2)
    plot(V(:,i),1:size(V,1),'ko-','markersize',8,'linewidth',2,'MarkerFaceColor',[.49 1 .63])
end
set(gca,'fontsize',14)
grid on
ylabel('Iteration k','fontsize',14)
axis([-1.1 1.1 0 ii+1])

hold on;plot(VXopt,ones(size(VXopt))*ii,'ro','markersize',8,'linewidth',2,'DisplayName','Optimal knots')
print(gcf,'-dpng','~/Desktop/IGA-DG/docs/multipatch/figs/knots.png')

figure
hold on
for i = 1:size(R,2)
    plot(R(:,i),1:size(R,1),'ko-','markersize',8,'linewidth',2,'MarkerFaceColor',[.49 1 .63])
end

set(gca,'fontsize',14)
grid on
ylabel('Iteration k','fontsize',14)
axis([-1.1 1.1 0 ii+1])


hold on;plot(ropt,ones(size(ropt))*(ii),'ro','markersize',8,'linewidth',2,'DisplayName','Optimal knots')
% legend show
print(gcf,'-dpng','~/Desktop/IGA-DG/docs/multipatch/figs/greville.png')


return
%%

clear
NB = 3;
Ksub = 8;
N = NB+Ksub-1;

VX = linspace(-1,1,Ksub+1);
t = [VX(1)*ones(1,NB) VX VX(end)*ones(1,NB)]; % open knot vec
VX0 = VX;
for i = 1:N+1
    r(i,1) = mean(t((i+1):(i+NB))); % greville
end

% smooth knots
for ii = 1:25
    re = linspace(-1,1,N+1)';
    reKsub = linspace(-1,1,Ksub+1)';
    Ve = bspline_basismatrix(NB+1, t, reKsub);
    VX = Ve*re; VX = VX(:)';
        
    if ii > 1
        rold = rnew;
    end
    t = [VX(1)*ones(1,NB) VX VX(end)*ones(1,NB)]; % open knot vec
    for i = 1:N+1
        rnew(i,1) = mean(t((i+1):(i+NB))); % greville
    end
    if ii > 1
        err(ii-1) = max(abs(rnew-rold));
    end
    if ii==1
        r1 = rnew;
        VX1 = VX;
    end
end

% semilogy(err,'o--');
% hold on
% semilogy(exp(-.5*(1:length(err))),'--');
% return

% rp = linspace(-1,1,500);
% Vp = bspline_basismatrix(NB+1, t, rp);

plot(VX0,VX0*0,'bo','markersize',15,'linewidth',2);
hold on
plot(VX1,VX1*0,'rx','markersize',15,'linewidth',2)
% plot(VX,VX*0,'ks','markersize',15,'linewidth',2)
set(gca,'fontsize',14)
grid on
set(gca,'yticklabel',[])
axis([-1.0 1. -1 1])
%legend('Original knots','Smoothed knots, k=1','Smoothed knots, k = 25')
legend('Original knots','Smoothed knots')
print(gcf,'-dpng','~/Desktop/IGA-DG/docs/multipatch/figs/knots.png')

clf
plot(r,r*0,'bo','markersize',15,'linewidth',2)
hold on
plot(r1,r1*0,'rx','markersize',15,'linewidth',2)
% plot(rnew,rnew*0,'ks','markersize',15,'linewidth',2)
%legend('Original Greville points','Greville points, k=1','Greville points, k = 25')
legend('Original Greville points','Smoothed Greville points')
set(gca,'fontsize',14)
grid on
set(gca,'yticklabel',[])
axis([-1.0 1. -1 1])
print(gcf,'-dpng','~/Desktop/IGA-DG/docs/multipatch/figs/greville.png')

% axis tight
% figure
% plot(VX0,VX0-VX,'o--') % plot knot displacement
% hold on