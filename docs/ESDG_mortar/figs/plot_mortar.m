N = 3;
r1D = JacobiGQ(0,0,N);
[rq sq] = meshgrid(r1D);

rq1D = JacobiGQ(0,0,N);
e = ones(size(r1D));
r1D2 = JacobiGQ(0,0,N);
r1D2 = [-1+(1+r1D2)/2; (1+r1D2)/2];
e2 = ones(size(r1D2));

rf = e;
sf = r1D;
rf2 = e2;
sf2 = r1D2;

% rf = [r1D; e; r1D; -e];
% sf = [e; r1D; -e; r1D];
% rf2 = [r1D2; e; r1D; -e];
% sf2 = [e2; r1D; -e; r1D];


hold on    
plot([-1 1 1 -1 -1],[-1 -1 1 1 -1],'k-','linewidth',2)
plot(rq,sq,'ko','linewidth',2,'markersize',12,'MarkerFaceColor',[.49 1 .63])
plot(rq,sq,'ko','linewidth',2,'markersize',16,'MarkerFaceColor',[.49 1 .63])

if 0
    % % nonmortar values
    b = .33;
    % non-mortar coupling
    rfi = rf2(2)+b;
    sfi = sf2(N+1+3);
    for i = 1:length(rq(:))
        plot([rq(i) .97*rfi],[sq(i) sfi],'r--','linewidth',2)
    end
    plot(rfi,sfi,'ro','linewidth',2,'markersize',20)
    
else
    % mortar
    a = .7;
    b = 1.4;
    plot(a+rf,sf,'k-','linewidth',2)
    plot(a+rf,sf,'rs','linewidth',2,'markersize',12,'MarkerFaceColor',[.49 1 .63])
    plot(b+rf2(1:N+1),sf2(1:N+1),'k-','linewidth',2)
    plot(b+rf2(N+2:end),sf2(N+2:end),'k-','linewidth',2)
    plot(b+rf2,sf2,'bs','linewidth',2,'markersize',12,'MarkerFaceColor',[.49 1 .63])
    
%     % draw mortar coupling
%     rfi = rf2(2)+b;
%     sfi = sf2(N+1+3);
%     for i = 1:N+1
%         plot([a+rf(i) .97*rfi],[sf(i) sfi],'r--','linewidth',2)
%     end
%     plot(rfi,sfi,'ro','linewidth',2,'markersize',20)
end


% plot(.5*([-1 1 1 -1 -1] + 1) + 1 + 2*b,.5*[-1 -1 1 1 -1]-.5,'k-','linewidth',2)
% plot(.5*([-1 1 1 -1 -1] + 1) + 1 + 2*b,.5*[-1 -1 1 1 -1]+.5,'k-','linewidth',2)
axis equal
axis off
% print(gcf,'-dpng','~/Desktop/nonmortar_coupling.png')
% print(gcf,'-dpng','~/Desktop/bbwadg/talks/RPI_2018/figs/mortar1.png')