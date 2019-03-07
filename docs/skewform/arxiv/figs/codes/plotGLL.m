%%
% useQuads = 1;mypath

[rq1D wq1D] = JacobiGQ(0,0,N);
[rq sq] = meshgrid(rq1D);

rfq = {}; sfq = {};
rfq{1} = rq(1,:); sfq{1} = sq(1,:);
rfq{2} = rq(N+1,:); sfq{2} = sq(N+1,:);
rfq{3} = rq(:,1)'; sfq{3} = sq(:,1)';
rfq{4} = rq(:,N+1)'; sfq{4} = sq(:,N+1)';
% rq = [rfq{:}]; sq = [sfq{:}];
rq = rq(:); sq = sq(:);

plot([-1 1 1 -1 -1],[-1 -1 1 1 -1],'k-','linewidth',2)
hold on;

plot(rq,sq,'bo','linewidth',2,'markersize',15,'MarkerFaceColor',[.49 1 .63]); 

a = 2.1;
plot([-1 1 1 -1 -1]+a,[-1 -1 1 1 -1],'k-','linewidth',2)
plot(rq+a,sq,'bo','linewidth',2,'markersize',15,'MarkerFaceColor',[.49 1 .63]); 
% plot(rfq,sfq,'s','linewidth',2,'markersize',15,'MarkerFaceColor',[.49 1 .63]); 
axis off
axis equal
% text(rq+.025,sq,num2str((1:(N+1)^2)'))
% text(rfq+.025,sfq,num2str((1:4*(N+1))'))
% quiver(rfq,sfq,nrhat,nshat)
% print(gcf,'-dpng','../gsbp.png')

%% staggered grid

[rq1D wq1D] = JacobiGQ(0,0,N);
[rq sq] = meshgrid(rq1D);
% rq = rq(:); sq = sq(:);

e = ones(size(rq1D));
rfq = [rq1D -e rq1D e];
sfq = [-e rq1D e rq1D];
[rq1D wq1D] = JacobiGL(0,0,N+1);
[rq2 sq2] = meshgrid(rq1D);
% rq2 = rq2(:); sq2 = sq2(:);

plot([-1 1 1 -1 -1],[-1 -1 1 1 -1],'k-','linewidth',2)
hold on;

plot(rq,sq,'bo','linewidth',2,'markersize',15,'MarkerFaceColor',[.49 1 .63]); 
plot(rq2,sq2,'rs','linewidth',2,'markersize',15,'MarkerFaceColor',[.49 1 .63]); 
for i = 1:N+1
%     plot(rq2(:,i),sq2(:,i),'k--','linewidth',1)
%     plot(rq2(i,:),sq2(i,:),'k--','linewidth',1)    
end
axis off
axis equal
print(gcf,'-dpng','~/Desktop/proposals/doe/figs/staggered.png')

clf
plot([-1 1 1 -1 -1],[-1 -1 1 1 -1],'k-','linewidth',2)
hold on;
plot(rq,sq,'bo','linewidth',2,'markersize',15,'MarkerFaceColor',[.49 1 .63]); 
plot(rfq,sfq,'rs','linewidth',2,'markersize',15,'MarkerFaceColor',[.49 1 .63]); 
axis off
axis equal
print(gcf,'-dpng','~/Desktop/proposals/doe/figs/gsbp.png')

%% coupled non-con

N = 3;
[rq1D wq1D] = JacobiGL(0,0,N);
[rq sq] = meshgrid(rq1D);

rfq{1} = rq(1,:); sfq{1} = sq(1,:);
rfq{2} = rq(N+1,:); sfq{2} = sq(N+1,:);
rfq{3} = rq(:,1)'; sfq{3} = sq(:,1)';
rfq{4} = rq(:,N+1)'; sfq{4} = sq(:,N+1)';

rq = [rfq{:}]; sq = [sfq{:}];
rq = rq(:); sq = sq(:);

hold on;

a = 1/3; b = 0;
plot([-1 1 1 -1 -1],a*[-1 -1 1 1 -1]+b,'k-','linewidth',2)
plot(rq,a*sq+b,'bo','linewidth',2,'markersize',15,'MarkerFaceColor',[.49 1 .63]); 

a = 1/3; b = -2/3;
plot([-1 1 1 -1 -1],a*[-1 -1 1 1 -1]+b,'k-','linewidth',2)
plot(rq,a*sq+b,'bo','linewidth',2,'markersize',15,'MarkerFaceColor',[.49 1 .63]); 

a = 1/3; b = 2/3;
plot([-1 1 1 -1 -1],a*[-1 -1 1 1 -1]+b,'k-','linewidth',2)
plot(rq,a*sq+b,'bo','linewidth',2,'markersize',15,'MarkerFaceColor',[.49 1 .63]); 

a = -2.25; 
plot(a+[-1 1 1 -1 -1],[-1 -1 1 1 -1],'k-','linewidth',2)
plot(a+rq,sq,'bo','linewidth',2,'markersize',15,'MarkerFaceColor',[.49 1 .63]); 
axis off
axis equal

print(gcf,'-dpng','~/Desktop/proposals/doe/figs/noncon.png')

%% decoupled non-con

N = 3;
[rq1D wq1D] = JacobiGL(0,0,N);
[rq sq] = meshgrid(rq1D);

rfq{1} = rq(1,:); sfq{1} = sq(1,:);
rfq{2} = rq(N+1,:); sfq{2} = sq(N+1,:);
rfq{3} = rq(:,1)'; sfq{3} = sq(:,1)';
rfq{4} = rq(:,N+1)'; sfq{4} = sq(:,N+1)';

rq = [rfq{:}]; sq = [sfq{:}];
rq = rq(:); sq = sq(:);

hold on;

a = 1/3; b = 0;
plot([-1 1 1 -1 -1],a*[-1 -1 1 1 -1]+b,'k-','linewidth',2)
plot(rq,a*sq+b,'bo','linewidth',2,'markersize',15,'MarkerFaceColor',[.49 1 .63]); 


a = 1/3; b = -2/3;
plot([-1 1 1 -1 -1],a*[-1 -1 1 1 -1]+b,'k-','linewidth',2)
plot(rq,a*sq+b,'bo','linewidth',2,'markersize',15,'MarkerFaceColor',[.49 1 .63]); 

a = 1/3; b = 2/3;
plot([-1 1 1 -1 -1],a*[-1 -1 1 1 -1]+b,'k-','linewidth',2)
plot(rq,a*sq+b,'bo','linewidth',2,'markersize',15,'MarkerFaceColor',[.49 1 .63]); 

a = -2.5; 
plot(a+[-1 1 1 -1 -1],[-1 -1 1 1 -1],'k-','linewidth',2)
plot(a+rq,sq,'bo','linewidth',2,'markersize',15,'MarkerFaceColor',[.49 1 .63]); 
axis off
axis equal

plot([-1 1]+a,[-1 -1]-.25,'k-','linewidth',2)
plot(rfq{1}+a,sfq{1}-.25,'rs','linewidth',2,'markersize',15,'MarkerFaceColor',[.49 1 .63]); 
plot([-1 1]+a,[1 1]+.25,'k-','linewidth',2)
plot(rfq{2}+a,sfq{2}+.25,'rs','linewidth',2,'markersize',15,'MarkerFaceColor',[.49 1 .63]); 
plot([-1 -1]+a-.25,[-1 1],'k-','linewidth',2)
plot(rfq{3}+a-.25,sfq{3},'rs','linewidth',2,'markersize',15,'MarkerFaceColor',[.49 1 .63]); 

plot([1 1]+a+.25,[-1 1],'k-','linewidth',2)
plot(rfq{4}+a+.25,sfq{4}/3,'rs','linewidth',2,'markersize',15,'MarkerFaceColor',[.49 1 .63]); 
plot(rfq{4}+a+.25,sfq{4}/3-2/3,'rs','linewidth',2,'markersize',15,'MarkerFaceColor',[.49 1 .63]); 
plot(rfq{4}+a+.25,sfq{4}/3+2/3,'rs','linewidth',2,'markersize',15,'MarkerFaceColor',[.49 1 .63]); 


print(gcf,'-dpng','~/Desktop/proposals/doe/figs/nonconQuad.png')