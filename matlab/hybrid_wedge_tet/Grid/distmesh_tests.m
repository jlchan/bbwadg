fd1=@(p) sqrt(sum((p).^2,2))-1;
fd2=@(p) sqrt(sum((p + 1).^2,2))-1;
fd = @(p) dunion(fd1(p),fd2(p));
%fd=inline('-0.3+abs(0.7-sqrt(sum(p.^2,2)))');
[p,t]=distmesh2d(fd,@huniform,.1,[-2,-2;2,2],[]);
axis on

%%

K = 10;
[x y] = meshgrid(linspace(-2,2,10*K));
% dd = .25*sin(4*x)-y;
dd = ones(size(x));
dd = dd.*(1-x).*(1+x);
dd = dd.*(1-y).*(1+y);
% dd = dd-1;
% dd = dd/max(abs(dd(:))); % normalize
% dd=max(sqrt(xx.^2+yy.^2),.96)-1;

% h = color_line3(x(:),y(:),dd,dd,'.');set(h,'markersize',32);
% ids = find(abs(dd)<5e-2); hold on;plot3(x(ids),y(ids),ones(size(ids)),'k.');return
fd1 = @(p) drectangle(p,-1,1,-1,1);
%fd2 = @(p) p(:,2)-.2*(sin(1*p(:,1)) + sin(5*p(:,1))/5);
fd2 = @(p) .2*(sin(1*p(:,1)) + sin(5*p(:,1))/5)+1-p(:,2);
%fd = @(p) dintersect(fd1(p),fd2(p));
fd = @(p) ddiff(fd1(p),fd2(p));
% [p,t]=distmesh2d(fd,@huniform,2/K,2*[-1,-1;1,1],[]);
[p,t]=distmesh2d(fd,@huniform,.1,[-1,-1;1,1],[-1,-1;-1,1;1,-1;1,1]);
%[p,t]=distmesh2d(@dmatrix,@huniform,2/K,[-1,-1;1,1],[],xx,yy,dd);
axis on

%%