function PAT_mesh

load PAT/PAT_breast_boundary.mat
load PAT/PAT_setup.mat
hmin = .125;

tb = linspace(min(t),max(t)-hmin,round((max(t)-min(t))/hmin));

xb = interp1(t,xb,tb);
yb = interp1(t,yb,tb);

tb = [tb tb(1)]; xb = [xb xb(1)]; yb = [yb yb(1)]; 
plot(xb,yb,'o--');return

node = [xb(:) yb(:)];
hdata.fun = @(x,y) hmin*ones(size(x));
[p,t] = mesh2d(node,[],hdata);

d1 = size(speed,1); d2 = size(speed,2);
a = d1*.2/100;
b = d2*.2/100;
[xm ym] = meshgrid(linspace(0,a,d1),linspace(0,b,d2));
pcolor(xm,ym,speed');shading interp
hold on

plot(p(:,1),p(:,2),'ro','markersize',8)
hold on;
triplot(t,p(:,1),p(:,2),'k','linewidth',2)



return

%% box mesh
global xm ym himg 

load PAT/PAT_setup.mat
d1 = size(speed,1); d2 = size(speed,2);
a = d1*.2/100;
b = d2*.2/100;
[xm ym] = meshgrid(linspace(0,a,d1),linspace(0,b,d2));

himg = speed/min(speed(:));
svals = sort(uniquetol(himg),'descend');
for i = 1:length(svals)
    ids = abs(himg-svals(i))<1e-8;
    if (i==1)        
        himg(ids) = 1;2*himg(ids);
    else
        himg(ids) = 1;
    end
end

for iter = 1:0
    savg = himg;
    for i = 2:d1-1
        for j = 2:d2-1
            sij = himg(i,j) + ...
                himg(i+1,j) + himg(i-1,j) + himg(i,j+1) + himg(i,j-1) + ...
                himg(i+1,j+1) + himg(i-1,j+1) + himg(i+1,j-1) + himg(i-1,j-1);
            savg(i,j) = sij/9;
        end
    end
    himg = savg;
end

node = [
    0   0
    a   0
    a   b
    0   b
    ];

hmin = .06;
hdata.fun = @hfun;
hdata.args = {0,a,0,b,hmin};
[p,t] = mesh2d(node,[],hdata);

h = pcolor(xm,ym,speed');shading interp;axis on
plot(p(:,1),p(:,2),'b.','markersize',1)
hold on;
triplot(t,p(:,1),p(:,2))
% for e = 1:size(t,1)
% %     plot(p([t(e,:) t(e,1)],1),p([t(e,:) t(e,1)],2),'k-')
% end
% patch('faces',edge,'vertices',node,'facecolor','none','edgecolor','k')
% h = hfun(xm,ym);
keyboard

function h = hfun(x,y,x1,x2,y1,y2,hmin)

global xm ym himg
c2 = interp2(xm,ym,himg',x,y);

h = hmin.*(1./c2);

in = (x>=x1)&(x<=x2)&(y>=y1)&(y<=y2);
h(~in) = hmin;


