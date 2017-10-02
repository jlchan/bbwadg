load PAT/PAT_setup.mat

d1 = size(speed,1); d2 = size(speed,2);
a = d1*.2*.001;
b = d2*.2*.001;
[xm ym] = meshgrid(linspace(0,a,d1),linspace(0,b,d2));

c = index;
c(index~=0) = 1;

% pcolor(xm,ym,c');shading interp

bw2 = bwperim(c);
[id] = find(bw2');

% imagesc(bw2')
% hold on
% grid on
% ve = linspace(-1,1,length(id));
% for i = 1:length(id)
%     plot(xm(id(i)),ym(id(i)),'o')
%     drawnow
% end
% colorbar

xb = xm(id); yb = ym(id);
[~,ia,ic] = uniquetol(xb);
x0 = mean(xb);
y0 = mean(yb);
t = atan2(yb-y0,xb-x0);
plot(t,xb,'x')
hold on
plot(t,yb,'.')

save PAT/PAT_breast_boundary.mat t xb yb