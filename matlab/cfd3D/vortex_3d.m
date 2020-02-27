x = linspace(0,15,60);
y = linspace(0,20,60);
z = linspace(0,5,20);
[x y z] = meshgrid(x,y,z);
bids = find(abs(x)<1e-10 | abs(y)<1e-10);

x0 = x; y0 = y; z0 = z;

aa = .125;
dx = aa*sin(1*pi*(x-7.5)/7.5);
dy = aa*sin(2*pi*(y-10)/10);
dz = aa*sin(pi*(z-2.5)/2.5);

gamma = 1.4;
c1 = 7.5;
c2 = 7.5;
t = 5;
p0 = 1/gamma;
PImax = .4;
PI = PImax.*exp(1-.5*((-(y-c2-t)).^2+(x-c1).^2));

rho = (1-(gamma-1)/2*PI.^2).^(1./(gamma-1));
u = PI.*(-(y-c2-t));
v = PI.*(x-c1);
w = 0*x;
E = p0/(gamma-1).*(1-(gamma-1)/2*PI.^2).^(gamma/(gamma-1))+rho/2.*(u.*u+v.*v+w.*w);

vv = rho;
% vv = dx.*dy.*dz;
color_line3(x,y,z,vv,'.')
 % ids = (y < 10)&(z < 4);
% color_line3(x(ids),y(ids),z(ids),vv(ids),'.')
colorbar
view(3)
xlabel('x')
ylabel('y')
axis equal
norm(u(bids))

