N = 20;

Lx = 15;
Ly = 20;
Lz = 5;
for a = 0:.01:1/5
    [x y z] = meshgrid(linspace(0,Lx,2*N),linspace(0,Ly,3*N),linspace(0,Lz,N));
    x = x(:);
    y = y(:);
    z = z(:);
    
    xx = (x-Lx/2)/Lx;
    yy = (y-Ly/2)/Ly;
    zz = 0*(z-Lz/2)/Lz;
    
    dy = cos(3*pi*xx).*cos(pi*yy).*cos(pi*zz);
    y = y + a*Ly*dy;
    
    yy = (y-Ly/2)/Ly;
    dx = cos(pi*xx).*sin(4*pi*yy).*cos(pi*zz);
    x = x + a*Lx*dx;
    
    xx = (x-Lx/2)/Lx;
    dz = cos(pi*xx).*cos(2*pi*yy).*cos(pi*zz);
%     z = z + a*Lz*dz;

    plot3(x,y,z,'o')
    axis equal
    view(90,-90)
%     view(2)
    
    drawnow
    
end

%% 2D vortex warping

N = 40;
Lx = 20;
Ly = 10;

[x y] = meshgrid(linspace(0,Lx,2*N),linspace(0,Ly,N));

a = 1/8;

x0 = Lx/2; y0 = Ly/2;

xx = (x-x0)/Lx;
yy = (y-y0)/Ly;
dx = cos(pi*xx).*cos(3*pi*yy);
x = x + Lx*a*dx;

xx = (x-x0)/Lx;
dy = sin(4*pi*xx).*cos(1*pi*yy);
y = y + Ly*a*dy;

figure
%mesh(x,y,0*x)
plot(x,y,'bo')
view(2)
axis equal
