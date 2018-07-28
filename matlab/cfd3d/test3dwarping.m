N = 15;

for a = 0:.025:1.0
    [x y z] = meshgrid(linspace(0,10,N),linspace(0,20,2*N),linspace(0,10,N));
    x = x(:);
    y = y(:);
    z = z(:);
    
    dx = cos(1/2*pi*(x-5)/5).*cos(3/2*pi*(y-10)/10).*cos(1/2*pi*(z-5)/5);
    dy = cos(3/2*pi*(x-5)/5).*cos(1/2*pi*(y-10)/10).*cos(3/2*pi*(z-5)/5);
    dz = cos(1/2*pi*(x-5)/5).*cos(3/2*pi*(y-10)/10).*cos(1/2*pi*(z-5)/5);
%     d = sin(pi*(x-5)/5).*sin(2*pi*(y-10)/10).*sin(pi*(z-5)/5);
    x = x + a*dx;
    y = y + a*dy;
    z = z + a*dz;
    plot3(x,y,z,'o')
    axis equal
    view(-120,10)
    
    drawnow
    
end