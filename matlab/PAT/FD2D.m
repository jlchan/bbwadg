% grid resolution
N = 500;
[x y] = meshgrid(linspace(-1,1,N));
h = 2/N; % grid spacing

e = ones(N,1);
Ax = diag(2*e) - diag(e(2:end),1) - diag(e(2:end),-1);
Ax = sparse(Ax);

Ay = Ax;

p = exp(-25^2*(x.^2+y.^2));

% previous p
pp = p;

c2 = 1 + 5*((x-.1).^2+y.^2 < .125);
% c2 = 1;

dt = .75*h;
dt = dt/max(c2(:));

time = 0;
FinalTime = 2;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;
for i = 1:Nsteps
    
    ptmp = 2*p-pp - c2.*((dt/h)^2*(Ax*p + p*Ay));
    pp = p;
    p = ptmp;
    
    % plot every 10th timestep
    if mod(i,10)==0
        clf
        surf(x,y,p);
        shading interp
        view(2)
        axis equal;   
        axis tight        
        drawnow
    end
end
