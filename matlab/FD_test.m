N = 200;
x = linspace(-1,1,N)';
h = 2/N;
dt = .5*h;

e = ones(N,1);
A = diag(2*e) - diag(e(2:end),1) - diag(e(2:end),-1);
A = sparse(A);
% A(1,1) = 1; A(end,end) = 1;

p0 = @(x) exp(-10^2*(x).^2);
% p0 = @(x) .5*(exp(-10^2*(x-.5).^2) + exp(-10^2*(x+.5).^2));
p = p0(x);
pp = p;

% c2 = 1 + (x > 0);
% c2 = 1+25*(x.^2+y.^2 < .125);
c2 = 1;

time = 0;
FinalTime = 4;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;
for i = 1:Nsteps

    %   p(k+1)-2*p(k)+p(k-1) + dt^2*A*p(k) = 0
    
    ptmp = 2*p-pp - c2.*((dt/h)^2*A*p);
    pp = p;
    p = ptmp;
    if mod(i,25)==0
        clf
        plot(x,p,'.-')               
        axis([-1 1 -1 1])
        drawnow        
    end
end
clf
plot(x,p,'.')               
hold on
plot(x,p0(x))



%%
N = 200;
[x y] = meshgrid(linspace(-1,1,N));
h = 2/N;
dt = .5*h;

e = ones(N,1);
A = diag(2*e) - diag(e(2:end),1) - diag(e(2:end),-1);
A = sparse(A);

p = exp(-25^2*((x).^2+y.^2));
% p = abs(x)<.5 & abs(y) <.5;

pp = p;

%c2 = 1 + (x > 0);
% c2 = 1+25*(x.^2+y.^2 < .125);
c2 = 1;

time = 0;
FinalTime = 2;
Nsteps = ceil(FinalTime/dt);
dt = FinalTime/Nsteps;
figure
for i = 1:Nsteps
    
    ptmp = 2*p-pp - c2.*((dt/h)^2*(A*p + p*A));
    pp = p;
    p = ptmp;
    if mod(i,10)==0
        clf
        %         color_line3(x,y,p,p,'.')
        surf(x,y,p);shading interp
%                 axis([-1 1 -1 1 -1 1])
                view(2)
                axis equal;   axis tight
%         [~,id] = min(abs(y(:,1)));
%         plot(x(id,:),p(id,:),'.-')
        
        drawnow
    end
end
max(p(:))
po = p;

%%
p = po;
A = (1/3)*(diag(e) + diag(e(2:end),1) + diag(e(2:end),-1));
for i = 1:25
    p = .5*(A*p + p*A);
end
figure
[~,id] = min(abs(y(:,1)));
        plot(x(id,:),p(id,:),'.-')
        return
surf(x,y,p) 
shading interp 
axis equal
axis tight

view(2)

